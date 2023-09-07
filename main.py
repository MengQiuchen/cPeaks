# Adapted from cPeaks https://github.com/MengQiuchen/cPeaks/blob/main/main.py

import os
import gzip
from tqdm import tqdm
from collections import defaultdict
from joblib import Parallel, delayed
import numpy as np
import pandas as pd

# import dask.dataframe as dd
import argparse
from time import time


def load_cpeaks(reference):
    """
    Format:
    chromosome  start   end
    """
    print("Info: Loading cPeaks...", end=" ", flush=True)
    path, filename = os.path.split(os.path.realpath(__file__))
    if reference == "hg38":
        cPeaks_path = os.path.join(path, "cPeaks_hg38.bed")
    else:
        cPeaks_path = os.path.join(path, "cPeaks_hg19.bed")

    try:
        cpeaks = pd.read_csv(cPeaks_path, delimiter="\t", names=["chr", "start", "end"])
    except:
        raise "There is something wrong with the cPeaks file!"
    print("done.", flush=True)

    # Order of chromosomes
    chrom_order = [f"chr{chrom}" for chrom in list(range(1, 23)) + ["X", "Y"]]

    cpeaks["chr"] = pd.Categorical(cpeaks["chr"], categories=chrom_order)
    cpeaks = cpeaks.sort_values(["chr", "start", "end"])

    print(f"Info: Loaded {len(cpeaks)} cPeak regions.", flush=True)
    return cpeaks.values


def prepare_cpeaks(cpeaks_starts, cpeaks_ends, min_bp_overlap):
    # Calculate half region length
    half_region_length = np.floor((cpeaks_ends - cpeaks_starts) * 0.5).astype(int)

    # cPeak correction value: min(half_region_length, min_bp_overlap)
    cpeak_correction = np.min(
        (half_region_length, np.full_like(half_region_length, min_bp_overlap)), axis=0
    )

    # For fragment mapping, the cPeak regions are shortened on both ends by a given number
    # of base pairs in order to ensure a minimal overlap of region and fragment
    #
    #            + correction value      - correction value
    #                |'''''''|               |'''''''|
    #
    #              cPeak                           cPeak
    #              start     |               |      end
    #    ------------|'''''''#################'''''''|-----------------
    #                    s--------e          |                            Overlap
    #           s---------e  |               |                            No overlap
    #                        |               |  s----------e              No overlap
    #               s-----------------------e|                            Overlap
    #     s------e           |               |                            No overlap
    #                        |               |           s-------e        No overlap

    # Starts + value ; Ends - value
    cpeaks_starts += cpeak_correction
    cpeaks_ends -= cpeak_correction

    return cpeaks_starts, cpeaks_ends


def load_fragments_df(
    file_path, chromosomes, min_frag_len, max_frag_len, chunksize=100000
):
    """
    Function to read the TSV file chunk by chunk and apply the progress bar
    """
    chunks = pd.read_csv(
        file_path,
        sep="\t",
        compression="gzip",
        header=None,
        names=["chr", "start", "end", "barcodes"],
        usecols=[0, 1, 2, 3],
        chunksize=chunksize,
        comment="#",
    )

    fragments_df = pd.concat(
        tqdm(
            chunks,
            desc="Info: Loading fragments",
            unit=" fragments",
            unit_scale=chunksize,
        ),
        ignore_index=True,
    )

    # Only keep fragments for the chromosomes of cPeaks (1-22 as well as X and Y)
    fragments_df["chr"] = pd.Categorical(
        fragments_df["chr"], categories=chromosomes, ordered=True
    )
    fragments_df.dropna(subset=["chr"], inplace=True)

    fragments_df["length"] = fragments_df["end"] - fragments_df["start"]
    quantiles_unfiltered = (fragments_df["length"]).quantile(q=np.arange(0, 1.1, 0.1))

    print(
        f"Info: Removing fragments shorter than {min_frag_len}bp or longer than {max_frag_len}bp."
    )
    fragments_df = fragments_df.loc[fragments_df["length"] >= min_frag_len]
    fragments_df = fragments_df.loc[fragments_df["length"] <= max_frag_len]

    quantiles_filtered = (fragments_df["length"]).quantile(q=np.arange(0, 1.1, 0.1))

    print(" Fragment length:")
    print(
        " Quantile:",
        "\t".join(np.round(quantiles_unfiltered.index, 2).astype(str)),
        sep="\t",
    )
    print(
        " Unfiltered:",
        "\t".join(np.round(quantiles_unfiltered, 0).astype(int).astype(str)),
        sep="\t",
    )
    print(
        " Filtered",
        "\t".join(np.round(quantiles_filtered, 0).astype(int).astype(str)),
        sep="\t",
        flush=True,
    )

    return fragments_df


def load_barcodes(barcode_path):
    print("Info: Loading barcodes...", end=" ", flush=True)
    if barcode_path[-3:] == ".gz":
        barcodes = [
            barcode.rstrip("\n")
            for barcode in gzip.open(barcode_path, "rt").readlines()
        ]
    else:
        barcodes = [
            barcode.rstrip("\n") for barcode in open(barcode_path, "r").readlines()
        ]
    print("done.", flush=True)
    # Create a set to remove duplicates while preserving order
    unique_barcodes = set()
    barcodes = [
        barcode
        for barcode in barcodes
        if not (barcode in unique_barcodes or unique_barcodes.add(barcode))
    ]
    return barcodes


def map_fragments_to_cpeaks(fragments, chr_indices, cpeaks_starts, cpeaks_ends):
    """
    Map fragment bed to cPeak bed.
    The output is a tuple of two numpy arrays: (indices, values)
    """
    overlap_counts = np.zeros(len(cpeaks_starts), dtype=int)

    for chr, s, e, barcode in fragments:
        #          chr     start   end          barcode
        # array([['chr1', 181500, 181531, 'AAACAGCCAACCCTAA-1'],
        #        ['chr1', 629913, 629992, 'AAACAGCCAACCCTAA-1'],
        #        ['chr1', 629914, 629985, 'AAACAGCCAACCCTAA-1'],
        #        ['chr1', 629914, 629986, 'AAACAGCCAACCCTAA-1'],
        #        ['chr1', 629939, 629986, 'AAACAGCCAACCCTAA-1'],
        #        ['chr1', 633935, 634098, 'AAACAGCCAACCCTAA-1'],
        #        ['chr1', 633957, 634087, 'AAACAGCCAACCCTAA-1'],
        #        ... ])

        # Get the relevant indices for the current chromosome
        chr_ind = chr_indices[chr]

        # Use binary search to find the start and end indices in cpeaks
        # (-> only check interesting regions)
        #
        # Example fragment on chromosome 1:
        # start = 14498 ; end = 14673

        start_ind = np.searchsorted(cpeaks_starts[chr_ind], e, side="left")
        #                                                  fragment ends after this point (start < e)
        #                                                  ↓
        # start: [ 9919, 11043, 12173, 13249, 14012, 14346, 15403, 16415, 16800, 17352, 19818, 20668, 21057, 23651, ...]
        # ends:  [10727, 11422, 12922, 13832, 14345, 14845, 16233, 16606, 17224, 17638, 20668, 21057, 22150, 24150, ...]
        #                                           ↑
        #                                           fragment starts before this point (s < end)
        end_ind = np.searchsorted(cpeaks_ends[chr_ind], s, side="right")

        #              cPeak          cPeak
        #              start           end
        #    ------------|##############|-----------------
        #                |  s--------e  |                   Overlap
        #              s---------e      |                   Overlap
        #                |      s----------e                Overlap
        #          s-----------------------e                Overlap
        #     s------e   |              |                   No overlap
        #                |              |    s-------e      No overlap
        #
        # Two possible cases:
        # start_ind > end_ind -> fragment overlaps with at least one region -> len(overlap_ind) > 0
        # start_ind == end_ind -> fragment lays between regions without overlap -> len(overlap_ind) == 0

        overlap_ind = chr_ind[end_ind:start_ind]

        if len(overlap_ind) > 0:
            overlap_counts[overlap_ind] += 1

    non_zero_indices = np.nonzero(overlap_counts)[0]
    non_zero_values = overlap_counts[non_zero_indices]

    return (non_zero_indices, non_zero_values)


# def function to trans fragment file to cPeaks referenced mtx
def frag2mm(fragments_df, cpeaks, num_cores):
    """
    fragmen_path: path of fragment.gz file
    savepath: path to save mtx file
    barcode_path: path of barcode file, ,txt format,each line is a barcode; if None, use all barcode in fragment.gz file
    """

    print("Info: Preparing fragments for mapping...", end=" ", flush=True)
    fragments_df.sort_values(["barcodes", "chr", "start", "end"], inplace=True)
    fragments_df.reset_index(drop=True, inplace=True)
    fragments_df = fragments_df.groupby("barcodes", sort=False)
    print("done.", flush=True)

    chr_indices, cpeaks_starts, cpeaks_ends = cpeaks

    results = Parallel(n_jobs=num_cores, backend="multiprocessing")(
        delayed(map_fragments_to_cpeaks)(
            barcode_subset_df.values[:, :4],
            chr_indices,
            cpeaks_starts,
            cpeaks_ends,
        )
        for _, barcode_subset_df in tqdm(
            fragments_df, desc="Mapping to cPeaks", unit=" barcodes"
        )
    )
    return results


def save_as_mm(savepath, results, n_features):
    n_cells = len(results)
    n_values = sum([len(result[0]) for result in results])

    print("Info: Start to write mtx file...", end=" ", flush=True)
    with gzip.open(os.path.join(savepath, "atac_matrix.mtx.gz"), "wb") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n".encode())
        f.write(f"{n_features} {n_cells} {n_values}\n".encode())
        for barcode_index, (feature_indices, values) in enumerate(results):
            for feature_index, value in zip(feature_indices, values):
                f.write(f"{feature_index + 1} {barcode_index + 1} {value}\n".encode())
    print("done.", flush=True)


def main(args):
    fragment_path = args.fragment_path
    barcode_path = args.barcode_path
    savepath = args.output
    save_type = args.type_saved
    num_cores = args.num_cores
    reference = args.reference
    min_bp_overlap = args.min_bp_overlap
    min_frag_len = args.min_frag_len
    max_frag_len = args.max_frag_len

    # Check if output directory exists
    if not os.path.exists(savepath):
        print("Info: Creating", savepath, flush=True)
        os.makedirs(savepath)

    if all(
        [
            os.path.isfile(os.path.join(savepath, file))
            for file in ["features.tsv.gz", "barcodes.tsv.gz", "atac_matrix.mtx.gz"]
        ]
    ):
        print("Warning: Sample has already been fully processed.")
        ignore = input("Do you want to continue? (yes/no) ")
        if ignore.lower() in ["yes", "y"]:
            print(f"Info: Overwriting non-empty directory {savepath}.")
        else:
            exit()

    ### cPeaks ###
    # Load cPeaks regions
    cpeaks = load_cpeaks(reference)

    # Save cPeaks regions to the output directory
    features_file = os.path.join(savepath, "features.tsv.gz")
    if not os.path.isfile(features_file):
        print("Info: Saving cPeaks features...", end=" ", flush=True)
        with gzip.open(features_file, "wt") as gzip_file:
            np.savetxt(gzip_file, cpeaks, delimiter="\t", fmt="%s")
        print("done.", flush=True)
    else:
        print("Warning: features.tsv.gz already exists in " + savepath, flush=True)

    # Create a dictionary to store the indices of each chromosome in cpeaks
    chr_indices = defaultdict(list)
    for i, (chr, _, _) in enumerate(cpeaks):
        chr_indices[chr].append(i)
    chr_indices = {chr: np.array(indices) for chr, indices in chr_indices.items()}
    cpeaks_starts = cpeaks[:, 1].astype(int)
    cpeaks_ends = cpeaks[:, 2].astype(int)
    cpeaks = (chr_indices, cpeaks_starts, cpeaks_ends)

    ### Fragments ###
    # Check fragment_path: does not end with 'tsv.gz' and is not None
    if fragment_path[-6:] != "tsv.gz":
        raise ("Expected a tsv.gz file but got " + fragment_path)

    # Read in the fragment file as a pandas dataframe
    fragments_df = load_fragments_df(
        fragment_path,
        chromosomes=chr_indices.keys(),
        min_frag_len=min_frag_len,
        max_frag_len=max_frag_len,
    )

    ### BARCODES ###
    # Get barcodes
    if barcode_path is None:
        print("Info: Start to get barcodes from fragment file", flush=True)
        barcodes = fragments_df["barcodes"].unique()
    else:
        barcodes = load_barcodes(barcode_path)
        # Only use barcodes for which fragments are available
        barcode_intersect = set(fragments_df["barcodes"]) & set(barcodes)

        barcodes = [barcode for barcode in barcodes if barcode in barcode_intersect]

    fragments_df["barcodes"] = pd.Categorical(
        fragments_df["barcodes"], categories=barcodes, ordered=True
    )
    fragments_df.dropna(subset=["barcodes"], inplace=True)

    print(f"Info: Loaded {len(barcodes)} barcodes.", flush=True)

    # Print number of peaks for the last three barcodes
    num_peaks = [
        len(fragments_df.loc[fragments_df["barcodes"] == barcode, :].index)
        for barcode in barcodes[-3:]
    ]
    for n_peaks in num_peaks:
        print(f"Info: Random barcode has {n_peaks} fragments.", flush=True)

    # Save barcodes to the output directory
    barcodes_file = os.path.join(savepath, "barcodes.tsv.gz")
    if not os.path.isfile(barcodes_file):
        print("Info: Saving barcodes...", end="", flush=True)
        with gzip.open(barcodes_file, "wt") as gzip_file:
            np.savetxt(gzip_file, np.array(barcodes), delimiter="\t", fmt="%s")
        print("done.", flush=True)
    else:
        print("Warning: barcodes.tsv.gz already exists in " + savepath, flush=True)

    ### PREPARE cPEAKS ###
    print(
        f"Info: Setting the minimal overlap of region and fragment to {min_bp_overlap}bp."
    )
    cpeaks_starts, cpeaks_ends = prepare_cpeaks(
        cpeaks_starts=cpeaks_starts,
        cpeaks_ends=cpeaks_ends,
        min_bp_overlap=min_bp_overlap,
    )

    ### MAP FRAGMENTS TO cPEAKS ###
    if save_type == "mtx":
        results = frag2mm(fragments_df, cpeaks, num_cores)
    else:
        raise ("save_type must be mtx, this function is not finished yet")

    ### SAVE MATRIX MARKET ###
    save_as_mm(savepath, results, n_features=len(cpeaks[1]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fragment_path", "-f", type=str, default=None, help="path of fragment file"
    )
    parser.add_argument(
        "--barcode_path", "-b", type=str, default=None, help="path of barcode file"
    )
    parser.add_argument(
        "--output", "-o", type=str, default="./res", help="path to save files"
    )
    parser.add_argument(
        "--type_saved", "-t", type=str, default="mtx", help="save type, h5ad or mtx"
    )
    parser.add_argument(
        "--num_cores", "-n", type=int, default=1, help="num of cores to use"
    )
    parser.add_argument(
        "--reference", type=str, default="hg38", help="cPeaks version: hg38 or hg19"
    )
    parser.add_argument(
        "--min_frag_len", type=int, default=20, help="Minimal allowed fragment length"
    )
    parser.add_argument(
        "--max_frag_len", type=int, default=1000, help="Maximal allowed fragment length"
    )
    parser.add_argument(
        "--min_bp_overlap",
        type=int,
        default=60,
        help="Minimal fragment overlap with region",
    )
    # "Select the region between 175-1000bp to determine the average fragment size of a Single Cell ATAC library. Lower molecular weight (≤ 150 bp) and/or a high molecular weight (~2,000 bp) product may be present. These do not affect sequencing and should not be included when determining average fragment size."
    # https://kb.10xgenomics.com/hc/en-us/articles/9817838634253-How-do-I-determine-the-average-fragment-size-of-a-Single-Cell-ATAC-or-Multiome-ATAC-final-library-

    args = parser.parse_args()

    start_time = time()

    main(args)

    runtime = time() - start_time
    hours = int(runtime // 3600)
    minutes = int((runtime % 3600) // 60)
    seconds = runtime % 60

    formatted_runtime = f"{hours}:{minutes:02d}:{seconds:.02f}"

    print("Info: Runtime:", formatted_runtime)
