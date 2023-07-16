# python 
# 2023/07/15
# version 1.1
# @auther : Xinze Wu

import os
import sys
import gzip
from tqdm import tqdm
from collections import OrderedDict
from joblib import Parallel, delayed
import numpy as np
import scanpy as sc
import anndata
import argparse
import pandas as pd


# map a bed to b bed, the output is a tuple of two numpy arrays: (indices, values)
def map_A_to_B(a_peaks):
    
    chromosomes = list(range(1, 23))
    chromosomes.extend(["X", "Y"])
    chrlist = ["chr" + str(chrom) for chrom in chromosomes]
    
    a_dict = {}
    for chrom in chrlist:
        a_dict[chrom] = []
    
    for peak in a_peaks:
        chrom = peak[0]
        start = peak[1]
        end = peak[2]
        if chrom in a_dict:
            a_dict[chrom].append((start, end))
    
    overlap_counts = np.zeros(len(b_peaks), dtype=int)
    for i, peak in enumerate(b_peaks):
        chrom = peak[0]
        start = peak[1]
        end = peak[2]
        if chrom not in a_dict:
            continue
        
        for a_peak in a_dict[chrom]:
            if a_peak[0] < end and start < a_peak[1]:
                overlap_counts[i] += 1
          
    non_zero_indices = np.nonzero(overlap_counts)[0]
    non_zero_values = overlap_counts[non_zero_indices]
    
    return (non_zero_indices+1, non_zero_values)


# def function to trans fragment file to cPeaks referenced mtx
def frag2mtx(fragment_path,savepath,barcode_path,num_cores):
    '''
    fragmen_path: path of fragment.gz file
    savepath: path to save mtx file
    barcode_path: path of barcode file, ,txt format,each line is a barcode; if None, use all barcode in fragment.gz file
    '''
    fragment_name = fragment_path.split('/')[-1]
    # read fragment.gz file by open
    fragment_data = gzip.open(fragment_path, 'rb')
    
    # get unique barcode
    if barcode_path is None:
        print('start to get barcode')

        barcodes = set()

        for line in fragment_data:
            line = line.decode()
            barcodes.add(line.split('\t')[3])
        print('num of barcode is {}'.format(len(barcodes)))
        # write barcode to file
        
        with open(os.path.join(savepath,'barcode.txt'), 'w') as f:
            for barcode in barcodes:
                f.write(barcode+'\n')

    else:
        barcodes = set(open(barcode_path, 'r').readlines())
        barcodes = [i.strip() for i in barcodes]
        print('num of barcode is {}'.format(len(barcodes)))
        
    barcodes = list(barcodes)
    print(barcodes[:10])
    # predict the time needed for trans

    # get the bed in each barcode
    print('start to get bed in each barcode')

    fragment_data = gzip.open(fragment_path, 'rb')
    
    bar2bed = OrderedDict()
    for barcode in barcodes:
        bar2bed[barcode] = []
    
    for line in tqdm(fragment_data):
        temp = line.decode().strip().split('\t')
        try:
            bar2bed[temp[3]].append((temp[0],int(temp[1]),int(temp[2])))
        except:
            continue
    print('num peak of random bar: ',len(bar2bed[barcodes[-1]]))
    print('num peak of random bar: ',len(bar2bed[barcodes[-2]]))
    print('num peak of random bar: ',len(bar2bed[barcodes[-3]]))


    print('mapping to cPeaks...')

    results = Parallel(n_jobs=num_cores, backend="multiprocessing")(
        delayed(map_A_to_B)(input_list) for input_list in tqdm(list(bar2bed.values()))
    )
    
    value_num = 0
    for res in results:
        value_num += len(res[0])
    
    print('the num of the value of mtx is: ',value_num)
    
    print('start to write mtx file')
    with open(os.path.join(savepath,fragment_name+'.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write(f'{len(b_peaks)} {len(barcodes)} {value_num}\n')
        
        for i,res in enumerate(results):
            for j in range(len(res[0])):
                f.write(f'{res[0][j]} {i+1} {res[1][j]}\n')
    print('done')


# mtx, cPeaks, barcode to h5ad
def mtx2h5ad(mtx_path, cPeaks_path, barcode_path, savepath):
    '''
    mtx_path: path of mtx file
    cPeaks_path: path of cPeaks file
    barcode_path: path of barcode file
    savepath: path to save h5ad file
    '''
    # read mtx file
    mtx = pd.read_csv(mtx_path, sep=' ', skiprows=2, header=None)
    mtx.columns = ['cPeaks', 'barcode', 'value']
    
    # read cPeaks file
    cPeaks = open(cPeaks_path, 'r').readlines()
    cPeaks = ['_'.join(i.strip().split(' ')) for i in cPeaks]

    # read barcode file
    barcode = open(barcode_path, 'r').readlines()
       
    # get the sparse matrix
    mtx_sparse = np.zeros((len(cPeaks), len(barcode)))
    for i in range(len(mtx)):
        mtx_sparse[mtx.iloc[i, 0]-1, mtx.iloc[i, 1]-1] = mtx.iloc[i, 2]
    
    # save to h5ad file
    adata = anndata.AnnData(X=mtx_sparse, obs=barcode, var=cPeaks)
    adata.write(savepath)
    
    print('done')

# map bed to bed 
def map_bed_to_bed(a_bed_path, b_peaks,savepath):

    """
    a_bed_path: path of a bed file
    b_bed: path of cpeak bed
    savepath: path to save result
    """ 

    a_peaks = open(a_bed_path, 'r').readlines()
    a_peaks = [i.strip().split('\t') for i in a_peaks]

    chromosomes = list(range(1, 23))
    chromosomes.extend(["X", "Y"])
    chrlist = ["chr" + str(chrom) for chrom in chromosomes]
    
    a_dict = {}
    for chrom in chrlist:
        a_dict[chrom] = []
    
    for peak in a_peaks:
        chrom = peak[0]
        start = peak[1]
        end = peak[2]
        if chrom in a_dict:
            a_dict[chrom].append((start, end))
    
    overlap_counts = np.zeros(len(b_peaks), dtype=int)
    for i, peak in tqdm(enumerate(b_peaks)):
        chrom = peak[0]
        start = peak[1]
        end = peak[2]
        if chrom not in a_dict:
            continue
        
        for a_peak in a_dict[chrom]:
            if int(a_peak[0]) < end and start < int(a_peak[1]):
                overlap_counts[i] += 1
    # get the index of non-zero value
    non_zero_indices = np.nonzero(overlap_counts)[0]
    # map to cPeaks and write to bed
    with open(savepath, 'w') as f:
        for i in non_zero_indices:
            f.write(b_peaks[i][0]+'\t'+str(b_peaks[i][1])+'\t'+str(b_peaks[i][2])+'\n')
    print('done')

# main 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragment_path", '-f',type=str, default=None,help="path of mtx file")
    parser.add_argument("--barcode_path", '-b',type=str, default=None,help="path of barcode file")
    parser.add_argument("--output", '-o',type=str, default='./res',help="path to save files")
    parser.add_argument("--type_saved", '-t',type=str, default='mtx',help="save type, h5ad or mtx")
    parser.add_argument("--num_cores",'-n', type=int, default= 1,help="num of cores to use")
    parser.add_argument("--bed_path", '-bed',type=str, default=None,help="bed file you called before")
    parser.add_argument("--reference",type=str, default='hg38',help="cPeaks version: hg38 or hg19")

    args = parser.parse_args()
   
    from time import time
    time_ = time()

    fragment_path = args.fragment_path
    barcode_path = args.barcode_path
    savepath = args.output
    save_type = args.type_saved
    num_cores = args.num_cores
    bed_path = args.bed_path
    reference = args.reference

    if reference == 'hg38':
        cPeaks_path = 'cPeaks_hg38.bed'
    else:
        cPeaks_path = 'cPeaks_hg19.bed'

    # get cPeaks
    try:
        with open(cPeaks_path, 'r') as b_file:
            b_lines = b_file.readlines()
        b_peaks = []
        for line in b_lines:
            fields = line.strip().split('\t')
            b_peaks.append((fields[0], int(fields[1]), int(fields[2])))  
    except:
        raise('there is sth wrong with cPeaks file')
    
    if not os.path.exists(savepath):
        os.makedirs(savepath)


    if bed_path is not None:
        print('attention: you use bed file you called before')
        map_bed_to_bed(bed_path, b_peaks,os.path.join(savepath,'res.bed'))
        sys.exit(0)



    # fragment_path not end with 'tsv.gz' and not None
    if fragment_path[-6:] != 'tsv.gz':
        print(fragment_path[-6:] )
        raise('fragment_path not end with tsv.gz')
    

    if save_type == 'mtx':
        frag2mtx(fragment_path,savepath,barcode_path,num_cores)
    else:
        raise('save_type must be mtx, this function is not finished yet')
        # frag2mtx(fragment_path,savepath,barcode_path,num_cores)
        # mtx2h5ad(os.path.join(savepath,fragment_path+'.mtx'), cPeaks_path, barcode_path, savepath)

    print('use time: ',time()-time_)
