Last updated: 8 October 2026

Source: [docs/Tutorials.md](https://github.com/MengQiuchen/cPeaks/blob/main/Tutorials/docs/Tutorials.md)

## 1. Overview of cPeaks

### Introduction

cPeaks serves as a unified reference for different cell types for ATAC-seq or scATAC-seq data, enhancing downstream analysis, particularly in cell annotation and rare cell type detection.

<img src=".\media\introduction.png" alt="1" width="500" style="zoom:100%;" />

### Download

#### Simple data

cPeaks reference files are available as ZIP archives from Zenodo: [hg19 file](https://zenodo.org/records/23241189/files/cpeaks_hg19_sorted.bed.zip?download=1) and [hg38 file](https://zenodo.org/records/23241189/files/cpeaks_hg38_sorted.bed.zip?download=1). You can download and extract them using the following commands:

```bash
wget -O YOUR_PATH/cpeaks_hg19_sorted.bed.zip 'https://zenodo.org/records/23241189/files/cpeaks_hg19_sorted.bed.zip?download=1'
wget -O YOUR_PATH/cpeaks_hg38_sorted.bed.zip 'https://zenodo.org/records/23241189/files/cpeaks_hg38_sorted.bed.zip?download=1'
unzip -p YOUR_PATH/cpeaks_hg19_sorted.bed.zip cpeaks_hg19_sorted.bed > YOUR_PATH/cpeaks_hg19.bed
unzip -p YOUR_PATH/cpeaks_hg38_sorted.bed.zip cpeaks_hg38_sorted.bed > YOUR_PATH/cpeaks_hg38.bed
```

We will utilize the two downloaded files in the upcoming tutorials.

#### Detailed data

For additional details, download [Detailed_cPeaks.tsv.zip from Zenodo](https://zenodo.org/records/23241189/files/Detailed_cPeaks.tsv.zip?download=1) ([Google Drive mirror](https://drive.google.com/file/d/1R_cY1ourPLEyohG65ZYhrUWOYQXtOPSd/view?usp=drivesdk)). The ZIP archive contains the corrected `Detailed_cPeaks.tsv` file with basic information, cPeak annotations, and integration with biological data. It replaces `Information of all cPeaks.tsv.zip` from earlier Zenodo versions. The [version 3 record](https://zenodo.org/records/23241189) includes a README with file descriptions and checksums. The detailed table contains 1,657,194 cPeaks and 14 columns (ZIP SHA-256: `9438096b09f4243a3e99f6d3752c8e9bc5788dfa71802208a4c71250169ed120`):

| Column Name | Column Description |
| ----------- | ------------------ |
| ID | A unique cPeak identifier beginning with “CPHS” (human), followed by nine digits. For example, the first cPeak is encoded as “CPHS000000001”. |
| source | Indicates the origin of the cPeak, either “observed” or “predicted”. |
| chr_hg38 | The chromosome where this cPeak locates in the hg38 reference genome. |
| start_hg38 | The start position of this cPeak in the hg38 reference genome. |
| end_hg38 | The end position of this cPeak in the hg38 reference genome. |
| housekeeping | Specifies whether the cPeak is accessible across nearly all datasets (“TRUE” or “FALSE”). |
| shape pattern | The shape pattern of the cPeak, categorized as “well-positioned”, “asymmetrically-positioned”, or “weakly-positioned”. |
| inferredElements | The inferred regulatory elements associated with the cPeak, such as “CTCF”, “TES”, “TSS”, “Enhancer” or “Promoter”. |
| chr_hg19 | The chromosome where this cPeak locates in the hg19 reference genome. |
| start_hg19 | The start position of this cPeak in the hg19 reference genome. |
| end_hg19 | The end position of this cPeak in the hg19 reference genome. |
| cDHS_ID | The ID of the overlapped cDHS region in the hg38 reference, formatted as “chr_start_end”. |
| CATLAS_ID | The ID of the overlapped CATLAS cCREs region in the hg38 reference, formatted as “chr_start_end”. |
| ReMap_ID | The ID of the overlapped ReMap region in the hg38 reference, formatted as “chr_start_end”. |

## 2. Quick Start Guide

<img src=".\media\methods.png" alt="1" style="zoom:100%;" />


cPeaks eliminates the need for the peak-calling step by providing a ready-to-use reference. If you are familiar with SnapATAC2, ArchR or Python, this section will help you quickly get started with cPeaks. However, if you have any questions, please refer to [detailed tutorials](#detail).

- SnapATAC2
    ```python
    # Read cPeaks file cpeaks_hg19.bed or cpeaks_hg38.bed
    cpeaks_path = 'YOUR_PATH/cpeaks_hg38.bed'
    with open(cpeaks_path) as cpeaks_file:
        cpeaks = cpeaks_file.read().strip().split('\n')
        cpeaks = [peak.split('\t')[0] + ':' + peak.split('\t')[1] + '-' + peak.split('\t')[2] for peak in cpeaks]
    # Set parameter use_rep to 'cpeaks' as mapping reference. 'data' is your AnnData object.
    data = snap.pp.make_peak_matrix(data, use_rep=cpeaks)
    ```
- ArchR
    ```r
    # Read cPeaks file cpeaks_hg19.bed or cpeaks_hg38.bed
    cpeaks <- read_table('YOUR_PATH/cpeaks_hg19.bed', col_names = F)
    cpeaks.gr <- GRanges(seqnames = cpeaks$X1, ranges = IRanges(cpeaks$X2, cpeaks$X3))
    # Set parameter 'features' to cpeaks.gr. 'proj' is your ArchRProject object.
    proj <- addFeatureMatrix(proj, features = cpeaks.gr, matrixName = 'FeatureMatrix')
    ```
- Run Python Script Manually
    ```bash
    git clone https://github.com/MengQiuchen/cPeaks.git
    cd cPeaks/map2cpeak
    python main.py --fragment_path PATH/to/YOUR_fragment.tsv.gz --output map2cpeaks_result
    ```
    After running `main.py` as above, you will get a matrix file `cell_cpeaks.mtx` and a text file `barcodes.txt` within a newly created folder named `map2cpeaks_result`.

    [Learn more about the arguments for main.py](#arguments).

## <a id="detail"></a>3. Comprehensive Guide

SnapATAC2 and ArchR stand out as two popular packages for scATAC-seq data analysis. Integrating cPeaks into the analysis workflow of these packages is straightforward and seamless. Additionally, we provide an user-friendly Python script for transforming fragment files into cell-by-peak matrices. In the following sections, we present detailed code examples and explanations for three scenarios corresponding to the aforementioned cases.

* [SnapATAC2](#method1): A Python/Rust package for single-cell epigenomics analysis. Click the [link](https://github.com/kaizhang/SnapATAC2) for detailed information.
* [ArchR](#method2): A full-featured R package for processing and analyzing single-cell ATAC-seq data. Click the [link](https://github.com/GreenleafLab/ArchR) for detailed information.
* [Run Python Script Manually](#method3): Run Python script [main.py](https://github.com/MengQiuchen/cPeaks/blob/main/map2cpeak/main.py) to obtain cPeaks-based data matrics from fragment files, facilitating downstream analysis steps.

### <a id="method1"></a>3.1 SnapATAC2

#### Install SnapATAC2

SnapATAC2 requires Python>=3.8. There have been changes in the functions and some function parameters between versions 2.4 and 2.5 of SnapATAC2. We recommend installing the 2.5 or higher versions, for compatibility and access to the most recent features and improvements.

```bash
pip install snapatac2==2.5
```

For more installation options, please refer to [SnapATAC2 installation instructions](https://snapatac2.scverse.org/install.html).


#### Integrating cPeaks with SnapATAC2

The example codes and descriptions in this section are adapted from the [SnapATAC2 standard pipeline](https://snapatac2.scverse.org/version/2.5/tutorials/pbmc.html). You can download the code file here: [cPeaks_SnapATAC2.ipynb](https://github.com/MengQiuchen/cPeaks/blob/main/Tutorials/docs/cPeaks_SnapATAC2.ipynb).

[cPeaks_SnapATAC2](docs/cPeaks_SnapATAC2.md ':include')

### <a id="method2"></a>3.2 ArchR

#### Install ArchR

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
```

If you encounter any installation issues, please refer to [ArchR installation instructions]( https://www.archrproject.com/).


#### Integrating cPeaks with ArchR

The example codes and descriptions in this section are adapted from [A Brief Tutorial of ArchR](https://www.archrproject.com/articles/Articles/tutorial.html). You can download the code file here: [cPeaks_ArchR.ipynb](https://github.com/MengQiuchen/cPeaks/blob/main/Tutorials/docs/cPeaks_ArchR.ipynb). 

[cPeaks_ArchR](docs/cPeaks_ArchR.md ':include')

### <a id="method3"></a>3.3 Run Python Script Manually 

In this section, we will illustrate how to map sequencing reads to cPeaks using Map2cPeak. Begin by navigating to the 'map2cpeak' directory. Once there, download and run the software. A demo is also available in the 'demo' folder for trial purposes.

#### Prerequisites

- Python version 3.7 or higher
- Required packages: numpy, gzip and tqdm. 
Ensure these packages are installed before proceeding.

<!-- 
To know time and memory consumption, we tested the code under the following conditions:
- System & version: 
- cpu memory and cores:
- dataset: 


<img src=".\media\time_consuming.png" alt="1" width=400 style="zoom:100%;" /> <img src=".\media\memory_consuming.png" alt="1" width=400 style="zoom:100%;" /> -->

#### Install

```bash
git clone https://github.com/MengQiuchen/cPeaks.git
```

#### Method 1: Map the sequencing reads (fragments.tsv.gz) in each sample/cell to generate cell-by-cPeak matrix (.mtx)

Generate a cell-by-cPeak matrix (.mtx) by mapping the sequencing reads (fragments.tsv.gz) for each sample or cell.

**Important Note**: A minimum of 20GB of memory is required to run this process (`--mode 'normal'`). For improved speed, with at least 30GB of memory, you may activate the default performance mode.

##### Usage Instructions:

1. Navigate to the `map2cpeak` directory:
```bash
cd map2cpeak
```

2. To map your sequencing reads:
```bash
python main.py -f path/to/your_fragment.tsv.gz
```

##### <a id="arguments"></a>Input Arguments

| Argument | Alternate Display Name | Default | Description |
| --------- | ---------------------- | ------- | ----------- | 
| fragment_path | -f | None | Path to the input fragment file (must be a .gz file). |
| barcode_path | -b |None | Path to a file containing barcodes to be used. If not specified, all barcodes in the fragment file will be used. |
| reference |  | hg38 | Specify the cPeaks version: 'hg38' or 'hg19'. |
| output | -o | map2cpeaks_result | Name of the output folder. |

##### Output

A new `output` folder will be created in the current directory, containing a matrix file `cell_cpeaks.mtx` and a text file `barcodes.txt`.

##### Example:

To run a demo mapping:
```bash
python main.py -f demo/test_fragment.tsv.gz
```

To use a provided barcode mapping (ensure 'barcodes.txt' is included in the fragments):
```bash
python main.py -f demo/test_fragment.tsv.gz -b demo/test_barcodes.txt --reference hg19
```

The resulting output includes a `barcode.txt` and a `.mtx` file housing the mapping matrix.

To use hg19 as a reference of mapping, you can run:
```bash
python main.py -f demo/test_fragment.tsv.gz --reference hg19
```

#### Method 2: Mapping Pre-Identified Features to cPeaks (Not Recommended)

**Caution**: This method can result in loss of genomic information as it only considers the pre-identified features. Moreover, the quantification of cPeaks may not be accurate for bulk ATAC-seq data.

##### Usage for Pre-Identified Feature Mapping:

```bash
python main.py --bed_path PATH/to/YOUR_feature.bed
```

##### Input Arguments

| Argument | Alternate Display Name | Default | Description |
| --------- | ---------------------- | ------- | ----------- | 
| bed_path | -bed | None | Path to the input .bed file of features (e.g., 'MACS2calledPeaks.bed'). |
| reference |  | hg38 | Specify the cPeaks version: 'hg38' or 'hg19'. |
| output | -o | map2cpeaks_result | Name of the output folder. |

##### Output

A new `output` folder will be created in the current directory, containing a bed file `map2cpeak.bed`.

Remember to pre-adjust your operational environment according to the system requirements and ensure you’ve properly understood the process to achieve optimal outcomes.

## 4. Reference

[1] Zhang, Kai, et al. "A fast, scalable and versatile tool for analysis of single-cell omics data." Nature methods 21.2 (2024): 217-227.

[2] Granja, Jeffrey M., et al. "ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis." Nature genetics 53.3 (2021): 403-411.

[3] Meng, Qiuchen, et al. "A generic reference defined by consensus peaks for single-cell ATAC-seq data analysis." Nature Communications 17, 2522 (2026). https://doi.org/10.1038/s41467-026-69461-6

## 5. Contact
Please reach out to Meng Qiuchen at [qiuchenmeng@outlook.com](mailto:qiuchenmeng@outlook.com) if you encounter any issues or have any recommendations.
