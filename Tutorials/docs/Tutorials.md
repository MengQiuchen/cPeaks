Compiled date: 11th Mar 2024

Source: [docs/Tutorials.md](https://github.com/MengQiuchen/cPeaks/blob/main/Tutorials/docs/Tutorials.md)

## 1. Overview of cPeaks

### Introduction

cPeaks is a unified reference across cell types for ATAC-seq or scATAC-seq data, enhancing downstream analysis, particularly for cell annotation and rare cell type detection.

<img src=".\media\Cover.png" alt="1" width="500" style="zoom:100%;" />

### Download

#### Simple data

cPeaks reference files can be obtained from the following links: [hg19 file](https://cloud.tsinghua.edu.cn/f/7b7664158dd7482c9a95/?dl=1) and [hg38 file](https://cloud.tsinghua.edu.cn/f/ff4591857f5d472d9401/?dl=1). Alternatively, you can download them using the following commands:

```bash
wget -O YOUR_PATH/cpeaks_hg19.bed https://cloud.tsinghua.edu.cn/f/7b7664158dd7482c9a95/?dl=1
wget -O YOUR_PATH/cpeaks_hg38.bed https://cloud.tsinghua.edu.cn/f/ff4591857f5d472d9401/?dl=1
```

We will utilize the two downloaded files in the upcoming tutorials.

#### Detailed data

For more information, you can download the cPeaks resource files from the following links: [cPeaks resource](https://cloud.tsinghua.edu.cn/f/6460b32917224d32aef1/?dl=1). It contained the basic information, cPeaks annotation, and integration with biological data. It includes 15 columns:

| Column Name | Column Description |
| ----------- | ------------------ |
| ID | The ID of this cPeak. We employed a 14-character string as the unique ID for each cPeak. cPeak ID begins with POR for potential open regions, then HS for human, and ends with nine-digital number. For example, the first cPeak was encoded as PORHS000000001. |
| chr_hg38 | The chromosome of this cPeak in hg38. |
| source | The source of this cPeak: “observed” or “predicted”. |
| start_hg38 | Start positions of this cPeak in hg38 reference. |
| end_hg38 | End positions of this cPeak in hg38 reference. |
| housekeeping | The housekeeping status, whether this cPeak is open in almost all datasets, "TRUE" or "FALSE". |
| openMode | The open mode of cPeaks, "well-positioned", "one-side", "weakly-positioned", or "else". |
| inferredElements | The inferred regulatory elements, "CTCF", "TES", "TSS", "Enhancer" or "Promoter". |
| chr_hg19 | The chromosome of this cPeak in hg19 reference. |
| start_hg19 | start and end positions of this cPeak in hg19 reference. |
| end_hg19 | Start and end positions of this cPeak in hg19 reference. |
| cDHS_ID | The overlapped cDHS ID of this cPeak, names with "chr_start_end" in hg38 reference. |
| rDHS_ID | The overlapped rDHS ID of this cPeak, names with "chr_start_end" in hg38 reference. |
| ReMap_ID | The overlapped ReMap ID of this cPeak, names with "chr_start_end" in hg38 reference. |
| CATLAS_ID | The overlapped CATLAS cCREs of this cPeak, names with "chr_start_end" in hg38 reference. |

## 2. Usages in a Nutshell

<img src=".\media\methods.png" alt="1" style="zoom:100%;" />


cPeaks simplifies the peak-calling step by providing a ready-to-use reference. If you are familiar with SnapATAC2, ArchR or python, this part we will help you get started with cPeaks very quickly. Otherwise if you get confused about anything in this part, please refer to [detailed tutorials](#detail).

- SnapATAC2
    ```python
    # read cPeaks file cpeaks_hg19.bed or cpeaks_hg38.bed
    cpeaks_path = 'YOUR_PATH/cpeaks_hg38.bed'
    with open(cpeaks_path) as cpeaks_file:
        cpeaks = cpeaks_file.read().strip().split('\n')
        cpeaks = [peak.split('\t')[0] + ':' + peak.split('\t')[1] + '-' + peak.split('\t')[2] for peak in cpeaks]
    # set parameter use_rep to user cPeaks as mapping reference
    data = snap.pp.make_peak_matrix(data, use_rep=cpeaks)
    ```
- ArchR
    ```r
    # read cPeaks file cpeaks_hg19.bed or cpeaks_hg38.bed
    cpeaks <- read_table('YOUR_PATH/cpeaks_hg19.bed', col_names = F)
    cpeaks.gr <- GRanges(seqnames = cpeaks$X1, ranges = IRanges(cpeaks$X2, cpeaks$X3))
    proj <- addFeatureMatrix(proj, features = cpeaks.gr, matrixName = 'FeatureMatrix')
    ```
- Run Python Script Manually
    ```bash
    git clone https://github.com/MengQiuchen/cPeaks.git
    cd cPeaks/map2cpeak
    python main.py --fragment_path PATH/to/YOUR_fragment.tsv.gz --output map2cpeaks_result
    ```
    After running `main.py` as above, you will get a matrix file `cell_cpeaks.mtx` and a text file `barcodes.txt` under a newly created folder `map2cpeaks_result`.

    [See more about paramters for main.py](#param).


## <a id="detail"></a>3. Usages in Detail

SnapATAC2 and ArchR are two popular packages for scATAC-seq data analysis. Integrating cPeaks into the analysis workflow of these packages is straightforward and seamless. Additionally, we provide an easy-to-use Python script for transforming fragment files into cell-by-peak matrices. In the following sections, we present detailed code examples and explanations for three scenarios corresponding to the aforementioned cases.

* [SnapATAC2](#method1): A Python/Rust package for single-cell epigenomics analysis. Click the [link](https://github.com/kaizhang/SnapATAC2) for detailed information.
* [ArchR](#method2): A full-featured R package for processing and analyzing single-cell ATAC-seq data. Click the [link](https://github.com/GreenleafLab/ArchR) for detailed information.
* [Run Python Script Manually](#method3): Run python script [main.py](https://github.com/MengQiuchen/cPeaks/blob/main/main.py) to transform fragment files to cPeaks-based data matrics, which can be use to downstream analysis steps.



### <a id="method1"></a>3.1 SnapATAC2

#### Install SnapATAC2

SnapATAC2 requires python>=3.8. There have been changes in the functions and some function parameters between versions 2.4 and 2.5 of SnapATAC2. We recommend installing the 2.5 or higher versions, for compatibility and access to the most recent features and improvements.

```bash
pip install snapatac2==2.5
```

For more installation options, please refer to [SnapATAC2 installation instructions](https://kzhang.org/SnapATAC2/install.html).


#### Integrating cPeaks with SnapATAC2

The example codes and descriptions in this section are adapted from [SnapATAC2 standard pipeline](https://kzhang.org/SnapATAC2/tutorials/pbmc.html). You can download the code file here: [cPeaks_SnapATAC2.ipynb](https://github.com/MengQiuchen/cPeaks/blob/main/Tutorials/docs/cPeaks_SnapATAC2.ipynb).

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

In this section, we will demonstrate how to map sequencing reads to cPeaks using Map2cPeak. Begin by navigating to the 'map2cpeak' directory. Once there, download and run the software. A demo is also available in the 'demo' folder for trial purposes.

#### Prerequisites

- Python version 3.7 or higher
- Required packages: numpy, gzip and tqdm. 
Ensure these packages are installed before proceeding.

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

##### <a id="param"></a>Input Arguments

| Parameter | Alternate display name | Default | Description |
| --------- | ---------------------- | ------- | ----------- | 
| fragment_path | -f | None | Input file path (must be a .gz file) |
| barcode_path | -b |None | Specify a file containing barcodes to be used. If not specified, all barcodes in the fragment file will be used. |
| reference |  | hg38 | Specify the cPeaks version, either 'hg38' or 'hg19'. |
| output | -o | map2cpeaks_result | Designate an output folder name. |

##### Output

You will get a matrix file `cell_cpeaks.mtx` and a text file `barcodes.txt` under the `output` folder.

##### Example:

To run a demo mapping:
```bash
python main.py -f demo/test_fragment.tsv.gz
```

To use a provided barcode mapping (ensure 'barcodes.txt' is included in the fragments):
```bash
python main.py -f demo/test_fragment.tsv.gz -b demo/test_barcodes.txt --reference hg19
```

The resulting output will include a `barcode.txt` and a `.mtx` file housing the mapping matrix.

To use hg19 as a reference of mapping, you can run:
```bash
python main.py -f demo/test_fragment.tsv.gz --reference hg19
```

#### Method 2: Mapping Pre-Identified Features to cPeaks (Not Recommended)

**Caution**: This method can result in loss of genomic information as it only considers the pre-identified features. Moreover, the quantification of cPeaks may not be accurate for bulk ATAC-seq data.

##### Usage for Pre-Identified Feature Mapping:

```bash
python main.py [--bed_path feature.bed]
```

##### Input Arguments

| Parameter | Alternate display name | Default | Description |
| --------- | ---------------------- | ------- | ----------- | 
| bed_path | -bed | None | Input the .bed file of features, such as 'MACS2calledPeaks.bed'. |
| reference |  | hg38 | Specify the cPeaks version, either 'hg38' or 'hg19'. |
| output | -o | map2cpeaks_result | Designate an output folder name. |

##### Output

You will get a matrix file `map2cpeak.bed` under the `output` folder.

Remember to priorly adjust your operational environment according to the system requirements and ensure you’ve properly understood the process to achieve optimal outcomes.

# Reference & Contact

[1]Zhang, K., Zemke, N. R., Armand, E. J. & Ren, B. (2024). A fast, scalable and versatile tool for analysis of single-cell omics data. Nature Methods, 1–11. https://doi.org/10.1038/s41592-023-02139-9

[2]Granja JM, Corces MR et al., ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nature Genetics (2021)

[3]Meng Q, Wu X, Li C, et al. The full set of potential open regions (PORs) in the human genome defined by consensus peaks of ATAC-seq data. doi:10.1101/2023.05.30.542889

Please reach out to Meng Qiuchen at mqc17@mails.tsinghua.edu.cn if you encounter any issues or have any recommendations.