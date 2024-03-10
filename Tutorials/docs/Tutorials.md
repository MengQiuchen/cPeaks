# cPeak Tutorials

Compiled date: 25th Feb

Source: [Tutorials/README.md](https://github.com/MengQiuchen/cPeaks/blob/dev-tutorials/Tutorials/README.md)

## 1. What is cPeak?

cPeak is xxxxxxx. (introduction and benefit) xxx
xxxxx
xxxxx
xxxxx

### Download Files

cPeak reference files can be obtained from the following links: [hg19 file](https://cloud.tsinghua.edu.cn/f/7b7664158dd7482c9a95/?dl=1) and [hg38 file](https://cloud.tsinghua.edu.cn/f/ff4591857f5d472d9401/?dl=1). Alternatively, you can download them using the following commands:

```bash
wget -O YOUR_PATH/cpeaks_hg19.bed https://cloud.tsinghua.edu.cn/f/7b7664158dd7482c9a95/?dl=1
wget -O YOUR_PATH/cpeaks_hg38.bed https://cloud.tsinghua.edu.cn/f/ff4591857f5d472d9401/?dl=1
```

We will utilize the downloaded files in the upcoming tutorials.

## 2. Usages in a Nutshell

<img src=".\media\methods.png" alt="1" style="zoom:100%;" />


Using cPeak is to replace the call peaking step by using cPeak reference. If you are familiar with SnapATAC2, ArchR or python, this part we will help you get started with cPeak very quickly. Otherwise if you get confused about anything in this part, please refer to [detailed tutorials](#detail).

- SnapATAC2
    ```python
    # read cPeak file
    cpeaks = open('cpeaks_hg38.bed').read().strip().split('\n')[1:]
    cpeaks = [i.split('\t')[0]+':'+i.split('\t')[1]+'-'+i.split('\t')[2] for i in cpeaks]
    # set parameter use_rep as
    data = snapatac2.pp.make_peak_matrix(data, use_rep=cpeaks)
    ```
- ArchR
    ```r
    cpeaks <- read_table('YOUR_PATH/cpeaks.v3.simple.bed',col_names=F)
    cpeaks.gr <- GRanges(seqnames=cpeaks$X1, ranges=IRanges(cpeaks$X2, cpeaks$X3))
    proj <- addFeatureMatrix(proj,features = gr,matrixName='FeatureMatrix')
    ```
- Run Python Script Manually
    ```bash
    git clone https://github.com/MengQiuchen/cPeaks.git
    cd cPeaks/map2cpeaks
    python main.py --fragment_path PATH/to/YOUR_fragment.tsv.gz --output map2cpeaks_result --output_name Cell_by_cPeak_Matrix --type_saved .mtx
    ```
    After manually running `main.py`, you will get a matrix file `Cell_by_cPeak_Matrix.mtx` under a newly created folder `map2cpeaks_result`.

    [See more about paramters for main.py](#param).


## <a id="detail"></a>3. Usages in Detail

SnapATAC2 and ArchR are two popular packages for scATAC-seq (single-cell epigenomics?) data analysis. Integrating cPeak into the analysis workflow of these packages is straightforward and seamless. Additionally, we provide an easy-to-use Python script for transforming fragment files into cell-by-peak matrices. In the following sections, we present detailed code examples and explanations for three scenarios corresponding to the aforementioned cases.

* [SnapATAC2](#method1): 用cpeak用这个包做了什么 Click the [link](https://github.com/kaizhang/SnapATAC2) for detailed information.
* [ArchR](#method2): 用cpeak用这个包做了什么 Click the [link](https://github.com/GreenleafLab/ArchR) for detailed information.
* [Run Python Script Manually](#method3): Use python script [main.py](https://github.com/MengQiuchen/cPeaks/blob/main/main.py) to transform fragment files to cPeak-based data matrics, which can be use to downstream analysis steps.



### <a id="method1"></a>3.1 SnapATAC2

#### Install SnapATAC2

SnapATAC2 requires python>=3.8. There have been changes in the functions and some function parameters between versions 2.4 and 2.5 of SnapATAC2. We recommend installing the 2.5 or higher versions, for compatibility and access to the most recent features and improvements.

```bash
pip install snapatac2==2.5
```

For more installation options, please refer to [SnapATAC2 installation instructions](https://kzhang.org/SnapATAC2/install.html).


#### Integrating cPeak with SnapATAC2

The example codes and descriptions in this section are adapted from [SnapATAC2 standard pipeline](https://kzhang.org/SnapATAC2/tutorials/pbmc.html). You can download the code file here: [cPeak_SnapATAC2.ipynb](https://TODO). TODO: add file link

[cPeak_SnapATAC2](media/cPeak_SnapATAC2.md ':include')

### <a id="method2"></a>3.2 ArchR

#### Install ArchR

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
```

If you encounter any installation issues, please refer to [ArchR installation instructions]( https://www.archrproject.com/).


#### Integrating cPeak with ArchR

The example codes and descriptions in this section are adapted from [A Brief Tutorial of ArchR](https://www.archrproject.com/articles/Articles/tutorial.html). You can download the code file here: [cPeaks_ArchR.ipynb](https://TODO). TODO: add file link

[cPeak_ArchR](media/cPeaks_ArchR.md ':include')

### <a id="method1"></a>3.3 Run Python Script Manually 

In this section, we will demonstrate how to map sequencing reads to cPeaks using Map2cpeak. Begin by navigating to the 'map2cpeaks' directory. Once there, download and run the software. A demo is also available in the 'demo' folder for trial purposes.

#### Prerequisites

- Python version 3.7 or higher
- Required Packages: numpy, gzip and tqdm. 
Ensure these packages are installed before proceeding.

#### Install

```bash
git clone https://github.com/MengQiuchen/cPeaks.git
```

#### Method 1: Map the sequencing reads (fragments.tsv.gz) in each sample/cell to generate cell-by-cPeak matrix (.mtx)

Generate a cell-by-cPeak matrix (.mtx) by mapping the sequencing reads (fragments.tsv.gz) for each sample or cell.

**Important Note**: A minimum of 20GB of memory is required to run this process (`--mode 'normal'`). For improved speed, with at least 30GB of memory, you may activate the default performance mode.

##### Usage Instructions:

1. Navigate to the `map2cpeaks` directory:
```bash
cd map2cpeaks
```

2. To map your sequencing reads:
```bash
python main.py -f path/to/your_fragment.tsv.gz
```
- `--fragment_path`, `-f`: Input file path (must be a .gz file).

##### <a id="param"></a>Arguments

| Parameter | Default | Description |
| ------ | ------ | ------ | 
| fragment_path | - | The input file must be *.tsv.gz file, which will be transformed to cPeak reference. |
| barcode_path | - | If barcode file is given, barcode in the file or all barcodes in the fragment will be used. |
| reference | hg38 | cPeak version, hg38 or hg19. |
| type_saved | .mtx | The type of output file, .mtx or .h5ad. |
| output | map2cpeaks_result | Output folder name. |
| output_name | cell-cpeaks | Name of output files. |
| num_cores | 10 |  Number of cores to use. |

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
- `--bed_path`, `-bed`: Input the `.bed` file of features, such as `MACS2calledPeaks.bed`.

Remember to priorly adjust your operational environment according to the system requirements and ensure you’ve properly understood the process to achieve optimal outcomes.


# Contact

Please reach out to Meng Qiuchen at mqc17@mails.tsinghua.edu.cn if you encounter any issues or have any recommendations.