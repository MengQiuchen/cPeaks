Compiled date: 25th Feb

Source: [Tutorials/README.md](https://github.com/MengQiuchen/cPeaks/blob/dev-tutorials/Tutorials/README.md)


# cPeak Tutorials

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


#### Utilize cPeak with SnapATAC2

The example codes and descriptions in this section are adapted from [SnapATAC2 standard pipeline](https://kzhang.org/SnapATAC2/tutorials/pbmc.html). You can download the code file here: [cPeak_SnapATAC2.ipynb](https://TODO). TODO: add file link


[testtest](media/cPeak_SnapATAC2.md ':include')


### <a id="method2"></a>3.2 ArchR

#### Install ArchR

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
```

If you encounter any installation issues, please refer to [ArchR installation instructions]( https://www.archrproject.com/).


#### Utilize cPeak with ArchR

Code references to https://www.archrproject.com/bookdown/index.html

```r
library(ArchR)
set.seed(1)
addArchRThreads(threads = 90) 
library(GenomicRanges)
library(tidyverse)
library(dplyr)
addArchRGenome(‘hg19')

cpeaks = read_table('/nfs/mqc/Consensus_peak/code/cpeaks/all/hg19/cpeaks.v3.simple.bed',col_names = F)
cpeaks.gr <- GRanges(seqnames = cpeaks$X1,
              ranges = IRanges(cpeaks$X2, cpeaks$X3))
get_cluster = function(fragement_path,output.name,num_cluster,gr=NULL,ref = 'hg19'){
    addArchRGenome(ref)
    if (length(list.files(fragement_path)) == 0 & (!endsWith(fragement_path,'tsv.gz'))){
        print("input should be tsv.gz or list of tsv.gz")
        return()
    }
    if (length(list.files(fragement_path))>0){
        file_list <- list.files(fragement_path, pattern = "fragments.tsv.gz$", full.names = TRUE)
        print(paste0('there are ',length(file_list),' fragments files'))
        file_name <- lapply(file_list, function(x) {
          tmp <- strsplit(x, '/')[[1]]
          tmp[length(tmp)]
        })
    }
	else{
        file_list = fragement_path
        print(paste0('there are 1 files: ',fragement_path))
        tmp <- strsplit(fragement_path, '/')[[1]]
        file_name = tmp[length(tmp)]   
    }
    ArrowFiles <- createArrowFiles(
      inputFiles = file_list,
      sampleNames = file_name,
      minTSS = 0, #Dont set this too high because you can always increase later
      minFrags = 0, 
      addTileMat = TRUE,
      addGeneScoreMat = TRUE
    )
    proj =  ArchRProject(
      ArrowFiles = ArrowFiles, 
      outputDirectory = output.name,
      copyArrows = F #This is recommened so that you maintain an unaltered copy for later usage.
    )
    #proj <- filterDoublets(ArchRProj = proj)
    if (!is.null(gr)){
        proj = addFeatureMatrix(proj,features = gr,matrixName='FeatureMatrix')
        proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "FeatureMatrix", name = "IterativeLSI")   
    }
    else{
        proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")    
    }
    proj <- addClusters(input = proj, reducedDims = "IterativeLSI",maxClusters=num_cluster)   
    proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI") 
    p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")    
    p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")  
    ggAlignPlots(p1, p2, type = "h")
    rname = row.names(getCellColData(proj))
    cluster = unlist(proj$Clusters%>%lapply(function(x){strsplit(x,'C')[[1]][2]%>%as.numeric}))
    pre = cbind(rname,cluster)%>%as.data.frame  
    pre$cluster = as.numeric(pre$cluster)   
    #saveArchRProject(ArchRProj = proj, outputDirectory = paste0(output.name,'_save'), load = FALSE)
    return(pre)
}
library(parallel)
res = get_cluster("sort.HSC.fragments.tsv.gz",'HSC_all',num_cluster = 10,gr=cpeaks.gr)

```

## <a id="method1"></a>Run Python Script Manually 

### <a id="param"></a>main.py Parameters

| Parameter | Default | Description |
| ------ | ------ | ------ | 
| fragment_path | - | The input file must be *.tsv.gz file, which will be transformed to cPeak reference. |
| barcode_path | - | If barcode file is given, barcode in the file or all barcodes in the fragment will be used. |
| reference | hg38 | cPeak version, hg38 or hg19. |
| type_saved | .mtx | The type of output file, .mtx or .h5ad. |
| output | map2cpeaks_result | Output folder name. |
| output_name | cell-cpeaks | Name of output files. |
| num_cores | 10 |  Number of cores to use. |


### Version
Python (>=3.7)

### Requirments

```
numpy
gzip
tqdm
```

### Method 1: Map the sequencing reads (fragments.tsv.gz) in each sample/cell to generate cell-by-cPeak matrix (.mtx/.h5ad)

Map the sequencing reads (fragments.tsv.gz) in each sample/cell to generate cell-by-cPeak matrix (.mtx/.h5ad)
```bash
usage:


The output file contains a barcode.txt and an mtx file that stores the matrix of map results.

### Method 2. Directly map the pre-identified features like peaks to cPeaks (NOT recommand)

**This is not a good idea.** It may lose information in the genomic regions which are not included in pre-identfied features. Also, for bulk ATAC-seq data, the quantification of each cPeak is inaccurate.

```bash
usage: python main.py [--bed_path feature.bed]

--bed_path, -bed: the input feature.bed file, for example, MACS2calledPeaks.bed.
```


