# How to use cPeak?
Compiled date: 25th Feb
Source: [Tutorials/xxx](https://breakdance.github.io/breakdance/)
# TODO
- [ ] 统一cpeaks文件
- [ ] 页内跳转


cPeak is xxxxxxx. (introduction and benefit) xxx
xxxxx
xxxxx
xxxxx

cPeak files can be downloaded from:  [hg19](https://cloud.tsinghua.edu.cn/f/d9f40ed01cf749478080/?dl=1) and [hg38](https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1). Or downloaded by:
```bash
wget -O cpeaks_hg19_features.txt 'https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1'
wget -O cpeak_hg38_features.txt 'https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1'
```
## cPeak Tutorials in a Nutshell
![Fig 1](https://github.com/MengQiuchen/cPeaks/blob/dev-tutorials/Tutorials/media/methods.png)
Using cPeak is to replace call peaking step by using cPeak reference. If you are familiar with snapATAC2 or ArchR or python, this part we will help you get started with cPeak very quickly. Otherwise if you get confused about anything in this part, please refer to [detailed tutorials](##2).

snapATAC2 is a widely used python package for ATAC data analysis. ArchR is a widely used R package for ATAC data analysis. Either can be used for ATAC data analysis.

- snapATAC2
    ```python
    # read cPeak file
    cpeaks = open('path_to_downloaded_files/cpeaks_hg19_features.txt').read().strip().split()
    # set parameter use_rep as
    data = snapatac2.pp.make_peak_matrix(data, use_rep=cpeaks)
    ```
- ArchR
    ```{r, message=FALSE, warning=FALSE}
    cpeaks <- read_table('path_to_downloaded_files/cpeaks.v3.simple.bed',col_names=F)
    cpeaks.gr <- GRanges(seqnames=cpeaks$X1, ranges=IRanges(cpeaks$X2, cpeaks$X3))
    proj <- addFeatureMatrix(proj,features = gr,matrixName='FeatureMatrix')
    ```
Otherwise, 

* [Direct Use](##1): Use python script [main.py](https://github.com/MengQiuchen/cPeaks/blob/main/main.py) to transform fragment file to cPeak-based data matrix. It can be use to downstream analysis steps.
* [snapATAC2](##2): snapATAC2 is a widely used python package for ATAC data analysis. Click the [link](https://github.com/kaizhang/SnapATAC2) for detailed information.
* [ArchR](##3): ArchR is a widely used R package for ATAC data analysis.

## Direct Use

in map2cpeaks folder, download it and use it, you can try to run a dome in the dome folder

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

1.


cd map2cpeaks


2. 

python main.py -f path/to/your_fragment.tsv.gz
               
--fragment_path, -f: the input file must be *.tsv.gz file. If barcode_path is None, tsv file need contain 4 columns:chr, start, end, barcode, sep with '\t'. If barcode_path is  not None, tsv file need contain 3 columns:chr, start, end, sep with '\t';
 
optional arguments:

 --help, -h:          show this help message
 --barcode_path, -b:  Each line is a barcode, the code will use the barcodes in the file, Default to use all barcodes in fragment
 --output, -o:        output folder, Default to ./map2cpeaks_result
 --output_name:       name of output files, Default to cell-cpeaks.
 --num_cores, -n:     number of cores to use, Default to 10.
 --reference:         cPeak version, hg38 or hg19, Default to hg38.
```

The output file contains a barcode.txt and an mtx file that stores the matrix of map results.

## Method 2. Directly map the pre-identified features like peaks to cPeaks (NOT recommand)

**This is not a good idea.** It may lose information in the genomic regions which are not included in pre-identfied features. Also, for bulk ATAC-seq data, the quantification of each cPeak is inaccurate.

```bash
usage: python main.py [--bed_path feature.bed]

--bed_path, -bed: the input feature.bed file, for example, MACS2calledPeaks.bed.
```

## snapATAC2

### 1.Download cPeak

```bash
wget -O cpeaks_hg19_features.txt 'https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1'
wget -O cpeak_hg38_features.txt 'https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1'
```

### 2.Install snapATAC2

https://kzhang.org/SnapATAC2/install.html

### 3.Analyze cPeak with snapATAC2

To start, we need to download a fragment file. This can be achieved from https://cloud.tsinghua.edu.cn/d/7f2ccf8067314a4b9a02/ 
```python
import snapatac2 as snap

data = snap.pp.import_data(
    fragment_paths,
    chrom_sizes=snap.genome.hg19,
    sorted_by_barcode=False,
    n_jobs=90
    )

snap.metrics.tsse(data, snap.genome.hg38)
data = snap.pp.make_peak_matrix(data,use_rep=gr)
snap.pp.select_features(data,n_features=len(gr)) # kan renbin 
snap.tl.spectral(data,n_comps=30)
snap.tl.umap(data)
snap.pp.knn(data)
snap.tl.leiden(data,resolution=0.5)
snap.pl.umap(data, color='leiden', interactive=False, height=500)
```

The resulting plot is as follows:



<img src="..\media\1.png" alt="1" style="zoom:70%;" />

Detailed methods refer to https://kzhang.org/SnapATAC2/tutorials/pbmc.html

## ArchR

### 1.Download cPeak

```bash
wget -O cpeaks_hg19_features.txt 'https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1'
wget -O cpeak_hg38_features.txt 'https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1'
```

### 2. Install ArchR

 https://www.archrproject.com/

### 3. Analyze cPeak with ArchR

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



