Compiled date: 25th Feb

Source: [Tutorials/README.md](https://github.com/MengQiuchen/cPeaks/blob/dev-tutorials2/Tutorials/README.md)


# cPeak Tutorials

## What is cPeak?

cPeak is xxxxxxx. (introduction and benefit) xxx
xxxxx
xxxxx
xxxxx

### Download Files

cPeak reference files are available from:  [hg19 file](https://cloud.tsinghua.edu.cn/f/d9f40ed01cf749478080/?dl=1) and [hg38 file](https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1). Or downloaded by:
```bash
wget -O YOUR_PATH/cpeaks_hg19_features.txt 'https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1'
wget -O YOUR_PATH/cpeak_hg38_features.txt 'https://cloud.tsinghua.edu.cn/f/dc1c89903e8744eea0aa/?dl=1'
```
## Usages in a Nutshell

<img src=".\media\methods.png" alt="1" style="zoom:70%;" />


Using cPeak is to replace call peaking step by using cPeak reference. If you are familiar with snapATAC2 or ArchR or python, this part we will help you get started with cPeak very quickly. Otherwise if you get confused about anything in this part, please refer to [detailed tutorials](#detail).

- snapATAC2
    ```python
    # read cPeak file
    cpeaks = open('YOUR_PATH/cpeaks_hg19_features.txt').read().strip().split()
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
    [See more about paramters for main.py](#param).


## <a id="detail"></a>Usages in Detail

You can use cPeak
snapATAC2 is a widely used python package for ATAC data analysis. ArchR is a widely used R package for ATAC data analysis. Either can be used for ATAC data analysis.

* [snapATAC2](#method2): snapATAC2 is a widely used python package for ATAC data analysis. Click the [link](https://github.com/kaizhang/SnapATAC2) for detailed information.
* [ArchR](#method3): ArchR is a widely used R package for ATAC data analysis.
* [Direct Use](#method1): Use python script [main.py](https://github.com/MengQiuchen/cPeaks/blob/main/main.py) to transform fragment file to cPeak-based data matrix. It can be use to downstream analysis steps.



## <a id="method2"></a>snapATAC2

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



<img src=".\media\umap1.png" alt="1" style="zoom:70%;" />

Detailed methods refer to https://kzhang.org/SnapATAC2/tutorials/pbmc.html

## <a id="method3"></a>ArchR

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


