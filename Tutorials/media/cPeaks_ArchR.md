##### Getting Set Up

First, we load the ArchR library. If loading ArchR fails, you have not properly installed ArchR and should [revisit the installation instructions](https://www.archrproject.com/index.html). We load other useful libraries too. Install them before you load them. We also recommend setting and remembering a known seed to facilitate replication of operations requiring randomization.

```R
library(ArchR)
library(GenomicRanges)
library(tidyverse)
library(dplyr)
library(parallel)
set.seed(1)
```

Next, we set the default number of threads for parallelized operations in ArchR functions. You should change the value passed to `threads` to match the specifications of your local machine.

```R
addArchRThreads(threads = 16)
```

The Hematopoeisis tutorial data can be downloaded using the [getTutorialData()](https://www.archrproject.com/reference/getTutorialData.html) function. The tutorial data is approximately 0.5 GB in size. If you have already downloaded the tutorial in the current working directory, ArchR will bypass downloading.


```R
inputFiles <- getTutorialData("Hematopoiesis")
```

Before we begin, we need add a reference genome annotation for ArchR to have access to chromosome and gene information. ArchR natively supports hg19, hg38, mm9, and mm10.


```R
addArchRGenome('hg19')
```

##### Creating Arrow Files

Now we will create our Arrow files which will take 10-15 minutes. For each sample, this step will:

1. Read accessible fragments from the provided input files.
2. Calculate quality control information for each cell (i.e. TSS enrichment scores and nucleosome info).
3. Filter cells based on quality control parameters.
4. Create a genome-wide TileMatrix using 500-bp bins.
5. Create a GeneScoreMatrix using the custom geneAnnotation that was defined when we called [addArchRGenome()](https://www.archrproject.com/reference/addArchRGenome.html).

```R
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
```

We can inspect the `ArrowFiles` object to see that it is actually just a character vector of Arrow file paths.

```R
ArrowFiles
# Output:
## scATAC_BMMC_R1
## “HemeFragments/scATAC_BMMC_R1.fragments.tsv.gz”
## scATAC_CD34_BMMC_R1
## “HemeFragments/scATAC_CD34_BMMC_R1.fragments.tsv.gz”
## scATAC_PBMC_R1
## “HemeFragments/scATAC_PBMC_R1.fragments.tsv.gz”
```

##### Inferring Doublets

After Arrow file creation, we can infer potential doublets (a single droplet containing multiple cells) that can confound downstream results. This is done using the [addDoubletScores()](https://www.archrproject.com/reference/addArchRGenome.html) function.

```R
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
```

##### Creating an ArchRProject

With our Arrow files in hand, we are now ready to create an `ArchRProject`. An `ArchRProject` is associated with a set of Arrow files and is the backbone of nearly all ArchR analyses.

```R
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
```

Now we can filter putative doublets based on the previously determined doublet scores using the [filterDoublets()](https://www.archrproject.com/reference/filterDoublets.html) function. This doesn’t physically remove data from the Arrow files but rather tells the `ArchRProject` to ignore these cells for downstream analysis.


```R
proj <- filterDoublets(ArchRProj = proj)
```

##### Dimensionality Reduction and Clustering

ArchR implements an iterative LSI dimensionality reduction via the [addIterativeLSI()](https://www.archrproject.com/reference/addIterativeLSI.html) function. Here, we use cPeaks as the Matrix features.

```R
cpeaks <- read_table('YOUR_PATH/cpeaks_hg19.bed', col_names = F)
cpeaks.gr <- GRanges(seqnames = cpeaks$X1, ranges = IRanges(cpeaks$X2, cpeaks$X3))
proj <- addFeatureMatrix(proj, features = cpeaks.gr, matrixName = 'FeatureMatrix')
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "FeatureMatrix", name = "IterativeLSI")
```

To call clusters in this reduced dimension sub-space, we use the [addClusters()](https://www.archrproject.com/reference/addClusters.html) function which uses [Seurat’s](https://satijalab.org/seurat/) graph clustering as the default clustering method.

```R
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
```

##### Visualizing in a 2D UMAP Embedding

We can visualize our scATAC-seq data using a 2-dimensional representation such as Uniform Manifold Approximation and Projection (UMAP). To do this, we add a UMAP embedding to our ArchRProject object with the [addUMAP()](https://www.archrproject.com/reference/addUMAP.html) function. This function uses the [uwot package](https://github.com/jlmelville/uwot) to perform UMAP.

```R
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")  
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")    
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")   
ggAlignPlots(p1, p2, type = "h")
```

<img src=".\media\archr_output1.png" alt="1" style="zoom:100%;" />