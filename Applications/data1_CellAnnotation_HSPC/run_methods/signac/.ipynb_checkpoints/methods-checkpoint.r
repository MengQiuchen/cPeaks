library(Signac)
library(Seurat)
library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v86)

library(SummarizedExperiment)
suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
    library(magrittr)
    library(viridis)
})

set.seed(1234)

fun_all <- function(datafr,dims=c(2:50)){
    chrom_assay <- CreateChromatinAssay(
                      counts = datafr,
                      sep = c("_", "_"),
                      genome = 'hg19',
                      # fragments = "/data1/lichen/code/ECA2.0/scATAC/data/annotation/PBMC//atac_pbmc_5k_nextgem_fragments.tsv.gz",
                      min.cells = 0,
                      min.features = 20
                    )
    pbmc <- CreateSeuratObject(
                      counts = chrom_assay,
                      assay = "peaks",
                      #meta.data = metadata
                    )
    
    pbmc <- RunTFIDF(pbmc)
    pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
    pbmc <- RunSVD(pbmc)
    
    DepthCor(pbmc)%>%print
    
    pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = dims)
    pbmc <- RunTSNE(object = pbmc, reduction = 'lsi', dims = dims)
    
    
    umap = pbmc@reductions$umap@cell.embeddings
    tsne = pbmc@reductions$tsne@cell.embeddings
    df_out <- pbmc@reductions$lsi@cell.embeddings %>%t%>%.[dims,]%>%as.data.frame
    
    return(list(df_out=df_out,
            tsne=tsne,
            umap=umap))
    
}