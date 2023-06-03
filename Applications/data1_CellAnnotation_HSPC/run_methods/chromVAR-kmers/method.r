library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library('JASPAR2016')
library(BSgenome.Hsapiens.UCSC.hg19)
library(umap)

register(MulticoreParam(60))
set.seed(2019)


run_umap <- function(fm_mat){
    umap_object = umap(t(fm_mat),random_state = 2019)
    df_umap = umap_object$layout
    return(df_umap)
}


fun_all <- function(peakfile){
     suppressMessages({   
        peaks <- getPeaks(peakfile, sort_peaks = TRUE)
        peaks <- resize(peaks, width = 500, fix = "center")

        seqinfo(peaks) <- Seqinfo(genome="hg19")
        peaks <- trim(peaks)

        cellnames <- sapply(strsplit(bamfile,'.',fixed = TRUE), "[[", 1)

        fragment_counts <- getCounts(paste0("../../../../../test_data/Buenrostro_2018/bam/files/sc-bams_nodup/",bamfile), 
                                     peaks, 
                                     paired =  TRUE, 
                                     by_rg = TRUE, 
                                     format = "bam", 
                                     colData = data.frame(celltype = cellnames))

        fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg19)

        counts_filtered <- filterPeaks(fragment_counts, non_overlapping = TRUE)

        bg <- getBackgroundPeaks(counts_filtered)
        # Potentially save the bg object
        # saveRDS(bg, file = "bulkPeaks_background_peaks_kmers.rds")

        kmer_ix <- matchKmers(6, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg19)

        dev <- computeDeviations(object = counts_filtered, annotations = kmer_ix,
                                 background_peaks = bg)

        df_zscores = dev@assays@data$z

        df_out <- df_zscores

        ############## UMAP ############

        ## subset according to variance, to make feature number less than observation to plot umap
        df = df_out

        df.vars=df%>%apply(1,var)
        df.means=df%>%apply(1,mean)

        df.plot <- df.vars%>%cbind(df.means)%>%as.data.frame%>%rename_with(~c('var','mean'))
        df.plot <- df.plot%>%mutate(rank.var= base::rank(plyr::desc(var)),
                                    rank.mean=base::rank(plyr::desc(mean)),
                                    labels=ifelse(rank.var<=1500,'variable','non-variable'))
        # psize()
        # df.plot%>%ggplot(aes(x=mean,y=var))+geom_point(aes(color=labels),cex=1,alpha=0.5)+theme_classic()#+xlim(c(0,20))+ylim(c(0,400))

        select.peaks <- df.plot%>%filter(labels=='variable')%>%rownames
        df_out.sub <- df_out[select.peaks,]

        umap <- run_umap(df_out.sub)

        ################# TSNE #################
        tsne <- deviationsTsne(dev, threshold = 1.5, perplexity = 10)

        # variability <- computeVariability(dev)
        # plotVariability(variability, use_plotly = FALSE)
    })
    return(list(df_out=df_out,
            tsne=tsne,
            umap=umap))
    
}