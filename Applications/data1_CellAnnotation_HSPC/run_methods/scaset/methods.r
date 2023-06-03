suppressPackageStartupMessages({
    library(prabclus)
    library(Matrix)
    library(Rtsne)
    library(ggplot2) 
    
})


run_umap <- function(fm_mat){
    umap_object = umap(t(fm_mat),random_state = 2019)
    df_umap = umap_object$layout
    return(df_umap)
}

filter_peaks <- function (datafr,cutoff = 0.01){
    binary_mat = as.matrix((datafr > 0) + 0)
    binary_mat = Matrix(binary_mat, sparse = TRUE) 
    num_cells_ncounted = Matrix::rowSums(binary_mat)
    ncounts = binary_mat[num_cells_ncounted >= dim(binary_mat)[2]*cutoff,]
    ncounts = ncounts[rowSums(ncounts) > 0,]    
    
    options(repr.plot.width=4, repr.plot.height=4)
    hist(log10(num_cells_ncounted),main="No. of Cells Each Site is Observed In",breaks=50)
    abline(v=log10(min(num_cells_ncounted[num_cells_ncounted >= dim(binary_mat)[2]*cutoff])),lwd=2,col="indianred")
#     hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
    datafr_filtered = datafr[rownames(ncounts),]
    return(datafr_filtered)
}


getJaccardDist <- function(cdBinary){
        
    if(colnames(cdBinary[,2:3])[1] == 'start' && colnames(cdBinary[,2:3])[2] == 'end'){
        SingleCell.Binary <- cdBinary[,4:(dim(cdBinary)[2])]
    }
    else
        SingleCell.Binary <- cdBinary
    
    
    SingleCell.Binary.Jaccard <- jaccard(as.matrix(SingleCell.Binary))
    
    return(SingleCell.Binary.Jaccard)
}


plot_tSNE <- function(fit, nDimToUSE,SingleCell.Binary.Jaccard,groups=NULL, cellName, perplexity_division, ret.val=FALSE , text.label=FALSE, title=""){


    if(missing(ret.val)){
        ret.val = FALSE
    }
    
    if(text.label==TRUE & missing(cellName)){
        stop("ERROR: Please give cellName")
    }
    
    if(missing(perplexity_division)){
        stop("ERROR: Please enter the number with which you want to divide the dimension of your data for perplexity setting")
    }
        
    

    
    #rtsne_pca_out <- Rtsne(as.matrix(pcaPRComp$x[,1:nPCAToUSE]), perplexity = dim(pcaPRComp$x)[1]/perplexity_division, check_duplicates = FALSE)
    rtsne_pca_out <- Rtsne(as.matrix(fit$points[,1:nDimToUSE]), perplexity = perplexity_division, 
                           check_duplicates = FALSE, pca=FALSE, theta=0.01, max_iter=3000)
    
    if(is.null(groups)){
        if(ret.val==TRUE)
            df<-data.frame(X=rtsne_pca_out$Y[,1],Y=rtsne_pca_out$Y[,2], cellName=cellName)
        else
            df<-data.frame(X=rtsne_pca_out$Y[,1],Y=rtsne_pca_out$Y[,2])
        p1 <- ggplot(df, aes_string(x="X",y ="Y"))
    }
    else{
        if(ret.val==TRUE)
            df<-data.frame(X=rtsne_pca_out$Y[,1],Y=rtsne_pca_out$Y[,2], Cell=colnames(SingleCell.Binary.Jaccard), Batch=groups, cellName=cellName)
        else
            df<-data.frame(X=rtsne_pca_out$Y[,1],Y=rtsne_pca_out$Y[,2], Cell=colnames(SingleCell.Binary.Jaccard), Batch=groups)
        # p1 <- ggplot(df, aes_string(x="X",y ="Y", color="Batch"))
    }
    
    return(df)
}

fun_all <- function(datafr,cut_dim=15,tsne_perplexity=30){
    datafr_filtered <- filter_peaks(datafr)
    cdBinary = as.matrix((datafr_filtered > 0) + 0)
    cdBinary = Matrix(cdBinary, sparse = TRUE) 
    
    ## use the similar number of dimensions as shown in tutorial https://github.com/ManchesterBioinference/Scasat/blob/master/ScAsAT_functions_Buenrostro_All_Bam_Together.ipynb
    SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)
    fit <- cmdscale(as.dist(SingleCell.Binary.Jaccard),eig=TRUE, k=cut_dim)
    fm_Scasat = t(fit$points)
    rownames(fm_Scasat) = paste('Dim',1:dim(fm_Scasat)[1])
    
    # umap
    df_out <- fm_Scasat
    umap <- run_umap(df_out)
    
    # tsne
    tsNE_out_2PCs <- plot_tSNE(fit,cut_dim,SingleCell.Binary.Jaccard,groups=metadata$label, perplexity=tsne_perplexity, 
                              cellName=rownames(metadata),ret.val=TRUE, text.label=F)
    tsne <- tsNE_out_2PCs[,1:2]


    return(list(df_out=df_out,
            tsne=tsne,
            umap=umap))
}