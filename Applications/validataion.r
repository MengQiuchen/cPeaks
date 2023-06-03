# pipeline
suppressPackageStartupMessages({   
    library(data.table)
    library(Matrix)
    library(proxy)
    library(Rtsne)
    library(densityClust)
    library(data.table)
    library(irlba)
    library(umap)
    library(ggplot2)
    library(parallel)
    library(tidyverse)    
})

# clustering
suppressPackageStartupMessages({   
    library(data.table)
    library(Matrix)
    library(proxy)
    library(Rtsne)
    library(densityClust)
    library(data.table)
    library(irlba)
    library(umap)
    library(ggplot2)
    library(RColorBrewer)
})

# umap
suppressPackageStartupMessages({  
    library(data.table)
    library(dplyr)
    library(Matrix)
    library(BuenColors)
    library(stringr)
    library(GenomicRanges)
    library(irlba)
    library(cicero)
    library(umap)
})    

# clustering
suppressPackageStartupMessages({  
library(data.table)
library(Matrix)
library(proxy)
library(Rtsne)
library(densityClust)
library(data.table)
library(irlba)
library(umap)
library(ggplot2)
library(RColorBrewer)
    library(ggrepel)
}) 

elbow_plot <- function(mat,num_pcs=50,scale=FALSE,center=FALSE,title='',width=3,height=3){
    set.seed(2019) 
    mat = data.matrix(mat)
    SVD = irlba(mat, num_pcs, num_pcs,scale=scale,center=center)
    options(repr.plot.width=width, repr.plot.height=height)
    df_plot = data.frame(PC=1:num_pcs, SD=SVD$d);
#     print(SVD$d[1:num_pcs])
    p <- ggplot(df_plot, aes(x = PC, y = SD)) +
      geom_point(col="#cd5c5c",size = 1) + 
      ggtitle(title)
    return(p)
}

plot.tsne <- function(x, labels,
         main="A tSNE visualization",n=20,
         pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
         cex.main=1, cex.legend=1) {
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  layout = x
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  }
  points(layout[,1], layout[,2], col=col_vector[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)

  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text,
         col=col_vector[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}


colormap = c(jdb_color_maps, "UNK" = "#333333" )

run_umap <- function(fm_mat){
    umap_object = umap(t(fm_mat),random_state = 2019)
    df_umap = umap_object$layout
    return(df_umap)
}


plot_umap <- function(df_umap,labels,title='UMAP',colormap=colormap,point_size=0.2){
    set.seed(2019) 
    df_umap = data.frame(cbind(df_umap,labels),stringsAsFactors = FALSE)
    colnames(df_umap) = c('umap_1','umap_2','celltype')
    df_umap$umap_1 = as.numeric(df_umap$umap_1)
    df_umap$umap_2 = as.numeric(df_umap$umap_2)
    #options(repr.plot.width=4, repr.plot.height=4)
    p <- ggplot(shuf(df_umap), aes(x = umap_1, y = umap_2, color = celltype)) +
      geom_point(size = point_size) + scale_color_manual(values = colormap) +
      ggtitle(title)+theme_classic()
    return(p)
}



fun_all <- function(datafr,filter_ratio=0.01,num_pcs=150,plot.hist=TRUE){
   
    start_time <- Sys.time()
    binary_mat = as.matrix((datafr > 0) + 0)
    binary_mat = Matrix(binary_mat, sparse = TRUE) 

    num_cells_ncounted = rowSums(binary_mat)
    ncounts = binary_mat[num_cells_ncounted >= dim(binary_mat)[2]*0.01,]
    new_counts = colSums(ncounts)
    ncounts = ncounts[rowSums(ncounts) > 0,]

    if(plot.hist){
        par(mfrow=c(1,2))
        hist(log10(num_cells_ncounted),main="No. of Cells Each Site is Observed In",breaks=50)
        abline(v=log10(min(num_cells_ncounted[num_cells_ncounted >= dim(binary_mat)[2]*filter_ratio])),lwd=2,col="indianred")
        hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
        # abline(v=log10(quantile(new_counts,probs=0.1)),lwd=2,col="indianred")
    }

    sexsites = c(grep("chrY",rownames(ncounts)),grep("chrX",rownames(ncounts)))
    if(length(sexsites)==0){
        ncounts.nosex=ncounts
    }else{
        ncounts.nosex = ncounts[-sexsites,]  
    }

    nfreqs = t(t(ncounts.nosex) / Matrix::colSums(ncounts.nosex))
    idf = as(log(1 + ncol(ncounts.nosex) / Matrix::rowSums(ncounts.nosex)), "sparseVector")
    tf_idf_counts = as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs


    p_elbow_LSI <- elbow_plot(tf_idf_counts,num_pcs = 200, title = 'PCA on LSI')
    p_elbow_LSI

    set.seed(2019)
    num_pcs = num_pcs
    SVDtsne = irlba(tf_idf_counts, num_pcs, num_pcs)
    d_diagtsne = matrix(0, nrow=num_pcs, ncol=num_pcs)
    diag(d_diagtsne) = SVDtsne$d
    SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v))



    df_out = t(SVDtsne_vd)
    colnames(df_out) = colnames(datafr)
    rownames(df_out) = paste('PC',1:dim(SVDtsne_vd)[2])

    end_time <- Sys.time()
    message(end_time - start_time)

    ### Downstream Analysis
    
    umap <- run_umap(df_out)

    #plot_umap(df_umap_Cusanovich2018,labels = labels,colormap = colormap,title='Cusanovich2018')
    
    set.seed(0)
    tsnetfidf = Rtsne(SVDtsne_vd,pca=F)
    library(RColorBrewer)
    #plot.tsne(tsnetfidf$Y,as.factor(metadata[,'label']))

    return(list(ncounts=ncounts,
                ncounts.nosex=ncounts.nosex,
                tf_idf_counts=tf_idf_counts,
                df_out=df_out,
                tsne=tsnetfidf$Y,
                umap=umap))
}


fun_densityClust_old <- function(res, labels=NULL, rho_ = 20, delta_ = 10){
    
    tsne <- res$tsne   
   #To identify clusters of cells, we use the density peak algorithm.
    tsnedist = dist(tsne)
    set.seed(0)
    dclust = densityClust(tsnedist,gaussian=T)
    dclust = findClusters(dclust, rho = rho_, delta = delta_)
    
    #number of clusters produced
    nClusters_produced = length(levels(as.factor(dclust$clusters)))
    
    #plot from tutorial
    #The density peak algorithm requires you to set two parameters - “delta” and “rho”. For each data point, the algorithm calculates a local density of other points within some set distance and the minimum distance to the next point that has a higher local density. On the basis of these two values, you can choose a set of points that are outliers both in local density and the distance to another point with a higher density, which become cluster “peaks”. Below, we show you the distribution of these two values in our data set and where we decided to draw the cutoff. You can read more about this algorithm here.
    options(repr.plot.width=6, repr.plot.height=6)
    plot(dclust$rho,dclust$delta,pch=20,cex=0.6)
    points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
    text(dclust$rho[dclust$peaks]-2,dclust$delta[dclust$peaks]+1.5,labels=dclust$clusters[dclust$peaks])
    abline(v=rho_)
    abline(h=delta_)
    
    #plot from tutorial
    tsnecols = c("#E31A1C","#FFD700","#771122","#777711","#1F78B4","#68228B","#AAAA44",
                     "#60CC52","#771155","#DDDD77","#774411","#AA7744","#AA4455","#117744",
                     "#000080","#44AA77","#AA4488","#DDAA77")
    plot(tsne,pch=20,col=tsnecols[as.factor(dclust$clusters)],main="Density Peak Clusters (TSNE)",cex=0.25)
    text(tsne[dclust$peaks,1],tsne[dclust$peaks,2],labels=dclust$clusters[dclust$peaks],cex=2.5)

    plot(res$umap,pch=20,col=tsnecols[as.factor(dclust$clusters)],main="Density Peak Clusters (UMAP)",cex=0.25)
    text(res$umap[dclust$peaks,1],res$umap[dclust$peaks,2],labels=dclust$clusters[dclust$peaks],cex=2.5)

    if(length(labels)!=0){
        plot.tsne(res$tsne,as.factor(metadata[,'label']),cex = 0.2) 
        plot_umap(res$umap,labels = metadata$label,colormap = colormap,title='Cusanovich2018')%>%print
    }
    return(dclust)
    
}

library(patchwork)

fun_densityClust <- function(res, labels=NULL, rho_ = 20, delta_ = 10,title='',plot=TRUE){
    
    tsne <- res$tsne    
    #To identify clusters of cells, we use the density peak algorithm.
    tsnedist = dist(tsne)
    set.seed(0)
    dclust = densityClust(tsnedist,gaussian=T)
    dclust = findClusters(dclust, rho = rho_, delta = delta_)
    
    #number of clusters produced
    nClusters_produced = length(levels(as.factor(dclust$clusters)))
    
    #plot from tutorial
    #The density peak algorithm requires you to set two parameters - “delta” and “rho”. For each data point, the algorithm calculates a local density of other points within some set distance and the minimum distance to the next point that has a higher local density. On the basis of these two values, you can choose a set of points that are outliers both in local density and the distance to another point with a higher density, which become cluster “peaks”. Below, we show you the distribution of these two values in our data set and where we decided to draw the cutoff. You can read more about this algorithm here.
    df.dclust <- cbind(rho=dclust$rho, delta=dclust$delta)%>%as.data.frame%>%
            mutate(id=row_number())%>%
            mutate(colors=ifelse(id%in%dclust$peaks,'red','black'))%>%
            mutate(a=dclust$clusters,
                   label=ifelse(id%in%(dclust$peaks),a,''))
    plot.clusterNum <- (ggplot(df.dclust)+geom_point(aes(x=rho,y=delta,color=colors),cex=0.4)+theme_classic()+
                        geom_text_repel(aes(x=rho,y=delta,label=label),size=2) +
                        geom_hline(yintercept = delta_,color='grey')+
                        geom_vline(xintercept = rho_,color='grey')+ggtitle(title)+
                        scale_color_manual(values=c('black',colors_[1])))%>%fun_ggtheme()+theme(legend.position = 'none')
    
    #plot from tutorial
    tsnecols = c("#E31A1C","#FFD700","#771122","#777711","#1F78B4","#68228B","#AAAA44",
                     "#60CC52","#771155","#DDDD77","#774411","#AA7744","#AA4455","#117744",
                     "#000080","#44AA77","#AA4488","#DDAA77")
    ####### plot tsne
    tmp=as.data.frame(cbind(x=tsne[dclust$peaks,1],y=tsne[dclust$peaks,2]))
    df.tsne <- tsne %>%as.data.frame%>%rename_with(~c('X','Y'))%>%mutate(color=as.character(dclust$clusters))
    plot.tsne.cluster <- (ggplot(df.tsne)+geom_point(aes(x=X,y=Y,color=color),cex=0.1,alpha=0.6)+
                geom_text(data=tmp,aes(x=x,y=y),
                          label=dclust$clusters[dclust$peaks],size=4)+
                scale_color_manual(values=tsnecols))%>%fun_ggtheme()+ggtitle('TSNE(cluster labels)')+
                theme(legend.position = 'none')


    ####### plot umap
    tmp=as.data.frame(cbind(x=res$umap[dclust$peaks,1],y=res$umap[dclust$peaks,2]))
    df.umap <- res$umap %>%as.data.frame%>%rename_with(~c('X','Y'))%>%mutate(color=as.character(dclust$clusters))
    plot.umap.cluster <-(ggplot(df.umap)+geom_point(aes(x=X,y=Y,color=color),cex=0.1,alpha=0.6)+
                geom_text(data=tmp,aes(x=x,y=y),
                          label=dclust$clusters[dclust$peaks],size=4)+
                scale_color_manual(values=tsnecols))%>%fun_ggtheme()+ggtitle('UMAP(cluster labels)')+
                theme(legend.position = 'none')
    
    if(length(labels)!=0){
  
        df.tsne2 <- tsne %>%as.data.frame%>%rename_with(~c('X','Y'))%>%mutate(color=labels)
        plot.tsne.label <- (ggplot(df.tsne2)+geom_point(aes(x=X,y=Y,color=color),cex=0.1,alpha=0.8)+
                        scale_color_manual(values=tsnecols))%>%fun_ggtheme()+ggtitle('TSNE(true labels)')+
                        theme(legend.position = 'right')

        df.umap2 <- res$umap %>%as.data.frame%>%rename_with(~c('X','Y'))%>%mutate(color=labels)
        plot.umap.label <- (ggplot(df.umap2)+geom_point(aes(x=X,y=Y,color=color),cex=0.1,alpha=0.8)+
                        scale_color_manual(values=tsnecols))%>%fun_ggtheme()+ggtitle('UMAP(true labels)')+
                        theme(legend.position = 'right')
        
    }else{
        plot.tsne.label=NULL
        plot.umap.label=NULL
    }
    if(plot){
            grid.arrange(plot.clusterNum,
                         plot.tsne.cluster,
                         plot.umap.cluster,
                         plot.tsne.label+theme(legend.position = 'none'),
                         plot.umap.label+theme(legend.position = 'none'),ncol=5)
    }

    
    return(list(clust=dclust,
                plot=list(plot.clusterNum=plot.clusterNum,
                          plot.tsne.cluster=plot.tsne.cluster,
                          plot.umap.cluster=plot.umap.cluster,
                          plot.tsne.label=plot.tsne.label,
                          plot.umap.label=plot.umap.label)))
    
}