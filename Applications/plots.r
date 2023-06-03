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
    library(patchwork)
}) 

colormap = c(jdb_color_maps, "UNK" = "#333333" )



run_umap <- function(fm_mat){
    umap_object = umap(t(fm_mat),random_state = 2019)
    df_umap = umap_object$layout
    return(df_umap)
}

fun_densityClust <- function(res, labels=NULL, rho_ = 20, delta_ = 10,title=''){
    
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
        
    }

    grid.arrange(plot.clusterNum,
                  plot.tsne.cluster,
                  plot.umap.cluster,
                  plot.tsne.label+theme(legend.position = 'none'),
                   plot.umap.label+theme(legend.position = 'none'),ncol=5)
    
    return(list(clust=dclust,
                plot=list(plot.clusterNum=plot.clusterNum,
                          plot.tsne.cluster=plot.tsne.cluster,
                          plot.umap.cluster=plot.umap.cluster,
                          plot.tsne.label=plot.tsne.label,
                          plot.umap.label=plot.umap.label)))
    
}