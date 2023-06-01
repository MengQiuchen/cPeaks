# 版本相对于v1
# 1. 修正fun_all最后res.all unimodal结果len=NULL 的错误
# 3. 在fun_hdbscan中hdbscan函数前，fun_perChr中res.point初次生成前，添加set.seed(0)保证每次结果的可重复性
# 4. 把fun_all里的tryCatch，warning去掉了，使得warning的时候不再返回NULL


########################               peak split                   ###################################
# HDBSCAN method
#minPts_=NULL
fun_hdbscan <- function(mat,minPts_len=50){
    
    # 假设minimum length of peak is 50 basepair,所以我们平均需要融合
    # 总reads/length * 50
    # 把50bp里面的所有位点做一个融合。
    
    #if(length(minPts_)==0) {minPts_=(nrow(mat)/length(unique(mat$data)))*50}
    minPts_=(nrow(mat)/length(unique(mat$data)))*minPts_len
    
    set.seed(0)
    cl <- hdbscan(mat%>%as.matrix, minPts = minPts_)    

    res.point <- cbind(data=mat$data,y=mat$y,clusterID=cl$cluster, prob=cl$membership_prob,idx=mat$idx)%>%as.data.frame
    
    return(list(res.point=res.point,
               cl=cl))
}


# dip_p_th=1
# min_split_len=800
# max_iteration=100
# max_peaksum=5e4
fun_perChr <- function(i,s,v,
                       max_iteration=100,min_split_len=800,
                       dip_p_th=0.1,max_peaksum=5e4){
    

#     # use res.chr
#     v <- res.chr[[i]]$num
#     s <- res.chr[[i]]$site
#     data=rep(s,times = v)
#     scale_=2 
    
    
#     for(ii in c(1:10)){
        
#         v <- (res.chr[[i]]$num/scale_)%>% floor
#         data= rep(s,times = v)
#         if(length(data)< max_peaksum) {break}
#         scale_ = scale_ +1
        
#         }
    set.seed(0)
    res.point <- cbind(data=rep(s,times = v), y= rnorm(length(rep(s,times = v)), mean=0, sd=0.1),
                    clusterID=0,prob=0,time=0)%>%as.data.frame%>%
                mutate(idx=1:nrow(.))
    
    focus.list <- list(res.point$idx)
    num.cluster = 0
    f = 1
    minPts_len_tmp = 50
    cl_list <- list()

    # t in case it turns infinit
    for(t in 1:max_iteration){

        focus <- focus.list[[f]]

        mat=res.point%>%filter(idx%in%focus)

        aa=fun_hdbscan(mat,minPts_len=minPts_len_tmp)
        cl_list[[f]]  <- aa$cl
        
        # add only 1 cluter condition
        # if only 1 cluster, reduce minPts
        # until reach 20, directly return
        if(length((aa$cl$cluster)%>%unique)==1){
            if(minPts_len_tmp>30){
                minPts_len_tmp = minPts_len_tmp - 5
                next 
            }else{
              res.df <- res.point%>%group_by(data)%>%
                        summarise(counts=n(),
                                time=Mode(time),
                                clusterID=Mode(clusterID),
                                clusterPorb=mean(prob))
                
            return(list(#p=NULL,
                       # cl_list= cl_list,
                       # res.point=res.point,
                       res.df=res.df,
                       #res.cluster=NULL,
                       minPts_len=minPts_len_tmp))
            }

        }

        # change all probability
        idx <- mat$idx
        res.point[idx,c('prob')] <- aa$res.point%>%.[,c('prob')] 

        # change all non-zero cluster ID by idx
        # to avoid same clusterID, so we add the max_clusterID to current ID
        res.point[idx,'clusterID'] <- 0
        idx.sub <- aa$res.point%>%filter(clusterID!=0)%>%pull(idx)
        res.point[idx.sub,'clusterID'] <- aa$res.point %>%filter(clusterID!=0)%>%pull(clusterID) + num.cluster
        # res.point[idx.sub,'clusterID'] <- aa$res.point %>%filter(clusterID!=0)%>%pull(prob)


        # update max_cluster num
        num.cluster = max(res.point$clusterID)
        # add time to record each point's iteration
        res.point[idx,'time'] <-   res.point[idx,'time']+1

        res.df <- res.point%>%group_by(data)%>%
                  summarise(counts=n(),
                            time=Mode(time),
                            clusterID=Mode(clusterID),
                            clusterPorb=mean(prob))
        res.cluster <- res.df%>%group_by(clusterID)%>%summarise(num=n(),start=min(data),end=max(data))%>%
                        mutate(length=end-start+1,check=length-num)#%>%filter(clusterID!=0)

        dip.test <- c()
        site.len <- c()
        id.list <- (res.cluster)%>%filter(clusterID!=0)%>%pull(clusterID)
        for(j in id.list){
            # res.cluster%>%filter(clusterID!=0)%>%pull(check)
            dip.test[[as.character(j)]] <- res.point%>%filter(clusterID==j)%>%pull(data)%>%dip.test%>%.$p.value
            site.len[[as.character(j)]] <- (res.df%>%filter(clusterID==j)%>%nrow)
        }

        focus <- lapply(id.list[intersect(which(dip.test < dip_p_th),
                                          which(site.len > min_split_len))],
                function(x){res.point%>%filter(clusterID==x)%>%pull(idx)})
        if(length(focus)){
            focus.list <- c(focus.list,focus)
        }
        

        f=f+1
        if(f > length(focus.list)) {break}


         message(f)
     }

    # res.df <- res.point%>%group_by(data)%>%
    #           summarise(counts=n(),
    #           clusterID=Mode(clusterID),
    #           clusterPorb=mean(prob))
    
    # plot distribution and density scatter
    res.df$clusterID <- as.character(res.df$clusterID)
    p1 <- (ggplot(res.df,aes(x = data,y=counts,group=clusterID))+
        geom_path(aes(color=(clusterID),alpha=clusterPorb),
                  lwd=2,group=1)+
        geom_ribbon(res.df%>%filter(clusterID!=0),
                    mapping = aes(x = data,y=counts,ymin=0,ymax=counts,fill=clusterID),
                    alpha=0.2)+
        scale_fill_manual(values=c(rep(colors_,20)))+
        scale_color_manual(values=c('black',rep(colors_,20)))+
        scale_y_continuous(expand = c(0, 0))+
        theme_classic()+
        xlab('')+
        ggtitle(paste0('chr',chr_,'_',i))+
        theme(aspect.ratio = 2/(nrow(res.df)/200),axis.text.x = element_blank(),
              legend.position = 'none'))%>%
        fun_ggtheme

    p2 <- (ggplot(res.point)+
            geom_bin2d(aes(x=data,y=y),bins=100)+
            scale_fill_continuous(type = "viridis")+

            scale_y_continuous(expand = c(0, 0))+
            theme_bw()+
            theme(aspect.ratio = 2/(nrow(res.df)/200),
                  legend.position = 'none'))%>%
            fun_ggtheme

    p=p1/p2
    p

    ggsave(paste0(dir_,'/chr',chr_,'/',i,'.pdf'),
          height = 4,width=15,dpi=80)


    # return plot, cluster results, matrix output
    return(list(# p=p1/p2,
               # cl_list= cl_list,
               # res.point=res.point,
               res.df=res.df,
               #res.cluster=res.cluster,
               minPts_len=minPts_len_tmp))

}

############## 运行函数 #########################
fun_all <- function(i,min_split_len=800,dip_p_th=0.1,max_peaksum=5e4){
    
    tryCatch({
        set.seed(0)
        
        samples <- rep(1:length(res.chr[[i]]$num),times = (res.chr[[i]]$num))

        v <- res.chr[[i]]$num
        s <- res.chr[[i]]$site
        data=rep(s,times = v)
        peak.sum <- sum(res.chr[[i]]$num) # peak sum/area
        scale_=1

        if(peak.sum > max_peaksum){  

            scale_=2 
            for(ii in c(1:10)){
                v <- (res.chr[[i]]$num/scale_)%>% floor
                data= rep(s,times = v)
                if(length(data)< max_peaksum) {break}
                scale_ = scale_ +1
                }
        }
        peak.multimodal.p <- dip.test(data)$p.value  

        if(peak.multimodal.p > dip_p_th | nrow(res.chr[[i]])<min_split_len){

            return(list(chr=chr_,
                        id=i,
                        scale_=scale_,
                        type='unimodal',
                        len=nrow(res.chr[[i]]),
                        num.reads=peak.sum,
                        dip.test=peak.multimodal.p,
                       res.peak=NULL))
        }else{

            res.peak <- fun_perChr(i,s,v,
                                  min_split_len=min_split_len,
                                  dip_p_th=dip_p_th)
                return(list(chr=chr_,
                        id=i,
                        scale_=scale_,
                        type='multimodal',
                        len=nrow(res.chr[[i]]),
                        num.reads=peak.sum,
                        dip.test=peak.multimodal.p,
                        res.peak=res.peak))
        }
    },#warning = function(w) {message("Warning!")},
      error = function(e) {message("Error!")}) 
}

############################################# peak split end ##########################################

fun_all_test <- function(i,min_split_len=800,dip_p_th=0.1,max_peaksum=5e4){
    
    #tryCatch({
        set.seed(0)
        
        samples <- rep(1:length(res.chr[[i]]$num),times = (res.chr[[i]]$num))

        v <- res.chr[[i]]$num
        s <- res.chr[[i]]$site
        data=rep(s,times = v)
        peak.sum <- sum(res.chr[[i]]$num) # peak sum/area
        scale_=1

        if(peak.sum > max_peaksum){  

            scale_=2 
            for(ii in c(1:10)){
                v <- (res.chr[[i]]$num/scale_)%>% floor
                data= rep(s,times = v)
                if(length(data)< max_peaksum) {break}
                scale_ = scale_ +1
                }
        }
        peak.multimodal.p <- dip.test(data)$p.value  

        if(peak.multimodal.p > dip_p_th | nrow(res.chr[[i]])<min_split_len){

            return(list(chr=chr_,
                        id=i,
                        scale_=scale_,
                        type='unimodal',
                        len=nrow(res.chr[[i]]),
                        num.reads=peak.sum,
                        dip.test=peak.multimodal.p,
                       res.peak=NULL))
        }else{

            res.peak <- fun_perChr(i,s,v,
                                  min_split_len=min_split_len,
                                  dip_p_th=dip_p_th)
                return(list(chr=chr_,
                        id=i,
                        scale_=scale_,
                        type='multimodal',
                        len=nrow(res.chr[[i]]),
                        num.reads=peak.sum,
                        dip.test=peak.multimodal.p,
                        res.peak=res.peak))
        }
#    },#warning = function(w) {message("Warning!")},
#      error = function(e) {message("Error!")}) 
}

############################################# peak split end ##########################################


########################################### bed generation ###############################################
cal_breaks <- function(res.df){
    

    
    res.cluster <- res.df%>%group_by(clusterID)%>%summarise(num=n(),start=min(data),end=max(data))%>%
                    mutate(length=end-start+1,check=length-num)#%>%filter(clusterID!=0)

    res.cluster.filter = res.cluster%>%filter(clusterID!=0)%>%
        mutate(center=(start+end)/2)%>%arrange(center)
    
    if(nrow(res.cluster.filter)<2){
        sites=NULL
    }else{
        sites <- c()
        for(i in 2:nrow(res.cluster.filter)){
            tmp = res.df%>%filter(data%in%c(round(unlist(res.cluster.filter[i-1,'center'])): 
                                            round(unlist(res.cluster.filter[i,'center']))))
            sites =c(tmp%>%filter(clusterPorb==min(tmp$clusterPorb))%>%
                     pull(data)%>%mean%>%round, sites) 
        }
        
    }


    return(sites)
}

cal_bed <- function(i){
    #tryCatch({
    
    # 这里有个bug修复，就是有些极端情况只返回了res.point而不是res.df
    # 所以需要写一下转换
    if(length(res.all[[i]]$res.peak) != 0){
        if(length(res.all[[i]]$res.peak$res.point)!=0)
            res.all[[i]]$res.peak$res.df <-  res.all[[i]]$res.peak$res.point%>%group_by(data)%>%
                                              summarise(counts=n(),
                                                        time=Mode(time),
                                                        clusterID=Mode(clusterID),
                                                        clusterPorb=mean(prob))
            res.all[[i]]$res.peak$res.point <- NULL
            res.all[[i]]$res.peak$cl_list <- NULL
            
    }
    if(length(res.all[[i]]$res.peak$res.df)==0){ # equal to 'type==unimodaal'

    #         sites=c(res.chr[[as.numeric(i)]]%>%pull(site)%>%min,
    #                 res.chr[[as.numeric(i)]]%>%pull(site)%>%max)

            return(list(cuts=NULL,
                        start=res.chr[[as.numeric(i)]]%>%pull(site)%>%min,
                        end=res.chr[[as.numeric(i)]]%>%pull(site)%>%max,
                        split_id=as.numeric(i),
                        chr=res.all[[(i)]]$chr,
                        type='unimodel'))
        }else{ # multimodal

            cuts=cal_breaks(res.all[[(i)]]$res.peak$res.df)
            # sites=c(res.chr[[as.numeric(i)]]%>%pull(site)%>%min,
            #         res.chr[[as.numeric(i)]]%>%pull(site)%>%max,
                    # cuts)

            if(length(cuts)==0){
                    return(list(cuts=NULL,
                    start=res.chr[[as.numeric(i)]]%>%pull(site)%>%min,
                    end=res.chr[[as.numeric(i)]]%>%pull(site)%>%max,
                    split_id=as.numeric(i),
                    chr=res.all[[(i)]]$chr,
                    type='unimodel'))
            }else{
                return(list(cuts=cuts,
                start=res.chr[[as.numeric(i)]]%>%pull(site)%>%min,
                end=res.chr[[as.numeric(i)]]%>%pull(site)%>%max,
                split_id=as.numeric(i),
                chr=res.all[[(i)]]$chr,
                type='multimodal'))
            }


        }
        
    #},warning = function(w) {message("Warning!")},
     #   error = function(e) {message("Error!")}) 

}

#################### bed generation end ##############################