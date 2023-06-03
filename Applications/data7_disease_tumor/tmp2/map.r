fun_map <- function(a,b){
    
    regions1 = a%>%as.data.frame%>%rename_with(~c('region'))%>%mutate(chr=str_extract(region,'^[^-]+'),
                                               start=str_extract(region,'(?<=-)[^-]*(?=-)'),
                                               end=str_extract(region,'(?<=-)[^-]*$'))%>%
                mutate(start=as.numeric(start), 
                       end=as.numeric(end),
                       id=row_number())
    regions2 = b%>%as.data.frame%>%rename_with(~c('region'))%>%mutate(chr=str_extract(region,'^[^-]+'),
                                                   start=str_extract(region,'(?<=-)[^-]*(?=-)'),
                                                   end=str_extract(region,'(?<=-)[^-]*$'))%>%
                    mutate(start=as.numeric(start), 
                           end=as.numeric(end),
                           id=row_number())
    
    source('../../map_to_cPeak_v1/map_function2.r')
    res = fun_map_bed_customRef(regions1,regions2)
    res.df=res$df.trans%>%mutate(regions1=paste0(chr,'-',start_q,'-',end_q),
                     regions2=paste0(chr,'-',start,'-',end))%>%dplyr::select(regions1,regions2)
    
    return(list(res.df=res.df, res=res))
}