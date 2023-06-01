# naccessary package
suppressPackageStartupMessages({
    library(tidyverse)
    library(dplyr)
    library(parallel)
    #library(progress)
    library(qs)
    library(diptest)
    library(ggplot2)
    library(ggplotify)
    library(patchwork)
    library(dbscan)
    library(HiClimR)
    library(pheatmap)
    library(mixtools)
    library(gridExtra)
    library("ape")
    library('dendextend')
})

# figure size
psize <- function (w = 6, h = 6) 
{
    options(repr.plot.width = w, repr.plot.height = h)
}

# npg colors
colors_ <- ggsci::pal_npg()(10)

# ggplot theme
fun_ggtheme <- function (p, aspect.ratio = 1, legend_point_size = 1, legend_size = 8,angle.x=0,angle.y=0) {
    pt <- p + guides(shape = guide_legend(override.aes = list(size = legend_point_size)), 
        color = guide_legend(override.aes = list(size = legend_point_size))) + 
        theme_classic() + theme(aspect.ratio = aspect.ratio, 
        legend.title = element_text(size = legend_size), legend.text = element_text(size = legend_size), 
        axis.text.x = element_text(size = 8, color = "black", face = "bold",angle = angle.x), 
        axis.text.y = element_text(size = 8, color = "black", face = "bold",angle = angle.y), 
        axis.title = element_text(size = 8, color = "black", 
            face = "bold"), plot.title = element_text(hjust = 0.5, 
            size = 10, color = "black", face = "bold"))
    return(pt)
}

# Mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}