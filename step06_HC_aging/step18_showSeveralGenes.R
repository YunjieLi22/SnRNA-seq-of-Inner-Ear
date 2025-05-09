

#/home/toolkit/tools/R4*/bin/R




library(Seurat)
library(dplyr)
library(patchwork)

library(ggplot2)
library(stringr)

library(DOSE)
library(GSEABase) 
library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)

library(org.Mm.eg.db)

#####################################################################
setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step18_showSeveralGenes')


pbmc=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/COM_SEURAT.rds')






library(Seurat)
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(),
               axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
               plot.margin = plot.margin )
       return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
       plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
            theme(axis.text.x=element_text(), axis.ticks.x = element_line())
       p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
       return(p)
}



my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175')


SELECT_GENE=c('Col1a2','Col2a1','Col11a1','Col11a2','Col9a1','Col23a1','Lama2','Tnf')


Idents(pbmc)=paste0(pbmc$batch,'.',pbmc$type)
StackedVlnPlot(pbmc, c(), pt.size=0, cols=rep(my36colors,10))


source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')

DATA=pbmc[['RNA']]@data
DATA=DATA[which(rownames(DATA) %in% SELECT_GENE),]

LABEL=paste0(pbmc$batch,'.',pbmc$type)

REF=.generate_mean(DATA,LABEL)

REF=REF[,order(colnames(REF))]


write.table(REF, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step18_showSeveralGenes/SelectGeneExp.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)

















