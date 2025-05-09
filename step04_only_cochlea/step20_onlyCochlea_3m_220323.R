
#/home/toolkit/tools/R4*/bin/R


setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS')



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

pbmc=readRDS( '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/COM_SEURAT_SUBCLST.rds')

DimPlot(pbmc, reduction='umap', group.by='new.clst.sub',pt.size=0.1,label=T,label.size=8)  +NoLegend()

Idents(pbmc)=pbmc$new.clst.sub


pbmc.3m=subset(pbmc, cells=colnames(pbmc)[which(pbmc$batch=='cochlea.3m')])
pbmc.12m=subset(pbmc, cells=colnames(pbmc)[which(pbmc$batch=='cochlea.12m')])
pbmc.24m=subset(pbmc, cells=colnames(pbmc)[which(pbmc$batch %in% c('cochlea.24m','cochlea.24mbu'))])


pbmc.3m.markers <- FindAllMarkers(pbmc.3m, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.12m.markers <- FindAllMarkers(pbmc.12m, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.24m.markers <- FindAllMarkers(pbmc.24m, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




saveRDS(pbmc.3m.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step20_onlyCochlea_3m/Cochlea_MAKER_3M.rds')
saveRDS(pbmc.12m.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step20_onlyCochlea_3m/Cochlea_MAKER_12M.rds')
saveRDS(pbmc.24m.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step20_onlyCochlea_3m/Cochlea_MAKER_24M.rds')

write.table(pbmc.3m.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step20_onlyCochlea_3m/Cochlea_MAKER_3M.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(pbmc.12m.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step20_onlyCochlea_3m/Cochlea_MAKER_12M.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(pbmc.24m.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step20_onlyCochlea_3m/Cochlea_MAKER_24M.tsv',col.names=T,row.names=F,sep='\t',quote=F)







