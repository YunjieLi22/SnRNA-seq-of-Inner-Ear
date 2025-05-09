
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

pbmc=readRDS( '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_20211025/UTRICLE_SEURAT_FINAL.rds')
PCUSE=readRDS( '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_PCUSE.rds')

DimPlot(pbmc, reduction='umap', group.by='final.clst.sub.label',pt.size=0.1,label=T,label.size=8)  +NoLegend()

Idents(pbmc)=pbmc$final.clst.sub.label


pbmc.3m=subset(pbmc, cells=colnames(pbmc)[which(pbmc$batch=='utricle.3m')])
pbmc.12m=subset(pbmc, cells=colnames(pbmc)[which(pbmc$batch=='utricle.12m')])
pbmc.22m=subset(pbmc, cells=colnames(pbmc)[which(pbmc$batch=='utricle.22m')])


pbmc.3m.markers <- FindAllMarkers(pbmc.3m, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.12m.markers <- FindAllMarkers(pbmc.12m, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.22m.markers <- FindAllMarkers(pbmc.22m, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




saveRDS(pbmc.3m.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step19_onlyUtricle_4type3m/UTRICLE_MAKER_3M.rds')
saveRDS(pbmc.12m.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step19_onlyUtricle_4type3m/UTRICLE_MAKER_12M.rds')
saveRDS(pbmc.22m.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step19_onlyUtricle_4type3m/UTRICLE_MAKER_22M.rds')

write.table(pbmc.3m.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step19_onlyUtricle_4type3m/UTRICLE_MAKER_3M.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(pbmc.12m.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step19_onlyUtricle_4type3m/UTRICLE_MAKER_12M.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(pbmc.22m.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step19_onlyUtricle_4type3m/UTRICLE_MAKER_22M.tsv',col.names=T,row.names=F,sep='\t',quote=F)








