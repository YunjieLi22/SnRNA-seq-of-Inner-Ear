


#/home/toolkit/tools/R4*/bin/R


setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')

library(Seurat)
library(dplyr)
library(patchwork)


utricle_all=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_20211025/UTRICLE_SEURAT_FINAL.rds')
cochlea_all=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/COM_SEURAT_SUBCLST.rds')


utricle_used=subset(utricle_all,cells=colnames(utricle_all)[which(!utricle_all$new.clst.sub %in% c('6.1','6.2'))])
cochlea_used=subset(cochlea_all,cells=colnames(cochlea_all)[which(!cochlea_all$new.clst.sub %in% c('3'))])


utricle_used$time=utricle_used$batch
utricle_used$time[which(utricle_used$time %in% c('utricle.3m'))]='3m'
utricle_used$time[which(utricle_used$time %in% c('utricle.12m'))]='12m'
utricle_used$time[which(utricle_used$time %in% c('utricle.22m'))]='22m'


cochlea_used$time=cochlea_used$batch
cochlea_used$time[which(cochlea_used$time %in% c('cochlea.3m'))]='3m'
cochlea_used$time[which(cochlea_used$time %in% c('cochlea.12m'))]='12m'
cochlea_used$time[which(cochlea_used$time %in% c('cochlea.24m','cochlea.24mbu'))]='24m'


DimPlot(utricle_used, reduction='umap', group.by='new.clst.sub',pt.size=0.1,label=T,label.size=8)  +NoLegend()
DimPlot(cochlea_used, reduction='umap', group.by='new.clst.sub',pt.size=0.1,label=T,label.size=8)  +NoLegend()

#########################################################################################################################

source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')

COUNT1=utricle_used@assays$RNA@counts
COUNT2=cochlea_used@assays$RNA@counts

BATCH=c(utricle_used$batch, cochlea_used$batch)
TYPE=c(paste0('U.',utricle_used$new.clst.sub),paste0('C.',cochlea_used$new.clst.sub))
TIME=c(paste0('U.',utricle_used$time),paste0('C.',cochlea_used$time))

COM=.simple_combine(COUNT1,COUNT2)$combine

mybeer=BEER(DATA=COM, BATCH=BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=3000, SEED=1, COMBAT=TRUE )

########################
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

###################################################
pbmc <- mybeer$seurat
PCUSE <- mybeer$select
pbmc <- RunUMAP(object = pbmc, reduction='pca',dims = PCUSE, check_duplicates=FALSE)

DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 

pbmc$clst=TYPE
pbmc$time=TIME

DimPlot(pbmc, reduction='umap', group.by='clst', pt.size=0.1, label=TRUE) + NoLegend()
DimPlot(pbmc, reduction='umap', group.by='time', pt.size=0.1, label=TRUE) + NoLegend()

################################################################################################################
saveRDS(mybeer, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_BEER.rds')
################################################################################################################
saveRDS(PCUSE, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_PCUSE.rds')
saveRDS(pbmc, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')
################################################################################################################
saveRDS(utricle_used, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
saveRDS(cochlea_used, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
################################################################################################################











