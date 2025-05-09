

#/home/toolkit/tools/R4*/bin/R


setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS')



library(Seurat)
library(dplyr)
library(patchwork)


pbmc_all=readRDS('COM_SEURAT.rds')

BATCH=pbmc_all@meta.data$batch
CLST=readRDS( 'COM_SUBCLST1.rds')


USED=which(BATCH %in% c('utricle.12m','utricle.22m','utricle.3m'))

DATA=as.matrix(pbmc_all@assays$RNA@counts)[,USED]
BATCH=BATCH[USED]
CLST=CLST[USED]


source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=3000, SEED=1, COMBAT=TRUE )

saveRDS(mybeer, 'UTRICLE_BEER.rds')


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



pbmc <- FindNeighbors(pbmc, dims = PCUSE)
#pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- FindClusters(pbmc, resolution = 0.7)


pbmc@meta.data$ori.clst=CLST
pbmc@meta.data$new.clst=Idents(pbmc)

saveRDS(pbmc, 'UTRICLE_SEURAT.rds')
saveRDS(PCUSE, 'UTRICLE_PCUSE.rds')

DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 
DimPlot(pbmc, reduction='umap', group.by='ori.clst', pt.size=0.1,label=T) 
DimPlot(pbmc, reduction='umap', group.by='new.clst', pt.size=0.1,label=T) 


pdf('UTRICLE_UMAP.pdf',width=12,height=12)
DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 
DimPlot(pbmc, reduction='umap', group.by='ori.clst', pt.size=0.1,label=T) 
DimPlot(pbmc, reduction='umap', group.by='new.clst', pt.size=0.1,label=T) 
dev.off()


TAB=table(pbmc@meta.data$ori.clst,pbmc@meta.data$new.clst)
write.table(TAB, file='UTRICLE_NEW_OLD_CLST_TABLE.tsv',col.names=T,row.names=T,sep='\t',quote=F)



pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pbmc.markers, 'UTRICLE_MAKER.rds')

write.table(pbmc.markers, file='UTRICLE_MAKER.tsv',col.names=T,row.names=F,sep='\t',quote=F)
