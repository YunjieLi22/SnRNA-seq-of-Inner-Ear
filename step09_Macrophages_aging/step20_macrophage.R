

# /home/toolkit/tools/R4.0.3/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')

library(Seurat)
library(dplyr)
library(patchwork)

##################################
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')


UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')


COMBINE$time[which(COMBINE$time %in% c('U.3m'))]='U.03m'
COMBINE$time[which(COMBINE$time %in% c('C.3m'))]='C.03m'
COUNT=COMBINE@assays$RNA@counts

#################################################

TAG=paste0(COMBINE$clst,'_',COMBINE$time)
TISSUE=c(rep('U',ncol(UTRICLE)),rep('C',ncol(COCHLEA)))






#################
# COCHLEA
#################

pbmc=COCHLEA
pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$new.clst.sub %in% c(14))])
DimPlot(pbmc)

VEC=pbmc@reductions$umap@cell.embeddings
pbmc=subset(pbmc, cells=colnames(pbmc)[which(VEC[,1]<0 & VEC[,2]>3)])
DimPlot(pbmc)


VEC=pbmc@reductions$umap@cell.embeddings
set.seed(123)
KM=kmeans(VEC,center=6)
Idents(pbmc)=factor(KM$cluster,levels=sort(unique(KM$cluster)))
DimPlot(pbmc, reduction = "umap",label=T)+NoLegend()


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


########################################################

write.table(pbmc.markers,
file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/OUT20220523/MACROPHAGE_COCHLEA_MARKER.txt',
quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/OUT20220523/MACROPHAGE_COCHLEA_PLOT.pdf',width=7,height=7)
DimPlot(pbmc, reduction = "umap",label=T)+NoLegend()
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()


################################################










#################
# UTRICLE
#################

pbmc=UTRICLE
pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$new.clst.sub %in% c(16.1))])
DimPlot(pbmc)

VEC=pbmc@reductions$umap@cell.embeddings
pbmc=subset(pbmc, cells=colnames(pbmc)[which(VEC[,1]< -5 & VEC[,2]< -4)])
DimPlot(pbmc)


VEC=pbmc@reductions$umap@cell.embeddings
set.seed(123)
KM=kmeans(VEC,center=5)
Idents(pbmc)=factor(KM$cluster,levels=sort(unique(KM$cluster)))
DimPlot(pbmc, reduction = "umap",label=T)+NoLegend()


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


########################################################

write.table(pbmc.markers,
file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/OUT20220523/MACROPHAGE_UTRICLE_MARKER.txt',
quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/OUT20220523/MACROPHAGE_UTRICLE_PLOT.pdf',width=7,height=7)
DimPlot(pbmc, reduction = "umap",label=T)+NoLegend()
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()


################################################
