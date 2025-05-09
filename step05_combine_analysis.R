
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
setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage')


pbmc.utricle=readRDS( '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_20211025/UTRICLE_SEURAT_FINAL.rds')
pbmc.cochlea=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/COM_SEURAT_SUBCLST.rds')
pbmc.cochlea.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/COM_MAKER.rds')


Idents(pbmc.utricle)=pbmc.utricle$final.clst.sub.label
DimPlot(pbmc.utricle,label=TRUE)+NoLegend()

Idents(pbmc.cochlea)=pbmc.cochlea$new.clst.sub
DimPlot(pbmc.cochlea,label=TRUE)+NoLegend()






pbmc.cochlea.markers=pbmc.cochlea.markers[order(as.numeric(as.character(pbmc.cochlea.markers$cluster))),]

pbmc.cochlea.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

Idents(pbmc.cochlea)=factor(pbmc.cochlea$new.clst.sub, levels=sort(unique(as.numeric(as.character(pbmc.cochlea$new.clst.sub)))))


all.genes <- rownames(pbmc.cochlea)
pbmc.cochlea <- ScaleData(pbmc.cochlea, features = all.genes)

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/cochlea_markers_heatmap.pdf',width=10,height=23)
DoHeatmap(pbmc.cochlea, features = top10$gene) + NoLegend()
dev.off()







############################################################



D1=pbmc.utricle[['RNA']]@counts
D2=pbmc.cochlea[['RNA']]@counts




source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')


DATA=.simple_combine(D1,D2,FILL=TRUE)$combine
BATCH=c(pbmc.utricle$batch, pbmc.cochlea$batch)

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=3000, SEED=1, COMBAT=TRUE )

saveRDS(mybeer, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/BEER_UTR_COC_20211221.rds')


################################





mybeer=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/BEER_UTR_COC_20211221.rds')

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

TYPE=c(pbmc.utricle$final.clst.sub.label, pbmc.cochlea$new.clst.sub )

pbmc$type=TYPE
Idents(pbmc)=TYPE


DimPlot(pbmc, reduction='umap', group.by='type', pt.size=0.1, label=TRUE) +NoLegend()



pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/combineUMAP.pdf',width=10,height=10)
DimPlot(pbmc, reduction='umap', group.by='type', pt.size=0.1, label=TRUE) +NoLegend()
dev.off()

saveRDS(pbmc, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/COM_SEURAT.rds')
saveRDS(PCUSE, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/COM_PCUSE.rds')


#####################################























