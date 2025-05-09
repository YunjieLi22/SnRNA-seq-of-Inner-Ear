
#/home/toolkit/tools/R4*/bin/R


setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS')

library(Seurat)
library(dplyr)
library(patchwork)


utricle_3m=readRDS('utricle_3m_SCT_UMAP.rds')
utricle_12m=readRDS('utricle_12m_SCT_UMAP.rds')
utricle_22m=readRDS('utricle_22m_SCT_UMAP.rds')

cochlea_3m=readRDS('cochlea_3m_SCT_UMAP.rds')
cochlea_12m=readRDS('cochlea_12m_SCT_UMAP.rds')


#################################################################



D1=as.matrix(utricle_3m@assays$RNA@counts)
D2=as.matrix(utricle_12m@assays$RNA@counts)
D3=as.matrix(utricle_22m@assays$RNA@counts)

D4=as.matrix(cochlea_3m@assays$RNA@counts)
D5=as.matrix(cochlea_12m@assays$RNA@counts)




source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')


D12=.simple_combine(D1,D2,FILL=TRUE)$combine
D123=.simple_combine(D12,D3,FILL=TRUE)$combine
D1234=.simple_combine(D123,D4,FILL=TRUE)$combine
DATA=.simple_combine(D1234,D5,FILL=TRUE)$combine


BATCH=c(rep('utricle.3m',ncol(utricle_3m)),
        rep('utricle.12m',ncol(utricle_12m)),
        rep('utricle.22m',ncol(utricle_22m)),
        rep('cochlea.3m',ncol(cochlea_3m)),
        rep('cochlea.12m',ncol(cochlea_12m))
        )


colnames(DATA)=paste0(BATCH,'_',colnames(DATA))

CLST=c(
    paste0('utricle.3m.',as.character(Idents(utricle_3m))),
    paste0('utricle.12m.',as.character(Idents(utricle_12m))),
    paste0('utricle.22m.',as.character(Idents(utricle_22m))),
    paste0('cochlea.3m.',as.character(Idents(cochlea_3m))),
    paste0('cochlea.12m.',as.character(Idents(cochlea_12m)))
    )


saveRDS(DATA,'COM_DATA.rds')
saveRDS(BATCH,'COM_BATCH.rds')
saveRDS(CLST,'COM_CLST.rds')

###########################################################

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=3000, SEED=1, COMBAT=TRUE )

saveRDS(mybeer, 'COM_BEER.rds')


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

pbmc@meta.data$cluster=CLST

saveRDS(pbmc, 'COM_SEURAT.rds')
saveRDS(PCUSE, 'COM_PCUSE.rds')

DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 



##############################################


pbmc <- FindNeighbors(pbmc, dims = PCUSE)
#pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- FindClusters(pbmc, resolution = 0.7)

all.genes=rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

DimPlot(pbmc, reduction='umap', group.by='seurat_clusters', pt.size=0.1, label=T)  + NoLegend()

saveRDS(pbmc, 'COM_SEURAT_CLST.rds')

saveRDS(pbmc$seurat_clusters, 'COM_SEURAT_CLST_ID.rds')
saveRDS(pbmc@meta.data, 'COM_SEURAT_CLST_META.rds')


pdf('COM_UMAP.pdf',width=7,height=7)
DimPlot(pbmc, reduction='umap', group.by='seurat_clusters', pt.size=0.1, label=T)  + NoLegend()
dev.off()



pdf('COM_BATCH_UMAP.pdf',width=9,height=7)
DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1, label=F)  
dev.off()

##################################################

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pbmc.markers, 'COM_MAKER.rds')


##################################################

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)



pdf('COM_HEAT.pdf',width=30,height=30)
this_plot=DoHeatmap(pbmc, features = top10 $gene) + NoLegend()
print(this_plot)
dev.off()















#######################################

#/home/toolkit/tools/R4*/bin/R


setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS')


library(Seurat)
library(dplyr)
library(patchwork)



pbmc=readRDS('COM_SEURAT_CLST.rds')

UMAP=pbmc@reductions$umap@cell.embeddings
CLST=as.character(pbmc$RNA_snn_res.0.7)

###########################
C11.INDEX=which(CLST=='11')

COL=rep('grey80',nrow(UMAP))
COL[C11.INDEX]='red'

plot(UMAP, cex=0.5, pch=16, col=COL)

C11.UMAP=UMAP[C11.INDEX,]

set.seed(321)
C11.CLST=kmeans(C11.UMAP,center=4)$cluster
plot(C11.UMAP, col=as.factor(C11.CLST))
##################################


C19.INDEX=which(CLST=='19')

COL=rep('grey80',nrow(UMAP))
COL[C19.INDEX]='red'

plot(UMAP, cex=0.5, pch=16, col=COL)

C19.UMAP=UMAP[C19.INDEX,]

set.seed(321)
C19.CLST=kmeans(C19.UMAP,center=2)$cluster
plot(C19.UMAP, col=as.factor(C19.CLST))


###########################################
NEW.CLST=CLST
NEW.CLST[C11.INDEX]=paste0(CLST[C11.INDEX],'.',C11.CLST)
NEW.CLST[C19.INDEX]=paste0(CLST[C19.INDEX],'.',C19.CLST)

pbmc@meta.data$subcluster1=NEW.CLST

saveRDS(pbmc@meta.data$subcluster1, 'COM_SUBCLST1.rds')
##################################################


pdf('COM_UMAP_SUBCLST1.pdf',width=7,height=7)
DimPlot(pbmc, reduction='umap', group.by='subcluster1', pt.size=0.1, label=T)  + NoLegend()
dev.off()

########################################################

Idents(pbmc)=pbmc@meta.data$subcluster1

C11.1=FindMarkers(pbmc, ident.1='11.1', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
C11.2=FindMarkers(pbmc, ident.1='11.2', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
C11.3=FindMarkers(pbmc, ident.1='11.3', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
C11.4=FindMarkers(pbmc, ident.1='11.4', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
C19.1=FindMarkers(pbmc, ident.1='19.1', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
C19.2=FindMarkers(pbmc, ident.1='19.2', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


OUT=C11.1
OUT=cbind(OUT,rep('C11.1',nrow(OUT)),rownames(OUT))
colnames(OUT)[(ncol(OUT)-1):ncol(OUT)]=c('cluster','gene')
OUT1=OUT

OUT=C11.2
OUT=cbind(OUT,rep('C11.2',nrow(OUT)),rownames(OUT))
colnames(OUT)[(ncol(OUT)-1):ncol(OUT)]=c('cluster','gene')
OUT2=OUT


OUT=C11.3
OUT=cbind(OUT,rep('C11.3',nrow(OUT)),rownames(OUT))
colnames(OUT)[(ncol(OUT)-1):ncol(OUT)]=c('cluster','gene')
OUT3=OUT

OUT=C11.4
OUT=cbind(OUT,rep('C11.4',nrow(OUT)),rownames(OUT))
colnames(OUT)[(ncol(OUT)-1):ncol(OUT)]=c('cluster','gene')
OUT4=OUT



OUT=C19.1
OUT=cbind(OUT,rep('C19.1',nrow(OUT)),rownames(OUT))
colnames(OUT)[(ncol(OUT)-1):ncol(OUT)]=c('cluster','gene')
OUT5=OUT


OUT=C19.2
OUT=cbind(OUT,rep('C19.2',nrow(OUT)),rownames(OUT))
colnames(OUT)[(ncol(OUT)-1):ncol(OUT)]=c('cluster','gene')
OUT6=OUT


OUT=rbind(OUT1,OUT2,OUT3,OUT4,OUT5,OUT6)

table(OUT[,6])

saveRDS(OUT, 'COM_SUB1_MAKER.rds')



write.table(OUT, file='COM_SUB1_MAKER.tsv',col.names=T,row.names=F,sep='\t',quote=F)

