
#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111')

library(Seurat)
library(dplyr)
library(patchwork)

cochlea_3m=readRDS('cochlea_3m_afterQC.rds')
cochlea_12m=readRDS('cochlea_12m_afterQC.rds')
cochlea_24m=readRDS('cochlea_24m_afterQC.rds')
cochlea_24mbu=readRDS('cochlea_24mbu_afterQC.rds')


D1=as.matrix(cochlea_3m@assays$RNA@counts)
D2=as.matrix(cochlea_12m@assays$RNA@counts)
D3=as.matrix(cochlea_24m@assays$RNA@counts)
D4=as.matrix(cochlea_24mbu@assays$RNA@counts)



source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')


D12=.simple_combine(D1,D2,FILL=TRUE)$combine
D123=.simple_combine(D12,D3,FILL=TRUE)$combine
DATA=.simple_combine(D123,D4,FILL=TRUE)$combine

BATCH=c(
        rep('cochlea.3m',ncol(cochlea_3m)),
        rep('cochlea.12m',ncol(cochlea_12m)),
        rep('cochlea.24m',ncol(cochlea_24m)),
        rep('cochlea.24mbu',ncol(cochlea_24mbu))
        )


colnames(DATA)=paste0(BATCH,'_',colnames(DATA))


saveRDS(DATA,'COM_DATA.rds')
saveRDS(BATCH,'COM_BATCH.rds')


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


saveRDS(pbmc, 'COM_SEURAT.rds')
saveRDS(PCUSE, 'COM_PCUSE.rds')

DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 


pbmc <- FindNeighbors(pbmc, dims = PCUSE)
pbmc <- FindClusters(pbmc, resolution = 0.7)


DimPlot(pbmc, reduction='umap', pt.size=0.1,label=TRUE) 


##################################
NEW.CLST=as.character(Idents(pbmc))
NEW.CLST.SUB=NEW.CLST
VEC=pbmc@reductions$umap@cell.embeddings

#################
USED=which(NEW.CLST==9)

USED_VEC=VEC[USED,]

    plot(USED_VEC, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(USED_VEC, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

NEW.CLST.SUB[USED]='9.1'
NEW.CLST.SUB[which(rownames(VEC) %in% SELECT)]='9.2'

#################
USED=which(NEW.CLST==1)

USED_VEC=VEC[USED,]

    plot(USED_VEC, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(USED_VEC, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

NEW.CLST.SUB[USED]='1.1'
NEW.CLST.SUB[which(rownames(VEC) %in% SELECT)]='1.2'



####################
pbmc$new.clst.sub=NEW.CLST.SUB

DimPlot(pbmc, reduction='umap',group.by='new.clst.sub', pt.size=0.1,label=TRUE) 


######################

pbmc$flag.3m=rep(0,ncol(pbmc))
pbmc$flag.3m[which(pbmc$batch=='cochlea.3m')]=1
pbmc$flag.12m=rep(0,ncol(pbmc))
pbmc$flag.12m[which(pbmc$batch=='cochlea.12m')]=1
pbmc$flag.24m=rep(0,ncol(pbmc))
pbmc$flag.24m[which(pbmc$batch=='cochlea.24m')]=1
pbmc$flag.24mbu=rep(0,ncol(pbmc))
pbmc$flag.24mbu[which(pbmc$batch=='cochlea.24mbu')]=1
##############

saveRDS(pbmc, 'COM_SEURAT_SUBCLST.rds')

###################################

pdf('Cochlea_UMAP.pdf',width=7,height=7)
DimPlot(pbmc, reduction='umap',group.by='new.clst.sub', pt.size=0.1,label=TRUE) +NoLegend()
FeaturePlot(pbmc,features=c('flag.3m','flag.12m','flag.24m','flag.24mbu'),order=TRUE)
dev.off()

########################################


Idents(pbmc)=pbmc$new.clst.sub
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pbmc.markers, 'COM_MAKER.rds')

write.table(pbmc.markers, file='cochlea_all_markers.tsv',col.names=T,row.names=F,sep='\t',quote=F)



