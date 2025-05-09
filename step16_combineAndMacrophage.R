
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

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')




##############
#14

############################
TYPE=pbmc.cochlea$new.clst.sub
this_pbmc=subset(pbmc.cochlea, cells=colnames(pbmc.cochlea)[which(TYPE=='14')])

    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)

saveRDS(this_pbmc, file='./pbmc.14.rds')


this_pbmc=readRDS(file='./pbmc.14.rds')


#########################################################


this_pbmc$time=this_pbmc$batch
this_pbmc$time[which(this_pbmc$batch=='cochlea.24mbu')]='cochlea.24m'

pdf('./TRA_14.pdf',width=3.5,height=4)


  DimPlot(this_pbmc, group.by='time',label=TRUE)

  DimPlot(this_pbmc, group.by='time',label=TRUE)+NoLegend()


this_pbmc$timev=rep(NA, length(ncol(this_pbmc)))
this_pbmc$timev[which(this_pbmc$time=='cochlea.3m')]=1 #3m
this_pbmc$timev[which(this_pbmc$time=='cochlea.12m')]=2 #12m
this_pbmc$timev[which(this_pbmc$time== 'cochlea.24m')]=3 #24m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')

    this_v=scale(this_pbmc$timev)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

    this_v=scale(pred.time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='fitted.time')



NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=F)

this_pbmc$pseudotime=OUT$P.PS
FeaturePlot(this_pbmc,features=c('pseudotime'))
######################################################
dev.off()







used=which(!is.na(this_pbmc$pseudotime))
this_cor=cor( t(as.matrix(this_pbmc@assays$RNA@data)[,used]),this_pbmc$pseudotime[used])
this_cor[which(is.na(this_cor))]=0
this_cor=this_cor[,1]
names(this_cor)=rownames(this_pbmc)
plot(sort(this_cor),type='p',pch=16)
write.table(sort(this_cor), paste0('./PseudotimeGeneCor_14.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(pbmc), 'ENTREZID', 'SYMBOL')

###############
POS.CUT=0.1
NEG.CUT= -0.1

###########################################
    this_sig_gene=names(this_cor)[which(this_cor>POS.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./PseudotimeGO_POS_14_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./PseudotimeGO_NEG_14_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################












