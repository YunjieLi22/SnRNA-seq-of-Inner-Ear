

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


library('ComplexHeatmap')
library('circlize')



source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')




#####################################################################
setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging')


pbmc=readRDS( '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_20211025/UTRICLE_SEURAT_FINAL.rds')
PCUSE=readRDS( '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_PCUSE.rds')

DimPlot(pbmc, reduction='umap', group.by='final.clst.sub.label',pt.size=0.1,label=T,label.size=8)  +NoLegend()


TIME=pbmc$batch
TYPE=pbmc$final.clst.sub.label


#########################################################################
#Fibrocyte.1'
#############################################################################


this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$final.clst.sub.label=='Fibrocyte.1')])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)


saveRDS(this_pbmc, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Utri.Fibrocyte.1.rds')





pdf('TRA_Utri_Fibrocyte.1.pdf',width=3.5,height=4)


  DimPlot(this_pbmc, group.by='batch',label=TRUE)

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
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
write.table(sort(this_cor), paste0('PseudotimeGeneCor_Utri_Fibrocyte.1.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('PseudotimeGO_POS_Utri_Fibrocyte.1_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('PseudotimeGO_NEG_Utri_Fibrocyte.1_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################







#########################################################################
#Fibrocyte.2
#############################################################################

this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$final.clst.sub.label=='Fibrocyte.2')])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)


saveRDS(this_pbmc, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Utri.Fibrocyte.2.rds')





pdf('TRA_Utri_Fibrocyte.2.pdf',width=3.5,height=4)


  DimPlot(this_pbmc, group.by='batch',label=TRUE)

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
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
write.table(sort(this_cor), paste0('PseudotimeGeneCor_Utri_Fibrocyte.2.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('PseudotimeGO_POS_Utri_Fibrocyte.2_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('PseudotimeGO_NEG_Utri_Fibrocyte.2_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################








#########################################################################
#Fibrocyte.3
#############################################################################

TAG='Fibrocyte.3'

this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$final.clst.sub.label==TAG)])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)


saveRDS(this_pbmc, file=paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Utri.',TAG,'.rds'))





pdf(paste0('TRA_Utri_',TAG,'.pdf'),width=3.5,height=4)


  DimPlot(this_pbmc, group.by='batch',label=TRUE)

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
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
write.table(sort(this_cor), paste0('PseudotimeGeneCor_Utri_',TAG,'.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('PseudotimeGO_POS_Utri_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('PseudotimeGO_NEG_Utri_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################






#########################################################################
#Stroma.cell
#############################################################################

TAG='Stroma.cell'

this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$final.clst.sub.label==TAG)])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)


saveRDS(this_pbmc, file=paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Utri.',TAG,'.rds'))





pdf(paste0('TRA_Utri_',TAG,'.pdf'),width=3.5,height=4)


  DimPlot(this_pbmc, group.by='batch',label=TRUE)

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
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
write.table(sort(this_cor), paste0('PseudotimeGeneCor_Utri_',TAG,'.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('PseudotimeGO_POS_Utri_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('PseudotimeGO_NEG_Utri_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################











##############################################


pbmc=readRDS( '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/COM_SEURAT_SUBCLST.rds')
ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(pbmc), 'ENTREZID', 'SYMBOL')


#########################
TYPE=pbmc$new.clst.sub
############################
TIME=pbmc$batch
TIME[which(pbmc$batch=='cochlea.24mbu')]='cochlea.24m'
LABEL=paste0(TIME,'.',TYPE)








#########################################################################
#0
#############################################################################

TAG='0'

this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$new.clst.sub==TAG)])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)


saveRDS(this_pbmc, file=paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Coch.',TAG,'.rds'))





pdf(paste0('TRA_Coch_',TAG,'.pdf'),width=3.5,height=4)


  DimPlot(this_pbmc, group.by='batch',label=TRUE)

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident %in% c('cochlea.24m','cochlea.24mbu'))]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
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
write.table(sort(this_cor), paste0('PseudotimeGeneCor_Coch_',TAG,'.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('PseudotimeGO_POS_Coch_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('PseudotimeGO_NEG_Coch_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################







#########################################################################
#2
#############################################################################

TAG='2'

this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$new.clst.sub==TAG)])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)


saveRDS(this_pbmc, file=paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Coch.',TAG,'.rds'))





pdf(paste0('TRA_Coch_',TAG,'.pdf'),width=3.5,height=4)


  DimPlot(this_pbmc, group.by='batch',label=TRUE)

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident %in% c('cochlea.24m','cochlea.24mbu'))]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
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
write.table(sort(this_cor), paste0('PseudotimeGeneCor_Coch_',TAG,'.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('PseudotimeGO_POS_Coch_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('PseudotimeGO_NEG_Coch_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################




#########################################################################
#9.1
#############################################################################

TAG='9.1'

this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$new.clst.sub==TAG)])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)


saveRDS(this_pbmc, file=paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Coch.',TAG,'.rds'))





pdf(paste0('TRA_Coch_',TAG,'.pdf'),width=3.5,height=4)


  DimPlot(this_pbmc, group.by='batch',label=TRUE)

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident %in% c('cochlea.24m','cochlea.24mbu'))]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
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
write.table(sort(this_cor), paste0('PseudotimeGeneCor_Coch_',TAG,'.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('PseudotimeGO_POS_Coch_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('PseudotimeGO_NEG_Coch_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################





#########################################################################
#9.2
#############################################################################

TAG='9.2'

this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$new.clst.sub==TAG)])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)


saveRDS(this_pbmc, file=paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Coch.',TAG,'.rds'))





pdf(paste0('TRA_Coch_',TAG,'.pdf'),width=3.5,height=4)


  DimPlot(this_pbmc, group.by='batch',label=TRUE)

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident %in% c('cochlea.24m','cochlea.24mbu'))]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
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
write.table(sort(this_cor), paste0('PseudotimeGeneCor_Coch_',TAG,'.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('PseudotimeGO_POS_Coch_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('PseudotimeGO_NEG_Coch_',TAG,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################