
#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111')

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


pbmc=readRDS('COM_SEURAT_SUBCLST.rds')
ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(pbmc), 'ENTREZID', 'SYMBOL')
#######################



#########################
TYPE=pbmc$new.clst.sub
############################
TIME=pbmc$batch
TIME[which(pbmc$batch=='cochlea.24mbu')]='cochlea.24m'
LABEL=paste0(TIME,'.',TYPE)




#############################################
Idents(pbmc)=TIME
pbmc$time=TIME


############################
this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(TYPE=='8')])

    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)

saveRDS(this_pbmc, file='./pbmc.8.rds')


this_pbmc=readRDS(file='./pbmc.8.rds')


#########################################################

pdf('./TRA_8.pdf',width=3.5,height=4)


  DimPlot(this_pbmc, group.by='time',label=TRUE)

  DimPlot(this_pbmc, group.by='time',label=TRUE)+NoLegend()


this_pbmc$timev=rep(NA, length(ncol(this_pbmc)))
this_pbmc$timev[which(this_pbmc$time=='cochlea.3m')]=1 #3m
this_pbmc$timev[which(this_pbmc$time=='cochlea.12m')]=2 #12m
this_pbmc$timev[which(this_pbmc$time=='cochlea.24m')]=3 #24m
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
write.table(sort(this_cor), paste0('./PseudotimeGeneCor_8.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('./PseudotimeGO_POS_8_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./PseudotimeGO_NEG_8_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################












##############
#11

############################
this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(TYPE=='11')])

    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)

saveRDS(this_pbmc, file='./pbmc.11.rds')


this_pbmc=readRDS(file='./pbmc.11.rds')


#########################################################

pdf('./TRA_11.pdf',width=3.5,height=4)


  DimPlot(this_pbmc, group.by='time',label=TRUE)

  DimPlot(this_pbmc, group.by='time',label=TRUE)+NoLegend()


this_pbmc$timev=rep(NA, length(ncol(this_pbmc)))
this_pbmc$timev[which(this_pbmc$time=='cochlea.3m')]=1 #3m
this_pbmc$timev[which(this_pbmc$time=='cochlea.12m')]=2 #12m
this_pbmc$timev[which(this_pbmc$time=='cochlea.24m')]=3 #24m
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
write.table(sort(this_cor), paste0('./PseudotimeGeneCor_11.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('./PseudotimeGO_POS_11_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./PseudotimeGO_NEG_11_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################


