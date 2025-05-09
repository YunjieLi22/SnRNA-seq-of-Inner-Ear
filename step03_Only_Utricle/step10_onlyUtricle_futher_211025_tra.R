
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

pbmc=readRDS( 'UTRICLE_20211025/UTRICLE_SEURAT_FINAL.rds')
PCUSE=readRDS( 'UTRICLE_PCUSE.rds')

DimPlot(pbmc, reduction='umap', group.by='final.clst.sub.label',pt.size=0.1,label=T,label.size=8)  +NoLegend()



TIME=pbmc$batch
TYPE=pbmc$final.clst.sub.label

############################################








this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$final.clst.sub.label=='Extrastriolar.SC')])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)

saveRDS(this_pbmc, file='UTRICLE_20211025/pbmc.Extrastriolar.SC.rds')


this_pbmc=readRDS(file='UTRICLE_20211025/pbmc.Extrastriolar.SC.rds')


pdf('UTRICLE_20211025/TRA_Extrastriolar.SC.pdf',width=3.5,height=4)


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
write.table(sort(this_cor), paste0('UTRICLE_20211025/PseudotimeGeneCor_Extrastriolar.SC.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('UTRICLE_20211025/PseudotimeGO_POS_Extrastriolar.SC_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_20211025/PseudotimeGO_NEG_Extrastriolar.SC_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################


















this_pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$final.clst.sub.label=='Striolar.SC')])


    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

this_pbmc=subset(this_pbmc, cells=SELECT)

saveRDS(this_pbmc, file='UTRICLE_20211025/pbmc.Striolar.SC.rds')


this_pbmc=readRDS( file='UTRICLE_20211025/pbmc.Striolar.SC.rds')






pdf('UTRICLE_20211025/TRA_Striolar.SC.pdf',width=3.5,height=4)


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
write.table(sort(this_cor), paste0('UTRICLE_20211025/PseudotimeGeneCor_Striolar.SC.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




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
    write.table(this_out, paste0('UTRICLE_20211025/PseudotimeGO_POS_Striolar.SC_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_20211025/PseudotimeGO_NEG_Striolar.SC_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################

























