
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
library(clusterProfiler)

library(org.Mm.eg.db)

#####################################################################

pbmc=readRDS( 'UTRICLE_SEURAT_NEW.rds')
PCUSE=readRDS( 'UTRICLE_PCUSE.rds')




pbmc.Macrophage=subset(pbmc, cells=colnames(pbmc)[which(pbmc$new.clst.sub.label=='Macrophage')])



this_pbmc=pbmc.Macrophage
    
    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

pbmc.Macrophage=subset(this_pbmc, cells=SELECT)


saveRDS(pbmc.Macrophage, file='UTRICLE_Macrophage_TIME/pbmc.Macrophage.rds')


################################




NNN=10


pdf('UTRICLE_Macrophage_TIME/Macrophage_TRA.pdf',width=3.5,height=4)


#################################
this_pbmc=pbmc.Macrophage


     DimPlot(this_pbmc, group.by='batch',label=TRUE)+NoLegend()

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

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=F)

pbmc.Macrophage$pseudotime=OUT$P.PS
FeaturePlot(pbmc.Macrophage,features=c('pseudotime'))
#######################################################################

dev.off()


saveRDS(pbmc.Macrophage, file='UTRICLE_Macrophage_TIME/pbmc.Macrophage.rds')





this_pbmc=pbmc.Macrophage

this_cor=cor( t(as.matrix(this_pbmc@assays$RNA@data)),this_pbmc$pseudotime)
this_cor[which(is.na(this_cor))]=0
this_cor=this_cor[,1]
names(this_cor)=rownames(this_pbmc)
plot(sort(this_cor),type='p',pch=16)
write.table(sort(this_cor), paste0('UTRICLE_Macrophage_TIME/PseudotimeGeneCor_Macrophage.tsv'),sep='\t',quote=F,col.names=F,row.names=T)





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
    write.table(this_out, paste0('UTRICLE_Macrophage_TIME/PseudotimeGO_POS_Macrophage_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_Macrophage_TIME/PseudotimeGO_NEG_Macrophage_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################









################################################

all_pbmc=readRDS('COM_SEURAT.rds')
CLST=readRDS( 'COM_SUBCLST1.rds')



Idents(all_pbmc)=CLST
DimPlot(all_pbmc,label=TRUE)+NoLegend()

Idents(all_pbmc)=paste0(all_pbmc$batch,'.',CLST)


table(Idents(all_pbmc))





TMP='22'

MK.RESULT=FindMarkers(all_pbmc, ident.1=paste0('cochlea.3m','.',TMP),ident.2=paste0('cochlea.12m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('COCHLEA_22/3m.over.12m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('COCHLEA_22/3m.over.12m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)



MK.RESULT=FindMarkers(all_pbmc, ident.1=paste0('cochlea.12m','.',TMP),ident.2=paste0('cochlea.3m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('COCHLEA_22/12m.over.3m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('COCHLEA_22/12m.over.3m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


