
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


#########################################

Idents(pbmc)=LABEL






###############################################


this_ident1=c('cochlea.3m.8')
this_ident2=c('cochlea.12m.8')

this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff1$gene=rownames(this_diff1)


this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff2$gene=rownames(this_diff2)

write.table(this_diff1,file='./DIFF.Clst8.3m.over.12m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(this_diff2,file='./DIFF.Clst8.12.over.3m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

####
this_diff=this_diff1

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst8.3m.over.12m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)
####
this_diff=this_diff2

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst8.12m.over.3m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)














###############################################

this_ident1=c('cochlea.12m.8')
this_ident2=c('cochlea.24m.8')

this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff1$gene=rownames(this_diff1)


this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff2$gene=rownames(this_diff2)

write.table(this_diff1,file='./DIFF.Clst8.12m.over.24m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(this_diff2,file='./DIFF.Clst8.24.over.12m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

####
this_diff=this_diff1

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst8.12m.over.24m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)
####
this_diff=this_diff2

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst8.24m.over.12m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)












################################################

this_ident1=c('cochlea.3m.8')
this_ident2=c('cochlea.24m.8')

this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff1$gene=rownames(this_diff1)


this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff2$gene=rownames(this_diff2)

write.table(this_diff1,file='./DIFF.Clst8.3m.over.24m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(this_diff2,file='./DIFF.Clst8.24.over.3m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)


####
this_diff=this_diff1

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst8.3m.over.24m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)
####
this_diff=this_diff2

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst8.24m.over.3m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)












###############################################


this_ident1=c('cochlea.3m.11')
this_ident2=c('cochlea.12m.11')

this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff1$gene=rownames(this_diff1)


this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff2$gene=rownames(this_diff2)

write.table(this_diff1,file='./DIFF.Clst11.3m.over.12m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(this_diff2,file='./DIFF.Clst11.12.over.3m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

####
this_diff=this_diff1

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst11.3m.over.12m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)
####
this_diff=this_diff2

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst11.12m.over.3m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)














###############################################

this_ident1=c('cochlea.12m.11')
this_ident2=c('cochlea.24m.11')

this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff1$gene=rownames(this_diff1)


this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff2$gene=rownames(this_diff2)

write.table(this_diff1,file='./DIFF.Clst11.12m.over.24m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(this_diff2,file='./DIFF.Clst11.24.over.12m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

####
this_diff=this_diff1

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst11.12m.over.24m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)
####
this_diff=this_diff2

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst11.24m.over.12m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)












################################################

this_ident1=c('cochlea.3m.11')
this_ident2=c('cochlea.24m.11')

this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff1$gene=rownames(this_diff1)


this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
this_diff2$gene=rownames(this_diff2)

write.table(this_diff1,file='./DIFF.Clst11.3m.over.24m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(this_diff2,file='./DIFF.Clst11.24.over.3m.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)


####
this_diff=this_diff1

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst11.3m.over.24m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)
####
this_diff=this_diff2

    this_sig_gene=rownames(this_diff)[which(this_diff$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./GO_DIFF.Clst11.24m.over.3m.tsv'),sep='\t',quote=F,col.names=T,row.names=F)






