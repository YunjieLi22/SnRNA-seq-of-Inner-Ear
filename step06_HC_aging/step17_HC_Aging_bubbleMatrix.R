


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




.simple_combine_NA <- function(exp_sc_mat1, exp_sc_mat2, FILL=FALSE){
    FILL=FILL
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    ##############################################
    if(FILL==TRUE){
        gene1=rownames(exp_sc_mat)
        gene2=rownames(exp_ref_mat)
        gene12=gene2[which(!gene2 %in% gene1)]
        gene21=gene1[which(!gene1 %in% gene2)]
        exp_sc_mat_add=matrix(NA,ncol=ncol(exp_sc_mat),nrow=length(gene12))
        rownames(exp_sc_mat_add)=gene12
        colnames(exp_sc_mat_add)=colnames(exp_sc_mat)
        exp_ref_mat_add=matrix(NA,ncol=ncol(exp_ref_mat),nrow=length(gene21))
        rownames(exp_ref_mat_add)=gene21
        colnames(exp_ref_mat_add)=colnames(exp_ref_mat)
        exp_sc_mat=rbind(exp_sc_mat, exp_sc_mat_add)
        exp_ref_mat=rbind(exp_ref_mat, exp_ref_mat_add)
    }
    ############################################
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }






###############
POS.CUT=0.1
NEG.CUT= -0.1
########################################

GO_UP=as.character(read.table('./GO_HC_UP.txt')[,1])
GO_DW=as.character(read.table('./GO_HC_DW.txt')[,1])









#############################################
# HC.type1.1
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/pbmc.HC.type1.1.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type1.1.tsv',row.names=1,header=F,sep='\t')


ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(this.pbmc), 'ENTREZID', 'SYMBOL')


#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]

#########################################
this_sig_gene=this.gene.pos
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_up=this_out[which(this_out$ID %in% GO_UP),]
this_up_info=cbind(this_up$pvalue,this_up$Count)
rownames(this_up_info)=this_up$ID
colnames(this_up_info)=c('Pvalue','Count')
########################################

#########################################
this_sig_gene=this.gene.neg
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_dw=this_out[which(this_out$ID %in% GO_DW),]
this_dw_info=cbind(this_dw$pvalue,this_dw$Count)
rownames(this_dw_info)=paste0(this_dw$ID)
colnames(this_dw_info)=c('Pvalue','Count')
########################################


COM_U=this_up_info
COM_D=this_dw_info







#############################################
# HC.type1.2
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/pbmc.HC.type1.2.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type1.2.tsv',row.names=1,header=F,sep='\t')



ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(this.pbmc), 'ENTREZID', 'SYMBOL')


#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]

#########################################
this_sig_gene=this.gene.pos
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_up=this_out[which(this_out$ID %in% GO_UP),]
this_up_info=cbind(this_up$pvalue,this_up$Count)
rownames(this_up_info)=this_up$ID
colnames(this_up_info)=c('Pvalue','Count')
########################################

#########################################
this_sig_gene=this.gene.neg
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_dw=this_out[which(this_out$ID %in% GO_DW),]
this_dw_info=cbind(this_dw$pvalue,this_dw$Count)
rownames(this_dw_info)=paste0(this_dw$ID)
colnames(this_dw_info)=c('Pvalue','Count')
########################################

COM_U=.simple_combine_NA(COM_U,this_up_info,FILL=TRUE)$combine
COM_D=.simple_combine_NA(COM_D,this_dw_info,FILL=TRUE)$combine






#############################################
# HC.type2.1
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/pbmc.HC.type2.1.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type2.1.tsv',row.names=1,header=F,sep='\t')



ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(this.pbmc), 'ENTREZID', 'SYMBOL')


#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]

#########################################
this_sig_gene=this.gene.pos
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_up=this_out[which(this_out$ID %in% GO_UP),]
this_up_info=cbind(this_up$pvalue,this_up$Count)
rownames(this_up_info)=this_up$ID
colnames(this_up_info)=c('Pvalue','Count')
########################################

#########################################
this_sig_gene=this.gene.neg
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_dw=this_out[which(this_out$ID %in% GO_DW),]
this_dw_info=cbind(this_dw$pvalue,this_dw$Count)
rownames(this_dw_info)=paste0(this_dw$ID)
colnames(this_dw_info)=c('Pvalue','Count')
########################################


COM_U=.simple_combine_NA(COM_U,this_up_info,FILL=TRUE)$combine
COM_D=.simple_combine_NA(COM_D,this_dw_info,FILL=TRUE)$combine









#############################################
# HC.type2.2
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/pbmc.HC.type2.2.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type2.2.tsv',row.names=1,header=F,sep='\t')



ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(this.pbmc), 'ENTREZID', 'SYMBOL')


#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]

#########################################
this_sig_gene=this.gene.pos
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_up=this_out[which(this_out$ID %in% GO_UP),]
this_up_info=cbind(this_up$pvalue,this_up$Count)
rownames(this_up_info)=this_up$ID
colnames(this_up_info)=c('Pvalue','Count')
########################################

#########################################
this_sig_gene=this.gene.neg
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_dw=this_out[which(this_out$ID %in% GO_DW),]
this_dw_info=cbind(this_dw$pvalue,this_dw$Count)
rownames(this_dw_info)=paste0(this_dw$ID)
colnames(this_dw_info)=c('Pvalue','Count')
########################################


COM_U=.simple_combine_NA(COM_U,this_up_info,FILL=TRUE)$combine
COM_D=.simple_combine_NA(COM_D,this_dw_info,FILL=TRUE)$combine





#############################################
# IHC
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/pbmc.8.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/PseudotimeGeneCor_8.tsv',row.names=1,header=F,sep='\t')


ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(this.pbmc), 'ENTREZID', 'SYMBOL')


#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]

#########################################
this_sig_gene=this.gene.pos
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_up=this_out[which(this_out$ID %in% GO_UP),]
this_up_info=cbind(this_up$pvalue,this_up$Count)
rownames(this_up_info)=this_up$ID
colnames(this_up_info)=c('Pvalue','Count')
########################################

#########################################
this_sig_gene=this.gene.neg
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_dw=this_out[which(this_out$ID %in% GO_DW),]
this_dw_info=cbind(this_dw$pvalue,this_dw$Count)
rownames(this_dw_info)=paste0(this_dw$ID)
colnames(this_dw_info)=c('Pvalue','Count')
########################################


COM_U=.simple_combine_NA(COM_U,this_up_info,FILL=TRUE)$combine
COM_D=.simple_combine_NA(COM_D,this_dw_info,FILL=TRUE)$combine






#############################################
# OHC
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/pbmc.11.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/PseudotimeGeneCor_11.tsv',row.names=1,header=F,sep='\t')


ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(this.pbmc), 'ENTREZID', 'SYMBOL')


#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]

#########################################
this_sig_gene=this.gene.pos
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_up=this_out[which(this_out$ID %in% GO_UP),]
this_up_info=cbind(this_up$pvalue,this_up$Count)
rownames(this_up_info)=this_up$ID
colnames(this_up_info)=c('Pvalue','Count')
########################################

#########################################
this_sig_gene=this.gene.neg
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_dw=this_out[which(this_out$ID %in% GO_DW),]
this_dw_info=cbind(this_dw$pvalue,this_dw$Count)
rownames(this_dw_info)=paste0(this_dw$ID)
colnames(this_dw_info)=c('Pvalue','Count')
########################################

COM_U=.simple_combine_NA(COM_U,this_up_info,FILL=TRUE)$combine
COM_D=.simple_combine_NA(COM_D,this_dw_info,FILL=TRUE)$combine
#############################################

saveRDS(COM_U,file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/GO_UP_COM_HC.rds')
saveRDS(COM_D,file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/GO_DW_COM_HC.rds')



COM_U_P=COM_U[,which(colnames(COM_U)=='Pvalue')]
COM_U_P[which(is.na(COM_U_P))]=1

COM_U_C=COM_U[,which(colnames(COM_U)=='Count')]
COM_U_C[which(is.na(COM_U_C))]=0


COM_D_P=COM_D[,which(colnames(COM_D)=='Pvalue')]
COM_D_P[which(is.na(COM_D_P))]=1

COM_D_C=COM_D[,which(colnames(COM_D)=='Count')]
COM_D_C[which(is.na(COM_D_C))]=0



colnames(COM_U_P)=c('HC1.1','HC1.2','HC2.1','HC2.2','IHC','OHC')
colnames(COM_U_C)=c('HC1.1','HC1.2','HC2.1','HC2.2','IHC','OHC')
colnames(COM_D_P)=c('HC1.1','HC1.2','HC2.1','HC2.2','IHC','OHC')
colnames(COM_D_C)=c('HC1.1','HC1.2','HC2.1','HC2.2','IHC','OHC')


library(reshape2)
library(ggplot2) 

data=COM_U_C
pdata=COM_U_P
data_melt<-melt (data)
names(data_melt) = c('GO', 'TYPE', 'COUNT')

data_melt$PV=rep(1,length(data_melt$GO))
i=1
while(i<=length(data_melt$PV)){
    data_melt$PV[i]=pdata[which(rownames(data)==data_melt$GO[i]),which(colnames(data)==data_melt$TYPE[i])]
    i=i+1 
    }
data_melt$LgPV=-log(data_melt$PV,10)

p1<-ggplot(data_melt, aes(x = TYPE, y = GO, size = COUNT, color=LgPV)) + geom_point() + 
scale_colour_gradient2(
  high = "red",
  mid = "grey80",
  low = "grey80",
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "colour"
)+theme(panel.background = element_rect(fill = 'white', colour = 'black'))
print(p1)


#########################


library(reshape2)
library(ggplot2) 

data=COM_D_C
pdata=COM_D_P
data_melt<-melt (data)
names(data_melt) = c('GO', 'TYPE', 'COUNT')

data_melt$PV=rep(1,length(data_melt$GO))
i=1
while(i<=length(data_melt$PV)){
    data_melt$PV[i]=pdata[which(rownames(data)==data_melt$GO[i]),which(colnames(data)==data_melt$TYPE[i])]
    i=i+1 
    }
data_melt$LgPV=-log(data_melt$PV,10)

p2<-ggplot(data_melt, aes(x = TYPE, y = GO, size = COUNT, color=LgPV)) + geom_point() + 
scale_colour_gradient2(
  high = "blue",
  mid = "grey80",
  low = "grey80",
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "colour"
)+theme(panel.background = element_rect(fill = 'white', colour = 'black'))
print(p2)



pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/BubbleMatrix_HC_UP.pdf',width=5,height=5)
print(p1)
dev.off()

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/BubbleMatrix_HC_DW.pdf',width=5,height=5)
print(p2)
dev.off()



################################################














###############
POS.CUT=0.1
NEG.CUT= -0.1
########################################

GO_UP=as.character(read.table('./GO_Mac_UP.txt')[,1])
GO_DW=as.character(read.table('./GO_Mac_DW.txt')[,1])









#############################################
# UTRICLE_Macrophage
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_Macrophage_TIME/pbmc.Macrophage.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_Macrophage_TIME/PseudotimeGeneCor_Macrophage.tsv',row.names=1,header=F,sep='\t')


ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(this.pbmc), 'ENTREZID', 'SYMBOL')


#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]

#########################################
this_sig_gene=this.gene.pos
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_up=this_out[which(this_out$ID %in% GO_UP),]
this_up_info=cbind(this_up$pvalue,this_up$Count)
rownames(this_up_info)=this_up$ID
colnames(this_up_info)=c('Pvalue','Count')
########################################

#########################################
this_sig_gene=this.gene.neg
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_dw=this_out[which(this_out$ID %in% GO_DW),]
this_dw_info=cbind(this_dw$pvalue,this_dw$Count)
rownames(this_dw_info)=paste0(this_dw$ID)
colnames(this_dw_info)=c('Pvalue','Count')
########################################


COM_U=this_up_info
COM_D=this_dw_info







#############################################
# Cochlea_Macrophage
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/pbmc.14.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/PseudotimeGeneCor_14.tsv',row.names=1,header=F,sep='\t')


ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(this.pbmc), 'ENTREZID', 'SYMBOL')


#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]

#########################################
this_sig_gene=this.gene.pos
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_up=this_out[which(this_out$ID %in% GO_UP),]
this_up_info=cbind(this_up$pvalue,this_up$Count)
rownames(this_up_info)=this_up$ID
colnames(this_up_info)=c('Pvalue','Count')
########################################

#########################################
this_sig_gene=this.gene.neg
this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)
this_out=this_ego@result

this_dw=this_out[which(this_out$ID %in% GO_DW),]
this_dw_info=cbind(this_dw$pvalue,this_dw$Count)
rownames(this_dw_info)=paste0(this_dw$ID)
colnames(this_dw_info)=c('Pvalue','Count')
########################################

COM_U=.simple_combine_NA(COM_U,this_up_info,FILL=TRUE)$combine
COM_D=.simple_combine_NA(COM_D,this_dw_info,FILL=TRUE)$combine




#############################################

saveRDS(COM_U,file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/GO_UP_COM_Mac.rds')
saveRDS(COM_D,file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/GO_DW_COM_Mac.rds')



COM_U_P=COM_U[,which(colnames(COM_U)=='Pvalue')]
COM_U_P[which(is.na(COM_U_P))]=1

COM_U_C=COM_U[,which(colnames(COM_U)=='Count')]
COM_U_C[which(is.na(COM_U_C))]=0


COM_D_P=COM_D[,which(colnames(COM_D)=='Pvalue')]
COM_D_P[which(is.na(COM_D_P))]=1

COM_D_C=COM_D[,which(colnames(COM_D)=='Count')]
COM_D_C[which(is.na(COM_D_C))]=0



colnames(COM_U_P)=c('Utri.Mac','Coch.Mac')
colnames(COM_U_C)=c('Utri.Mac','Coch.Mac')
colnames(COM_D_P)=c('Utri.Mac','Coch.Mac')
colnames(COM_D_C)=c('Utri.Mac','Coch.Mac')


library(reshape2)
library(ggplot2) 

data=COM_U_C
pdata=COM_U_P
data_melt<-melt (data)
names(data_melt) = c('GO', 'TYPE', 'COUNT')

data_melt$PV=rep(1,length(data_melt$GO))
i=1
while(i<=length(data_melt$PV)){
    data_melt$PV[i]=pdata[which(rownames(data)==data_melt$GO[i]),which(colnames(data)==data_melt$TYPE[i])]
    i=i+1 
    }
data_melt$LgPV=-log(data_melt$PV,10)

p1<-ggplot(data_melt, aes(x = TYPE, y = GO, size = COUNT, color=LgPV)) + geom_point() + 
scale_colour_gradient2(
  high = "red",
  mid = "grey80",
  low = "grey80",
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "colour"
)+theme(panel.background = element_rect(fill = 'white', colour = 'black'))
print(p1)


#########################


library(reshape2)
library(ggplot2) 

data=COM_D_C
pdata=COM_D_P
data_melt<-melt (data)
names(data_melt) = c('GO', 'TYPE', 'COUNT')

data_melt$PV=rep(1,length(data_melt$GO))
i=1
while(i<=length(data_melt$PV)){
    data_melt$PV[i]=pdata[which(rownames(data)==data_melt$GO[i]),which(colnames(data)==data_melt$TYPE[i])]
    i=i+1 
    }
data_melt$LgPV=-log(data_melt$PV,10)

p2<-ggplot(data_melt, aes(x = TYPE, y = GO, size = COUNT, color=LgPV)) + geom_point() + 
scale_colour_gradient2(
  high = "blue",
  mid = "grey80",
  low = "grey80",
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "colour"
)+theme(panel.background = element_rect(fill = 'white', colour = 'black'))
print(p2)



pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/BubbleMatrix_Mac_UP.pdf',width=4,height=5)
print(p1)
dev.off()

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/BubbleMatrix_Mac_DW.pdf',width=4,height=5)
print(p2)
dev.off()













#############################################
# UTRICLE_Macrophage
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_Macrophage_TIME/pbmc.Macrophage.rds')


this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$batch=='utricle.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$batch=='utricle.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$batch=='utricle.22m')]=3 #24m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=F)

this.pbmc$pseudotime=OUT$P.PS



#this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))

this.exp=as.matrix(this.pbmc[['RNA']]@data)



pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/Mki67_Top2a_Mac_Utricle.pdf',width=5,height=3)


OOO=order(this.pbmc$pseudotime)
show.gene='Mki67'
plot(this.exp[which(rownames(this.exp)==show.gene),][OOO],ylab=show.gene,xlab='pseudotime',pch=16,type='p')


OOO=order(this.pbmc$pseudotime)
show.gene='Top2a'
plot(this.exp[which(rownames(this.exp)==show.gene),][OOO],ylab=show.gene,xlab='pseudotime',pch=16,type='p')

Idents(this.pbmc)=factor(this.pbmc$batch,levels=c('utricle.3m','utricle.12m','utricle.22m'))

VlnPlot(this.pbmc, features=c('Mki67','Top2a'),pt.size=2)

dev.off()




#############################################
# Cochlea_Macrophage
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/pbmc.14.rds')


this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$batch=='cochlea.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$batch=='cochlea.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$batch %in% c('cochlea.24m','cochlea.24mbu'))]=3 #24m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=F)

this.pbmc$pseudotime=OUT$P.PS



#this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))

this.exp=as.matrix(this.pbmc[['RNA']]@data)



pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/Mki67_Top2a_Mac_Cochlea.pdf',width=5,height=3)


OOO=order(this.pbmc$pseudotime)
show.gene='Mki67'
plot(this.exp[which(rownames(this.exp)==show.gene),][OOO],ylab=show.gene,xlab='pseudotime',pch=16,type='p')


OOO=order(this.pbmc$pseudotime)
show.gene='Top2a'
plot(this.exp[which(rownames(this.exp)==show.gene),][OOO],ylab=show.gene,xlab='pseudotime',pch=16,type='p')

TMP=this.pbmc$batch
TMP[which(TMP=='cochlea.24mbu')]='cochlea.24m'

Idents(this.pbmc)=factor(TMP,levels=c('cochlea.3m','cochlea.12m','cochlea.24m'))

VlnPlot(this.pbmc, features=c('Mki67','Top2a'),pt.size=2)

dev.off()


