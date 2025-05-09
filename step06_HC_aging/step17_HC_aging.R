#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')

library(Seurat)
library(dplyr)
library(patchwork)

UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')

##################################
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
#################################

.smooth<-function(x,df=10){
    y=smooth.spline(x,df=df)$y
    return(y)

}


#####################################
# Part4, HC aging
#####################################

##################
#UTRICLE
##################

pbmc=UTRICLE
DIFF=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new//DIFF_UTRICLE_ALLTIME.rds') 
MTGENE=rownames(pbmc)[grep('mt-',rownames(pbmc))]


#########################
TYPE=pbmc$new.clst.sub
############################
TIME=pbmc$batch
TIME[which(pbmc$batch=='cochlea.24mbu')]='cochlea.24m'
LABEL=paste0(TIME,'.',TYPE)

IDENT=LABEL
Idents(pbmc)=IDENT

ORDER=rep(1,length(TIME))
ORDER[which(TIME=='utricle.12m')]=2
ORDER[which(TIME=='utricle.22m')]=3
####################################################

TIME_LIST=unique(TIME)

################################################
# 需要做 UTRICLE： 5, 13, 1, 2 和 COCHLEA： 8, 11 
###################################################
this_clst='5'
################################################
DIFF_GENE=c()

     this_type=this_clst
     j1=1
     while(j1<=length(TIME_LIST)){
         j2=1
         while(j2<=length(TIME_LIST)){
              if(j1!=j2){
                  this_time1=TIME_LIST[j1]
                  this_time2=TIME_LIST[j2]
                 this_time_label=c(this_time1,this_time2)
 
                  this_time_tag=paste0('ident1.',this_time_label[1],'.vs.ident2.',this_time_label[2])
                  #############################
                 this_ident1=paste0(this_time_label[1],'.',this_type)
                 this_ident2=paste0(this_time_label[2],'.',this_type)
                  this_tag=paste0(this_time_tag,'.',this_type)
                  #############################
                  this_diff=DIFF[[this_tag]]
                  this_diff_gene=rownames(this_diff)[which( this_diff$p_val_adj<0.1 & this_diff$p_val <0.05 )]
                   DIFF_GENE=c(DIFF_GENE, this_diff_gene)
                  }
              j2=j2+1}
         j1=j1+1}

 DIFF_GENE=unique(DIFF_GENE)
DIFF_GENE=DIFF_GENE[which(!DIFF_GENE %in% MTGENE)]

 ##############
 MAT=pbmc[['RNA']]@scale.data
 #MAT=pbmc[['RNA']]@data
 #############

 USED_COL=which(TYPE %in%  this_clst )
 USED_ROW=which(rownames(MAT) %in% DIFF_GENE)

 USED_MAT=MAT[USED_ROW,USED_COL]

 USED_ORDER=ORDER[USED_COL]
 USED_TIME=TIME[USED_COL]

 USED_MAT=USED_MAT[,order(USED_ORDER)]
 USED_TIME=USED_TIME[order(USED_ORDER)]
 USED_ORDER=USED_ORDER[order(USED_ORDER)]

 USED_MAT=as.matrix(USED_MAT)
 COR=cor(t(USED_MAT),USED_ORDER)

USED_MAT=USED_MAT[which(!is.na(COR)),]
 COR=COR[which(!is.na(COR))]

USED_MAT=USED_MAT[order(COR),]

 SMAT=USED_MAT
 SMAT=t(apply(SMAT,1,.smooth,df=200))
 SMAT=t(apply(SMAT,1,scale))

 #SMAT=apply(SMAT,2,.smooth,df=100)
 #SMAT=t(apply(SMAT,1,scale))

 rownames(SMAT)=rownames(USED_MAT)
 colnames(SMAT)=colnames(USED_MAT)


 library('ComplexHeatmap')
 library('circlize')

 ha_bot = HeatmapAnnotation(
      TIME = USED_ORDER,
      col = list(TIME=c('1'='grey70','2'='grey55','3'='grey30') )
      )


 o.mat=SMAT



color_fun_3 =colorRamp2(c(-1, -0.2, 0, 0.2, 1 ), c('royalblue1','grey80','grey95','grey80','indianred1'))
ht=Heatmap(o.mat,row_title='',name="exp",
        cluster_columns=FALSE, cluster_rows=FALSE,
      show_column_dend = FALSE, show_row_dend = FALSE,
      show_column_names=FALSE, show_row_names=TRUE,
      col=color_fun_3, border = TRUE,
            top_annotation = ha_bot
       )
   
 print(ht)

    
ht_small=Heatmap(o.mat,row_title='',name="exp",
        cluster_columns=FALSE, cluster_rows=FALSE,
      show_column_dend = FALSE, show_row_dend = FALSE,
      show_column_names=FALSE, show_row_names=FALSE,
      col=color_fun_3, border = TRUE,
            top_annotation = ha_bot
       )



pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/H01_UTRICLE_C5_HEAT_BIG.pdf',width=60,height=100)
print(ht)
dev.off()


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/H01_UTRICLE_C5_HEAT_SMALL.pdf',width=6,height=10)
print(ht_small)
dev.off()

#####################################


























####################################################
# 绘制 Bubble Plot
####################################################

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



ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(COMBINE), 'ENTREZID', 'SYMBOL')




source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')





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


#Downregulated：
#GO:0007155~cell adhesion、 GO:0030010~establishment of cell polarity、 GO:0007409~axonogenesis、 GO:0097484~dendrite extension、 GO:0003341~cilium movement、 GO:0060088~auditory receptor cell stereocilium organization 
#Upregulated：
#GO:0008380~RNA splicing、 GO:0006457~protein folding、 GO:0006974~cellular response to DNA damage stimulus、 GO:0034063~stress granule assembly、 GO:2001233~regulation of apoptotic signaling pathway


GO_UP=c('GO:0008380','GO:0006457','GO:0006974','GO:0034063','GO:2001233')
GO_DW=c('GO:0007155','GO:0030010','GO:0007409','GO:0097484','GO:0003341','GO:0060088')



UTRICLE_DIFF=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new//DIFF_UTRICLE_ALLTIME.rds')
COCHLEA_DIFF=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new//DIFF_COCHLEA_ALLTIME.rds')




#############################################
# HC.type1.1,  3m v.s  12m
#############################################


#########################################
this.diff.pos=UTRICLE_DIFF[["ident1.utricle.12m.vs.ident2.utricle.3m.5"]]
this.diff.neg=UTRICLE_DIFF[["ident1.utricle.3m.vs.ident2.utricle.12m.5"]]

this.gene.pos=rownames(this.diff.pos)[which( this.diff.pos$p_val_adj<0.1 & this.diff.pos$p_val <0.05 )]
this.gene.neg=rownames(this.diff.neg)[which( this.diff.neg$p_val_adj<0.1 & this.diff.neg$p_val <0.05 )]


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
# HC.type1.1,  3m v.s  22m
#############################################


#########################################
this.diff.pos=UTRICLE_DIFF[["ident1.utricle.22m.vs.ident2.utricle.3m.5"]]
this.diff.neg=UTRICLE_DIFF[["ident1.utricle.3m.vs.ident2.utricle.22m.5"]]

this.gene.pos=rownames(this.diff.pos)[which( this.diff.pos$p_val_adj<0.1 & this.diff.pos$p_val <0.05 )]
this.gene.neg=rownames(this.diff.neg)[which( this.diff.neg$p_val_adj<0.1 & this.diff.neg$p_val <0.05 )]


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
# HC.type1.2,  3m v.s  12m
#############################################


#########################################
this.diff.pos=UTRICLE_DIFF[["ident1.utricle.12m.vs.ident2.utricle.3m.13"]]
this.diff.neg=UTRICLE_DIFF[["ident1.utricle.3m.vs.ident2.utricle.12m.13"]]

this.gene.pos=rownames(this.diff.pos)[which( this.diff.pos$p_val_adj<0.1 & this.diff.pos$p_val <0.05 )]
this.gene.neg=rownames(this.diff.neg)[which( this.diff.neg$p_val_adj<0.1 & this.diff.neg$p_val <0.05 )]


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
# HC.type1.2,  3m v.s  22m
#############################################


#########################################
this.diff.pos=UTRICLE_DIFF[["ident1.utricle.22m.vs.ident2.utricle.3m.13"]]
this.diff.neg=UTRICLE_DIFF[["ident1.utricle.3m.vs.ident2.utricle.22m.13"]]

this.gene.pos=rownames(this.diff.pos)[which( this.diff.pos$p_val_adj<0.1 & this.diff.pos$p_val <0.05 )]
this.gene.neg=rownames(this.diff.neg)[which( this.diff.neg$p_val_adj<0.1 & this.diff.neg$p_val <0.05 )]


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
# HC.type2.1,  3m v.s  12m
#############################################


#########################################
this.diff.pos=UTRICLE_DIFF[["ident1.utricle.12m.vs.ident2.utricle.3m.1"]]
this.diff.neg=UTRICLE_DIFF[["ident1.utricle.3m.vs.ident2.utricle.12m.1"]]

this.gene.pos=rownames(this.diff.pos)[which( this.diff.pos$p_val_adj<0.1 & this.diff.pos$p_val <0.05 )]
this.gene.neg=rownames(this.diff.neg)[which( this.diff.neg$p_val_adj<0.1 & this.diff.neg$p_val <0.05 )]


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
# HC.type2.1,  3m v.s  22m
#############################################


#########################################
this.diff.pos=UTRICLE_DIFF[["ident1.utricle.22m.vs.ident2.utricle.3m.1"]]
this.diff.neg=UTRICLE_DIFF[["ident1.utricle.3m.vs.ident2.utricle.22m.1"]]

this.gene.pos=rownames(this.diff.pos)[which( this.diff.pos$p_val_adj<0.1 & this.diff.pos$p_val <0.05 )]
this.gene.neg=rownames(this.diff.neg)[which( this.diff.neg$p_val_adj<0.1 & this.diff.neg$p_val <0.05 )]


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
# HC.type2.2,  3m v.s  12m
#############################################


#########################################
this.diff.pos=UTRICLE_DIFF[["ident1.utricle.12m.vs.ident2.utricle.3m.2"]]
this.diff.neg=UTRICLE_DIFF[["ident1.utricle.3m.vs.ident2.utricle.12m.2"]]

this.gene.pos=rownames(this.diff.pos)[which( this.diff.pos$p_val_adj<0.1 & this.diff.pos$p_val <0.05 )]
this.gene.neg=rownames(this.diff.neg)[which( this.diff.neg$p_val_adj<0.1 & this.diff.neg$p_val <0.05 )]


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
# HC.type2.2,  3m v.s  22m
#############################################


#########################################
this.diff.pos=UTRICLE_DIFF[["ident1.utricle.22m.vs.ident2.utricle.3m.2"]]
this.diff.neg=UTRICLE_DIFF[["ident1.utricle.3m.vs.ident2.utricle.22m.2"]]

this.gene.pos=rownames(this.diff.pos)[which( this.diff.pos$p_val_adj<0.1 & this.diff.pos$p_val <0.05 )]
this.gene.neg=rownames(this.diff.neg)[which( this.diff.neg$p_val_adj<0.1 & this.diff.neg$p_val <0.05 )]


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




###########################################



COM_U_P=COM_U[,which(colnames(COM_U)=='Pvalue')]
COM_U_P[which(is.na(COM_U_P))]=1

COM_U_C=COM_U[,which(colnames(COM_U)=='Count')]
COM_U_C[which(is.na(COM_U_C))]=0


COM_D_P=COM_D[,which(colnames(COM_D)=='Pvalue')]
COM_D_P[which(is.na(COM_D_P))]=1

COM_D_C=COM_D[,which(colnames(COM_D)=='Count')]
COM_D_C[which(is.na(COM_D_C))]=0



colnames(COM_U_P)=c('HC1.1_1','HC1.1_2','HC1.2_1','HC1.2_2','HC2.1_1','HC2.1_2','HC2.2_1','HC2.2_2')
colnames(COM_U_C)=c('HC1.1_1','HC1.1_2','HC1.2_1','HC1.2_2','HC2.1_1','HC2.1_2','HC2.2_1','HC2.2_2')
colnames(COM_D_P)=c('HC1.1_1','HC1.1_2','HC1.2_1','HC1.2_2','HC2.1_1','HC2.1_2','HC2.2_1','HC2.2_2')
colnames(COM_D_C)=c('HC1.1_1','HC1.1_2','HC1.2_1','HC1.2_2','HC2.1_1','HC2.1_2','HC2.2_1','HC2.2_2')





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










