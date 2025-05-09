
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

#####################################################################
setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging')

.smooth<-function(x,df=10){
    x=x
    df=df
    y=smooth.spline(x=x,df=df)$y
    return(y)
    }

###############
POS.CUT=0.1
NEG.CUT= -0.1
########################################


RANDOM_NUMBER=100
FONT_SIZE=10

#############################################
# HC.type1.1
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/pbmc.HC.type1.1.rds')
print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type1.1.tsv',row.names=1,header=F,sep='\t')

#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]


this.exp=as.matrix(this.pbmc[['RNA']]@scale.data)

this.exp.pos=this.exp[which(rownames(this.exp)%in% this.gene.pos),order(this.pbmc$pseudotime)] 
this.exp.pos.smooth=t(apply(this.exp.pos,1,.smooth,10))

this.exp.neg=this.exp[which(rownames(this.exp)%in% this.gene.neg),order(this.pbmc$pseudotime)] 
this.exp.neg.smooth=t(apply(this.exp.neg,1,.smooth,10))

##############
set.seed(123)
show.index=sort(sample(1:ncol(this.pbmc),RANDOM_NUMBER))
####################
o.mat=this.exp.pos.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht1=Heatmap(o.mat,name='HC1.1pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht2=Heatmap(o.mat,name='HC1.1neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)




#############################################
# HC.type1.2
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/pbmc.HC.type1.2.rds')
print(dim(this.pbmc))

this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type1.2.tsv',row.names=1,header=F,sep='\t')

#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]
this.exp=as.matrix(this.pbmc[['RNA']]@scale.data)

this.exp.pos=this.exp[which(rownames(this.exp)%in% this.gene.pos),order(this.pbmc$pseudotime)] 
this.exp.pos.smooth=t(apply(this.exp.pos,1,.smooth,10))

this.exp.neg=this.exp[which(rownames(this.exp)%in% this.gene.neg),order(this.pbmc$pseudotime)] 
this.exp.neg.smooth=t(apply(this.exp.neg,1,.smooth,10))



##############
set.seed(123)
show.index=sort(sample(1:ncol(this.pbmc),RANDOM_NUMBER))
####################

o.mat=this.exp.pos.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht3=Heatmap(o.mat,name='HC1.2pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht4=Heatmap(o.mat,name='HC1.2neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)




#############################################
# HC.type2.1
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/pbmc.HC.type2.1.rds')
print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type2.1.tsv',row.names=1,header=F,sep='\t')

#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]
this.exp=as.matrix(this.pbmc[['RNA']]@scale.data)

this.exp.pos=this.exp[which(rownames(this.exp)%in% this.gene.pos),order(this.pbmc$pseudotime)] 
this.exp.pos.smooth=t(apply(this.exp.pos,1,.smooth,10))

this.exp.neg=this.exp[which(rownames(this.exp)%in% this.gene.neg),order(this.pbmc$pseudotime)] 
this.exp.neg.smooth=t(apply(this.exp.neg,1,.smooth,10))


##############
set.seed(123)
show.index=sort(sample(1:ncol(this.pbmc),RANDOM_NUMBER))
####################

o.mat=this.exp.pos.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht5=Heatmap(o.mat,name='HC2.1pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht6=Heatmap(o.mat,name='HC2.1neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)




#############################################
# HC.type2.2
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/pbmc.HC.type2.2.rds')
print(dim(this.pbmc))




this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type2.2.tsv',row.names=1,header=F,sep='\t')

#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]
this.exp=as.matrix(this.pbmc[['RNA']]@scale.data)

this.exp.pos=this.exp[which(rownames(this.exp)%in% this.gene.pos),order(this.pbmc$pseudotime)] 
this.exp.pos.smooth=t(apply(this.exp.pos,1,.smooth,10))

this.exp.neg=this.exp[which(rownames(this.exp)%in% this.gene.neg),order(this.pbmc$pseudotime)] 
this.exp.neg.smooth=t(apply(this.exp.neg,1,.smooth,10))


##############
set.seed(123)
show.index=sort(sample(1:ncol(this.pbmc),RANDOM_NUMBER))
####################

o.mat=this.exp.pos.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht7=Heatmap(o.mat,name='HC2.2pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht8=Heatmap(o.mat,name='HC2.2neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)






#############################################
# IHC
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/pbmc.8.rds')
print(dim(this.pbmc))

this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/PseudotimeGeneCor_8.tsv',row.names=1,header=F,sep='\t')

################################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$time=='cochlea.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$time=='cochlea.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$time=='cochlea.24m')]=3 #24m
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


#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]
this.exp=as.matrix(this.pbmc[['RNA']]@scale.data)

this.exp.pos=this.exp[which(rownames(this.exp)%in% this.gene.pos),order(this.pbmc$pseudotime)] 
this.exp.pos.smooth=t(apply(this.exp.pos,1,.smooth,10))

this.exp.neg=this.exp[which(rownames(this.exp)%in% this.gene.neg),order(this.pbmc$pseudotime)] 
this.exp.neg.smooth=t(apply(this.exp.neg,1,.smooth,10))


##############
set.seed(123)
show.index=sort(sample(1:ncol(this.pbmc),RANDOM_NUMBER))
####################

o.mat=this.exp.pos.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht9=Heatmap(o.mat,name='IHCpos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht10=Heatmap(o.mat,name='IHCneg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)







#############################################
# OHC
#############################################

this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/pbmc.11.rds')
print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111/PseudotimeGeneCor_11.tsv',row.names=1,header=F,sep='\t')
###############################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$time=='cochlea.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$time=='cochlea.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$time=='cochlea.24m')]=3 #24m
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


########################################

#########################################
this.gene.pos=rownames(this.gene.cor)[which(this.gene.cor>POS.CUT)]
this.gene.neg=rownames(this.gene.cor)[which(this.gene.cor<NEG.CUT)]
this.exp=as.matrix(this.pbmc[['RNA']]@scale.data)

this.exp.pos=this.exp[which(rownames(this.exp)%in% this.gene.pos),order(this.pbmc$pseudotime)] 
this.exp.pos.smooth=t(apply(this.exp.pos,1,.smooth,10))

this.exp.neg=this.exp[which(rownames(this.exp)%in% this.gene.neg),order(this.pbmc$pseudotime)] 
this.exp.neg.smooth=t(apply(this.exp.neg,1,.smooth,10))

##############
set.seed(123)
show.index=sort(sample(1:ncol(this.pbmc),RANDOM_NUMBER))
####################


o.mat=this.exp.pos.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht11=Heatmap(o.mat,name='OHCpos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht12=Heatmap(o.mat,name='OHCneg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)























pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/HEAT_POS.pdf',width=4,height=8)
print(ht1 %v% ht3 %v% ht5 %v% ht7 %v% ht9 %v% ht11)  
dev.off()


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/HEAT_NEG.pdf',width=4,height=8)
print(ht2 %v% ht4 %v% ht6 %v% ht8 %v% ht10 %v% ht12)  
dev.off()







