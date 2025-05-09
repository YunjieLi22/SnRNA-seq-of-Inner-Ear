

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
# Fibrocyte.1
#############################################



this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Utri.Fibrocyte.1.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/PseudotimeGeneCor_Utri_Fibrocyte.1.tsv',row.names=1,header=F,sep='\t')

print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))





################################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=F, COL=OUT$COL, SHOW.SUMMIT=F)

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
ht1=Heatmap(o.mat,name='F1pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht2=Heatmap(o.mat,name='F1neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)








#############################################
# Fibrocyte.2
#############################################



this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Utri.Fibrocyte.2.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/PseudotimeGeneCor_Utri_Fibrocyte.2.tsv',row.names=1,header=F,sep='\t')

print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))





################################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=F, COL=OUT$COL, SHOW.SUMMIT=F)

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
ht3=Heatmap(o.mat,name='F2pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht4=Heatmap(o.mat,name='F2neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)





#############################################
# Fibrocyte.3
#############################################



this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Utri.Fibrocyte.3.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/PseudotimeGeneCor_Utri_Fibrocyte.3.tsv',row.names=1,header=F,sep='\t')

print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))





################################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=F, COL=OUT$COL, SHOW.SUMMIT=F)

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
ht5=Heatmap(o.mat,name='F3pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht6=Heatmap(o.mat,name='F3neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)








#############################################
# Stroma.cell
#############################################



this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Utri.Stroma.cell.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/PseudotimeGeneCor_Utri_Stroma.cell.tsv',row.names=1,header=F,sep='\t')

print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))





################################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=F, COL=OUT$COL, SHOW.SUMMIT=F)

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
ht7=Heatmap(o.mat,name='Spos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht8=Heatmap(o.mat,name='Sneg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)



pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/HEAT_POS_Utri.pdf',width=4,height=6)
print(ht1 %v% ht3 %v% ht5 %v% ht7 )  
dev.off()


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/HEAT_NEG_Utri.pdf',width=4,height=6)
print(ht2 %v% ht4 %v% ht6 %v% ht8 )  
dev.off()



#########################################################################################################






#############################################
# 0
#############################################



this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Coch.0.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/PseudotimeGeneCor_Coch_0.tsv',row.names=1,header=F,sep='\t')

print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))





################################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$orig.ident=='cochlea.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$orig.ident=='cochlea.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$orig.ident %in% c('cochlea.24m','cochlea.24mbu'))]=3 #22m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=F, COL=OUT$COL, SHOW.SUMMIT=F)

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
ht1=Heatmap(o.mat,name='0pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht2=Heatmap(o.mat,name='0neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


#############################################
# 2
#############################################



this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Coch.2.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/PseudotimeGeneCor_Coch_2.tsv',row.names=1,header=F,sep='\t')

print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))





################################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$orig.ident=='cochlea.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$orig.ident=='cochlea.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$orig.ident %in% c('cochlea.24m','cochlea.24mbu'))]=3 #22m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=F, COL=OUT$COL, SHOW.SUMMIT=F)

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
ht3=Heatmap(o.mat,name='2pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht4=Heatmap(o.mat,name='2neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)



#############################################
# 9.1
#############################################



this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Coch.9.1.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/PseudotimeGeneCor_Coch_9.1.tsv',row.names=1,header=F,sep='\t')

print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))





################################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$orig.ident=='cochlea.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$orig.ident=='cochlea.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$orig.ident %in% c('cochlea.24m','cochlea.24mbu'))]=3 #22m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=F, COL=OUT$COL, SHOW.SUMMIT=F)

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
ht5=Heatmap(o.mat,name='9.1pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht6=Heatmap(o.mat,name='9.1neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)



#############################################
# 9.2
#############################################



this.pbmc=readRDS(file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/pbmc.Coch.9.2.rds')
this.gene.cor=read.table('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/PseudotimeGeneCor_Coch_9.2.tsv',row.names=1,header=F,sep='\t')

print(dim(this.pbmc))


this.pbmc <- ScaleData(object = this.pbmc, features = rownames(this.pbmc[['RNA']]@data))





################################################

this.pbmc$timev=rep(NA, length(ncol(this.pbmc)))
this.pbmc$timev[which(this.pbmc$orig.ident=='cochlea.3m')]=1 #3m
this.pbmc$timev[which(this.pbmc$orig.ident=='cochlea.12m')]=2 #12m
this.pbmc$timev[which(this.pbmc$orig.ident %in% c('cochlea.24m','cochlea.24mbu'))]=3 #22m
this_vec=this.pbmc@reductions$umap@cell.embeddings

this_fit=lm(this.pbmc$timev~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

NNN=10

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=F)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=F, COL=OUT$COL, SHOW.SUMMIT=F)

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
ht7=Heatmap(o.mat,name='9.2pos',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)


o.mat=this.exp.neg.smooth[,show.index]
col_fun =colorRamp2(c(-1,0,1,2 ), c('royalblue3','grey95','red1','gold1'))
ht8=Heatmap(o.mat,name='9.2neg',cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,row_title=paste0(as.character(nrow(o.mat))),row_title_gp = gpar(fontsize = FONT_SIZE),
    row_names_side='left'
	)





pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/HEAT_POS_Coch.pdf',width=4,height=6)
print(ht1 %v% ht3 %v% ht5 %v% ht7 )  
dev.off()


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/HEAT_NEG_Coch.pdf',width=4,height=6)
print(ht2 %v% ht4 %v% ht6 %v% ht8 )  
dev.off()

