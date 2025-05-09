
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

pbmc=readRDS( 'UTRICLE_SEURAT_NEW.rds')
PCUSE=readRDS( 'UTRICLE_PCUSE.rds')


#DimPlot(pbmc, reduction='umap', group.by='new.clst.sub',pt.size=0.1,label=T,label.size=8)  +NoLegend()

pbmc$final.clst.sub.label=pbmc$new.clst.sub.label

pbmc$final.clst.sub.label[pbmc$new.clst.sub=='0']='Extrastriolar.SC'
pbmc$final.clst.sub.label[pbmc$new.clst.sub=='11']='Striolar.SC'   
pbmc$final.clst.sub.label[pbmc$new.clst.sub=='12']='Stroma.cell'

DimPlot(pbmc, reduction='umap', group.by='final.clst.sub.label',pt.size=0.1,label=T,label.size=8)  +NoLegend()



pdf('UTRICLE_20211025/UTRICLE_UMAP_LABEL_FINAL.pdf',width=12,height=12)
DimPlot(pbmc, reduction='umap', group.by='final.clst.sub.label', 
         pt.size=0.1,label=T,label.size=5) +NoLegend()
dev.off()


Idents(pbmc)=pbmc$final.clst.sub.label


#saveRDS(pbmc, file='UTRICLE_20211025/UTRICLE_SEURAT_FINAL.rds')

#######################################################



colnames(pbmc@meta.data)

table(pbmc$batch)
table(pbmc$final.clst.sub.label)

TIME=pbmc$batch
TYPE=pbmc$final.clst.sub.label


pbmc$final.clst.sub.label.time=paste0(TIME,'.',TYPE)
Idents(pbmc)=pbmc$final.clst.sub.label.time


write.table(table(pbmc$final.clst.sub.label.time),row.names=T,col.names=F,file='UTRICLE_HC_TIME/TIME_CELL_NUMBER.txt',sep='\t',quote=F)

#saveRDS(pbmc, file='UTRICLE_20211025/UTRICLE_SEURAT_FINAL.rds')


###################################################
IDENT=pbmc$final.clst.sub.label.time
TYPE_LIST=unique(TYPE)
TIME_LIST=unique(TIME)

DIFF=list()
i=1 #TYPE
while(i<=length(TYPE_LIST)){
    this_type=TYPE_LIST[i]
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
                 if(length(which(IDENT == this_ident1))>=3 & length(which(IDENT == this_ident2))>=3){
                     this_diff=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
                     DIFF[[this_tag]]=this_diff
                     }
                 }
             j2=j2+1}
        j1=j1+1}
    print(i)
    i=i+1
    }

#saveRDS(DIFF,file='UTRICLE_20211025/UTRICLE_DIFF_ALLTIME.rds' )

DIFF=readRDS('UTRICLE_20211025/UTRICLE_DIFF_ALLTIME.rds')

DIFF_NAME=names(DIFF)
###########################################################################



TIME=pbmc$batch
TYPE=pbmc$final.clst.sub.label
IDENT=pbmc$final.clst.sub.label.time
TYPE_LIST=unique(TYPE)
TIME_LIST=unique(TIME)



OUT=matrix(0,nrow=length(TYPE_LIST),ncol=3)
colnames(OUT)=c(    'ident1.utricle.12m.vs.ident2.utricle.3m',
                    'ident1.utricle.22m.vs.ident2.utricle.12m',
                    'ident1.utricle.22m.vs.ident2.utricle.3m'
                    )
rownames(OUT)=TYPE_LIST

PCUT=0.05

i=1
while(i<=nrow(OUT)){
    this_row_name=rownames(OUT)[i]
    j=1
    while(j<=ncol(OUT)){
        this_col_name=colnames(OUT)[j]
        this_tag=paste0(this_col_name,'.',this_row_name)
        if(this_tag %in% names(DIFF)){
            this_diff=DIFF[[this_tag]]
            this_n=length(which(this_diff$p_val_adj<PCUT))
            #this_n=length(which(this_diff$p_val<PCUT))
            OUT[i,j]=this_n
            }
        j=j+1}
    i=i+1}


colnames(OUT)=c('12m.3m','22m.12m','22m.3m')


OUT1=OUT


library('ComplexHeatmap')
library('circlize')

o.mat=OUT
color_fun_3 =colorRamp2(c(0,300,500 ), c('grey70','red1','yellow1'))
ht=Heatmap(o.mat,row_title='',name="exp",
        cluster_columns=FALSE, cluster_rows=TRUE,
	      show_column_dend = FALSE, show_row_dend = TRUE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun_3, border = TRUE
	      )

print(ht)


pdf('UTRICLE_20211025/DIFF_UP_HEAT.pdf',width=6,height=5)
print(ht)
dev.off()
#############################




OUT=matrix(0,nrow=length(TYPE_LIST),ncol=3)
colnames(OUT)=c(    'ident1.utricle.3m.vs.ident2.utricle.12m',
                    'ident1.utricle.12m.vs.ident2.utricle.22m',
                    'ident1.utricle.3m.vs.ident2.utricle.22m'
                    )
rownames(OUT)=TYPE_LIST

PCUT=0.05

i=1
while(i<=nrow(OUT)){
    this_row_name=rownames(OUT)[i]
    j=1
    while(j<=ncol(OUT)){
        this_col_name=colnames(OUT)[j]
        this_tag=paste0(this_col_name,'.',this_row_name)
        if(this_tag %in% names(DIFF)){
            this_diff=DIFF[[this_tag]]
            this_n=length(which(this_diff$p_val_adj<PCUT))
            #this_n=length(which(this_diff$p_val<PCUT))
            OUT[i,j]=this_n
            }
        j=j+1}
    i=i+1}


colnames(OUT)=c('3m.12m','12m.22m','3m.22m')
OUT2=OUT

library('ComplexHeatmap')
library('circlize')

o.mat=OUT
color_fun_3 =colorRamp2(c(0,300,500 ), c('grey70','royalblue3','cyan'))
ht=Heatmap(o.mat,row_title='',name="exp",
        cluster_columns=FALSE, cluster_rows=TRUE,
	      show_column_dend = FALSE, show_row_dend = TRUE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun_3, border = TRUE
	      )

print(ht)

pdf('UTRICLE_20211025/DIFF_DW_HEAT.pdf',width=6,height=5)
print(ht)
dev.off()

#############################





COM=cbind(-OUT2,OUT1)




library('ComplexHeatmap')
library('circlize')

o.mat=COM
color_fun_3 =colorRamp2(c(-500,-300,0,300,500 ), c('cyan','royalblue3','grey70','red1','yellow1'))
ht=Heatmap(o.mat,row_title='',name="exp",
        cluster_columns=FALSE, cluster_rows=TRUE,
	      show_column_dend = FALSE, show_row_dend = TRUE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun_3, border = TRUE
	      )

print(ht)

pdf('UTRICLE_20211025/DIFF_UPandDW_HEAT.pdf',width=6,height=5)
print(ht)
dev.off()



#############
write.table(table(pbmc$final.clst.sub.label.time),row.names=T,col.names=F,file='UTRICLE_20211025/TIME_CELL_NUMBER_FINAL.txt',sep='\t',quote=F)
##############






TAB=table(TYPE,TIME)
TAB.MIN=apply(TAB,1,min)

#TAB.MIN[which(TAB.MIN>100)]=100


###############
#  转录噪音
###########################
set.seed(123)
DS_CELL_INDEX=c()

i=1
while(i<=length(TAB.MIN)){
    this_type=names(TAB.MIN)[i]
    this_index=which(TYPE==this_type)
    this_min=TAB.MIN[i]
    this_time=TIME[this_index]
    if(this_min>=3){
        this_index1=sample(x=this_index[which(this_time=='utricle.3m')],this_min)
        this_index2=sample(x=this_index[which(this_time=='utricle.12m')],this_min)
        this_index3=sample(x=this_index[which(this_time=='utricle.22m')],this_min)
        }
    DS_CELL_INDEX=c(DS_CELL_INDEX, this_index1,this_index2,this_index3)
    i=i+1
}




#########################################
COUNT_NUMBER=1000

EXP_MAT=as.matrix(pbmc@assays$RNA@counts[,DS_CELL_INDEX])

#EXP_DS_MAT=EXP_MAT

PROP=COUNT_NUMBER/colSums(EXP_MAT)
PROP[which(PROP>1)]=1

set.seed(123)
EXP_DS_MAT=scuttle::downsampleMatrix(x=EXP_MAT, prop=PROP, bycol = TRUE, sink = NULL)
EXP_DS_MAT=as.matrix(EXP_DS_MAT)


########################################################################
#########################################


GENE_SUM_EXP=rowSums(EXP_DS_MAT)
EXP_DS_MAT=EXP_DS_MAT[order(GENE_SUM_EXP),]
GENE_SUM_EXP=GENE_SUM_EXP[order(GENE_SUM_EXP)]


##########################
EACH_BIN_SIZE=nrow(EXP_DS_MAT) %/% 10 +1
BIN=rep(1:10,each=EACH_BIN_SIZE)
BIN=BIN[1:nrow(EXP_DS_MAT)]
###################
GENE_VAR_EXP=apply(EXP_DS_MAT,1,var)

USED_INDEX=c()
i=2
while(i<max(BIN)){
    this_index=which(BIN==i)
    this_var=GENE_VAR_EXP[this_index]
    this_used_index=this_index[which(this_var<quantile(this_var,0.1))]
    USED_INDEX=c(USED_INDEX,this_used_index)
    i=i+1
    } 

EXP_DS_MAT=EXP_DS_MAT[USED_INDEX,]
EXP_DS_MAT=sqrt(EXP_DS_MAT)
############################################



EXP_DS_TAG=as.character(Idents(pbmc)[DS_CELL_INDEX])
EXP_DS_TAG_LIST=unique(EXP_DS_TAG)

NOISE=list()

i=1
while(i<=length(EXP_DS_TAG_LIST)){
    this_tag=EXP_DS_TAG_LIST[i]
    this_index=which(EXP_DS_TAG==this_tag)
    if(length(this_index)>3){
        this_exp=EXP_DS_MAT[,this_index]
        this_mean= apply(this_exp,2,mean)
        this_noise= sqrt(apply((this_exp-this_mean)**2,2,sum))
        }else{
        this_noise=0    
        }
    NOISE[[this_tag]]=this_noise
    print(i)
    i=i+1
    }

par(las=2)
par(mar=c(15,5,5,5))
boxplot(NOISE,outline=FALSE,col=c('grey70','red1','gold1'))


#saveRDS(NOISE, file='UTRICLE_20211025/NOISE.rds')


pdf('UTRICLE_20211025/NOISE.pdf',width=20,height=10)
par(las=2)
par(mar=c(20,5,5,5))
boxplot(NOISE,outline=FALSE,col=c('grey70','red1','gold1'))
dev.off()


################################################################



























