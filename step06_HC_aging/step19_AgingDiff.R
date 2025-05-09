
#######################################################################################################
# Part3-Organ/cell-specific aging
#######################################################################################################


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


#########################
#UTRICLE
###############

pbmc=UTRICLE

#########################
TYPE=pbmc$new.clst.sub
############################
TIME=pbmc$batch
TIME[which(TIME %in% c('utricle.3m'))]='utricle.03m'

#####################################
#TIME[which(pbmc$batch=='cochlea.24mbu')]='cochlea.24m'
LABEL=paste0(TYPE,'.',TIME)

########################
DATA=pbmc[['RNA']]@data

REF=.generate_mean(DATA,LABEL)

REF=REF[,order(colnames(REF))]
OUT=cbind(rownames(REF),REF)
colnames(OUT)[1]='gene'

write.table(OUT,'/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/TIME_CLUSTER_UTRICLE_EXP.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)



#########################
#COCHLEA
###############

pbmc=COCHLEA

#########################
TYPE=pbmc$new.clst.sub
############################
TIME=pbmc$batch
TIME[which(TIME %in% c('cochlea.3m'))]='cochlea.03m'
TIME[which(pbmc$batch=='cochlea.24mbu')]='cochlea.24m'

#####################################
LABEL=paste0(TYPE,'.',TIME)

########################
DATA=pbmc[['RNA']]@data

REF=.generate_mean(DATA,LABEL)

REF=REF[,order(colnames(REF))]
OUT=cbind(rownames(REF),REF)
colnames(OUT)[1]='gene'

write.table(OUT,'/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/TIME_CLUSTER_COCHLEA_EXP.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)






###################
## 计算 overdispersion
#############
library(scran)

.cal_overdisp <-function(counts){
    means <- rowMeans(counts)
    cv2 <- apply(counts, 1, var)/means^2
    dm.stat <- DM(means, cv2)
    return(dm.stat)
    }
##################################






COMBINE$time[which(COMBINE$time %in% c('U.3m'))]='U.03m'
COMBINE$time[which(COMBINE$time %in% c('C.3m'))]='C.03m'

COUNT=COMBINE@assays$RNA@counts



#################################################

TAG=paste0(COMBINE$clst,'_',COMBINE$time)
TISSUE=c(rep('U',ncol(UTRICLE)),rep('C',ncol(COCHLEA)))



###############################################################################################################
# UTRICLE
##########

USED=which(TISSUE=='U')

COUNT_USED=COUNT[which(rownames(COUNT) %in% VariableFeatures(COMBINE)),USED]

TAG_USED=TAG[USED]
TAG_USED_LIST=sort(unique(TAG_USED))


DISP=list()

i=1
while(i<=length(TAG_USED_LIST)){
    this_tag=TAG_USED_LIST[i]
    this_index=which(TAG_USED==this_tag)
    if(length(this_index)>1){
        this_count=COUNT_USED[,this_index]
        this_count=this_count[which(rowSums(this_count)>0),]
        this_disp=.cal_overdisp(this_count)
    }else{
        this_disp=NA    
        }    
    DISP[[this_tag]]=this_disp 
    print(i)
    print(this_tag)
    i=i+1
    }


MEAN=unlist(lapply(DISP,mean))
SD=unlist(lapply(DISP,sd))
N=unlist(lapply(DISP,length))
SEM=SD/sqrt(N)
UP=MEAN+SEM
LW=MEAN-SEM


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G01_UTRICLE_OVERDISP.pdf',width=18,height=6)
par(mar=c(12,6,6,6))
bp=barplot(MEAN,las=2,ylim=c(0,max(UP[which(!is.na(UP))])),col=c('royalblue1','grey70','indianred1'),
           space=c(0, c(rep(c(0,0,0.5),length(MEAN)/3)))[1:length(MEAN)]
           )

LWD=1
segments(x0=bp,x1=bp,y0=UP,y1=LW,lwd=LWD)
epsilon=0.2
segments(x0=bp-epsilon,y0=UP,x1=bp+epsilon,y1=UP,lwd=LWD)
segments(x0=bp-epsilon,y0=LW,x1=bp+epsilon,y1=LW,lwd=LWD)
dev.off()
##################################################################








###############################################################################################################
# COCHLEA
##########

USED=which(TISSUE=='C')

COUNT_USED=COUNT[which(rownames(COUNT) %in% VariableFeatures(COMBINE)),USED]


TAG_USED=TAG[USED]
TAG_USED_LIST=sort(c(unique(TAG_USED),'C.18_C.03m','C.19_C.24m'))


DISP=list()

i=1
while(i<=length(TAG_USED_LIST)){
    this_tag=TAG_USED_LIST[i]
    this_index=which(TAG_USED==this_tag)
    if(length(this_index)>1){
        this_count=COUNT_USED[,this_index]
        this_count=this_count[which(rowSums(this_count)>0),]
        this_disp=.cal_overdisp(this_count)
    }else{
        this_disp=NA    
        }    
    DISP[[this_tag]]=this_disp 
    print(i)
    print(this_tag)
    i=i+1
    }


MEAN=unlist(lapply(DISP,mean))
SD=unlist(lapply(DISP,sd))
N=unlist(lapply(DISP,length))
SEM=SD/sqrt(N)
UP=MEAN+SEM
LW=MEAN-SEM



pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G02_COCHLEA_OVERDISP.pdf',width=18,height=6)
par(mar=c(12,6,6,6))
bp=barplot(MEAN,las=2,ylim=c(0,max(UP[which(!is.na(UP))])),col=c('royalblue1','grey70','indianred1'),
           space=c(0, c(rep(c(0,0,0.5),length(MEAN)/3)))[1:length(MEAN)]
           )

LWD=1
segments(x0=bp,x1=bp,y0=UP,y1=LW,lwd=LWD)
epsilon=0.2
segments(x0=bp-epsilon,y0=UP,x1=bp+epsilon,y1=UP,lwd=LWD)
segments(x0=bp-epsilon,y0=LW,x1=bp+epsilon,y1=LW,lwd=LWD)
dev.off()

##################################################################














###################
## 计算 heterogeneity
#############
library(scran)

.cal_het <-function(counts){
    means <- rowMeans(counts)
    diff=counts-means
    distance=sqrt(colSums(diff**2))
    return(distance)
    }
##################################




COMBINE$time[which(COMBINE$time %in% c('U.3m'))]='U.03m'
COMBINE$time[which(COMBINE$time %in% c('C.3m'))]='C.03m'

COUNT=COMBINE@assays$RNA@counts



#################################################

TAG=paste0(COMBINE$clst,'_',COMBINE$time)
TISSUE=c(rep('U',ncol(UTRICLE)),rep('C',ncol(COCHLEA)))



###############################################################################################################
# UTRICLE
##########

USED=which(TISSUE=='U')

COUNT_USED=COUNT[,USED]

TAG_USED=TAG[USED]
TAG_USED_LIST=sort(unique(TAG_USED))


HET=list()

i=1
while(i<=length(TAG_USED_LIST)){
    this_tag=TAG_USED_LIST[i]
    this_index=which(TAG_USED==this_tag)
    if(length(this_index)>1){
        this_count=COUNT_USED[,this_index]
        this_count=this_count[which(rowSums(this_count)>0),]
        this_het=.cal_het(this_count)
    }else{
        this_het=NA    
        }    
    HET[[this_tag]]=this_het 
    print(i)
    print(this_tag)
    i=i+1
    }


MEAN=unlist(lapply(HET,mean))
SD=unlist(lapply(HET,sd))
N=unlist(lapply(HET,length))
SEM=SD/sqrt(N)
UP=MEAN+SEM
LW=MEAN-SEM




pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G03_UTRICLE_DISTANCE.pdf',width=18,height=6)
par(mar=c(12,6,6,6))
bp=barplot(MEAN,las=2,ylim=c(0,max(UP[which(!is.na(UP))])),col=c('royalblue1','grey70','indianred1'),
           space=c(0, c(rep(c(0,0,0.5),length(MEAN)/3)))[1:length(MEAN)]
           )

LWD=1
segments(x0=bp,x1=bp,y0=UP,y1=LW,lwd=LWD)
epsilon=0.2
segments(x0=bp-epsilon,y0=UP,x1=bp+epsilon,y1=UP,lwd=LWD)
segments(x0=bp-epsilon,y0=LW,x1=bp+epsilon,y1=LW,lwd=LWD)




par(mar=c(12,6,6,6))
boxplot(HET,las=2,outline=FALSE,col=c('royalblue1','grey70','indianred1'))

dev.off()
##################################################################









###############################################################################################################
# COCHLEA
##########

USED=which(TISSUE=='C')

COUNT_USED=COUNT[,USED]

TAG_USED=TAG[USED]
TAG_USED_LIST=sort(c(unique(TAG_USED),'C.18_C.03m','C.19_C.24m'))


HET=list()

i=1
while(i<=length(TAG_USED_LIST)){
    this_tag=TAG_USED_LIST[i]
    this_index=which(TAG_USED==this_tag)
    if(length(this_index)>1){
        this_count=COUNT_USED[,this_index]
        this_count=this_count[which(rowSums(this_count)>0),]
        this_het=.cal_het(this_count)
    }else{
        this_het=NA    
        }    
    HET[[this_tag]]=this_het 
    print(i)
    print(this_tag)
    i=i+1
    }


MEAN=unlist(lapply(HET,mean))
SD=unlist(lapply(HET,sd))
N=unlist(lapply(HET,length))
SEM=SD/sqrt(N)
UP=MEAN+SEM
LW=MEAN-SEM


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G04_COCHLEA_DISTANCE.pdf',width=18,height=6)
par(mar=c(12,6,6,6))
bp=barplot(MEAN,las=2,ylim=c(0,max(UP[which(!is.na(UP))])),col=c('royalblue1','grey70','indianred1'),
           space=c(0, c(rep(c(0,0,0.5),length(MEAN)/3)))[1:length(MEAN)]
           )

LWD=1
segments(x0=bp,x1=bp,y0=UP,y1=LW,lwd=LWD)
epsilon=0.2
segments(x0=bp-epsilon,y0=UP,x1=bp+epsilon,y1=UP,lwd=LWD)
segments(x0=bp-epsilon,y0=LW,x1=bp+epsilon,y1=LW,lwd=LWD)


par(mar=c(12,6,6,6))
boxplot(HET,las=2,outline=FALSE,col=c('royalblue1','grey70','indianred1'))

dev.off()
##################################################################
















###################
## 计算 DEG 数量，绘制HEATMAP
#############



#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')

library(Seurat)
library(dplyr)
library(patchwork)

UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')



#########################
#UTRICLE
###############

pbmc=UTRICLE

#########################
TYPE=pbmc$new.clst.sub
############################
TIME=pbmc$batch
#TIME[which(pbmc$batch=='cochlea.24mbu')]='cochlea.24m'
LABEL=paste0(TIME,'.',TYPE)

########################

IDENT=LABEL
Idents(pbmc)=IDENT
#########################
#########################


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


saveRDS(DIFF,file='./DIFF_UTRICLE_ALLTIME.rds' )

DIFF=readRDS('./DIFF_UTRICLE_ALLTIME.rds')

DIFF_NAME=names(DIFF)





##############################################


TYPE_LIST=unique(TYPE)
TIME_LIST=unique(TIME)



OUT=matrix(0,nrow=length(TYPE_LIST),ncol=3)
colnames(OUT)=c(    'ident1.utricle.12m.vs.ident2.utricle.3m',
                    'ident1.utricle.22m.vs.ident2.utricle.12m',
                    'ident1.utricle.22m.vs.ident2.utricle.3m'
                    )
rownames(OUT)=TYPE_LIST


i=1
while(i<=nrow(OUT)){
    this_row_name=rownames(OUT)[i]
    j=1
    while(j<=ncol(OUT)){
        this_col_name=colnames(OUT)[j]
        this_tag=paste0(this_col_name,'.',this_row_name)
        if(this_tag %in% names(DIFF)){
            this_diff=DIFF[[this_tag]]

            #this_n=length(which(this_diff$p_val_adj<0.05  ))
            this_n=length(which( this_diff$p_val_adj<0.1 & this_diff$p_val <0.05 & this_diff$avg_log2FC >0.25 ))
            

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


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G05_DIFF_UTRICLE_UP_HEAT.pdf',width=6,height=5)
print(ht)
dev.off()
############################




OUT=matrix(0,nrow=length(TYPE_LIST),ncol=3)
colnames(OUT)=c(    'ident1.utricle.3m.vs.ident2.utricle.12m',
                    'ident1.utricle.12m.vs.ident2.utricle.22m',
                    'ident1.utricle.3m.vs.ident2.utricle.22m'
                    )
rownames(OUT)=TYPE_LIST


i=1
while(i<=nrow(OUT)){
    this_row_name=rownames(OUT)[i]
    j=1
    while(j<=ncol(OUT)){
        this_col_name=colnames(OUT)[j]
        this_tag=paste0(this_col_name,'.',this_row_name)
        if(this_tag %in% names(DIFF)){
            this_diff=DIFF[[this_tag]]
            #this_n=length(which(this_diff$p_val_adj<0.05))
            this_n=length(which( this_diff$p_val_adj<0.1 & this_diff$p_val <0.05 & this_diff$avg_log2FC >0.25 ))
            OUT[i,j]=this_n
            }
        j=j+1}
    i=i+1}


colnames(OUT)=c('3m.12m','12m.22m','3m.24m')
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

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G06_DIFF_UTRICLE_DW_HEAT.pdf',width=6,height=5)
print(ht)
dev.off()

#############################






COM=cbind(-OUT2,OUT1)




library('ComplexHeatmap')
library('circlize')

o.mat=COM
color_fun_3 =colorRamp2(c(-500,-300,0,300,500 ), c('cyan','royalblue3','grey80','red1','yellow1'))
ht=Heatmap(o.mat,row_title='',name="exp",
        cluster_columns=FALSE, cluster_rows=TRUE,
	      show_column_dend = FALSE, show_row_dend = TRUE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun_3, border = TRUE
	      )

print(ht)

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G07_DIFF_UTRICLE_UPandDW_HEAT.pdf',width=4,height=5)
print(ht)
dev.off()


#####################################################






#########################
#COCHLEA
###############

pbmc=COCHLEA

#########################
TYPE=pbmc$new.clst.sub
############################
TIME=pbmc$batch
TIME[which(pbmc$batch=='cochlea.24mbu')]='cochlea.24m'
LABEL=paste0(TIME,'.',TYPE)

########################

IDENT=LABEL
Idents(pbmc)=IDENT
#########################
#########################


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


saveRDS(DIFF,file='./DIFF_COCHLEA_ALLTIME.rds' )

DIFF=readRDS('./DIFF_COCHLEA_ALLTIME.rds')

DIFF_NAME=names(DIFF)





##############################################


TYPE_LIST=unique(TYPE)
TIME_LIST=unique(TIME)



OUT=matrix(0,nrow=length(TYPE_LIST),ncol=3)
colnames(OUT)=c(    'ident1.cochlea.12m.vs.ident2.cochlea.3m',
                    'ident1.cochlea.24m.vs.ident2.cochlea.12m',
                    'ident1.cochlea.24m.vs.ident2.cochlea.3m'
                    )
rownames(OUT)=TYPE_LIST


i=1
while(i<=nrow(OUT)){
    this_row_name=rownames(OUT)[i]
    j=1
    while(j<=ncol(OUT)){
        this_col_name=colnames(OUT)[j]
        this_tag=paste0(this_col_name,'.',this_row_name)
        if(this_tag %in% names(DIFF)){
            this_diff=DIFF[[this_tag]]

            #this_n=length(which(this_diff$p_val_adj<0.05  ))
            this_n=length(which( this_diff$p_val_adj<0.1 & this_diff$p_val <0.05 & this_diff$avg_log2FC >0.25 ))
            

            OUT[i,j]=this_n
            }
        j=j+1}
    i=i+1}


colnames(OUT)=c('12m.3m','24m.12m','24m.3m')


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


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G08_DIFF_COCHLEA_UP_HEAT.pdf',width=6,height=5)
print(ht)
dev.off()
############################




OUT=matrix(0,nrow=length(TYPE_LIST),ncol=3)
colnames(OUT)=c(    'ident1.cochlea.3m.vs.ident2.cochlea.12m',
                    'ident1.cochlea.12m.vs.ident2.cochlea.24m',
                    'ident1.cochlea.3m.vs.ident2.cochlea.24m'
                    )
rownames(OUT)=TYPE_LIST


i=1
while(i<=nrow(OUT)){
    this_row_name=rownames(OUT)[i]
    j=1
    while(j<=ncol(OUT)){
        this_col_name=colnames(OUT)[j]
        this_tag=paste0(this_col_name,'.',this_row_name)
        if(this_tag %in% names(DIFF)){
            this_diff=DIFF[[this_tag]]
            #this_n=length(which(this_diff$p_val_adj<0.05))
            this_n=length(which( this_diff$p_val_adj<0.1 & this_diff$p_val <0.05 & this_diff$avg_log2FC >0.25 ))
            OUT[i,j]=this_n
            }
        j=j+1}
    i=i+1}


colnames(OUT)=c('3m.12m','12m.24m','3m.24m')
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

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G09_DIFF_COCHLEA_DW_HEAT.pdf',width=6,height=5)
print(ht)
dev.off()

#############################






COM=cbind(-OUT2,OUT1)




library('ComplexHeatmap')
library('circlize')

o.mat=COM
color_fun_3 =colorRamp2(c(-500,-300,0,300,500 ), c('cyan','royalblue3','grey80','red1','yellow1'))
ht=Heatmap(o.mat,row_title='',name="exp",
        cluster_columns=FALSE, cluster_rows=TRUE,
	      show_column_dend = FALSE, show_row_dend = TRUE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun_3, border = TRUE
	      )

print(ht)

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G10_DIFF_COCHLEA_UPandDW_HEAT.pdf',width=4,height=5)
print(ht)
dev.off()


#####################################################






#########################################################
# DEG overlap 分析
#######################################################


###
# UTRICLE
###

pbmc=UTRICLE
DIFF=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new//DIFF_UTRICLE_ALLTIME.rds')

#########################
TYPE=pbmc$new.clst.sub
############################

TYPE.LIST=unique(TYPE)


SPLIT=strsplit(names(DIFF),'m.')

MAT=matrix(0,nrow=nrow(pbmc),ncol=length(TYPE.LIST))
rownames(MAT)=rownames(pbmc)
colnames(MAT)=TYPE.LIST

#MAT=MAT[grep(rownames(MAT),pattern='mt-',invert=TRUE),]


i=1
while(i<=length(DIFF)){
    if(!'ident1.utricle.12' %in% SPLIT[[i]]){
        this_diff=DIFF[[i]]
        this_cluster=SPLIT[[i]][3]
        this_tag.list_index=which(TYPE.LIST==this_cluster)
        #######################################################
        this_gene=rownames(this_diff)[which( this_diff$p_val_adj<0.1 & this_diff$p_val <0.05 & this_diff$avg_log2FC >0.25 )]
        ########################################################
        this_gene_index=which(rownames(MAT) %in% this_gene)
        MAT[this_gene_index,this_tag.list_index]=MAT[this_gene_index,this_tag.list_index]+1
        #####################
        }
    print(i)
    i=i+1
    }


saveRDS(MAT, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/DIFF_UTRICLE_ALLTIME_MATno12m.rds')

BMAT=MAT
BMAT[which(MAT>0)]=1

BMAT=BMAT[which(rowSums(BMAT)>0),]

NUM=rowSums(BMAT)
TNUM=NUM
TNUM[which(TNUM>10)]=10
TAB=table(TNUM)

barplot(TAB)

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G11_DIFF_UTRICLE_ALLTIME_MATno12m.pdf')
barplot(TAB)
dev.off()

SNUM=sort(NUM, decreasing=T)


write.table(SNUM, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/DIFF_UTRICLE_ALLTIME_MATno12m_SORT.txt',sep='\t',
   row.names=TRUE,col.names=FALSE,quote=FALSE)



library(ggplot2)
library(stringr)

library(DOSE)
library(GSEABase) 
library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)

library(org.Mm.eg.db)

ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(COMBINE), 'ENTREZID', 'SYMBOL')

this_sig_entrez=mapIds(org.Mm.eg.db, names(SNUM[1:100]), 'ENTREZID', 'SYMBOL')

this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    ################################################################
    write.table(this_out, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/GO_UTRICLE_ALLTIME_MATno12m_SORT_TOP100.txt',sep='\t',quote=F,col.names=T,row.names=F)
   



this_sig_entrez=mapIds(org.Mm.eg.db, names(SNUM[1:50]), 'ENTREZID', 'SYMBOL')

this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    ################################################################
    write.table(this_out, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/GO_UTRICLE_ALLTIME_MATno12m_SORT_TOP50.txt',sep='\t',quote=F,col.names=T,row.names=F)
   

##########################




###
# COCHLEA
###

pbmc=COCHLEA
DIFF=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new//DIFF_COCHLEA_ALLTIME.rds')

#########################
TYPE=pbmc$new.clst.sub
############################

TYPE.LIST=unique(TYPE)


SPLIT=strsplit(names(DIFF),'m.')

MAT=matrix(0,nrow=nrow(pbmc),ncol=length(TYPE.LIST))
rownames(MAT)=rownames(pbmc)
colnames(MAT)=TYPE.LIST

#MAT=MAT[grep(rownames(MAT),pattern='mt-',invert=TRUE),]


i=1
while(i<=length(DIFF)){
    if(!'ident1.cochlea.12' %in% SPLIT[[i]]){
        this_diff=DIFF[[i]]
        this_cluster=SPLIT[[i]][3]
        this_tag.list_index=which(TYPE.LIST==this_cluster)
        #######################################################
        this_gene=rownames(this_diff)[which( this_diff$p_val_adj<0.1 & this_diff$p_val <0.05 & this_diff$avg_log2FC >0.25 )]
        ########################################################
        this_gene_index=which(rownames(MAT) %in% this_gene)
        MAT[this_gene_index,this_tag.list_index]=MAT[this_gene_index,this_tag.list_index]+1
        #####################
        }
    print(i)
    i=i+1
    }


saveRDS(MAT, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/DIFF_COCHLEA_ALLTIME_MATno12m.rds')

BMAT=MAT
BMAT[which(MAT>0)]=1

BMAT=BMAT[which(rowSums(BMAT)>0),]

NUM=rowSums(BMAT)
TNUM=NUM
TNUM[which(TNUM>10)]=10
TAB=table(TNUM)

barplot(TAB)

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/G12_DIFF_COCHLEA_ALLTIME_MATno12m.pdf')
barplot(TAB)
dev.off()

SNUM=sort(NUM, decreasing=T)


write.table(SNUM, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/DIFF_COCHLEA_ALLTIME_MATno12m_SORT.txt',sep='\t',
   row.names=TRUE,col.names=FALSE,quote=FALSE)



library(ggplot2)
library(stringr)

library(DOSE)
library(GSEABase) 
library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)

library(org.Mm.eg.db)

ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(COMBINE), 'ENTREZID', 'SYMBOL')

this_sig_entrez=mapIds(org.Mm.eg.db, names(SNUM[1:100]), 'ENTREZID', 'SYMBOL')

this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    ################################################################
    write.table(this_out, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/GO_COCHLEA_ALLTIME_MATno12m_SORT_TOP100.txt',sep='\t',quote=F,col.names=T,row.names=F)
   



this_sig_entrez=mapIds(org.Mm.eg.db, names(SNUM[1:50]), 'ENTREZID', 'SYMBOL')

this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    ################################################################
    write.table(this_out, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/GO_COCHLEA_ALLTIME_MATno12m_SORT_TOP50.txt',sep='\t',quote=F,col.names=T,row.names=F)
   

##########################















###############################################################################################################
# COCHLEA
##########

USED=which(TISSUE=='C')

COUNT_USED=COUNT[which(rownames(COUNT) %in% VariableFeatures(COMBINE)),USED]


TAG_USED=TAG[USED]
TAG_USED_LIST=sort(c(unique(TAG_USED),'C.18_C.03m','C.19_C.24m'))


DISP=list()

i=1
while(i<=length(TAG_USED_LIST)){
    this_tag=TAG_USED_LIST[i]
    this_index=which(TAG_USED==this_tag)
    if(length(this_index)>1){
        this_count=COUNT_USED[,this_index]
        this_count=this_count[which(rowSums(this_count)>0),]
        this_disp=.cal_overdisp(this_count)
    }else{
        this_disp=NA    
        }    
    DISP[[this_tag]]=this_disp 
    print(i)
    print(this_tag)
    i=i+1
    }


MEAN=unlist(lapply(DISP,mean))
SD=unlist(lapply(DISP,sd))
N=unlist(lapply(DISP,length))
SEM=SD/sqrt(N)
UP=MEAN+SEM
LW=MEAN-SEM



par(mar=c(12,6,6,6))
bp=barplot(MEAN,las=2,ylim=c(0,max(UP[which(!is.na(UP))])),col=c('royalblue1','grey70','indianred1'),
           space=c(0, c(rep(c(0,0,0.5),length(MEAN)/3)))[1:length(MEAN)]
           )

LWD=2
segments(x0=bp,x1=bp,y0=UP,y1=LW,lwd=LWD)
epsilon=0.2
segments(x0=bp-epsilon,y0=UP,x1=bp+epsilon,y1=UP,lwd=LWD)
segments(x0=bp-epsilon,y0=LW,x1=bp+epsilon,y1=LW,lwd=LWD)


##################################################################

















#######################################################################################################
# Part4, 5, 6 - Aging 
#######################################################################################################


#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')

library(Seurat)
library(dplyr)
library(patchwork)

UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')


library(ggplot2)
library(stringr)

library(DOSE)
library(GSEABase) 
library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)

library(org.Mm.eg.db)

ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(COMBINE), 'ENTREZID', 'SYMBOL')


.get_diffAndGo = function(pbmc, ident1=ident1, ident2=ident2, folder=folder){
    ################################################################
    this_ident1=ident1
    this_ident2=ident2
    this_folder=folder
    ################################################################
    this_diff=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)
    this_diff$gene=rownames(this_diff)
    ################################################################
    write.table(this_diff, file=paste0(this_folder,'/DIFF.',paste(this_ident1,collapse='.'),'.over.',paste(this_ident2,collapse='.'),'.txt'),col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
    saveRDS(this_diff, file=paste0(this_folder,'/DIFF.',paste(this_ident1,collapse='.'),'.over.',paste(this_ident2,collapse='.'),'.rds'))
    ################################################################

    this_sig_gene=this_diff$gene[which( this_diff$p_val_adj<0.1 & this_diff$p_val <0.05 & this_diff$avg_log2FC >0.25 )]

    #########################################################
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    ################################################################
    write.table(this_out, paste0(this_folder,'/GO.',paste(this_ident1,collapse='.'),'.over.',paste(this_ident2,collapse='.'),'.txt'),sep='\t',quote=F,col.names=T,row.names=F)
    saveRDS(this_out, file=paste0(this_folder,'/GO.',paste(this_ident1,collapse='.'),'.over.',paste(this_ident2,collapse='.'),'.rds'))
    ################################################################
    }


FOLDER='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/diffAndGo_20220509_p0.05adp0.1lfd0.25'


####################################################
pbmc=UTRICLE
LABEL='U'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(5,13,1,2)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)




########################################



####################################################
pbmc=COCHLEA
LABEL='C'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(8,11)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)


########################################








####################################################
pbmc=UTRICLE
LABEL='U'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(16.1)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)


########################################







####################################################
pbmc=COCHLEA
LABEL='C'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(14)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)


########################################







####################################################
pbmc=COCHLEA
LABEL='C'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(4, 5, 9.1, 9.2, 13)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)


########################################




####################################################
pbmc=COCHLEA
LABEL='C'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(0, 2)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)


########################################




####################################################
pbmc=COCHLEA
LABEL='C'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(1.1, 1.2, 5, 6, 7, 10, 15, 16)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_24m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_24m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)


########################################










####################################################
pbmc=UTRICLE
LABEL='U'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(3.1, 3.2, 12, 14)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)


########################################





####################################################
pbmc=UTRICLE
LABEL='U'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(7)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)


########################################







####################################################
pbmc=UTRICLE
LABEL='U'

TAG=paste0(LABEL,'.',pbmc$new.clst.sub,'_',pbmc$time)
table(TAG)
Idents(pbmc)=TAG

CLST=c(0, 11)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_3m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_3m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_12m')
ident2=paste0(LABEL,'.',CLST,'_22m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)

ident1=paste0(LABEL,'.',CLST,'_22m')
ident2=paste0(LABEL,'.',CLST,'_12m')
.get_diffAndGo(pbmc, ident1=ident1, ident2=ident2, folder=FOLDER)


########################################















