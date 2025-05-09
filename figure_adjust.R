#######################################################################################################
# Loading Packages
#######################################################################################################
library(devtools)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidydr)
library(tidyverse)
library(tidyr)
library(ggrepel)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(stringr)
library(DOSE)
library(GSEABase)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)

#######################################################################################################
# UMAP 
#######################################################################################################


#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')

UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F01_COCHLEA_ALL_UMAP.pdf',width=10,height=10)
DimPlot(COCHLEA, reduction='umap', group.by='new.clst.sub',pt.size=0.1,label=T,label.size=8) + NoLegend()
dev.off()

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F02_COCHLEA_3M_UMAP.pdf',width=10,height=10)
USED_CELL=colnames(COCHLEA)[which(COCHLEA$time=='3m')]
DimPlot(COCHLEA, cells=USED_CELL, reduction='umap', group.by='new.clst.sub',pt.size=0.5,label=T,label.size=8) + NoLegend()
dev.off()

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F03_COCHLEA_12M_UMAP.pdf',width=10,height=10)
USED_CELL=colnames(COCHLEA)[which(COCHLEA$time=='12m')]
DimPlot(COCHLEA, cells=USED_CELL, reduction='umap', group.by='new.clst.sub',pt.size=0.5,label=T,label.size=8) + NoLegend()
dev.off()

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F04_COCHLEA_24M_UMAP.pdf',width=10,height=10)
USED_CELL=colnames(COCHLEA)[which(COCHLEA$time=='24m')]
DimPlot(COCHLEA, cells=USED_CELL, reduction='umap', group.by='new.clst.sub',pt.size=0.5,label=T,label.size=8) + NoLegend()
dev.off()


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F05_UTRICLE_ALL_UMAP.pdf',width=10,height=10)
DimPlot(UTRICLE, reduction='umap', group.by='new.clst.sub',pt.size=0.1,label=T,label.size=8) + NoLegend()
dev.off()

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F06_UTRICLE_3M_UMAP.pdf',width=10,height=10)
USED_CELL=colnames(UTRICLE)[which(UTRICLE$time=='3m')]
DimPlot(UTRICLE, cells=USED_CELL, reduction='umap', group.by='new.clst.sub',pt.size=0.5,label=T,label.size=8) + NoLegend()
dev.off()

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F07_UTRICLE_12M_UMAP.pdf',width=10,height=10)
USED_CELL=colnames(UTRICLE)[which(UTRICLE$time=='12m')]
DimPlot(UTRICLE, cells=USED_CELL, reduction='umap', group.by='new.clst.sub',pt.size=0.5,label=T,label.size=8) + NoLegend()
dev.off()

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F08_UTRICLE_22M_UMAP.pdf',width=10,height=10)
USED_CELL=colnames(UTRICLE)[which(UTRICLE$time=='22m')]
DimPlot(UTRICLE, cells=USED_CELL, reduction='umap', group.by='new.clst.sub',pt.size=0.5,label=T,label.size=8) + NoLegend()
dev.off()


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F09_COMBINE_ALL_UMAP.pdf',width=12,height=10)
DimPlot(COMBINE, reduction='umap', group.by='clst',pt.size=0.1,label=T,label.size=8) + NoLegend()
dev.off()


#######################################################################################################
# DotPlot
#######################################################################################################


####################################################################
pbmc=UTRICLE
LEVEL=c('15','18','0','3.1','14','3.2','16.1','19','16.2','17','4','8','7',
        '10','11','12','9','5','13','1','2')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

#pbmc@assays$RNA@data=pbmc@assays$RNA@counts
#pbmc=NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.gene=rownames(pbmc)
pbmc=ScaleData(object = pbmc, features =all.gene)

################################################

MK1=c('Drp2','Plp1')
MK2=c('Cyyr1','Abcb1a')
MK3=c('Bmp6','Adam12')
MK4=c('Map2k6','Cacna1e')
MK5=c('Actn2','Nrp1')
MK6=c('Ucma','Cntn6')
MK7=c('Ptprc','Ly86')
MK8=c('Dct','Tyr')
MK9=c('C1qc','Tyrobp')
MK10=c('Nrg1','Nrg3')
MK11=c('Nell1','Slc26a4')
MK12=c('Kcnj3','Unc5cl')
MK13=c('Frem2','Ntn1')
MK14=c('Esrrb','Rspo3')
MK15=c('Isl1','Tectb')
MK16=c('Creb5','Dcn')
MK17=c('Rnf149','Tbx3')
MK18=c('Ptprq','Spp1')
MK19=c('Xirp2','Tnc')
MK20=c('Anxa4','Myo7a')
MK21=c('Mapt','Myo15b')



MKS=c(MK1,MK2,MK3,MK4,MK5,MK6,MK7,MK8,MK9,MK10,MK11,MK12,MK13,MK14,MK15,MK16,MK17,MK18,MK19,MK20,MK21)

p=DotPlot(pbmc, features=MKS,cols = c("lightgrey", "indianred3"))+ theme(axis.text.x = element_text(angle = 90))
custom_order <- rev(LEVEL)
p$data$id <- factor(p$data$id, levels = custom_order)
print(p)


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F10_UTRICLE_DOT_PLOT.pdf',width=14,height=10)

print(p)

dev.off()





####################################################################
pbmc=COCHLEA
LEVEL=c('0','1.1','1.2','2','4','5','6','7','8','9.1','9.2','10','11',
        '12','13','14','15','16','17','18','19')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

#pbmc@assays$RNA@data=pbmc@assays$RNA@counts
#pbmc=NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.gene=rownames(pbmc)
pbmc=ScaleData(object = pbmc, features =all.gene)

################################################


MK1=c('Emcn','Emilin2')
MK2=c('Tgfb3','Slc26a7')
MK3=c('Epyc','Lgr5')
MK4=c('Efemp1','Kcnb2')
MK5=c('Dcn','Scn7a')
MK6=c('Cd44') # Epyc
MK7=c('Otx2','Cnmd')
MK8=c('Rspo2','Ltbp2')
MK9=c('Dnm1','Slc17a8')
MK10=c('Actn2','Mapk4')
MK11=c('Sorcs3','Car3')
MK12=c('Gas2','Cep41')
MK13=c('Strip2','Slc26a5')
MK14=c('Nrg3','Snap25')
MK15=c('Col24a1','Grm7')
MK16=c('Ptprc','Ly86')
MK17=c('Cpa6','Gjb6')
MK18=c('Mdm1','Fgfr3')
MK19=c('Drp2','Prx')
MK20=c('Dct','Tyr')
MK21=c('Ptprb','Cyyr1')




MKS=c(MK1,MK2,MK3,MK4,MK5,MK6,MK7,MK8,MK9,MK10,MK11,MK12,MK13,MK14,MK15,MK16,MK17,MK18,MK19,MK20,MK21)

p=DotPlot(pbmc, features=MKS, cols = c("lightgrey", "indianred3"))+ theme(axis.text.x = element_text(angle = 90))
custom_order <- rev(LEVEL)
p$data$id <- factor(p$data$id, levels = custom_order)
print(p)


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F11_COCHLEA_DOT_PLOT.pdf',width=14,height=10)

print(p)

dev.off()



#################################################################














#######################################################################################################
# Finding Markers， Heatmap
#######################################################################################################


#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')


UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')


#########################################
pbmc=UTRICLE 
LEVEL=c('15','18','0','3.1','14','3.2','16.1','19','16.2','17','4','8','7','10','11','12','9','5','13','1','2')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pbmc.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_MAKER.rds')
write.table(pbmc.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_MAKER.tsv',col.names=T,row.names=F,sep='\t',quote=F)
############################



#########################################
pbmc=COCHLEA
LEVEL=c('0','1.1','1.2','2','4','5','6','7','8','9.1','9.2','10','11','12','13','14','15','16','17','18','19')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pbmc.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_MAKER.rds')
write.table(pbmc.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_MAKER.tsv',col.names=T,row.names=F,sep='\t',quote=F)
############################




#########################################
pbmc=UTRICLE 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_MAKER.rds')
LEVEL=c('15','18','0','3.1','14','3.2','16.1','19','16.2','17','4','8','7','10','11','12','9','5','13','1','2')
###############
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

DATA=pbmc@assays$RNA@data
CLST=as.character(Idents(pbmc))


##################################
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
MEAN=.generate_mean(DATA,CLST)
MEAN.BAC=MEAN
colnames(MEAN)=LEVEL
i=1
while(i<=ncol(MEAN)){
  MEAN[,i]=MEAN.BAC[,which(colnames(MEAN.BAC)==colnames(MEAN)[i])]
  i=i+1
}

##################################
library(dplyr)
N=10
topN <- pbmc.markers %>% group_by(cluster) %>% top_n(n = N, wt = avg_log2FC)


MAT=matrix(0,ncol=ncol(MEAN),nrow=nrow(topN ))
i=1
while(i<=nrow(topN )){
  this_gene=topN$gene[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN$gene
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F12_UTRICLE_HEAT.pdf',width=30,height=30)
color_fun_3=colorRamp2(c(-4,-2,0.5,0,0.5,2,4 ), c('grey95','grey95','grey90','grey90','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
        show_column_dend = FALSE, show_row_dend = FALSE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()






#########################################
pbmc=COCHLEA
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_MAKER.rds')
LEVEL=c('0','1.1','1.2','2','4','5','6','7','8','9.1','9.2','10','11','12','13','14','15','16','17','18','19')
######################

Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

DATA=pbmc@assays$RNA@data
CLST=as.character(Idents(pbmc))


##################################
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
MEAN=.generate_mean(DATA,CLST)
MEAN.BAC=MEAN
colnames(MEAN)=LEVEL
i=1
while(i<=ncol(MEAN)){
  MEAN[,i]=MEAN.BAC[,which(colnames(MEAN.BAC)==colnames(MEAN)[i])]
  i=i+1
}

##################################

N=10
topN <- pbmc.markers %>% group_by(cluster) %>% top_n(n = N, wt = avg_log2FC)


MAT=matrix(0,ncol=ncol(MEAN),nrow=nrow(topN ))
i=1
while(i<=nrow(topN )){
  this_gene=topN$gene[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN$gene
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F13_COCHLEA_HEAT.pdf',width=30,height=30)
color_fun_3=colorRamp2(c(-4,-2,0.5,0,0.5,2,4 ), c('grey95','grey95','grey90','grey90','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
        show_column_dend = FALSE, show_row_dend = FALSE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()


#################################################################################################

#######################################################################################################
# Finding Top Markers， Dotplot
#######################################################################################################

cochlea_markers <- c("Emcn", "Emilin2", "Tgfb3", "Slc26a7", "Ерус", "Lgr5", "Efemp1", "Kcnb2", 
                     "Dcn",  "Scn7a", "Cd44","Otx2", "Cnmd", "Rspo2", "Ltbp2", "Dnm1", "Slc17a8",
                     "Actn2","Mapk4", "Sorcs3", "Car3", "Gas2", "Cep41", "Strip2", "Slc26a5","Nrg3",
                     "Snap25", "Col24a1", "Grm7","Ptprc", "Ly86", "Cpa6", "Gjb6", "Mdm1", "Fgfr3", "Drp2", "Prx",
                     "Dct", "Tyr", "Ptprb", "Cyyr1")

utricle_markers <- c("Drp2","Plp1", "Cyyr1", "Abcb1a","Bmp6","Adam12", "Map2k6", "Cacna1e","Actn2", "Nrp1", "Ucma", "Cntn6" , "Ptprc" , "Ly86","Dct", "Tyr", "C1qc", "Tyrobp", "Nrg1", "Nrg3", "Nell1","Slc26a4", "Kcnj3", "Unc5cl", "Frem2", "Ntn1", "Esrrb", "Rspo3",
                     "Isl1", "Tectb", "Creb5", "Dcn", "Rnf149", "Tbx3", "Ptprq", "Spp1",
                     "Xirp2", "Tnc", "Anxa4", "Myo7a", "Mapt", "Myo15b")
common_markers <- c("Ptgds", "Apoe", "Atp1a2", "Ddx5", "Otos", "Apod",  "Rplp1", "Fth1", "Hsp90aa1",  "Rpl13", "Hsp90ab1", "Tcf25", "Tpt1",
                    "Msi2", "Gpc6", "Grip1", "Lsamp", "Pde4d", "Ptprt", "Txnip",
                    "Atp1a2", "Slc1a3", "Zdhhc", "Camk1d", "Ddx60", "Cdk8", "Cntfr", "Gabra2", "Rasgrf1",
                    "Auts2", "Gpc6", "Cdh18", "Msi2", "Lsamp", "Rora", "Adgrl3", "Maml3", "Ptprt", "Fgfr2", "Foxp1",  "Grip1", "Kirrel3", "Mmp16", "Syn3", "Ust", "Wdr17", "Fbn2", "H3f3b", "Luc7l2", "Mtmr3","Tenm3", "Timp3", "Txnip")

this_pbmc = UTRICLE
markers = utricle_markers

this_pbmc = COCHLEA
markers = cochlea_markers 

this_pbmc = COMBINE
markers = common_markers

DotPlot(this_pbmc, features = markers) +theme(axis.text.x = element_text(angle = 90))+theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())+scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+theme(panel.border = element_rect(color = "black"), panel.spacing = unit(1, "mm"))



#######################################################################################################
# GO analysis
#######################################################################################################


library(ggplot2)
library(stringr)

library(DOSE)
library(GSEABase) 
library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Mm.eg.db)

#########################################
pbmc=UTRICLE 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_MAKER.rds')
LEVEL=c('15','18','0','3.1','14','3.2','16.1','19','16.2','17','4','8','7','10','11','12','9','5','13','1','2')
FOLDER='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_GO/'
################

Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

UNIQ_CLST=unique(Idents(pbmc))
ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(pbmc), 'ENTREZID', 'SYMBOL')


i=1
while(i<=length(UNIQ_CLST)){
  this_clst=UNIQ_CLST[i]
  this_sig_gene=pbmc.markers$gene[which(pbmc.markers$cluster==this_clst & pbmc.markers$p_val_adj<0.05)]
  this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
  
  this_ego <- enrichGO(gene = this_sig_entrez, 
                       universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                       pAdjustMethod = "BH",pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1,readable = TRUE)
  
  this_out=this_ego@result
  write.table(this_out, paste0(FOLDER,this_clst,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)
  
  #############################               
  print(paste0(i,' / ',length(UNIQ_CLST)))
  print(this_clst)
  i=i+1
}


#########################################
pbmc=COCHLEA 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_MAKER.rds')
LEVEL=c('0','1.1','1.2','2','4','5','6','7','8','9.1','9.2','10','11','12','13','14','15','16','17','18','19')
FOLDER='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_GO/'
################

Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

UNIQ_CLST=unique(Idents(pbmc))
ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(pbmc), 'ENTREZID', 'SYMBOL')


i=1
while(i<=length(UNIQ_CLST)){
  this_clst=UNIQ_CLST[i]
  this_sig_gene=pbmc.markers$gene[which(pbmc.markers$cluster==this_clst & pbmc.markers$p_val_adj<0.05)]
  this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
  
  this_ego <- enrichGO(gene = this_sig_entrez, 
                       universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                       pAdjustMethod = "BH",pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1,readable = TRUE)
  
  this_out=this_ego@result
  write.table(this_out, paste0(FOLDER,this_clst,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)
  
  #############################               
  print(paste0(i,' / ',length(UNIQ_CLST)))
  print(this_clst)
  i=i+1
}


#################################################################




#######################################################################################################
# FeaturePlot & VlnPlot
#######################################################################################################

library(ggplot2)
library(stringr)

#########################################
pbmc=UTRICLE 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_MAKER.rds')
LEVEL=c('15','18','0','3.1','14','3.2','16.1','19','16.2','17','4','8','7','10','11','12','9','5','13','1','2')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)
################

show_genes=c('Mapt','Isl1','Bmp6','Kcnj3','Nell1','Esrrb','Frem2','Ntn1','Rspo2')
################

OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F14_UTRICLE_FEATUREPLOT.pdf'
pdf(OUTPUT_PATH,width=16,height=14)

N1=3
N2=2

p1 <- FeaturePlot(pbmc, features=show_genes,ncol=3,order=TRUE,pt.size=0.1,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c(rep('grey80',N1),rep('red1',N2)))
p2 <- lapply(p1, function (x) x + fix.sc)
print(CombinePlots(p2))

dev.off()

###########################

OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F15_UTRICLE_VLNPLOT.pdf'

pdf(OUTPUT_PATH,width=10,height=20)

VlnPlot(pbmc, features=show_genes, ncol=2,pt.size=0)

dev.off()









#########################################
pbmc=COCHLEA 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_MAKER.rds')
LEVEL=c('0','1.1','1.2','2','4','5','6','7','8','9.1','9.2','10','11','12','13','14','15','16','17','18','19')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)
################

show_genes_1=c('Klhl14','Liph','Lockd','Epyc','Gpc2','Lgr5','Isl1','Bmp6','Esrrb')
show_genes_2=c('Epyc','Gpc2','Lgr5','Isl1','Bmp6','Esrrb','Frem2','Ntn1','Rspo2')
################

OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F16_COCHLEA_FEATUREPLOT.pdf'

pdf(OUTPUT_PATH,width=16,height=14)

N1=3
N2=2

p1 <- FeaturePlot(pbmc, features=show_genes_1, ncol=3, order=TRUE,pt.size=0.1,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c(rep('grey80',N1),rep('red1',N2)))
p2 <- lapply(p1, function (x) x + fix.sc)
print(CombinePlots(p2))

p1 <- FeaturePlot(pbmc, features=show_genes_2, ncol=3, order=TRUE,pt.size=0.1,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c(rep('grey80',N1),rep('red1',N2)))
p2 <- lapply(p1, function (x) x + fix.sc)
print(CombinePlots(p2))

dev.off()



OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F17_COCHLEA_VLNPLOT.pdf'

pdf(OUTPUT_PATH,width=10,height=20)

VlnPlot(pbmc, features=show_genes_1, ncol=2,pt.size=0)
VlnPlot(pbmc, features=show_genes_2, ncol=2,pt.size=0)

dev.off()













#######################################################################################################
# FeaturePlot & VlnPlot
#######################################################################################################

library(ggplot2)
library(stringr)

#########################################


pbmc=UTRICLE 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_MAKER.rds')
LEVEL=c('15','18','0','3.1','14','3.2','16.1','19','16.2','17','4','8','7','10','11','12','9','5','13','1','2')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)
################

show_genes=c('Kcnh8','Chrm3','Dnm1','Dnm3','Ankfn1','Dgkg','Mapt')
################

OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F18_UTRICLE_FEATUREPLOT.pdf'
pdf(OUTPUT_PATH,width=16,height=14)

N1=3
N2=2

p1 <- FeaturePlot(pbmc, features=show_genes,ncol=3,order=TRUE,pt.size=0.1,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c(rep('grey80',N1),rep('red1',N2)))
p2 <- lapply(p1, function (x) x + fix.sc)
print(CombinePlots(p2))

dev.off()

###########################

OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F19_UTRICLE_VLNPLOT.pdf'

pdf(OUTPUT_PATH,width=10,height=20)

VlnPlot(pbmc, features=show_genes, ncol=2,pt.size=0)

dev.off()









#########################################
pbmc=COCHLEA 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_MAKER.rds')
LEVEL=c('0','1.1','1.2','2','4','5','6','7','8','9.1','9.2','10','11','12','13','14','15','16','17','18','19')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)
################

show_genes_1=c('Ap3b2','Atp2a3','Kcnh8','Chrm3','Dnm1','Dnm3','Ankfn1','Dgkg','Slc26a5')
show_genes_2=c('Atp2a3','Kcnh8','Chrm3','Dnm1','Dnm3','Ankfn1','Dgkg','Slc26a5','Slc17a8')
################

OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F20_COCHLEA_FEATUREPLOT.pdf'

pdf(OUTPUT_PATH,width=16,height=14)

N1=3
N2=2

p1 <- FeaturePlot(pbmc, features=show_genes_1, ncol=3, order=TRUE,pt.size=0.1,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c(rep('grey80',N1),rep('red1',N2)))
p2 <- lapply(p1, function (x) x + fix.sc)
print(CombinePlots(p2))

p1 <- FeaturePlot(pbmc, features=show_genes_2, ncol=3, order=TRUE,pt.size=0.1,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c(rep('grey80',N1),rep('red1',N2)))
p2 <- lapply(p1, function (x) x + fix.sc)
print(CombinePlots(p2))

dev.off()



OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F21_COCHLEA_VLNPLOT.pdf'

pdf(OUTPUT_PATH,width=10,height=20)

VlnPlot(pbmc, features=show_genes_1, ncol=2,pt.size=0)
VlnPlot(pbmc, features=show_genes_2, ncol=2,pt.size=0)

dev.off()







#######################################################################################################
# Time-related VlnPlot
#######################################################################################################


library(ggplot2)
library(stringr)



#########################################
pbmc=UTRICLE 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_MAKER.rds')
LEVEL=c('15','18','0','3.1','14','3.2','16.1','19','16.2','17','4','8','7','10','11','12','9','5','13','1','2')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)
################

# TYPE1.1，5
# TYPE1.2, 13
# TYPE2.1, 1
# TYPE2.2, 2

#####################
pbmc.sub=subset(pbmc, cells=colnames(pbmc)[which(Idents(pbmc) %in% c('5','13','1','2'))])
Idents(pbmc.sub)=factor(pbmc.sub$time,levels=c('3m','12m','22m'))


OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F22_UTRICLE_VLNPLOT_HC_TIME.pdf'
pdf(OUTPUT_PATH,width=10,height=10)
VlnPlot(pbmc.sub, features=c('Kcnh8','Chrm3','Ankfn1','Dgkg'), ncol=2,pt.size=0)
dev.off()


#####################
pbmc.sub=subset(pbmc, cells=colnames(pbmc)[which(Idents(pbmc) %in% c('5','13'))])
Idents(pbmc.sub)=factor(pbmc.sub$time,levels=c('3m','12m','22m'))


OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F23_UTRICLE_VLNPLOT_HC.TYPE1_TIME.pdf'
pdf(OUTPUT_PATH,width=5.5,height=4)
VlnPlot(pbmc.sub, features=c('Dnm1','Dnm3'), ncol=2,pt.size=0)
dev.off()




#####################
pbmc.sub=subset(pbmc, cells=colnames(pbmc)[which(Idents(pbmc) %in% c('1','2'))])
Idents(pbmc.sub)=factor(pbmc.sub$time,levels=c('3m','12m','22m'))


OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F24_UTRICLE_VLNPLOT_HC.TYPE2_TIME.pdf'
pdf(OUTPUT_PATH,width=3.5,height=4)
VlnPlot(pbmc.sub, features=c('Mapt'), ncol=1,pt.size=0)
dev.off()









#########################################
pbmc=COCHLEA 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_MAKER.rds')
LEVEL=c('0','1.1','1.2','2','4','5','6','7','8','9.1','9.2','10','11','12','13','14','15','16','17','18','19')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)
################

# IHC, 8
# OHC, 11


#####################
pbmc.sub=subset(pbmc, cells=colnames(pbmc)[which(Idents(pbmc) %in% c('8','11'))])
Idents(pbmc.sub)=factor(pbmc.sub$time,levels=c('3m','12m','24m'))


OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F25_COCHLEA_VLNPLOT_HC_TIME.pdf'
pdf(OUTPUT_PATH,width=10,height=10)
VlnPlot(pbmc.sub, features=c('Ap3b2','Kcnh8','Ankfn1'), ncol=2,pt.size=0)
dev.off()




#####################
pbmc.sub=subset(pbmc, cells=colnames(pbmc)[which(Idents(pbmc) %in% c('8'))])
Idents(pbmc.sub)=factor(pbmc.sub$time,levels=c('3m','12m','24m'))


OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F26_COCHLEA_VLNPLOT_IHC_TIME.pdf'
pdf(OUTPUT_PATH,width=10,height=13)
VlnPlot(pbmc.sub, features=c('Atp2a3','Chrm3','Dnm1','Dgkg','Slc17a8'), ncol=2,pt.size=0)
dev.off()




#####################
pbmc.sub=subset(pbmc, cells=colnames(pbmc)[which(Idents(pbmc) %in% c('11'))])
Idents(pbmc.sub)=factor(pbmc.sub$time,levels=c('3m','12m','24m'))


OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F27_COCHLEA_VLNPLOT_OHC_TIME.pdf'
pdf(OUTPUT_PATH,width=10,height=4)
VlnPlot(pbmc.sub, features=c('Dnm3','Slc26a5'), ncol=2,pt.size=0)
dev.off()










#######################################################################################################
# UTRICLE HC
#######################################################################################################



#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')

library(Seurat)
library(dplyr)
library(patchwork)

UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')


#########################################
pbmc=UTRICLE 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_MAKER.rds')
LEVEL=c('15','18','0','3.1','14','3.2','16.1','19','16.2','17','4','8','7','10','11','12','9','5','13','1','2')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)
################
# TYPE1.1，5
# TYPE1.2, 13
# TYPE2.1, 1
# TYPE2.2, 2
#####################
pbmc.sub=subset(pbmc, cells=colnames(pbmc)[which(Idents(pbmc) %in% c('5','13','1','2'))])

pbmc=pbmc.sub

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc=FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
all.genes=rownames(pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)
print('Calculating PCs ...')
pbmc <- RunPCA(object = pbmc, seed.use=123, npcs=50, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
pbmc <- RunUMAP(pbmc, dims = 1:30,seed.use = 123,n.components=2)

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F28_UTRICLE.HC_UMAP.pdf',width=10,height=10)
DimPlot(pbmc,label=TRUE,label.size=12)+NoLegend()
dev.off()
###########
saveRDS(pbmc, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.HC_SEURAT.rds')
###########

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(pbmc.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.HC_MAKER.rds')
write.table(pbmc.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.HC_MAKER.tsv',col.names=T,row.names=F,sep='\t',quote=F)
#################




#########################################
pbmc=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.HC_SEURAT.rds')
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.HC_MAKER.rds')
LEVEL=c('5','13','1','2')
###############
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

DATA=pbmc@assays$RNA@data
CLST=as.character(Idents(pbmc))


##################################
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
MEAN=.generate_mean(DATA,CLST)
MEAN.BAC=MEAN
colnames(MEAN)=LEVEL
i=1
while(i<=ncol(MEAN)){
  MEAN[,i]=MEAN.BAC[,which(colnames(MEAN.BAC)==colnames(MEAN)[i])]
  i=i+1
}

##################################
library(dplyr)
N=10
topN <- pbmc.markers %>% group_by(cluster) %>% top_n(n = N, wt = avg_log2FC)


MAT=matrix(0,ncol=ncol(MEAN),nrow=nrow(topN ))
i=1
while(i<=nrow(topN )){
  this_gene=topN$gene[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN$gene
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F29_UTRICLE.HC_HEAT.pdf',width=15,height=15)
color_fun_3=colorRamp2(c(-2,-1,0.5,0,0.5,1,2 ), c('grey95','grey95','grey90','grey90','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
        show_column_dend = FALSE, show_row_dend = FALSE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()


library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

this_plot=DoHeatmap(pbmc, features = topN$gene) + scale_fill_gradientn(colors = c("grey95", "grey90","red"))
pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F29_UTRICLE.HC_HEAT_CELL.pdf',width=15,height=15)
print(this_plot)
dev.off()

#############################################




##################################################################
pbmc.sub=pbmc
Idents(pbmc.sub)=factor(paste0(pbmc.sub$new.clst.sub,'_', pbmc.sub$time),levels=c('5_3m','5_12m','5_22m','13_3m','13_12m','13_22m',
                                                                                  '1_3m','1_12m','1_22m','2_3m','2_12m','2_22m'))


OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F30_UTRICLE.HC_VLNPLOT_TIME.pdf'
pdf(OUTPUT_PATH,width=5,height=13)
VlnPlot(pbmc.sub, features=c('Agbl1','Cdk15','Kcnj16','Kcnj6','Gpc5'), ncol=1,pt.size=0)
dev.off()

######################################################################








#########################################
pbmc=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.HC_SEURAT.rds')
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.HC_MAKER.rds')
LEVEL=c('5','13','1','2')
###############
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

DATA=pbmc@assays$RNA@data
CLST=as.character(Idents(pbmc))

###########################

GS1=c('Nalcn','Scn1a','Scn2a','Scn3a','Scn4a','Scn5a','Scn7a','Scn8a','Scn9a')
GS2=c('Hcn1','Kcna10','Kcnab1','Kcnb2','Kcnc2','Kcnd3','Kcnh7','Kcnh8','Kcnip4','Kcnj3','Kcnj6','Kcnj16','Kcnma1','Kcnq1','Kcnq1ot1','Kcnq3')
GS3=c('Atp1a1','Atp1a2','Atp1b1','Atp2a2','Atp2b1','Atp2b2','Atp2b3','Atp2b4','Atp2c1','Atp6v0a1','Atp8a1','Atp8a2','Atp11a','Atp11c','Atp13a3')
GS4=c('Camk2d','Camsap2','Camta1','Cacna1d','Cacna1e','Cacna2d1','Cacna2d3','Cacna2d4')
GS5=c('Slc1a3','Slc1a7','Slc17a8','Slc25a12')
GS6=c('Gria3','Grid1','Grid2','Grip1')


###############################

##################################
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
MEAN=.generate_mean(DATA,CLST)
MEAN.BAC=MEAN
colnames(MEAN)=LEVEL
i=1
while(i<=ncol(MEAN)){
  MEAN[,i]=MEAN.BAC[,which(colnames(MEAN.BAC)==colnames(MEAN)[i])]
  i=i+1
}

#################################
topN <- GS1


MAT=matrix(0,ncol=ncol(MEAN),nrow=length(topN ))
i=1
while(i<=length(topN )){
  this_gene=topN[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F31_UTRICLE.HC.GS1_HEAT.pdf',width=15,height=15)
color_fun_3=colorRamp2(c(-2,-1,0,1,2 ), c('grey95','grey95','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=TRUE,
        show_column_dend = FALSE, show_row_dend = TRUE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()
#############################################




#################################
topN <- GS2


MAT=matrix(0,ncol=ncol(MEAN),nrow=length(topN ))
i=1
while(i<=length(topN )){
  this_gene=topN[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F31_UTRICLE.HC.GS2_HEAT.pdf',width=15,height=15)
color_fun_3=colorRamp2(c(-2,-1,0,1,2 ), c('grey95','grey95','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=TRUE,
        show_column_dend = FALSE, show_row_dend = TRUE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()
#############################################





#################################
topN <- GS3


MAT=matrix(0,ncol=ncol(MEAN),nrow=length(topN ))
i=1
while(i<=length(topN )){
  this_gene=topN[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F31_UTRICLE.HC.GS3_HEAT.pdf',width=15,height=15)
color_fun_3=colorRamp2(c(-2,-1,0,1,2 ), c('grey95','grey95','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=TRUE,
        show_column_dend = FALSE, show_row_dend = TRUE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()
#############################################


#################################
topN <- GS4


MAT=matrix(0,ncol=ncol(MEAN),nrow=length(topN ))
i=1
while(i<=length(topN )){
  this_gene=topN[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F31_UTRICLE.HC.GS4_HEAT.pdf',width=15,height=15)
color_fun_3=colorRamp2(c(-2,-1,0,1,2 ), c('grey95','grey95','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=TRUE,
        show_column_dend = FALSE, show_row_dend = TRUE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()
#############################################


#################################
topN <- GS5


MAT=matrix(0,ncol=ncol(MEAN),nrow=length(topN ))
i=1
while(i<=length(topN )){
  this_gene=topN[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F31_UTRICLE.HC.GS5_HEAT.pdf',width=15,height=15)
color_fun_3=colorRamp2(c(-2,-1,0,1,2 ), c('grey95','grey95','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=TRUE,
        show_column_dend = FALSE, show_row_dend = TRUE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()
#############################################


#################################
topN <- GS6


MAT=matrix(0,ncol=ncol(MEAN),nrow=length(topN ))
i=1
while(i<=length(topN )){
  this_gene=topN[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F31_UTRICLE.HC.GS6_HEAT.pdf',width=15,height=15)
color_fun_3=colorRamp2(c(-2,-1,0,1,2 ), c('grey95','grey95','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=TRUE,
        show_column_dend = FALSE, show_row_dend = TRUE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()
#############################################

























#######################################################################################################
# UTRICLE supporting cell
#######################################################################################################

#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')

library(Seurat)
library(dplyr)
library(patchwork)

UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')


#########################################
pbmc=UTRICLE 
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_MAKER.rds')
LEVEL=c('15','18','0','3.1','14','3.2','16.1','19','16.2','17','4','8','7','10','11','12','9','5','13','1','2')
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)
################
# ESC，0
# SSC, 11
#####################
pbmc.sub=subset(pbmc, cells=colnames(pbmc)[which(Idents(pbmc) %in% c('0','11'))])

pbmc=pbmc.sub

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc=FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
all.genes=rownames(pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)
print('Calculating PCs ...')
pbmc <- RunPCA(object = pbmc, seed.use=123, npcs=50, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
pbmc <- RunUMAP(pbmc, dims = 1:30,seed.use = 123,n.components=2)

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F32_UTRICLE.SC_UMAP.pdf',width=6,height=6)
DimPlot(pbmc,label=TRUE,label.size=10)+NoLegend()
dev.off()
###########
saveRDS(pbmc, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.SC_SEURAT.rds')
###########

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(pbmc.markers, '/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.SC_MAKER.rds')
write.table(pbmc.markers, file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.SC_MAKER.tsv',col.names=T,row.names=F,sep='\t',quote=F)
#################




#########################################
pbmc=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.SC_SEURAT.rds')
pbmc.markers=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE.SC_MAKER.rds')
LEVEL=c('0','11')
###############
Idents(pbmc)=factor(pbmc$new.clst.sub, levels=LEVEL)

DATA=pbmc@assays$RNA@data
CLST=as.character(Idents(pbmc))


##################################
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
MEAN=.generate_mean(DATA,CLST)
MEAN.BAC=MEAN
colnames(MEAN)=LEVEL
i=1
while(i<=ncol(MEAN)){
  MEAN[,i]=MEAN.BAC[,which(colnames(MEAN.BAC)==colnames(MEAN)[i])]
  i=i+1
}

##################################
library(dplyr)
N=10
topN <- pbmc.markers %>% group_by(cluster) %>% top_n(n = N, wt = avg_log2FC)


MAT=matrix(0,ncol=ncol(MEAN),nrow=nrow(topN ))
i=1
while(i<=nrow(topN )){
  this_gene=topN$gene[i]
  MAT[i,]=MEAN[which(rownames(MEAN)==this_gene),]
  i=i+1
}

colnames(MAT)=colnames(MEAN)
rownames(MAT)=topN$gene
##################################


mat=t(apply(MAT,1,scale))
colnames(mat)=colnames(MAT)
rownames(mat)=rownames(MAT)

library('ComplexHeatmap')
library('circlize')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F33_UTRICLE.SC_HEAT.pdf',width=15,height=15)
color_fun_3=colorRamp2(c(-2,-1,0.5,0,0.5,1,2 ), c('grey95','grey95','grey90','grey90','grey90','indianred1','red1'))
Heatmap(mat,row_title='',name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
        show_column_dend = FALSE, show_row_dend = FALSE, 
        show_column_names=TRUE, show_row_names=TRUE,
        col=color_fun_3, border = TRUE
)
dev.off()


library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)


this_plot=DoHeatmap(pbmc, features = topN$gene) + scale_fill_gradientn(colors = c("grey95", "grey90","red"))

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F33_UTRICLE.SC_HEAT_CELL.pdf',width=10,height=10)
print(this_plot)
dev.off()



#############################################




##################################################################
pbmc.sub=pbmc
Idents(pbmc.sub)=factor(paste0(pbmc.sub$new.clst.sub,'_', pbmc.sub$time),levels=c('0_3m','0_12m','0_22m','11_3m','11_12m','11_22m'))


OUTPUT_PATH='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/figure/F34_UTRICLE.SC_VLNPLOT_TIME.pdf'
pdf(OUTPUT_PATH,width=5,height=10)
VlnPlot(pbmc.sub, features=c('Adam12','Bmp6','Isl1','Tectb'), ncol=1,pt.size=0)
dev.off()


GENE=c('Rbm25','Prpf4b','Luc7l3','Srrm1','Rpl13','Rps24','Rps8','Rpsa','Rps23')
pdf("pictures/vlnplot.pdf", width = 10, height = 7)
i=1
while (i <=  length(GENE)) {plot_cochlea1 <- VlnPlot(cochlea1, features = GENE[i], pt.size = 0)+NoLegend()
print(plot_cochlea1)
print(i)
i=i+1
}
dev.off()


#################################################################################################################
# Drawing Adjustment
#################################################################################################################

#Dotplot
this_pbmc= xxx
this_markers= xxx
DotPlot(this_pbmc, features = this_markers) +theme(axis.text.x = element_text(angle = 90))+theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())+scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+theme(panel.border = element_rect(color = "black"), panel.spacing = unit(1, "mm"))
##using complexheatmap
col_fun <- colorRamp2(c(-2, 0, 1,2),c('#330066','#336699','#66CC66','#FFCC33'))
p <- Dotppot()
df <- p$data

exp_mat<-df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()
head(exp_mat)

percent_mat<-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 
row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()
head(percent_mat)

cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
Heatmap(exp_mat,
        heatmap_legend_param=list(title="expression"),
        column_title = "clustered dotplot", 
        col=col_fun,
        cluster_rows = T,
        cluster_columns = F,
        rect_gp = gpar(type = "none"),
        cell_fun = cell_fun,
        row_names_gp = gpar(fontsize = 5),
        border = "black")

layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= pindex(percent_mat, i, j)/100 * unit(2, "mm"),
              gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}
lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                             gp = gpar(fill = "black")))
  ))
set.seed(123)
hp<- Heatmap(exp_mat,
             heatmap_legend_param=list(title="expression"),
             column_title = "clustered dotplot", 
             col=col_fun,
             rect_gp = gpar(type = "none"),
             layer_fun = layer_fun,
             row_names_gp = gpar(fontsize = 5),
             border = "black",
)
draw( hp, annotation_legend_list = lgd_list)


#Umap
FeaturePlot(this_pbmc, features = this_markers, reduction = "umap", label = T)+scale_color_gradientn(colours = c('#E0E0E0','#FFFBBB','#00DDDD','#005AB5','#000079'))+ theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +theme(panel.grid = element_blank(), axis.title = element_text(face = 2,hjust = 0.03))


#Volcano plot
EnhancedVolcano(this_pbmc,lab = this_markers, x="avg_log2FC", y= "p_val", selectLab = c(down, up),
                 xlab = bquote(~Log[2]~'fold change'),
                title = "OHC VS IHC",  
                 labSize = 3.0, 
                 labCol = 'black', 
                 labFace = 'bold', 
                 boxedLabels = TRUE, 
                 colAlpha = 1, 
                 legendPosition = 'right',
                 legendLabSize = 14, 
                 legendIconSize = 4.0, 
                 drawConnectors = TRUE, 
                 widthConnectors = 1.0, 
                 colConnectors = "black")

#Heatmap
.createMatrix <- function(pbmc){
    scale_data <- pbmc@assays[["RNA"]]@scale.data
    pbmc_matrix <- scale_data[rownames(scale_data) %in% RNA_splicing$`RNA splicing` , ]
    write.csv(pbmc_matrix, paste0(pbmc, "_RNA_splicing_matrix.csv"))
    return(pbmc)
}

Heatmap(utricle_matrix,
        name = "Gene Expression",
        col = color1,
        color_space = "LAB",
        cluster_rows = T,
        cluster_columns = F,
        show_row_names = TRUE,
        show_column_names = FALSE,
        column_title = "TypeII HC1",
        column_title_side = "top",
        row_names_gp = gpar(fontsize = 1),
        border = NA, 
        column_split = c(rep("3m", 732), rep("12m", 642), rep("22m", 467)),
        top_annotation= HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:6), labels = c("3m","12m","22m") )))











