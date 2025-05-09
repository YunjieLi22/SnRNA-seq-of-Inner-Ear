


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

###############################
# HC Type I vs. Type II

Idents(pbmc)=TYPE

this_ident1=c('Type.I.HC.1','Type.I.HC.2')
this_ident2=c('Type.II.HC.1','Type.II.HC.2')

this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)

this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)

#this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, test.use='t', 
#                      min.pct = 0.25, logfc.threshold = 0)

#this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, test.use='t',
#                     min.pct = 0.25, logfc.threshold = 0)


tmp=this_diff2
tmp[,3]=this_diff2[,4]
tmp[,4]=this_diff2[,3]
tmp[,2]=-this_diff2[,2]

this_diff=rbind(this_diff1,tmp)
this_diff=as.data.frame(this_diff)

plot(this_diff$avg_log2FC, -log10(this_diff$p_val_adj))

##########################################################################
#saveRDS(this_diff1,file='UTRICLE_20211025/DIFF.HC.TypeI.vs.TypeII_TypeI.rds')
#saveRDS(this_diff2,file='UTRICLE_20211025/DIFF.HC.TypeI.vs.TypeII_TypeII.rds')
#saveRDS(this_diff,file='UTRICLE_20211025/DIFF.HC.TypeI.vs.TypeII.rds')

this_diff$gene=rownames(this_diff)

write.table(this_diff,file='UTRICLE_20211025/DIFF.HC.TypeI.vs.TypeII.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

GENE1_INDEX=which(this_diff$avg_log2FC>0 & this_diff$p_val_adj==0)
GENE2_INDEX=which(this_diff$avg_log2FC<0 & this_diff$p_val_adj==0)
###############
GENE1_INDEX=GENE1_INDEX[order(-this_diff$avg_log2FC[GENE1_INDEX])[1:10]]
GENE2_INDEX=GENE2_INDEX[order(this_diff$avg_log2FC[GENE2_INDEX])[1:10]]
###############
GENE1_NAME=rownames(this_diff)[GENE1_INDEX]
GENE2_NAME=rownames(this_diff)[GENE2_INDEX]

########################################################
adjPvalue=this_diff$p_val_adj + min(this_diff$p_val_adj[which(this_diff$p_val_adj>0)])/10
log2FC=this_diff$avg_log2FC





####################################
plot(log2FC, -log10(adjPvalue),
     xlim=c(-4,4),ylim=c(0,300), pch=16, col='grey80',cex=0.5)
X=log2FC
Y=-log10(adjPvalue)

##############################
X1=X[GENE1_INDEX]
Y1=Y[GENE1_INDEX]
T1=GENE1_NAME

#########
OOO=order(X1)
X1=X1[OOO]
Y1=Y1[OOO]
T1=T1[OOO]

##########
points(X1,Y1,pch=16,col='red')
#####################
X11=X1 +(c(1:length(X1))/length(X1) - 0.5) *1
Y11=Y1-( 1 - c(1:length(X1))/length(X1) ) * 200
segments(x0=X1, y0=Y1, x1=X11, y1=Y11)
text(x=X11,y=Y11,labels=T1,srt=0,pos=c(4))



##############################
X2=X[GENE2_INDEX]
Y2=Y[GENE2_INDEX]
T2=GENE2_NAME

#########
OOO=order(X2)
X2=X2[OOO]
Y2=Y2[OOO]
T2=T2[OOO]

##########
points(X2,Y2,pch=16,col='red')
#####################
X22=X2 +(c(1:length(X2))/length(X2) - 0.5) *1
Y22=Y2-( c(1:length(X2))/length(X2) ) * 200
########
X22[1]=X2[1]+0.1
X22[2]=X2[2]-0.1
#####
segments(x0=X2, y0=Y2, x1=X22, y1=Y22)
text(x=X22,y=Y22,labels=T2,srt=0,pos=c(2))

#######################################################






#########################


pdf('UTRICLE_20211025/DIFF.HC.TypeI.vs.TypeII.pdf',width=7,height=6)

#########################


####################################
plot(log2FC, -log10(adjPvalue),
     xlim=c(-4,4),ylim=c(0,300), pch=16, col='grey80',cex=0.5)
X=log2FC
Y=-log10(adjPvalue)

##############################
X1=X[GENE1_INDEX]
Y1=Y[GENE1_INDEX]
T1=GENE1_NAME

#########
OOO=order(X1)
X1=X1[OOO]
Y1=Y1[OOO]
T1=T1[OOO]

##########
points(X1,Y1,pch=16,col='red')
#####################
X11=X1 +(c(1:length(X1))/length(X1) - 0.5) *1
Y11=Y1-( 1 - c(1:length(X1))/length(X1) ) * 200
segments(x0=X1, y0=Y1, x1=X11, y1=Y11)
text(x=X11,y=Y11,labels=T1,srt=0,pos=c(4))



##############################
X2=X[GENE2_INDEX]
Y2=Y[GENE2_INDEX]
T2=GENE2_NAME

#########
OOO=order(X2)
X2=X2[OOO]
Y2=Y2[OOO]
T2=T2[OOO]

##########
points(X2,Y2,pch=16,col='blue')
#####################
X22=X2 +(c(1:length(X2))/length(X2) - 0.5) *1
Y22=Y2-( c(1:length(X2))/length(X2) ) * 200
########
X22[1]=X2[1]+0.2
X22[2]=X2[2]
#####
segments(x0=X2, y0=Y2, x1=X22, y1=Y22)
text(x=X22,y=Y22,labels=T2,srt=0,pos=c(2))

#########################

dev.off()

#########################


#########################
# GO Enrichment

ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(pbmc), 'ENTREZID', 'SYMBOL')

###################

    this_sig_gene=GENE1_NAME
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_20211025/DIFF.HC.TypeI.vs.TypeII_TypeI_GO.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

    this_out1=this_out


###############################


    this_sig_gene=GENE2_NAME
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_20211025/DIFF.HC.TypeI.vs.TypeII_TypeII_GO.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

    this_out2=this_out

####################################

this_out=this_out1

par(las=2)
par(mar=c(35,5,5,5))
#adjPvalue=this_out$p.adjust
Pvalue=this_out$pvalue
NAME=paste0(this_out$ID, ', ',this_out$ONTOLOGY,', ', this_out$Description)
names(Pvalue)=NAME

Pvalue=sort(Pvalue)

barplot(-log(Pvalue[1:20]),ylab='-log10(Pvalue)')

pdf('UTRICLE_20211025/DIFF.HC.TypeI.vs.TypeII_TypeI_GO.pdf',height=10,width=7)
par(las=2)
par(mar=c(35,5,5,5))
barplot(-log(Pvalue[1:20]),ylab='-log10(Pvalue)')
dev.off()


####################################

this_out=this_out2

par(las=2)
par(mar=c(45,5,5,5))
#adjPvalue=this_out$p.adjust
Pvalue=this_out$pvalue
NAME=paste0(this_out$ID, ', ',this_out$ONTOLOGY,', ', this_out$Description)
names(Pvalue)=NAME

Pvalue=sort(Pvalue)

barplot(-log(Pvalue[1:20]),ylab='-log10(Pvalue)')

pdf('UTRICLE_20211025/DIFF.HC.TypeI.vs.TypeII_TypeII_GO.pdf',height=11,width=7)
par(las=2)
par(mar=c(45,5,5,5))
barplot(-log(Pvalue[1:20]),ylab='-log10(Pvalue)')
dev.off()

####################################



###############################################

# HC 4 subtypes

##############################################



###########################
Idents(pbmc)=pbmc$final.clst.sub.label

this_pbmc=subset(pbmc,idents=c('Type.I.HC.1','Type.I.HC.2','Type.II.HC.1','Type.II.HC.2'))

this_pbmc= ScaleData(this_pbmc, features = rownames(this_pbmc))


Idents(this_pbmc)=factor(Idents(this_pbmc), levels = c('Type.I.HC.1','Type.I.HC.2','Type.II.HC.1','Type.II.HC.2'))


#######################
MK=FindAllMarkers(this_pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MK=as.data.frame(MK)

#####################
TOP <- MK %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
##############

this_plot=DoHeatmap(this_pbmc, features = TOP$gene) + NoLegend()



print(this_plot)


write.table(MK, paste0('UTRICLE_20211025/4HC_Marker.tsv'),sep='\t',quote=F,col.names=T,row.names=F)



pdf('UTRICLE_20211025/4HC_Marker_Heatmap.pdf',width=7,height=7)
print(this_plot)
dev.off()


##################################################################################

GS1="Adgrl3/Erbb4/Gpc6/Il1rapl1/Lrrc4c/Palld/Unc13b/Ephb1/Kirrel3/Sipa1l1/Farp1/Efna5/Cadm1/Lmx1a/Col4a5/Malat1/Epha7/Tiam1/Gpm6a/Prkca/Actn1/Clstn2/Myo5b/Nfia/Pdlim5/Fgf13/Ppfibp1/Zdhhc2/Utrn/Adgrl2/Filip1/Snta1/Ndrg2/Nedd4/Efnb2/Epha5/Cttnbp2/Ank3/Cyfip1/Sorbs2/Adam10/Nrp2/Sparc/Wasf2/Sparcl1/Cntn5/Glrb/Rapgef4/Ptn/St8sia2/Igf1r/App/Sdk1/Nrcam/Fgfr2/Grid2"
GS2="Nfib/Robo1/Ephb1/Efna5/Bmpr1b/Lmx1a/Prtg/Epha7/Chl1/Boc/Tgfb2/Notch1/Apbb2/Ext1/Alcam/Efnb2/Cntn4/Epha5/Ank3/Cyfip1/Lrp1/Slit2/Nrp2/Sema6a/Enah/Cntn5/Arhgap35/Zswim6/Kif5b/App/Nrcam/Lama2"
GS3="Rab3ip/Syne1/Armc4/Cfap43/Spag16/Ttll3/Kif27/Ttc29/Pkd2/Drc1/Spef2/Dync2h1/Cfap65/Spag6l/Rfx3/Ift74/Cfap69/Unc119b/Ttc21a/Lrguk/Ahi1/Eno4/Ccdc65/Ablim3/Armc2/Cyld/Rpgr/Ccdc40/Kif3a/Bbs9"
GS4="Atp2b2/Kalrn/Kcnma1/Ankfn1/Dnm1/Espn/Myo15/Apba1/Oxr1/Fign/Agtpbp1/App/Mapk10/Eps8/Btbd9/Tnr"
GS5="Pcsk5/Armc4/Pkd2/Drc1/Dync2h1/Rfx3/Ift74/Ahi1/Ccdc40/Kif3a"
GS6="Shank2/Kalrn/Shisa6/Nlgn1/Plcb1/Grid2ip"
GS7="Bdnf/Nlgn1/Dnm1/Gabrb3/Car7"
GS8="Sorbs1/Chrm3/Plcb1/Ly6h/Atp2b4"

GS1=strsplit(GS1, '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
GS2=strsplit(GS2, '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
GS3=strsplit(GS3, '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
GS4=strsplit(GS4, '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
GS5=strsplit(GS5, '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
GS6=strsplit(GS6, '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
GS7=strsplit(GS7, '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
GS8=strsplit(GS8, '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]



pdf('UTRICLE_20211025/4HC_Selected_GO.pdf',width=7,height=7)

DoHeatmap(this_pbmc, features = GS1) + NoLegend()
DoHeatmap(this_pbmc, features = GS2) + NoLegend()
DoHeatmap(this_pbmc, features = GS3) + NoLegend()
DoHeatmap(this_pbmc, features = GS4) + NoLegend()
DoHeatmap(this_pbmc, features = GS5) + NoLegend()
DoHeatmap(this_pbmc, features = GS6) + NoLegend()
DoHeatmap(this_pbmc, features = GS7) + NoLegend()
DoHeatmap(this_pbmc, features = GS8) + NoLegend()
dev.off()





Idents(this_pbmc)=factor(this_pbmc$final.clst.sub.label.time, levels = c('utricle.3m.Type.I.HC.1','utricle.12m.Type.I.HC.1','utricle.22m.Type.I.HC.1',
                                                                         'utricle.3m.Type.I.HC.2','utricle.12m.Type.I.HC.2','utricle.22m.Type.I.HC.2',
                                                                         'utricle.3m.Type.II.HC.1','utricle.12m.Type.II.HC.1','utricle.22m.Type.II.HC.1',
                                                                         'utricle.3m.Type.II.HC.2','utricle.12m.Type.II.HC.2','utricle.22m.Type.II.HC.2'
                                                                         ))


GS1=c('Nalcn','Scn1a','Scn2a','Scn3a','Scn4a','Scn5a','Scn6a','Scn7a','Scn8a','Scn9a')
GS2=c('Hcn1','Kcna10','Kcnab1','Kcnb2','Kcnc2','Kcnd3','Kcnh7','Kcnh8','Kcnip4','Kcnj3','Kcnj6','Kcnj16','Kcnma1','Kcnq1','Kcnq1ot1','Kcnq3')
GS3=c('Atp1a1','Atp1a2','Atp1b1','Atp2a2','Atp2b1','Atp2b2','Atp2b3','Atp2b4','Atp2c1','Atp6v0a1','Atp8a1','Atp8a2','Atp11a','Atp11c','Atp13a3')
GS4=c('Camk2d','Camsap2','Camta1','Cacna1d','Cacna1e','Cacna2d1','Cacna2d3','Cacna2d4')
GS5=c('Slc1a3','Slc1a7','Slc17a8','Slc25a12')
GS6=c('Gria3','Grid1','Grid2','Grip1')
GS7=c('Chrm3','Chrna10')



pdf('UTRICLE_20211025/4HC_Selected_Genes.pdf',width=7,height=7)

DoHeatmap(this_pbmc, features = GS1) + NoLegend() #+scale_fill_gradientn(colors = c("grey80","grey80", 'red1'))
DoHeatmap(this_pbmc, features = GS2) + NoLegend() #+scale_fill_gradientn(colors = c("grey80","grey80", 'red1'))
DoHeatmap(this_pbmc, features = GS3) + NoLegend() #+scale_fill_gradientn(colors = c("grey80","grey80", 'red1'))
DoHeatmap(this_pbmc, features = GS4) + NoLegend() #+scale_fill_gradientn(colors = c("grey80","grey80", 'red1'))
DoHeatmap(this_pbmc, features = GS5) + NoLegend() #+scale_fill_gradientn(colors = c("grey80","grey80", 'red1'))
DoHeatmap(this_pbmc, features = GS6) + NoLegend() #+scale_fill_gradientn(colors = c("grey80","grey80", 'red1'))
DoHeatmap(this_pbmc, features = GS7) + NoLegend() #+scale_fill_gradientn(colors = c("grey80","grey80", 'red1'))
dev.off()


##########################################################









###########################
Idents(pbmc)=pbmc$final.clst.sub.label

this_pbmc=subset(pbmc,idents=c('Extrastriolar.SC','Striolar.SC'))

this_pbmc= ScaleData(this_pbmc, features = rownames(this_pbmc))

Idents(this_pbmc)=factor(Idents(this_pbmc), levels = c('Extrastriolar.SC','Striolar.SC'))




#######################
MK=FindAllMarkers(this_pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MK=as.data.frame(MK)

#####################
TOP <- MK %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
##############

this_plot=DoHeatmap(this_pbmc, features = TOP$gene) + NoLegend()



print(this_plot)


write.table(MK, paste0('UTRICLE_20211025/2SC_Marker.tsv'),sep='\t',quote=F,col.names=T,row.names=F)



pdf('UTRICLE_20211025/2SC_Marker_Heatmap.pdf',width=7,height=10)
print(this_plot)
dev.off()


GS=c('Tectb','Isl1','Tspan8','Cyp26b1','Prkar2b','Dab1','Dync1i1','Bmp6','Adam12')

VlnPlot(pbmc,features=GS,ncol=3,pt.size=0 )



pdf('UTRICLE_20211025/2SC_Selected_Genes_Vln.pdf',width=14,height=12)
VlnPlot(pbmc,features=GS,ncol=3,pt.size=0 )
dev.off()



pdf('UTRICLE_20211025/2SC_Selected_Genes_UMAP.pdf',width=14,height=12)
FeaturePlot(pbmc,features=GS,ncol=3 )
dev.off()















###########################
Idents(pbmc)=pbmc$final.clst.sub.label.time

this_pbmc=subset(pbmc,idents=c('utricle.3m.Extrastriolar.SC','utricle.3m.Striolar.SC'))

this_pbmc= ScaleData(this_pbmc, features = rownames(this_pbmc))

Idents(this_pbmc)=factor(Idents(this_pbmc), levels = c('utricle.3m.Extrastriolar.SC','utricle.3m.Striolar.SC'))




#######################
MK=FindAllMarkers(this_pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MK=as.data.frame(MK)

#####################
TOP <- MK %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
##############

this_plot=DoHeatmap(this_pbmc, features = TOP$gene) + NoLegend()



print(this_plot)


write.table(MK, paste0('UTRICLE_20211025/2SC_3M_Marker.tsv'),sep='\t',quote=F,col.names=T,row.names=F)



pdf('UTRICLE_20211025/2SC_3M_Marker_Heatmap.pdf',width=7,height=10)
print(this_plot)
dev.off()










###########################
Idents(pbmc)=pbmc$final.clst.sub.label.time

this_pbmc=subset(pbmc,idents=c('utricle.22m.Extrastriolar.SC','utricle.12m.Extrastriolar.SC','utricle.3m.Extrastriolar.SC',
                               'utricle.3m.Striolar.SC','utricle.12m.Striolar.SC','utricle.22m.Striolar.SC'))

this_pbmc= ScaleData(this_pbmc, features = rownames(this_pbmc))

Idents(this_pbmc)=factor(Idents(this_pbmc), levels = c('utricle.22m.Extrastriolar.SC','utricle.12m.Extrastriolar.SC','utricle.3m.Extrastriolar.SC',
                               'utricle.3m.Striolar.SC','utricle.12m.Striolar.SC','utricle.22m.Striolar.SC'))




#######################
MK=FindAllMarkers(this_pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MK=as.data.frame(MK)

#####################
TOP <- MK %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
##############

this_plot=DoHeatmap(this_pbmc, features = TOP$gene) + NoLegend()



print(this_plot)


write.table(MK, paste0('UTRICLE_20211025/2SC_3M.12M.22M_Marker.tsv'),sep='\t',quote=F,col.names=T,row.names=F)



pdf('UTRICLE_20211025/2SC_3M.12M.22M_Marker_Heatmap.pdf',width=7,height=10)
print(this_plot)
dev.off()

################################################





