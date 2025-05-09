
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

#####################################################################
setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging')


pbmc=readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step16_combineAndMacrophage/COM_SEURAT.rds')

pbmc$type1=pbmc$type
pbmc$type1[which(pbmc$type=='0')]='TBC.1'
pbmc$type1[which(pbmc$type=='1.1')]='CC'
pbmc$type1[which(pbmc$type=='1.2')]='DC'
pbmc$type1[which(pbmc$type=='2')]='TBC.1'
pbmc$type1[which(pbmc$type=='4')]='SGC'
pbmc$type1[which(pbmc$type=='6')]='OSC'
pbmc$type1[which(pbmc$type=='7')]='Dentate.interdental.cell'
pbmc$type1[which(pbmc$type=='8')]='IHC'
pbmc$type1[which(pbmc$type=='10')]='OPC'
pbmc$type1[which(pbmc$type=='11')]='OHC'
pbmc$type1[which(pbmc$type=='12')]='SGN'
pbmc$type1[which(pbmc$type=='14')]='Macrophage'
pbmc$type1[which(pbmc$type=='16')]='IPC'
pbmc$type1[which(pbmc$type=='17')]='Astrocytes'
pbmc$type1[which(pbmc$type=='18')]='Melanocyte'
pbmc$type1[which(pbmc$type=='19')]='Endothelial.cell'


##########################################
TYPE=pbmc$type1

############################################
Idents(pbmc)=TYPE
#################################################


############################################################
# Cochlea HC v.s.Utricle HC
###########################################################

this_ident1=c('IHC','OHC')
this_ident2=c('Type.I.HC.1','Type.I.HC.2','Type.II.HC.1','Type.II.HC.2')

this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)

this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)

tmp=this_diff2
tmp[,3]=this_diff2[,4]
tmp[,4]=this_diff2[,3]
tmp[,2]=-this_diff2[,2]

this_diff=rbind(this_diff1,tmp)
this_diff=as.data.frame(this_diff)


plot(this_diff$avg_log2FC, -log10(this_diff$p_val_adj))


this_diff$gene=rownames(this_diff)

write.table(this_diff,file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF.Cochlea.HC.v.s.Utricle.HC.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)




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
     xlim=c(-6.5,6.5),ylim=c(0,320), pch=16, col='grey80',cex=0.5)
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
X22[1]=X2[1]+0.1
X22[2]=X2[2]-0.1
#####
segments(x0=X2, y0=Y2, x1=X22, y1=Y22)
text(x=X22,y=Y22,labels=T2,srt=0,pos=c(2))

#######################################################


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF.Cochlea.HC.v.s.Utricle.HC.pdf',width=7,height=6)


####################################
plot(log2FC, -log10(adjPvalue),
     xlim=c(-6.5,6.5),ylim=c(0,320), pch=16, col='grey80',cex=0.5)
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
X22[1]=X2[1]+0.1
X22[2]=X2[2]-0.1
#####
segments(x0=X2, y0=Y2, x1=X22, y1=Y22)
text(x=X22,y=Y22,labels=T2,srt=0,pos=c(2))

#######################################################
dev.off()



#####################
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
    write.table(this_out, paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF_Cochlea.HC.HighTop10_GO.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

    this_out1=this_out


###############################


    this_sig_gene=GENE2_NAME
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF_Utricle.HC.HighTop10_GO.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

    this_out2=this_out




####################################

this_out=this_out1

Pvalue=this_out$pvalue
NAME=paste0(this_out$ID, ', ',this_out$ONTOLOGY,', ', this_out$Description)
names(Pvalue)=NAME
Pvalue=sort(Pvalue)


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF_Cochlea.HC.HighTop10_GO.pdf',height=10,width=7)
par(las=2)
par(mar=c(40,5,5,5))
barplot(-log(Pvalue[1:20]),ylab='-log10(Pvalue)')
dev.off()


####################################

this_out=this_out2

Pvalue=this_out$pvalue
NAME=paste0(this_out$ID, ', ',this_out$ONTOLOGY,', ', this_out$Description)
names(Pvalue)=NAME
Pvalue=sort(Pvalue)


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF_Utricle.HC.HighTop10_GO.pdf',height=10,width=7)
par(las=2)
par(mar=c(25,5,5,5))
barplot(-log(Pvalue[1:20]),ylab='-log10(Pvalue)')
dev.off()

####################################






############################################################
# Cochlea IHC v.s. Cochlea OHC
###########################################################

this_ident1=c('IHC')
this_ident2=c('OHC')

this_diff1=FindMarkers(pbmc, ident.1=this_ident1,ident.2=this_ident2,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)

this_diff2=FindMarkers(pbmc, ident.1=this_ident2,ident.2=this_ident1,only.pos = TRUE, 
                      min.pct = 0.25, logfc.threshold = 0)

tmp=this_diff2
tmp[,3]=this_diff2[,4]
tmp[,4]=this_diff2[,3]
tmp[,2]=-this_diff2[,2]

this_diff=rbind(this_diff1,tmp)
this_diff=as.data.frame(this_diff)


plot(this_diff$avg_log2FC, -log10(this_diff$p_val_adj))


this_diff$gene=rownames(this_diff)

write.table(this_diff,file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF.Cochlea.IHC.v.s.Cochlea.OHC.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)




#GENE1_INDEX=which(this_diff$avg_log2FC>0 & this_diff$p_val_adj==0)
#GENE2_INDEX=which(this_diff$avg_log2FC<0 & this_diff$p_val_adj==0)

GENE1_INDEX=which(this_diff$avg_log2FC>0)
GENE2_INDEX=which(this_diff$avg_log2FC<0 )


###############
GENE1_INDEX=GENE1_INDEX[order(this_diff$p_val_adj[GENE1_INDEX])[1:10]]
GENE2_INDEX=GENE2_INDEX[order(this_diff$p_val_adj[GENE2_INDEX])[1:10]]
###############
GENE1_NAME=rownames(this_diff)[GENE1_INDEX]
GENE2_NAME=rownames(this_diff)[GENE2_INDEX]

########################################################
adjPvalue=this_diff$p_val_adj + min(this_diff$p_val_adj[which(this_diff$p_val_adj>0)])/10
log2FC=this_diff$avg_log2FC

####################################
plot(log2FC, -log10(adjPvalue),
     xlim=c(-6.5,6.5),ylim=c(0,80), pch=16, col='grey80',cex=0.5)
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
Y11=Y1-( 1 - c(1:length(X1))/length(X1) ) * 30
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
Y22=Y2-( c(1:length(X2))/length(X2) ) * 20
########
X22[1]=X2[1]+0.1
X22[2]=X2[2]-0.1
#####
segments(x0=X2, y0=Y2, x1=X22, y1=Y22)
text(x=X22,y=Y22,labels=T2,srt=0,pos=c(2))

#######################################################


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF.Cochlea.IHC.v.s.Cochlea.OHC.pdf',width=8,height=7)


####################################
plot(log2FC, -log10(adjPvalue),
     xlim=c(-6.5,6.5),ylim=c(0,80), pch=16, col='grey80',cex=0.5)
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
Y11=Y1-( 1 - c(1:length(X1))/length(X1) ) * 30
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
Y22=Y2-( c(1:length(X2))/length(X2) ) * 20
########
X22[1]=X2[1]+0.1
X22[2]=X2[2]-0.1
#####
segments(x0=X2, y0=Y2, x1=X22, y1=Y22)
text(x=X22,y=Y22,labels=T2,srt=0,pos=c(2))


#######################################################
dev.off()



#####################
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
    write.table(this_out, paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF_Cochlea.IHC.HighTop10_GO.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

    this_out1=this_out


###############################


    this_sig_gene=GENE2_NAME
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 1,
                   qvalueCutoff = 1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF_Cochlea.OHC.HighTop10_GO.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

    this_out2=this_out




####################################

this_out=this_out1

Pvalue=this_out$pvalue
NAME=paste0(this_out$ID, ', ',this_out$ONTOLOGY,', ', this_out$Description)
names(Pvalue)=NAME
Pvalue=sort(Pvalue)


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF_Cochlea.IHC.HighTop10_GO.pdf',height=10,width=7)
par(las=2)
par(mar=c(40,5,5,5))
barplot(-log(Pvalue[1:20]),ylab='-log10(Pvalue)')
dev.off()


####################################

this_out=this_out2

Pvalue=this_out$pvalue
NAME=paste0(this_out$ID, ', ',this_out$ONTOLOGY,', ', this_out$Description)
names(Pvalue)=NAME
Pvalue=sort(Pvalue)


pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/step17_HC_Aging/DIFF_Cochlea.OHC.HighTop10_GO.pdf',height=10,width=7)
par(las=2)
par(mar=c(25,5,5,5))
barplot(-log(Pvalue[1:20]),ylab='-log10(Pvalue)')
dev.off()

####################################



