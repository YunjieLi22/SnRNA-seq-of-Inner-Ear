

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
library(clusterProfiler)

library(org.Mm.eg.db)

#####################################################################

pbmc=readRDS( 'UTRICLE_SEURAT.rds')
PCUSE=readRDS( 'UTRICLE_PCUSE.rds')

pbmc <- FindClusters(pbmc, resolution = 0.5)
DimPlot(pbmc, reduction='umap', pt.size=0.1,label=T,label.size=8)  +NoLegend()

pbmc@meta.data$new.clst=Idents(pbmc)



NEW.CLST=as.character(Idents(pbmc))
NEW.CLST.SUB=NEW.CLST

VEC=pbmc@reductions$umap@cell.embeddings

###################################
USED=which(NEW.CLST==6)

USED_VEC=VEC[USED,]

    plot(USED_VEC, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(USED_VEC, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

NEW.CLST.SUB[USED]='6.1'
NEW.CLST.SUB[which(rownames(VEC) %in% SELECT)]='6.2'
#####################################################

###################################
USED=which(NEW.CLST==16)

USED_VEC=VEC[USED,]

    plot(USED_VEC, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(USED_VEC, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

NEW.CLST.SUB[USED]='16.1'
NEW.CLST.SUB[which(rownames(VEC) %in% SELECT)]='16.2'
#####################################################

###################################
USED=which(NEW.CLST==3)

USED_VEC=VEC[USED,]

    plot(USED_VEC, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(USED_VEC, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

NEW.CLST.SUB[USED]='3.1'
NEW.CLST.SUB[which(rownames(VEC) %in% SELECT)]='3.2'
#####################################################



pbmc$new.clst.sub=NEW.CLST.SUB
DimPlot(pbmc, reduction='umap', group.by='new.clst.sub', pt.size=0.1,label=T,label.size=12) +NoLegend()




#saveRDS(pbmc, 'UTRICLE_SEURAT_NEW.rds')


pdf('UTRICLE_UMAP.pdf',width=12,height=12)
DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 
DimPlot(pbmc, reduction='umap', group.by='ori.clst', pt.size=0.1,label=T,label.size=12) +NoLegend()
DimPlot(pbmc, reduction='umap', group.by='new.clst', pt.size=0.1,label=T,label.size=12) +NoLegend()
DimPlot(pbmc, reduction='umap', group.by='new.clst.sub', pt.size=0.1,label=T,label.size=12) +NoLegend()
dev.off()




NEW.CLST.SUB.LABEL=as.character(NEW.CLST.SUB)


NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='0')]='Utricle.SC'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='11')]='Utricle.SC'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='1')]='Type.II.HC.1'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='2')]='Type.II.HC.2'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='4')]='Nonsensory.epithelial.cell.1'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='6.1')]='Camk1d.Cmss1.cells'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='5')]='Type.I.HC.1'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='3.1')]='Fibrocyte.1'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='7')]='Ntn1.epithelial.cell'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='8')]='Nonsensory.epithelial.cell.2'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='14')]='Fibrocyte.2'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='9')]='Transitional.cell'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='10')]='Roof.epithelial.cell '
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='17')]='Neuron'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='15')]='Astrocytes'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='3.2')]='Fibrocyte.3'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='16.1')]='Macrophage'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='13')]='Type.I.HC.2'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='18')]='Endothelial.cell'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='6.2')]='Lymphocyte'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='16.2')]='Microglia'
NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='19')]='Melanocyte'

NEW.CLST.SUB.LABEL[which(NEW.CLST.SUB=='12')]='UnDefined'



pbmc$new.clst.sub.label=NEW.CLST.SUB.LABEL

pdf('UTRICLE_UMAP_LABEL.pdf',width=12,height=12)
DimPlot(pbmc, reduction='umap', group.by='new.clst.sub.label', 
         pt.size=0.1,label=T,label.size=5) +NoLegend()
dev.off()

#saveRDS(pbmc, 'UTRICLE_SEURAT_NEW.rds')

#########################################################

Idents(pbmc)=pbmc$new.clst.sub

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pbmc.markers, 'UTRICLE_MAKER_NEW.rds')
write.table(pbmc.markers, file='UTRICLE_MAKER_NEW.tsv',col.names=T,row.names=F,sep='\t',quote=F)


##########################################

all.gene=rownames(pbmc)
pbmc=ScaleData(object = pbmc, features =all.gene)


################################################
MK1=c('Nckap5','Col11a1','Chl1')
MK2=c('Anxa4','Thsd7b','Elfn1')
MK3=c('Mapt','Myo15b','Ush2a')
MK4=c('Nell1','Nox3','Slc26a4')
MK5=c('Camk1d','Cmss1','Otos')
MK6=c('Spp1','Mgat4c','Agbl1')
MK7=c('Lama1','Cacna1e','Map2k6')
MK8=c('Frem2','Ntn1','Cndp1')
MK9=c('Kcnj3','Unc5cl','Il17rb')
MK10=c('Actn2','Mapk4','Nrp1')
MK11=c('Rnf149','Kdr','Tbx3')
MK12=c('Esrrb','Kcne1','Rspo3')
MK13=c('Nrg1','Nrg3','Csmd2')
MK14=c('Drp2','Prx','Plp1')
#MK15=c('Cacna1e','Ucma','Cntn6')
MK15=c('Ucma','Cntn6')
MK16=c('Ptprc','Ly86','P2ry6')
MK17=c('Kcnj16','Tnc','Wdr49')
MK18=c('Cyyr1','Abcb1a','Adgrl4')
MK19=c('Xlr5a','Mc1r','Pde6a')
MK20=c('C1qc','Cd52','Tyrobp')
MK21=c('Dct','Tyr','Gpnmb')

#MKS=c(MK1,MK2,MK3,MK4,MK5,MK6,MK7,MK8,MK9,MK10,MK11,MK12,MK13,MK14,MK15,MK16,MK17,MK18,MK19,MK20,MK21)

MKS=c(MK1,MK3,MK2,MK17,MK6,MK11,MK12,MK8,MK9,MK4,MK13,MK20,MK21,MK16,MK19,MK15,MK10,MK7,MK18,MK5,MK14)


Idents(pbmc)=pbmc$new.clst.sub.label

p=DotPlot(pbmc, features=MKS)+ theme(axis.text.x = element_text(angle = 90))

custom_order <- sort(names(table(pbmc$new.clst.sub.label)))
p$data$id <- factor(p$data$id, levels = custom_order)

print(p)


pdf('UTRICLE_UMAP_LABEL_DOT.pdf',width=16,height=10)

print(p)

dev.off()


#saveRDS(pbmc, 'UTRICLE_SEURAT_NEW.rds')

#############################################################


UNIQ_CLST=unique(pbmc$new.clst.sub)

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
    write.table(this_out, paste0('UTRICLE_GO/',this_clst,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)
    
    #############################               
    print(paste0(i,' / ',length(UNIQ_CLST)))
    print(this_clst)
    i=i+1
    }



#############################################################
#############################################################



MKS=c('Anxa4','Mapt','Tnc','Spp1','Thsd7b','Elfn1')


FeaturePlot(pbmc,features=MKS,ncol=3 )
VlnPlot(pbmc,features=MKS,ncol=3,pt.size=0 )



pdf('UTRICLE_UMAP_LABEL_FEATURE1.pdf',width=16,height=10)
FeaturePlot(pbmc,features=MKS,ncol=3 )
VlnPlot(pbmc,features=MKS,ncol=3,pt.size=0 )
dev.off()

############################################################



Idents(pbmc)=pbmc$new.clst.sub.label

Type.I.over.Type.II=FindMarkers(pbmc, ident.1=c('Type.I.HC.1','Type.I.HC.2'),ident.2=c('Type.II.HC.1','Type.II.HC.2'),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(Type.I.over.Type.II)
Type.I.over.Type.II$gene=gene

Type.II.over.Type.I=FindMarkers(pbmc, ident.1=c('Type.II.HC.1','Type.II.HC.2'),ident.2=c('Type.I.HC.1','Type.I.HC.2'),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(Type.II.over.Type.I)
Type.II.over.Type.I$gene=gene



write.table(Type.I.over.Type.II, file='UTRICLE_HC/Type.I.over.Type.II.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(Type.II.over.Type.I, file='UTRICLE_HC/Type.II.over.Type.I.tsv',col.names=T,row.names=F,sep='\t',quote=F)


############################
    this_sig_gene=rownames(Type.I.over.Type.II)[which(Type.I.over.Type.II$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC/Type.I.over.Type.II_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


############################
    this_sig_gene=rownames(Type.II.over.Type.I)[which(Type.II.over.Type.I$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC/Type.II.over.Type.I_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

##################################




MK.RESULT=FindMarkers(pbmc, ident.1=c('Type.I.HC.1'),ident.2=c('Type.II.HC.1','Type.II.HC.2','Type.I.HC.2'),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC/Type.I.HC.1.over.OtherThree.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC/Type.I.HC.1.over.OtherThree_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


##################################


MK.RESULT=FindMarkers(pbmc, ident.1=c('Type.I.HC.2'),ident.2=c('Type.II.HC.1','Type.II.HC.2','Type.I.HC.1'),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC/Type.I.HC.2.over.OtherThree.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC/Type.I.HC.2.over.OtherThree_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


##################################


MK.RESULT=FindMarkers(pbmc, ident.1=c('Type.II.HC.1'),ident.2=c('Type.II.HC.2','Type.I.HC.1','Type.I.HC.2'),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC/Type.II.HC.1.over.OtherThree.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC/Type.II.HC.1.over.OtherThree_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


##################################



MK.RESULT=FindMarkers(pbmc, ident.1=c('Type.II.HC.2'),ident.2=c('Type.II.HC.1','Type.I.HC.1','Type.I.HC.2'),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC/Type.II.HC.2.over.OtherThree.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC/Type.II.HC.2.over.OtherThree_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

##################################





MKS=c('Dok6','Mgat4c','Agbl1','Wdr49','Kcnj16')


FeaturePlot(pbmc,features=MKS,ncol=3 )
VlnPlot(pbmc,features=MKS,ncol=3,pt.size=0 )



pdf('UTRICLE_HC/UTRICLE_UMAP_LABEL_FEATURE2.pdf',width=16,height=10)
FeaturePlot(pbmc,features=MKS,ncol=3 )
VlnPlot(pbmc,features=MKS,ncol=3,pt.size=0 )
dev.off()


###############################################################


pbmc$new.clst.sub.label.time=paste0(pbmc@meta.data$orig.ident,'.',pbmc$new.clst.sub.label)

Idents(pbmc)=pbmc$new.clst.sub.label.time



###############################################################

TMP=c('Type.II.HC.1','Type.II.HC.2','Type.I.HC.1','Type.I.HC.2')

MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.12m','.',TMP),ident.2=paste0('utricle.3m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC_TIME/12m.over.3m.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/12m.over.3m_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################


###############################################################

TMP=c('Type.II.HC.1','Type.II.HC.2','Type.I.HC.1','Type.I.HC.2')

MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.3m','.',TMP),ident.2=paste0('utricle.12m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC_TIME/3m.over.12m.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/3m.over.12m_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################



###############################################################

TMP=c('Type.II.HC.1','Type.II.HC.2','Type.I.HC.1','Type.I.HC.2')

MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.22m','.',TMP),ident.2=paste0('utricle.12m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC_TIME/22m.over.12m.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/22m.over.12m_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################



TMP=c('Type.II.HC.1','Type.II.HC.2','Type.I.HC.1','Type.I.HC.2')

MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.12m','.',TMP),ident.2=paste0('utricle.22m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC_TIME/12m.over.22m.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/12m.over.22m_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################


TMP=c('Type.II.HC.1','Type.II.HC.2','Type.I.HC.1','Type.I.HC.2')

MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.22m','.',TMP),ident.2=paste0('utricle.3m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC_TIME/22m.over.3m.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/22m.over.3m_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################



TMP=c('Type.II.HC.1','Type.II.HC.2','Type.I.HC.1','Type.I.HC.2')

MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.3m','.',TMP),ident.2=paste0('utricle.22m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file='UTRICLE_HC_TIME/3m.over.22m.tsv',col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/3m.over.22m_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################








pbmc$new.clst.sub.label.time=paste0(pbmc@meta.data$orig.ident,'.',pbmc$new.clst.sub.label)

Idents(pbmc)=pbmc$new.clst.sub.label.time


###############################################################

GetMkGo<-function(TMP){



MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.12m','.',TMP),ident.2=paste0('utricle.3m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('UTRICLE_HC_TIME/12m.over.3m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/12m.over.3m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################


###############################################################

MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.3m','.',TMP),ident.2=paste0('utricle.12m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('UTRICLE_HC_TIME/3m.over.12m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/3m.over.12m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################

###############################################################

MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.22m','.',TMP),ident.2=paste0('utricle.12m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('UTRICLE_HC_TIME/22m.over.12m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/22m.over.12m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################



MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.12m','.',TMP),ident.2=paste0('utricle.22m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('UTRICLE_HC_TIME/12m.over.22m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/12m.over.22m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################

MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.22m','.',TMP),ident.2=paste0('utricle.3m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('UTRICLE_HC_TIME/22m.over.3m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/22m.over.3m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

###############################################################



MK.RESULT=FindMarkers(pbmc, ident.1=paste0('utricle.3m','.',TMP),ident.2=paste0('utricle.22m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('UTRICLE_HC_TIME/3m.over.22m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('UTRICLE_HC_TIME/3m.over.22m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

}
###############################################################


###############################################################

GetMkGo(TMP='Type.I.HC.1')
GetMkGo(TMP='Type.I.HC.2')
GetMkGo(TMP='Type.II.HC.1')
GetMkGo(TMP='Type.II.HC.2')


######################################################################



























###############################################################

write.table(table(pbmc$new.clst.sub.label.time),row.names=T,col.names=F,file='UTRICLE_HC_TIME/TIME_CELL_NUMBER.txt',sep='\t',quote=F)

##############################################################

# Aging Trajectory



Idents(pbmc)=pbmc$new.clst.sub.label
USED=which(Idents(pbmc) %in% TMP)

VEC=pbmc@reductions$umap@cell.embeddings

USED_VEC=VEC[USED,]

plot(USED_VEC)

USED_VEC=VEC[USED,]
USED_TIME=rep(NA, length(USED))
USED_TIME[which(pbmc$orig.ident[USED]=='utricle.3m')]=1 #3m
USED_TIME[which(pbmc$orig.ident[USED]=='utricle.12m')]=2 #12m
USED_TIME[which(pbmc$orig.ident[USED]=='utricle.22m')]=3 #22m
USED_TYPE=pbmc$new.clst.sub.label[USED]

USED_DATA=as.matrix(pbmc@assays$RNA@data)[,USED]
USED_COUNT=as.matrix(pbmc@assays$RNA@counts)[,USED]
USED_PCA=pbmc@reductions$pca@cell.embeddings[USED,]
USED_META=pbmc@meta.data[USED,]

########################################



        pbmc.HC=CreateSeuratObject(counts = USED_COUNT, min.cells = 0, min.features = 0, project = "ALL")
        pbmc.HC <- NormalizeData(object = pbmc.HC, normalization.method = "LogNormalize", scale.factor = 10000)
        all.genes=rownames(pbmc.HC)
        pbmc.HC <- ScaleData(object = pbmc.HC, features = all.genes)

    pbmc.HC=FindVariableFeatures(object = pbmc.HC, selection.method = "vst", nfeatures = 2000)
    pbmc.HC <- RunPCA(object = pbmc.HC, seed.use=123, npcs=150, features = VariableFeatures(object = pbmc.HC), ndims.print=1,nfeatures.print=1)
  
  
    pbmc.HC <- RunUMAP(pbmc.HC, dims = 1:150,n.components=2)

    pbmc.HC@meta.data=USED_META


#saveRDS(pbmc.HC, file='UTRICLE_HC_TIME/pbmc.HC.rds')

pdf('UTRICLE_HC_TIME/OnlyHC_UMAP.pdf',width=10,height=10)
DimPlot(pbmc.HC,group.by='batch')
DimPlot(pbmc.HC,group.by='new.clst.sub.label',label=TRUE, label.size=10)+NoLegend()
dev.off()



pbmc.HC.type1.1=subset(pbmc.HC, cells=colnames(pbmc.HC)[which(pbmc.HC$new.clst.sub.label=='Type.I.HC.1')])
pbmc.HC.type1.2=subset(pbmc.HC, cells=colnames(pbmc.HC)[which(pbmc.HC$new.clst.sub.label=='Type.I.HC.2')])
pbmc.HC.type2.1=subset(pbmc.HC, cells=colnames(pbmc.HC)[which(pbmc.HC$new.clst.sub.label=='Type.II.HC.1')])
pbmc.HC.type2.2=subset(pbmc.HC, cells=colnames(pbmc.HC)[which(pbmc.HC$new.clst.sub.label=='Type.II.HC.2')])

#############

this_pbmc=pbmc.HC.type1.1
    
    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

pbmc.HC.type1.1=subset(this_pbmc, cells=SELECT)

#############


this_pbmc=pbmc.HC.type1.2
    
    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

pbmc.HC.type1.2=subset(this_pbmc, cells=SELECT)

#############


this_pbmc=pbmc.HC.type2.1
    
    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

pbmc.HC.type2.1=subset(this_pbmc, cells=SELECT)

#############


this_pbmc=pbmc.HC.type2.2
    
    this_vec=this_pbmc@reductions$umap@cell.embeddings
    plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
    library(gatepoints)
    selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
    SELECT=selectedPoints

pbmc.HC.type2.2=subset(this_pbmc, cells=SELECT)

#####################



#saveRDS(pbmc.HC.type1.1, file='UTRICLE_HC_TIME/pbmc.HC.type1.1.rds')
#saveRDS(pbmc.HC.type1.2, file='UTRICLE_HC_TIME/pbmc.HC.type1.2.rds')
#saveRDS(pbmc.HC.type2.1, file='UTRICLE_HC_TIME/pbmc.HC.type2.1.rds')
#saveRDS(pbmc.HC.type2.2, file='UTRICLE_HC_TIME/pbmc.HC.type2.2.rds')




pdf('UTRICLE_HC_TIME/SubHC_UMAP.pdf',width=6,height=6)
DimPlot(pbmc.HC.type1.1,group.by='new.clst.sub.label',label=TRUE, label.size=10)+NoLegend()
DimPlot(pbmc.HC.type1.2,group.by='new.clst.sub.label',label=TRUE, label.size=10)+NoLegend()
DimPlot(pbmc.HC.type2.1,group.by='new.clst.sub.label',label=TRUE, label.size=10)+NoLegend()
DimPlot(pbmc.HC.type2.2,group.by='new.clst.sub.label',label=TRUE, label.size=10)+NoLegend()
DimPlot(pbmc.HC.type1.1,group.by='batch',label=TRUE, label.size=10)+NoLegend()
DimPlot(pbmc.HC.type1.2,group.by='batch',label=TRUE, label.size=10)+NoLegend()
DimPlot(pbmc.HC.type2.1,group.by='batch',label=TRUE, label.size=10)+NoLegend()
DimPlot(pbmc.HC.type2.2,group.by='batch',label=TRUE, label.size=10)+NoLegend()
dev.off()






NNN=10


pdf('UTRICLE_HC_TIME/SubHC_TRA.pdf',width=3.5,height=4)
#################################
this_pbmc=pbmc.HC.type1.1

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

    this_v=scale(pred.time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='fitted.time')

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=F)

pbmc.HC.type1.1$pseudotime=OUT$P.PS
FeaturePlot(pbmc.HC.type1.1,features=c('pseudotime'))
#######################################################################

#######################################################################
this_pbmc=pbmc.HC.type1.2

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

    this_v=scale(pred.time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='fitted.time')

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=F)

pbmc.HC.type1.2$pseudotime=OUT$P.PS
FeaturePlot(pbmc.HC.type1.2,features=c('pseudotime'))
#######################################################################


#######################################################################
this_pbmc=pbmc.HC.type2.1

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

    this_v=scale(pred.time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='fitted.time')

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=F)

pbmc.HC.type2.1$pseudotime=OUT$P.PS
FeaturePlot(pbmc.HC.type2.1,features=c('pseudotime'))
#######################################################################


#######################################################################
this_pbmc=pbmc.HC.type2.2

this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.3m')]=1 #3m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.12m')]=2 #12m
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.22m')]=3 #22m
this_vec=this_pbmc@reductions$umap@cell.embeddings

    this_v=scale(this_pbmc$time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')

this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)

    this_v=scale(pred.time)
    this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
    plot(this_vec, cex=0.5,pch=16,col=this_col,main='fitted.time')

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=F)

pbmc.HC.type2.2$pseudotime=OUT$P.PS
FeaturePlot(pbmc.HC.type2.2,features=c('pseudotime'))
#######################################################################

dev.off()







this_pbmc=pbmc.HC.type1.1

this_cor=cor( t(as.matrix(this_pbmc@assays$RNA@data)),this_pbmc$pseudotime)
this_cor[which(is.na(this_cor))]=0
this_cor=this_cor[,1]
names(this_cor)=rownames(this_pbmc)
plot(sort(this_cor),type='p',pch=16)
write.table(sort(this_cor), paste0('UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type1.1.tsv'),sep='\t',quote=F,col.names=F,row.names=T)


###############
POS.CUT=0.1
NEG.CUT= -0.1

###########################################
    this_sig_gene=names(this_cor)[which(this_cor>POS.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC_TIME/PseudotimeGO_POS_HC.type1.1_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC_TIME/PseudotimeGO_NEG_HC.type1.1_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)








this_pbmc=pbmc.HC.type1.2

used=which(!is.na(this_pbmc$pseudotime))
this_cor=cor( t(as.matrix(this_pbmc@assays$RNA@data[,used])),this_pbmc$pseudotime[used])
this_cor[which(is.na(this_cor))]=0
this_cor=this_cor[,1]
names(this_cor)=rownames(this_pbmc)
plot(sort(this_cor),type='p',pch=16)
write.table(sort(this_cor), paste0('UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type1.2.tsv'),sep='\t',quote=F,col.names=F,row.names=T)



###############
POS.CUT=0.1
NEG.CUT= -0.1

###########################################
    this_sig_gene=names(this_cor)[which(this_cor>POS.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC_TIME/PseudotimeGO_POS_HC.type1.2_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC_TIME/PseudotimeGO_NEG_HC.type1.2_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

################################################################







this_pbmc=pbmc.HC.type2.1

used=which(!is.na(this_pbmc$pseudotime))
this_cor=cor( t(as.matrix(this_pbmc@assays$RNA@data[,used])),this_pbmc$pseudotime[used])
this_cor[which(is.na(this_cor))]=0
this_cor=this_cor[,1]
names(this_cor)=rownames(this_pbmc)
plot(sort(this_cor),type='p',pch=16)
write.table(sort(this_cor), paste0('UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type2.1.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




###############
POS.CUT=0.1
NEG.CUT= -0.1

###########################################
    this_sig_gene=names(this_cor)[which(this_cor>POS.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC_TIME/PseudotimeGO_POS_HC.type2.1_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC_TIME/PseudotimeGO_NEG_HC.type2.1_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

############################################











this_pbmc=pbmc.HC.type2.2

used=which(!is.na(this_pbmc$pseudotime))
this_cor=cor( t(as.matrix(this_pbmc@assays$RNA@data[,used])),this_pbmc$pseudotime[used])
this_cor[which(is.na(this_cor))]=0
this_cor=this_cor[,1]
names(this_cor)=rownames(this_pbmc)
plot(sort(this_cor),type='p',pch=16)
write.table(sort(this_cor), paste0('UTRICLE_HC_TIME/PseudotimeGeneCor_HC.type2.2.tsv'),sep='\t',quote=F,col.names=F,row.names=T)




###########################################
    this_sig_gene=names(this_cor)[which(this_cor>POS.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC_TIME/PseudotimeGO_POS_HC.type2.2_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('UTRICLE_HC_TIME/PseudotimeGO_NEG_HC.type2.2_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

############################################



saveRDS(pbmc.HC.type1.1, file='UTRICLE_HC_TIME/pbmc.HC.type1.1.rds')
saveRDS(pbmc.HC.type1.2, file='UTRICLE_HC_TIME/pbmc.HC.type1.2.rds')
saveRDS(pbmc.HC.type2.1, file='UTRICLE_HC_TIME/pbmc.HC.type2.1.rds')
saveRDS(pbmc.HC.type2.2, file='UTRICLE_HC_TIME/pbmc.HC.type2.2.rds')











USED.GENE=c('Adamts12','Aldh2','Apod','Apoe','Appl2','Cdh5','Ets1','Fgfr1','Hgf','Lbp','Lpl','Lrfn5','Lrrk2','Nr1h3','Pla2g4a','Pmp22','Prkca',
'Rora','Slc7a2','Slit2','Stat5b','Wnt5a')



this_pbmc=pbmc.HC.type2.1



used=which(!is.na(this_pbmc$pseudotime))
this_cor=cor( t(as.matrix(this_pbmc@assays$RNA@data[,used])),this_pbmc$pseudotime[used])
this_cor[which(is.na(this_cor))]=0
this_cor=this_cor[,1]
names(this_cor)=rownames(this_pbmc)


this_pseudotime=this_pbmc$pseudotime[used]
this_mat=t(apply(as.matrix(this_pbmc@assays$RNA@data[,used]),1,scale))
this_mat=this_mat[order(this_cor),]
this_mat[which(is.na(this_mat))]=0
this_mat=this_mat[which(rownames(this_mat) %in% USED.GENE),]


this_mat=this_mat[,order(this_pseudotime)]



.smooth<-function(x,df=50){
   y=smooth.spline(x,df=df)$y
   return(y)
}


this_smat=t(apply(this_mat,1,.smooth,df=500))
rownames(this_smat)=rownames(this_mat)
library('ComplexHeatmap')
library('circlize')

o.mat=this_smat
color_fun_3 =colorRamp2(c(-1,0,1 ), c('royalblue3','white','indianred3'))
ht=Heatmap(o.mat,row_title='',name="exp",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=FALSE, show_row_names=TRUE,
	      col=color_fun_3, border = TRUE
	      )


print(ht)


pdf('UTRICLE_HC_TIME/Infla_HEAT_TypeII.HC.1.pdf',width=6,height=5)
print(ht)
dev.off()

############################################################################################################












all_pbmc=readRDS('COM_SEURAT.rds')
CLST=readRDS( 'COM_SUBCLST1.rds')


Idents(all_pbmc)=CLST
DimPlot(all_pbmc,label=TRUE)+NoLegend()

Idents(all_pbmc)=paste0(all_pbmc$batch,'.',CLST)


table(Idents(all_pbmc))



TMP='25'

MK.RESULT=FindMarkers(all_pbmc, ident.1=paste0('cochlea.3m','.',TMP),ident.2=paste0('cochlea.12m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('COCHLEA_2526/3m.over.12m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('COCHLEA_2526/3m.over.12m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)



MK.RESULT=FindMarkers(all_pbmc, ident.1=paste0('cochlea.12m','.',TMP),ident.2=paste0('cochlea.3m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('COCHLEA_2526/12m.over.3m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('COCHLEA_2526/12m.over.3m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)




TMP='26'

MK.RESULT=FindMarkers(all_pbmc, ident.1=paste0('cochlea.3m','.',TMP),ident.2=paste0('cochlea.12m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('COCHLEA_2526/3m.over.12m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('COCHLEA_2526/3m.over.12m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)



MK.RESULT=FindMarkers(all_pbmc, ident.1=paste0('cochlea.12m','.',TMP),ident.2=paste0('cochlea.3m','.',TMP),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gene=rownames(MK.RESULT)
MK.RESULT$gene=gene


    this_sig_gene=rownames(MK.RESULT)[which(MK.RESULT$p_val_adj<0.05)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result

    write.table(MK.RESULT, file=paste0('COCHLEA_2526/12m.over.3m.',TMP,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
    write.table(this_out, paste0('COCHLEA_2526/12m.over.3m.',TMP,'_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)



pdf('COCHLEA_2526/VlnUMAP_Apod.pdf',width=7,height=7)
VlnPlot(all_pbmc, features=c('Apod'),pt.size=0,idents=c('cochlea.3m.26','cochlea.12m.26'))+NoLegend()
FeaturePlot(all_pbmc, features=c('Apod'),order=TRUE)
dev.off()








































































