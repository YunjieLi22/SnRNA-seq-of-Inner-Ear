

# /home/toolkit/tools/R4.0.3/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new')

library(Seurat)
library(dplyr)
library(patchwork)

##################################
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')


UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')


COMBINE$time[which(COMBINE$time %in% c('U.3m'))]='U.03m'
COMBINE$time[which(COMBINE$time %in% c('C.3m'))]='C.03m'
COUNT=COMBINE@assays$RNA@counts

#################################################

TAG=paste0(COMBINE$clst,'_',COMBINE$time)
TISSUE=c(rep('U',ncol(UTRICLE)),rep('C',ncol(COCHLEA)))






#################
# COCHLEA
#################

pbmc=COCHLEA
pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$new.clst.sub %in% c(14))])
DimPlot(pbmc)

VEC=pbmc@reductions$umap@cell.embeddings
pbmc=subset(pbmc, cells=colnames(pbmc)[which(VEC[,1]<0 & VEC[,2]>3)])
DimPlot(pbmc)


VEC=pbmc@reductions$umap@cell.embeddings
set.seed(123)
KM=kmeans(VEC,center=6)
Idents(pbmc)=factor(KM$cluster,levels=sort(unique(KM$cluster)))
DimPlot(pbmc, reduction = "umap",label=T)+NoLegend()


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


########################################################

write.table(pbmc.markers,
file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/OUT20220523/MACROPHAGE_COCHLEA_MARKER.txt',
quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/OUT20220523/MACROPHAGE_COCHLEA_PLOT.pdf',width=7,height=7)
DimPlot(pbmc, reduction = "umap",label=T)+NoLegend()
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()


################################################










#################
# UTRICLE
#################

pbmc=UTRICLE
pbmc=subset(pbmc, cells=colnames(pbmc)[which(pbmc$new.clst.sub %in% c(16.1))])
DimPlot(pbmc)

VEC=pbmc@reductions$umap@cell.embeddings
pbmc=subset(pbmc, cells=colnames(pbmc)[which(VEC[,1]< -5 & VEC[,2]< -4)])
DimPlot(pbmc)


VEC=pbmc@reductions$umap@cell.embeddings
set.seed(123)
KM=kmeans(VEC,center=5)
Idents(pbmc)=factor(KM$cluster,levels=sort(unique(KM$cluster)))
DimPlot(pbmc, reduction = "umap",label=T)+NoLegend()


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


########################################################

write.table(pbmc.markers,
file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/OUT20220523/MACROPHAGE_UTRICLE_MARKER.txt',
quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')

pdf('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/OUT20220523/MACROPHAGE_UTRICLE_PLOT.pdf',width=7,height=7)
DimPlot(pbmc, reduction = "umap",label=T)+NoLegend()
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()


################################################

this_pbmc=pbmc.Macrophage

this_vec=this_pbmc@reductions$umap@cell.embeddings
plot(this_vec, pch=16,col='grey50',cex=0.5,mark = TRUE)
library(gatepoints)
selectedPoints <- fhs(this_vec, pch=16,col='red3',cex=0.5,mark = TRUE)
SELECT=selectedPoints

pbmc.Macrophage=subset(this_pbmc, cells=SELECT)

saveRDS(pbmc.Macrophage, file='data/cochlea_Macrophage.rds')

NNN=10
this_pbmc=pbmc.Macrophage

DimPlot(this_pbmc, group.by='batch',label=TRUE)+NoLegend()

*utricle
this_pbmc$time=rep(NA, length(ncol(this_pbmc)))
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.3m')]=1  
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.12m')]=2   
this_pbmc$time[which(this_pbmc$orig.ident=='utricle.22m')]=3  
*cochlea
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.3m')]=1
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.12m')]=2
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.24m')]=3
this_pbmc$time[which(this_pbmc$orig.ident=='cochlea.24mbu')]=3

this_pbmc$time[which(this_pbmc$time == "3m")]=1 
this_pbmc$time[which(this_pbmc$time == "12m")]=2
this_pbmc$time[which(this_pbmc$time == "24m")]=3

this_vec=this_pbmc@reductions$umap@cell.embeddings

source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')


this_v=scale(as.numeric(this_pbmc$time))
this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
plot(this_vec, cex=0.5,pch=16,col=this_col,main='orig.time')


this_fit=lm(this_pbmc$time~.:.,data=as.data.frame(this_vec))
pred.time=predict(this_fit)


this_v=scale(pred.time)
this_col=vector.vcol(this_v,CN=c('indianred1','gold1','royalblue1'),CV=c(-1,0,1))
plot(this_vec, cex=0.5,pch=16,col=this_col,main='fitted.time')

"#BFC6FF", "#FFC1DF", "#8EFFB5"

OUT=vector.buildGrid(this_vec, N=NNN,SHOW=F)
OUT=vector.buildNet(OUT, CUT=1, SHOW=F)
OUT$VALUE=max(pred.time)-pred.time
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=F)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=F)
pbmc.Macrophage$pseudotime=OUT$P.PS
FeaturePlot(pbmc.Macrophage,features=c('pseudotime'))








ALL_ENTREZ=mapIds(org.Mm.eg.db, rownames(pbmc), 'ENTREZID', 'SYMBOL')

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
    write.table(this_out, paste0('./PseudotimeGO_POS_14_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)


#######################################
    this_sig_gene=names(this_cor)[which(this_cor< NEG.CUT)]
    this_sig_entrez=mapIds(org.Mm.eg.db, this_sig_gene, 'ENTREZID', 'SYMBOL')
    this_ego <- enrichGO(gene = this_sig_entrez, 
                   universe = ALL_ENTREZ,OrgDb = org.Mm.eg.db,ont = "ALL", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,readable = TRUE)

    this_out=this_ego@result
    write.table(this_out, paste0('./PseudotimeGO_NEG_14_go.tsv'),sep='\t',quote=F,col.names=T,row.names=F)

#########################################################


co_expression_matrix <- GetAssayData(this_pbmc, slot = "data")[rownames(this_pbmc) %in% co_genes, ]
cochclea_matrix <- as.data.frame(co_expression_matrix)
co_pseudotime <- pbmc.Macrophage$pseudotime
co_pseudotime <- t(as.data.frame(co_pseudotime ))
new_co_matrix <- rbind(cochlea_matrix, co_pseudotime) 
new_co_matrix <- as.matrix(new_co_matrix)
new_co_matrix1 <- new_co_matrix[ ,order(new_co_matrix["co_pseudotime", ])] 
  
row_sums <- apply(utri_all_matrix, 1, sum)
utri_all_matrix<- utri_all_matrix[row_sums != 0, ]
  
y_list <- list()
for (i in 1:nrow(utricle_matrix)) {
value <- as.numeric(utricle_matrix[i, ])
y <- smooth.spline(value, df = 8)$y
y_list[[i]] <- y
}
row_names <- row.names(utricle_matrix)
new_utri_matrix <- as.data.frame(do.call(rbind, y_list), row.names = row_names)


color1 <- colorRamp2(c(-1, 0, 10), c("blue", "white", "red"))
color2 <- colorRamp2(c(-0.1, 0, 0.1,0.25,0.5, 1), c("white",  "beige", "lightgoldenrod1", "gold1", "orange", "orangered"))
color3 <- colorRamp2(c(-0.1, 0, 0.25,0.5, 0.75, 1), c("white", "aliceblue","lightskyblue","skyblue2", "royalblue", "blue"))
color4 <- colorRamp2(c(-0.1, 0, 0.25,0.5, 0.75, 1), c("white", "pink", "lightcoral", "coral", "tomato", "red"))
colorRampPalette(c("blue", "white", "red"))
color5 <- colorRamp2(c(-2, 0, 2), c("#F16D6C", "white", "#4F94CD"))

Heatmap(new_utri_matrix,
name = "Gene Expression",
col = color1,
color_space = "LAB",
cluster_rows = T,
cluster_columns = F,
show_row_names = TRUE,
show_column_names = FALSE,
column_title = "Pseudotime",
column_title_side = "top",
row_names_gp = gpar(fontsize = 8),
border = NA )


Heatmap(utricle_matrix,
        name = "Gene Expression",
        col = color_fun,
        cluster_rows = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        show_row_names = TRUE,
        show_column_names = FALSE,
        column_title = "Pseudotime",
        column_title_side = "top",
        column_order = pseudotime_data, 
        row_names_gp = gpar(fontsize = 8), 
        border = NA )



