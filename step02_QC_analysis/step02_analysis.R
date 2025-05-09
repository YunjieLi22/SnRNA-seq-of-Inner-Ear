
#/home/toolkit/tools/R4*/bin/R


setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS')

library(Seurat)
library(dplyr)
library(patchwork)


utricle_3m=readRDS('utricle_3m_afterQC.rds')
utricle_12m=readRDS('utricle_12m_afterQC.rds')
utricle_22m=readRDS('utricle_22m_afterQC.rds')

cochlea_3m=readRDS('cochlea_3m_afterQC.rds')
cochlea_12m=readRDS('cochlea_12m_afterQC.rds')



####################################

.runSCTandUMAP<-function(pbmc){
    pbmc <- SCTransform(pbmc, variable.features.n =3000)
    pbmc <- RunPCA(pbmc,npcs=30, ndims.print=1,nfeatures.print=1)
    pbmc <- FindNeighbors(pbmc, dims = 1:10)
    pbmc <- FindClusters(pbmc, resolution = 0.5)
    pbmc <- RunUMAP(object = pbmc, reduction = "pca", dims = 1:30, n.components=2)  
    ################
    return(pbmc)
    }

utricle_3m=.runSCTandUMAP(utricle_3m)
utricle_12m=.runSCTandUMAP(utricle_12m)
utricle_22m=.runSCTandUMAP(utricle_22m)

cochlea_3m=.runSCTandUMAP(cochlea_3m)
cochlea_12m=.runSCTandUMAP(cochlea_12m)



#######################################

.findMarkers <-function(pbmc){
     pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
     saveRDS(pbmc.markers, paste0(pbmc@project.name,'_markers.rds'))
     return(pbmc)
     }

utricle_3m=.findMarkers(utricle_3m)
utricle_12m=.findMarkers(utricle_12m)
utricle_22m=.findMarkers(utricle_22m)

cochlea_3m=.findMarkers(cochlea_3m)
cochlea_12m=.findMarkers(cochlea_12m)







#######################################


.viewUMAP=function(pbmc){
    pdf(paste0(pbmc@project.name,'_UMAP.pdf'),width=7,height=6)
    this_plot=DimPlot(pbmc, label=TRUE )
    print(this_plot)
    dev.off()
    return(pbmc)
    }



utricle_3m=.viewUMAP(utricle_3m)
utricle_12m=.viewUMAP(utricle_12m)
utricle_22m=.viewUMAP(utricle_22m)

cochlea_3m=.viewUMAP(cochlea_3m)
cochlea_12m=.viewUMAP(cochlea_12m)





#######################################

.viewHeat=function(pbmc){
    pbmc.markers=readRDS(paste0(pbmc@project.name,'_markers.rds'))
    top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    pdf(paste0(pbmc@project.name,'_HEAT.pdf'),width=12,height=16)
    this_plot=DoHeatmap(pbmc, features = top10 $gene) + NoLegend()
    print(this_plot)
    dev.off()
    return(pbmc)
    }



utricle_3m=.viewHeat(utricle_3m)
utricle_12m=.viewHeat(utricle_12m)
utricle_22m=.viewHeat(utricle_22m)

cochlea_3m=.viewHeat(cochlea_3m)
cochlea_12m=.viewHeat(cochlea_12m)







###################################################################################



.saveRDS<-function(pbmc){
    saveRDS(pbmc, paste0(pbmc@project.name,'_SCT_UMAP.rds'))
    return(pbmc)
    }

utricle_3m=.saveRDS(utricle_3m)
utricle_12m=.saveRDS(utricle_12m)
utricle_22m=.saveRDS(utricle_22m)

cochlea_3m=.saveRDS(cochlea_3m)
cochlea_12m=.saveRDS(cochlea_12m)
