
#/home/toolkit/tools/R4*/bin/R


setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS')

library(Seurat)
library(dplyr)
library(patchwork)

####################################
# 处理 Cochlea 24m 和 24mbu
###################################

D6P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_24M_Cochlea/outs/filtered_feature_bc_matrix'
D7P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_24M_Cochlea_BU/outs/filtered_feature_bc_matrix'

cochlea_24m=CreateSeuratObject(Read10X(D6P),project='cochlea_24m')
cochlea_24mbu=CreateSeuratObject(Read10X(D7P),project='cochlea_24mbu')


.MT=function(pbmc){
    pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^mt-")
    return(pbmc)
    }


cochlea_24m=.MT(cochlea_24m)
cochlea_24mbu=.MT(cochlea_24mbu)

#################################################################

.viewQC=function(pbmc){
    pdf(paste0(pbmc@project.name,'_QC.pdf'),width=6,height=4)
    this_plot=VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(this_plot)
    dev.off()
    return(pbmc)
    }


cochlea_24m=.viewQC(cochlea_24m)
cochlea_24mbu=.viewQC(cochlea_24mbu)


###############################################################
.subSet=function(pbmc){
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
    return(pbmc)
    }

cochlea_24m=.subSet(cochlea_24m)
cochlea_24mbu=.subSet(cochlea_24mbu)

##############################################




.saveRDS<-function(pbmc){
    saveRDS(pbmc, paste0(pbmc@project.name,'_afterQC.rds'))
    return(pbmc)
    }

cochlea_24m=.saveRDS(cochlea_24m)
cochlea_24mbu=.saveRDS(cochlea_24mbu)


####################################################

#######################################################
# 开始细分
######################################################
utricle_3m=readRDS('utricle_3m_afterQC.rds')
utricle_12m=readRDS('utricle_12m_afterQC.rds')
utricle_22m=readRDS('utricle_22m_afterQC.rds')

cochlea_3m=readRDS('cochlea_3m_afterQC.rds')
cochlea_12m=readRDS('cochlea_12m_afterQC.rds')
cochlea_24m=readRDS('cochlea_24m_afterQC.rds')
cochlea_24mbu=readRDS('cochlea_24mbu_afterQC.rds')
########################################





#####################################
# 各种参数
####################################
# 以下的参数只选一组
#####################################

PC_calculated=50
Cluster_usedPC=1:30
Cluster_resolution=1
UMAP_usedPC=1:30
TAG='C1' # 这种参数组合的名字

#################################################

PC_calculated=50
Cluster_usedPC=1:30
Cluster_resolution=1.5
UMAP_usedPC=1:30
TAG='C2' # 这种参数组合的名字
####################################################





####################################

.runSCTandUMAP<-function(pbmc){
    pbmc <- SCTransform(pbmc, variable.features.n =3000)
    pbmc <- RunPCA(pbmc,npcs=PC_calculated, ndims.print=1,nfeatures.print=1)
    pbmc <- FindNeighbors(pbmc, dims = Cluster_usedPC)
    pbmc <- FindClusters(pbmc, resolution = Cluster_resolution)
    pbmc <- RunUMAP(object = pbmc, reduction = "pca", dims = UMAP_usedPC, n.components=2)  
    ################
    return(pbmc)
    }

utricle_3m=.runSCTandUMAP(utricle_3m)
utricle_12m=.runSCTandUMAP(utricle_12m)
utricle_22m=.runSCTandUMAP(utricle_22m)

cochlea_3m=.runSCTandUMAP(cochlea_3m)
cochlea_12m=.runSCTandUMAP(cochlea_12m)
cochlea_24m=.runSCTandUMAP(cochlea_24m)
cochlea_24mbu=.runSCTandUMAP(cochlea_24mbu)



#######################################

.findMarkers <-function(pbmc){
     pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
     saveRDS(pbmc.markers, paste0(pbmc@project.name,'_markers','_',TAG,'.rds'))
     return(pbmc)
     }

utricle_3m=.findMarkers(utricle_3m)
utricle_12m=.findMarkers(utricle_12m)
utricle_22m=.findMarkers(utricle_22m)

cochlea_3m=.findMarkers(cochlea_3m)
cochlea_12m=.findMarkers(cochlea_12m)
cochlea_24m=.findMarkers(cochlea_24m)
cochlea_24mbu=.findMarkers(cochlea_24mbu)







#######################################


.viewUMAP=function(pbmc){
    pdf(paste0(pbmc@project.name,'_UMAP','_',TAG,'.pdf'),width=7,height=6)
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
cochlea_24m=.viewUMAP(cochlea_24m)
cochlea_24mbu=.viewUMAP(cochlea_24mbu)





#######################################

.viewHeat=function(pbmc){
    pbmc.markers=readRDS(paste0(pbmc@project.name,'_markers','_',TAG,'.rds'))
    top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    pdf(paste0(pbmc@project.name,'_HEAT','_',TAG,'.pdf'),width=12,height=19)
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
cochlea_24m=.viewHeat(cochlea_24m)
cochlea_24mbu=.viewHeat(cochlea_24mbu)







###################################################################################



.saveRDS<-function(pbmc){
    saveRDS(pbmc, paste0(pbmc@project.name,'_SCT_UMAP','_',TAG,'.rds'))
    return(pbmc)
    }

utricle_3m=.saveRDS(utricle_3m)
utricle_12m=.saveRDS(utricle_12m)
utricle_22m=.saveRDS(utricle_22m)

cochlea_3m=.saveRDS(cochlea_3m)
cochlea_12m=.saveRDS(cochlea_12m)
cochlea_24m=.saveRDS(cochlea_24m)
cochlea_24mbu=.saveRDS(cochlea_24mbu)



##########################################################
# 输出 tsv
############################################

m1=readRDS(paste0('utricle_3m_markers','_',TAG,'.rds'))
m2=readRDS(paste0('utricle_12m_markers','_',TAG,'.rds'))
m3=readRDS(paste0('utricle_22m_markers','_',TAG,'.rds'))

m4=readRDS(paste0('cochlea_3m_markers','_',TAG,'.rds'))
m5=readRDS(paste0('cochlea_12m_markers','_',TAG,'.rds'))
m6=readRDS(paste0('cochlea_24m_markers','_',TAG,'.rds'))
m7=readRDS(paste0('cochlea_24mbu_markers','_',TAG,'.rds'))


write.table(m1, file=paste0('utricle_3m_markers','_',TAG,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
write.table(m2, file=paste0('utricle_12m_markers','_',TAG,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
write.table(m3, file=paste0('utricle_22m_markers','_',TAG,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)

write.table(m4, file=paste0('cochlea_3m_markers','_',TAG,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
write.table(m5, file=paste0('cochlea_12m_markers','_',TAG,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
write.table(m6, file=paste0('cochlea_24m_markers','_',TAG,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)
write.table(m7, file=paste0('cochlea_24mbu_markers','_',TAG,'.tsv'),col.names=T,row.names=F,sep='\t',quote=F)

