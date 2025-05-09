#/home/toolkit/tools/R4*/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS')

D1P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_3M_Utricle/outs/filtered_feature_bc_matrix'
D2P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_12M_Utricle/outs/filtered_feature_bc_matrix'
D3P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_22M_Utricle/outs/filtered_feature_bc_matrix'

D4P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_3M_Cochlea/outs/filtered_feature_bc_matrix'
D5P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_12M_Cochlea/outs/filtered_feature_bc_matrix'

library(Seurat)

utricle_3m_data=Read10X(D1P)
utricle_12m_data=Read10X(D2P)
utricle_22m_data=Read10X(D3P)

cochlea_3m_data=Read10X(D4P)
cochlea_12m_data=Read10X(D5P)


utricle_3m=CreateSeuratObject(utricle_3m_data,project='utricle_3m')
utricle_12m=CreateSeuratObject(utricle_12m_data,project='utricle_12m')
utricle_22m=CreateSeuratObject(utricle_22m_data,project='utricle_22m')

cochlea_3m=CreateSeuratObject(cochlea_3m_data,project='cochlea_3m')
cochlea_12m=CreateSeuratObject(cochlea_12m_data,project='cochlea_12m')




.MT=function(pbmc){
    pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^mt-")
    return(pbmc)
    }



utricle_3m=.MT(utricle_3m)
utricle_12m=.MT(utricle_12m)
utricle_22m=.MT(utricle_22m)

cochlea_3m=.MT(cochlea_3m)
cochlea_12m=.MT(cochlea_12m)

###########################################################################################

.viewQC=function(pbmc){
    pdf(paste0(pbmc@project.name,'_QC.pdf'),width=6,height=4)
    this_plot=VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(this_plot)
    dev.off()
    return(pbmc)
    }



utricle_3m=.viewQC(utricle_3m)
utricle_12m=.viewQC(utricle_12m)
utricle_22m=.viewQC(utricle_22m)

cochlea_3m=.viewQC(cochlea_3m)
cochlea_12m=.viewQC(cochlea_12m)

############################################################



.subSet=function(pbmc){
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
    return(pbmc)
    }



utricle_3m=.subSet(utricle_3m)
utricle_12m=.subSet(utricle_12m)
utricle_22m=.subSet(utricle_22m)

cochlea_3m=.subSet(cochlea_3m)
cochlea_12m=.subSet(cochlea_12m)

###################################################################################



.saveRDS<-function(pbmc){
    saveRDS(pbmc, paste0(pbmc@project.name,'_afterQC.rds'))
    return(pbmc)
    }



utricle_3m=.saveRDS(utricle_3m)
utricle_12m=.saveRDS(utricle_12m)
utricle_22m=.saveRDS(utricle_22m)

cochlea_3m=.saveRDS(cochlea_3m)
cochlea_12m=.saveRDS(cochlea_12m)

