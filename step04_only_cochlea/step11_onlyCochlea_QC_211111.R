
#/home/toolkit/tools/R4*/bin/R


setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS/Cochlea_211111')



D1P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_3M_Cochlea/outs/filtered_feature_bc_matrix'
D2P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_12M_Cochlea/outs/filtered_feature_bc_matrix'
D3P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_24M_Cochlea/outs/filtered_feature_bc_matrix'
D4P='/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2021-11028-snRNA_24M_Cochlea_BU/outs/filtered_feature_bc_matrix'




library(Seurat)


cochlea_3m_data=Read10X(D1P)
cochlea_12m_data=Read10X(D2P)
cochlea_24m_data=Read10X(D3P)
cochlea_24mbu_data=Read10X(D4P)


cochlea_3m=CreateSeuratObject(cochlea_3m_data,project='cochlea_3m')
cochlea_12m=CreateSeuratObject(cochlea_12m_data,project='cochlea_12m')
cochlea_24m=CreateSeuratObject(cochlea_24m_data,project='cochlea_24m')
cochlea_24mbu=CreateSeuratObject(cochlea_24mbu_data,project='cochlea_24mbu')



.MT=function(pbmc){
    pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^mt-")
    return(pbmc)
    }


cochlea_3m=.MT(cochlea_3m)
cochlea_12m=.MT(cochlea_12m)
cochlea_24m=.MT(cochlea_24m)
cochlea_24mbu=.MT(cochlea_24mbu)


.viewQC=function(pbmc){
    pdf(paste0(pbmc@project.name,'_QC.pdf'),width=6,height=4)
    this_plot=VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(this_plot)
    dev.off()
    return(pbmc)
    }



cochlea_3m=.viewQC(cochlea_3m)
cochlea_12m=.viewQC(cochlea_12m)
cochlea_24m=.viewQC(cochlea_24m)
cochlea_24mbu=.viewQC(cochlea_24mbu)




############################################################



.subSet=function(pbmc){
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
    return(pbmc)
    }


cochlea_3m=.subSet(cochlea_3m)
cochlea_12m=.subSet(cochlea_12m)
cochlea_24m=.subSet(cochlea_24m)
cochlea_24mbu=.subSet(cochlea_24mbu)


###################################################################################



.saveRDS<-function(pbmc){
    saveRDS(pbmc, paste0(pbmc@project.name,'_afterQC.rds'))
    return(pbmc)
    }


cochlea_3m=.saveRDS(cochlea_3m)
cochlea_12m=.saveRDS(cochlea_12m)
cochlea_24m=.saveRDS(cochlea_24m)
cochlea_24mbu=.saveRDS(cochlea_24mbu)




