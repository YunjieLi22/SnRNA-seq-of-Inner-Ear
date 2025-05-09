
#/home/toolkit/tools/R4.2.0/bin/R

setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_20221213/')

library(Seurat)
library(dplyr)
library(patchwork)

UTRICLE = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/UTRICLE_SEURAT.rds')
COCHLEA = readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COCHLEA_SEURAT.rds')
COMBINE= readRDS('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_202204new/COM_SEURAT.rds')

COCHLEA$final.clst.sub.label = COCHLEA$new.clst.sub
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('0','2'))]='Tympanic.border.cell.1'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('1.1'))]='Claudius.cell'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('1.2'))]='Deiters.cell'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('4'))]='Satellite.glial.cells'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('5'))]='Inner.sulcus.cell'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('6'))]='Outer.sulcus.cell'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('7'))]='Dentate.interdental.cell'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('8'))]='IHC'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('9.1'))]='Cochlea.Fibrocyte.1'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('9.2'))]='Cochlea.Fibrocyte.2'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('10'))]='Outer.pillar.cell'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('11'))]='OHC'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('12'))]='Spiral.ganglion.neuron'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('13'))]='Cochlea.Fibrocyte.3'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('14'))]='Macrophage'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('15'))]='Inner.phalangeal.cells.Inner.border'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('16'))]='Inner.pillar.cell'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('17'))]='Astrocytes'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('18'))]='Melanocyte'
COCHLEA$final.clst.sub.label[which(COCHLEA$new.clst.sub %in% c('19'))]='Endothelial.cell'

UTRICLE$final.clst.sub.label = UTRICLE$new.clst.sub
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in% c('0'))]='Extrastriolar.SC'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in% c('1'))]='Type.II.HC.1'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('2'))]='Type.II.HC.2'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in% c('3.1'))]='Fibrocyte.1'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in% c('3.2'))]='Fibrocyte.3'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('4'))]='Roof.epithelial.cell.1'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('5'))]='Type.I.HC.1'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('7'))]='Translational epithelial cell'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('8'))]='Roof.epithelial.cell.2'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('9'))]='Roof.epithelial.cell.2'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('10'))]='Roof.epithelial.cell 2'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('11'))]='Striolar.SC'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('12'))]='Fibrocyte.2'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('13'))]='Type.I.HC.2'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('14'))]='Fibrocyte.4'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('15'))]='Astrocytes'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('16.1'))]='Macrophage'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('16.2'))]='Microglia'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in% c('17'))]='Vestibular ganglion neuron'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in%  c('18'))]='Endothelial.cell'
UTRICLE$final.clst.sub.label[which(UTRICLE$new.clst.sub %in% c('19'))]='Melanocyte'

UL=cbind(paste0('Utricle.',UTRICLE$final.clst.sub.label), paste0('U.',UTRICLE$new.clst.sub))
CL=cbind(paste0('Cochlea.',COCHLEA$final.clst.sub.label), paste0('C.',COCHLEA$new.clst.sub))
LABEL=unique(rbind(UL,CL))
TYPE=LABEL[match(COMBINE$clst,LABEL[,2]),1]
COMBINE$final.clst.sub.label=TYPE

COMBINE$final.type=paste0(COMBINE$final.clst.sub.label, '.',COMBINE$time)


USED_TYPE=c( 'Cochlea.IHC.C.3m','Cochlea.OHC.C.3m','Utricle.Type.I.HC.1.U.3m','Utricle.Type.I.HC.2.U.3m','Utricle.Type.II.HC.1.U.3m','Utricle.Type.II.HC.2.U.3m',
             'Cochlea.IHC.C.12m','Cochlea.OHC.C.12m','Utricle.Type.I.HC.1.U.12m','Utricle.Type.I.HC.2.U.12m','Utricle.Type.II.HC.1.U.12m','Utricle.Type.II.HC.2.U.12m',
             'Cochlea.IHC.C.24m','Cochlea.OHC.C.24m','Utricle.Type.I.HC.1.U.22m','Utricle.Type.I.HC.2.U.22m','Utricle.Type.II.HC.1.U.22m','Utricle.Type.II.HC.2.U.22m')

this_pbmc=subset(COMBINE, cells=colnames(COMBINE)[which(COMBINE$final.type %in% USED_TYPE)])

Idents(this_pbmc)=factor( this_pbmc$final.type,levels=USED_TYPE )




saveRDS(COMBINE@meta.data,file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_20221213/COMBINE_META.rds')
saveRDS(UTRICLE@meta.data,file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_20221213/UTRICLE_META.rds')
saveRDS(COCHLEA@meta.data,file='/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS_20221213/COCHLEA_META.rds')

