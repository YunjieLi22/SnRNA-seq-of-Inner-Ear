


#/home/toolkit/tools/R4*/bin/R


setwd('/home/disk/database/data/SingleCellAuditory_Xia/analysis/RDS')


m1=readRDS('utricle_3m_markers.rds')
m2=readRDS('utricle_12m_markers.rds')
m3=readRDS('utricle_22m_markers.rds')
m4=readRDS('cochlea_3m_markers.rds')
m5=readRDS('cochlea_12m_markers.rds')


m6=readRDS('COM_MAKER.rds')


write.table(m1, file='utricle_3m_markers.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(m2, file='utricle_12m_markers.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(m3, file='utricle_22m_markers.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(m4, file='cochlea_3m_markers.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(m5, file='cochlea_12m_markers.tsv',col.names=T,row.names=F,sep='\t',quote=F)
write.table(m6, file='COM_MAKER.tsv',col.names=T,row.names=F,sep='\t',quote=F)

