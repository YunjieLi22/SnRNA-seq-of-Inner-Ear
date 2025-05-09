
CELLRANGER=/home/toolkit/tools/cellranger-6.1.0/cellranger

CPU=40
MEM=400

REFERENCE=/home/database/CellRanger/refdata-gex-mm10-2020-A


TAG=HT2021-11028-snRNA_22M_Utricle
SAMPLE=R21016303-OES210230150-OES210230150
FASTQ=/home/disk/database/data/SingleCellAuditory_Xia/HT2021-11028-snRNA/HT2021-11028-2-22M/22M_Utricle
$CELLRANGER count --include-introns --fastqs $FASTQ --localcores $CPU --sample $SAMPLE --localmem $MEM --id $TAG --transcriptome $REFERENCE


TAG=HT2021-11028-snRNA_3M_Cochlea
SAMPLE=3M_Cochlea
FASTQ=/home/disk/database/data/SingleCellAuditory_Xia/HT2021-11028-snRNA/HT2021-11028-4-3M/3M_Cochlea
$CELLRANGER count --include-introns --fastqs $FASTQ --localcores $CPU --sample $SAMPLE --localmem $MEM --id $TAG --transcriptome $REFERENCE


TAG=HT2021-11028-snRNA_3M_Utricle
SAMPLE=R21073209-OES211322150-OES211322150
FASTQ=/home/disk/database/data/SingleCellAuditory_Xia/HT2021-11028-snRNA/HT2021-11028-4-3M/3M_Utricle
$CELLRANGER count --include-introns --fastqs $FASTQ --localcores $CPU --sample $SAMPLE --localmem $MEM --id $TAG --transcriptome $REFERENCE


TAG=HT2021-11028-snRNA_12M_Cochlea
SAMPLE=R21046880-OES210895150-OES210895150
FASTQ=/home/disk/database/data/SingleCellAuditory_Xia/HT2021-11028-snRNA/HT2021-11028-3-12M/12M_Cochlea
$CELLRANGER count --include-introns --fastqs $FASTQ --localcores $CPU --sample $SAMPLE --localmem $MEM --id $TAG --transcriptome $REFERENCE


TAG=HT2021-11028-snRNA_12M_Utricle
SAMPLE=R21046879-OES210894150-OES210894150
FASTQ=/home/disk/database/data/SingleCellAuditory_Xia/HT2021-11028-snRNA/HT2021-11028-3-12M/12M_Utricle
$CELLRANGER count --include-introns --fastqs $FASTQ --localcores $CPU --sample $SAMPLE --localmem $MEM --id $TAG --transcriptome $REFERENCE



