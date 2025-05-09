CELLRANGER_ARC=/home/toolkit/tools/cellranger-arc-2.0.0/cellranger-arc

CPU=40
MEM=400

REFERENCE=/home/database/CellRanger/refdata-cellranger-arc-mm10-2020-A-2.0.0


TAG=HT2020-21864-MultiOmics_24M_AC
CSV=/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2020-21864-MultiOmics_24M_AC.csv
$CELLRANGER_ARC count --id $TAG --reference $REFERENCE --libraries $CSV --localcores $CPU --localmem $MEM


TAG=HT2020-21864-MultiOmics_24M_AC_BU
CSV=/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2020-21864-MultiOmics_24M_AC_BU.csv
$CELLRANGER_ARC count --id $TAG --reference $REFERENCE --libraries $CSV --localcores $CPU --localmem $MEM


TAG=HT2020-21864-MultiOmics_3M_AC
CSV=/home/disk/database/data/SingleCellAuditory_Xia/analysis/HT2020-21864-MultiOmics_3M_AC.csv
$CELLRANGER_ARC count --id $TAG --reference $REFERENCE --libraries $CSV --localcores $CPU --localmem $MEM


