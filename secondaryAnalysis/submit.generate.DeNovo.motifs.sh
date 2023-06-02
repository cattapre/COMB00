#!/bin/bash

#SBATCH --job-name=mot.gen
#SBATCH --time=10:00


SOURCEBED_DIR="/group/nordlab/users/rinaldo/COMB00/ChIPseq/data/bedfiles"
echo
echo "generate.motifs.sh"
echo $1
PV_THSHLD=$1
echo

MYF=`ls ${SOURCEBED_DIR}/*.motifs.bed`
MYF="${MYF//${SOURCEBED_DIR}\//}"

for MYFILE in $MYF
do
    echo 
    echo $MYFILE
    sbatch /group/nordlab/users/rinaldo/COMB00/ChIPseq/codes/COMB00.generate.DeNovo.motifs.sh $SOURCEBED_DIR $MYFILE $PV_THSHLD
done
