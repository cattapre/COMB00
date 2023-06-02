#!/bin/bash
# Homer.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=Homer

echo
echo 
echo Homer.sbatch version 1.0.0
echo ... started running on $(date)
echo  
echo  

RANDOM=$$
MYTEMP="env.${RANDOM}.tmp"
sinfo > $MYTEMP
MY=`grep 'low' $MYTEMP`
rm $MYTEMP
if [[ "$MY" == "" ]]
then
	MYENV="barbera"
else
	MYENV="crick"
fi

echo
echo This job is running on
/bin/hostname
echo $MYENV
echo
echo

# Load Modules, and Run Programs
module load homer


GENOME=$1
HM_DIR=$2
F_TYPE=$3
SIZE_UP=$4
SIZE_DWN=$5

if [[ "$6" == "void" ]];
then
	BKGRND=""
else
	BKGRND="$6"
fi

if [[ "$7" == "void" ]];
then
	MASK=""
else
	MASK="$7"
fi

if [[ "$8" == "void" ]];
then
	RNA=""
else
	RNA="$8"
fi

if [[ "${10}" == "void" ]];
then
	NORM=""
else
	NORM="${10}"
fi

THREADS=${11}
NARROW_DIR=${12}

if [[ "${13}" == "binomial" ]];
then
	SCORING=""          # negative binomial distribution
else
	SCORING="-h"        # hypergeometric distribution
fi

BAYES_TYPE_FOR_DUAL=${14}
BAYES_DIR=${15}

t=$((16 + 5 * $SLURM_ARRAY_TASK_ID));
c=$((17 + 5 * $SLURM_ARRAY_TASK_ID));
o=$((18 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((19 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((20 + 5 * $SLURM_ARRAY_TASK_ID));

trt=${!t}
ctl=${!c}
otp=${!o}

if [[ "$F_TYPE" == ".bayesPeaks.bed" ]];
then
	SOURCE_DIR=$BAYES_DIR
else
	SOURCE_DIR=$NARROW_DIR
fi

homerfile="${SOURCE_DIR}${otp}${F_TYPE}"
SIZE_ARG=""
OUTDIR="${HM_DIR}$otp"
mkdir -p $OUTDIR

mkdir -p "${OUTDIR}/enrich"
touch "${OUTDIR}/enrich/selected.motifs"

if [[ "$9" == "void" ]];
then
	KNOWN=""
else
	if [ -f "${OUTDIR}/enrich/$9" ]
	then
		OUTDIR="${OUTDIR}/enrich"
		KNOWN="-nomotif -mknown ${OUTDIR}/$9"
	fi
fi

echo  
echo
echo homer name is $homerfile
echo 
echo 
echo known_motif_file is $KNOWN
echo genome is $GENOME
echo hm_dir is $HM_DIR
echo filetype is $F_TYPE
echo sizeup is $SIZE_UP
echo sizedown is $SIZE_DWN
echo bkgrnd_flag is $BKGRND
echo mask_flag is $MASK
echo rna_flag is $RNA
echo normalization is $NORM
echo no_threads is $THREADS
echo scoring is $SCORING
echo src_dir is $SOURCE_DIR
echo 

if [[ "$SIZE_UP" == "given" ]];
then
	SIZE_ARG="-size given"
else
	if [[ "$SIZE_DWN" == "given" ]];
	then
		SIZE_ARG="-size given"
	else
		SIZE_ARG="-size -${SIZE_DWN},${SIZE_UP}"
	fi
fi

if [[ "$MYENV" == "crick" ]]
then
	PREPARSE="/group/nordlab/libraries/externalDatasets/$GENOME"
else
	PREPARSE="/share/nordlab/libraries/externalDatasets/$GENOME"
fi

# Perform motif discovery on peak file
findMotifsGenome.pl $homerfile $GENOME $OUTDIR $SIZE_ARG $BKGRND $MASK $RNA $NORM $SCORING \
					-p $THREADS $KNOWN -preparsedDir $PREPARSE

if [[ ! "$BAYES_TYPE_FOR_DUAL" == "null" ]];
then
	mkdir -p "${OUTDIR}.tmp"
	cd "${OUTDIR}.tmp"
	homerfile="${BAYES_DIR}${otp}${BAYES_TYPE_FOR_DUAL}"
	findMotifsGenome.pl $homerfile $1 $OUTDIR $SIZE_ARG $BKGRND $MASK $RNA $NORM $SCORING \
						-p $THREADS $KNOWN -preparsedDir $PREPARSE
	cd ..
	rm -rf "${OUTDIR}.tmp"
fi

echo  
echo
echo 
echo Completed execution on $(date)
echo 
echo 
echo 

# changelog
# 1.0.0	first working edition