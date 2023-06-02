#!/bin/bash
# Homer.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=Homer_annot

echo
echo 
echo Homer_annot.sbatch version 1.0.0
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

# PREFIX="-"

GENOME=$1
HM_DIR=$2
F_TYPE=$3
SIZE_UP=$4
SIZE_DWN=$5

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

HM_BIN=$6
NARROW_DIR=$7
KNOWN=$8
THREADS=$9

SOURCE_DIR=$NARROW_DIR

t=$((10 + 5 * $SLURM_ARRAY_TASK_ID));
c=$((11 + 5 * $SLURM_ARRAY_TASK_ID));
o=$((12 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((13 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((14 + 5 * $SLURM_ARRAY_TASK_ID));

trt=${!t}
ctl=${!c}
otp=${!o}

if [[ "$KNOWN" == "void" ]];
then
	OUTDIR="${HM_DIR}$otp"
	allmotiffs="${OUTDIR}/homerMotifs.all.motifs"
else
	OUTDIR="${HM_DIR}${otp}/enrich"
	mkdir -p $OUTDIR
	allmotiffs="${OUTDIR}/$KNOWN"
fi

homerfile="${SOURCE_DIR}${otp}$F_TYPE"
homerbed="${OUTDIR}/${otp}.motif.instances.bed"

ANNOT_HOMER="${OUTDIR}/${otp}.motif.counts"


echo  
echo
echo homer name is $homerfile
echo allmottifs is $allmotiffs
echo homerbed is $homerbed 
echo Annotationfile is $ANNOT_HOMER
echo 
echo known is $KNOWN
echo genome is $GENOME
echo hm_dir is $HM_DIR
echo filetype is $F_TYPE
echo sizeup is $SIZE_UP
echo sizedown is $SIZE_DWN
echo binsize is $HM_BIN
echo outdir is $OUTDIR
echo no_threads is $THREADS
echo src_dir is $SOURCE_DIR
echo 

# Perform known motif analysis on all subset bed files
# (not implemented yet)


# Map individual motifs of interest to peak bedfiles
# motif counts for all motifs
mkdir -p "${NARROW_DIR}${otp}${F_TYPE}.tmp"
cd "${NARROW_DIR}${otp}${F_TYPE}.tmp"
annotatePeaks.pl $homerfile $GENOME -m $allmotiffs -nmotifs \
				 -cpu $THREADS > $ANNOT_HOMER


# motif co-occurence for all motifs
MATRIX_HOMER="${ANNOT_HOMER}.cooccur"

annotatePeaks.pl $homerfile $GENOME -m $allmotiffs -matrix $MATRIX_HOMER -cpu $THREADS > $MATRIX_HOMER.out


# motif histogram for primary motifs
annotatePeaks.pl $homerfile $GENOME -m $allmotiffs $SIZE_ARG -hist $HM_BIN \
				 -cpu $THREADS > $ANNOT_HOMER.hist.txt
annotatePeaks.pl $homerfile $GENOME -m $allmotiffs $SIZE_ARG -hist $HM_BIN -ghist \
				 -cpu $THREADS > $ANNOT_HOMER.heatmap.hist.txt
annotatePeaks.pl $homerfile $GENOME -m $allmotiffs $SIZE_ARG -hist $HM_BIN -CpG \
				 -cpu $THREADS > $ANNOT_HOMER.CpG.hist.txt
annotatePeaks.pl $homerfile $GENOME -m $allmotiffs -mbed $homerbed \
				 -cpu $THREADS > "${homerbed}.annotated"
tail -n +2 $homerbed | sort -k1,1 -k2,2n > "${homerbed}.sorted"
# /usr/bin/Rscript formatHomerBed.R "${homerbed}.sorted"

annotatePeaks.pl tss $GENOME -m $allmotiffs $SIZE_ARG -hist $HM_BIN \
				 -cpu $THREADS > $ANNOT_HOMER.tss.hist.txt


# GO analysis

GO_dir="${OUTDIR}"

mkdir -p $GO_dir

annotatePeaks.pl $homerfile $GENOME -m $allmotiffs -nmotifs -go $GO_dir \
				 -genomeOntology $GO_dir -cpu $THREADS > $MATRIX_HOMER.GO

cd ..
rm -rf "${SOURCE_DIR}${otp}${F_TYPE}.tmp"
echo  
echo
echo 
echo Completed execution on $(date)
echo 
echo 
echo 

# changelog
# 1.0.0	first working edition