#!/bin/bash

#SBATCH --job-name=Bayes

echo
echo BayesPeak.sh version 1.0.0
echo
echo ... started running on $(date)
echo
echo

echo
echo This job is running on
/bin/hostname
echo
echo

t=$((5 + 5 * $SLURM_ARRAY_TASK_ID))
c=$((6 + 5 * $SLURM_ARRAY_TASK_ID))
o=$((7 + 5 * $SLURM_ARRAY_TASK_ID))
n=$((8 + 5 * $SLURM_ARRAY_TASK_ID))
m=$((9 + 5 * $SLURM_ARRAY_TASK_ID))

chip=${!t}

PAIR="_PE"

if [ "${chip/$PAIR}" = "$chip" ];
then
    SUFF="srt.bam"
else
    SUFF="bam"
fi

suff="bed"

chip=${chip/$SUFF/$suff}
echo $chip

ctrl=${!c}
ctrl=${ctrl/$SUFF/$suff}
echo $ctrl

name=${!o}
name="$name.bayesPeaks"
echo $name

OUTDIR=${!o}

threads="$1"
bayes_dir="$2$OUTDIR"
echo $bayes_dir

lambda="$3"
echo "lambda is $lambda"
calls="$4"

mkdir -p $bayes_dir

sleep $SLURM_ARRAY_TASK_ID

/usr/bin/Rscript BayesPeak.R $threads $bayes_dir $lambda $calls $chip $ctrl $name

perl Bayes2Bed.pl $bayes_dir $name 

cd $bayes_dir

shopt -s extglob

mv !(clustering.txt) ../

echo
echo
echo
echo Completed execution on $(date)
echo
