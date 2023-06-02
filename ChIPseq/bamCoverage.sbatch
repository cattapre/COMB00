#!/bin/bash
# bamCoverage.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=bamCoverage

echo
echo 
echo bamCoverage.sbatch version 1.0.1
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
if [[ "$MYENV" == "crick" ]]
then
	module load bio
else
	module load deeptools/2.5.4
	module load samtools
	export PYTHON_EGG_CACHE='/share/nordlab/users/rinaldo/chipseq/code/cache'
fi


t=$((10 + $SLURM_ARRAY_TASK_ID));

MACS_t=${!t}

BROAD=$1
NARROW=$2
DEEPDIR=$3
THREADS=$4
EXT=$5
IGN_DUP=$6
NORMALIZ=$7
GENOME=$8
BLKLST=$9

echo
echo $BROAD
echo $NARROW
echo $DEEPDIR
echo $THREADS
echo $EXT
echo $IGN_DUP
echo $NORMALIZ
echo $GENOME
echo

echo
echo "sample is $MACS_t"
echo 

if [[ "$NORMALIZ" == "normalizeTo1x" ]]
then
	NORMAL_SUFF="normalized1x"
elif [[ "$NORMALIZ" == "normalizeUsingRPKM" ]]
then
	NORMAL_SUFF="normalizedRPKM"
fi

if [[ "$GENOME" == "mm9" ]]
then
	GENFORNORM="2150570000"
elif [[ "$GENOME" == "hg19" ]]
then
	GENFORNORM="2451960000"
elif [[ "$GENOME" == "rn5" ]]
then
	GENFORNORM="2572853723"
	# http://genomewiki.ucsc.edu/index.php/Mm10_Genome_size_statistics
else
	echo
	echo
	echo
	echo " ... Curently supported genomes are mm9, hg19 and rn5"
	echo
	echo
fi

if [[ "$BLKLST" == "null" ]]
then
	BLKLST=""
else
	BLKLST="-bl $BLKLST"
fi

if [[ "$NORMALIZ" == "normalizeTo1x" ]]
then
	NORMALIZ="${NORMALIZ} $GENFORNORM"
fi

if [[ "$EXT" > 0 ]]
then
	EXT="--extendReads $EXT"
else
	EXT=''
fi

INDEX='.bai'

echo
echo

if [ "${MACS_t}" != "" ]
then
	if [[ ! -f "${MACS_t}$INDEX" ]]
	then
		samtools index $MACS_t
	else
		echo
		echo "... indexing $MACS_t not necessary"
		echo
	fi

	cd $DEEPDIR

	echo "... working on $MACS_t coverage file"
	echo 

	if [[ "$NORMALIZ" == "none" ]]
	then
		OUTFILE="${MACS_t}.coverage.bw"
		if [[ "$IGN_DUP" == "no" ]]
		then
			bamCoverage -b "${MACS_t}" -o "$OUTFILE" -of bigwig -p "$THREADS" -bs 10 $EXT $BLKLST
		else
			bamCoverage -b "${MACS_t}" -o "$OUTFILE" -of bigwig -p "$THREADS" -bs 10 $EXT "$IGN_DUP" $BLKLST
		fi
	else
		OUTFILE="${MACS_t}.${NORMAL_SUFF}.coverage.bw"
		if [[ "$IGN_DUP" == "no" ]]
		then
			bamCoverage -b "${MACS_t}" -o "$OUTFILE" -of bigwig -p "$THREADS" -bs 10 $EXT --"$NORMALIZ" $BLKLST
		else
			echo "bamCoverage -b ${MACS_t} -o $OUTFILE -of bigwig -p $THREADS -bs 10 $EXT $IGN_DUP --$NORMALIZ $BLKLST"
			bamCoverage -b "${MACS_t}" -o "$OUTFILE" -of bigwig -p "$THREADS" -bs 10 $EXT "$IGN_DUP" --"$NORMALIZ" $BLKLST
		fi
	fi

	mv $OUTFILE ../deeptools/
fi

echo "... done with ${OUTFILE}"
echo 
echo 
          


echo  
echo
echo 
echo Completed execution on $(date)
echo 
echo 
echo 

# changelog
# 1.0.1	corrected some stability issues
# 1.0.0	first working edition