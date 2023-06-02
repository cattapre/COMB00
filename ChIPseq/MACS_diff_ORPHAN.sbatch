#!/bin/bash
# MACS_diff.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=MACS2diff

echo
echo 
echo MACS_diff.sbatch version 1.0.0
echo ... started running on $(date)
echo  
echo  

echo
echo This job is running on
/bin/hostname
echo
echo


# Load Modules, and Run Programs
# module load macs2/2.1.0.20140616
module load macs2
export PYTHON_EGG_CACHE='/share/nordlab/users/rinaldo/chipseq/code/cache'

echo 

BROAD_DIR=$1
NARROW_DIR=${2}
EXTSIZE=$3
DIFFCALL_DIR=$4
MEMORY=$5
CUTOFF=$6
MIN_LENGTH=$7
MAX_GAP=$8
DEPTH1=$9
DEPTH2=${10}
LOG_DIR=${11}
cold=${12}
hot=${13}

echo
echo
echo $BROAD_DIR
echo $NARROW_DIR
echo $EXTSIZE 
echo $DIFFCALL_DIR
echo $MEMORY
echo $CUTOFF
echo $MIN_LENGTH
echo $MAX_GAP
echo $DEPTH1
echo $DEPTH2
echo $LOG_DIR
echo $cold
echo $hot
echo

origterm=""
pile=".no_model_treat_pileup.bdg"
lambda=".no_model_control_lambda.bdg"
xls=".no_model_peaks.xls"
TYPES="$NARROW_DIR $BROAD_DIR"
# TYPES="$NARROW_DIR"
# TYPES="$BROAD_DIR"

SUFF="\-vs\.Input"
suff=""

for k in $TYPES
do
	cd $k

	hotprefix="${hot/$SUFF/$suff}"
	coldprefix="${cold/$SUFF/$suff}"
	
	echo 
	echo 

	TRT=(`egrep "tags after filtering in treatment" "$k$cold$xls"`)
	CTL=(`egrep "tags after filtering in treatment" "$k$cold$xls"`)
	trt=${TRT[-1]}
	ctl=${CTL[-1]}

	D1=$( [ $trt -le $ctl ] && echo "$trt" || echo "$ctl" )

	TRT=(`egrep "tags after filtering in treatment" "$k$hot$xls"`)
	CTL=(`egrep "tags after filtering in treatment" "$k$hot$xls"`)
	trt=${TRT[-1]}
	ctl=${CTL[-1]}

	D2=$( [ $trt -le $ctl ] && echo "$trt" || echo "$ctl" )

	echo "D1: $D1   D2: $D2"
	echo
	echo

	cd $DIFFCALL_DIR

	KOchecker="KO"
	Mutchecker="Mut"

	if [ "${cold/$KOchecker}" = "$cold" ]
	then
		if [ "${cold/$Mutchecker}" = "$cold" ]
		then
			if [ "$k" == "$BROAD_DIR"  ]
			then
				oprefix="${coldprefix}_vs_${hotprefix}_diff_BROAD"
			else
				oprefix="${coldprefix}_vs_${hotprefix}_diff_NARROW"
			fi

			echo "$oprefix"
		
			macs2 bdgdiff --t1 "$k$cold$pile" --c1 "$k$cold$lambda" \
						  --t2 "$k$hot$pile" --c2 "$k$hot$lambda" \
						  --d1 "$D1" --d2 "$D2" -g $MAX_GAP --o-prefix "$oprefix" \
						  --outdir "$DIFFCALL_DIR" -C $CUTOFF -l $MIN_LENGTH 
		else
			if [ "$k" == "$BROAD_DIR"  ]
			then
				oprefix="${hotprefix}_vs_${coldprefix}_diff_BROAD"
			else
				oprefix="${hotprefix}_vs_${coldprefix}_diff_NARROW"
			fi

			echo "$oprefix"
		
			macs2 bdgdiff --t1 "$k$hot$pile" --c1 "$k$hot$lambda" \
						  --t2 "$k$cold$pile" --c2 "$k$cold$lambda" \
						  --d1 "$D2" --d2 "$D1" -g $MAX_GAP --o-prefix "$oprefix" \
						  --outdir "$DIFFCALL_DIR" -C $CUTOFF -l $MIN_LENGTH 
		fi
	else
		if [ "$k" == "$BROAD_DIR"  ]
		then
			oprefix="${hotprefix}_vs_${coldprefix}_diff_BROAD"
		else
			oprefix="${hotprefix}_vs_${coldprefix}_diff_NARROW"
		fi

		echo "$oprefix"
		
		macs2 bdgdiff --t1 "$k$hot$pile" --c1 "$k$hot$lambda" \
					  --t2 "$k$cold$pile" --c2 "$k$cold$lambda" \
					  --d1 "$D2" --d2 "$D1" -g $MAX_GAP --o-prefix "$oprefix" \
					  --outdir "$DIFFCALL_DIR" -C $CUTOFF -l $MIN_LENGTH 
	fi

done

wait

echo  
echo
echo 
echo Completed execution with peaks on $(date)
echo 
echo 
echo 



# changelog
# 1.2.0 changed parameter uptake
# 1.1.0 added narrow peak differential calls
# 1.0.0	first working edition