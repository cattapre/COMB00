#!/bin/bash
# MACS_exp.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=MACS2exp

echo
echo 
echo MACS_exp.sbatch version 1.4.0
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
	MYENV="cabernet"
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
	#module load python-libs/2.7.6-ubuntu
	module load macs2
	export PYTHON_EGG_CACHE='/share/nordlab/users/rinaldo/chipseq/code/cache'
fi


t=$((12 + 5 * $SLURM_ARRAY_TASK_ID));
c=$((13 + 5 * $SLURM_ARRAY_TASK_ID));
o=$((14 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((15 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((16 + 5 * $SLURM_ARRAY_TASK_ID));

MACS_t=${!t}
MACS_c=${!c}
MACS_o=${!o}

MACS_sple=${!n}
MACS_ctrl=${!m}

MEM=${10}
LOG=${11}

PAIR="_PE"

if [ "${MACS_o/$PAIR}" = "$MACS_o" ];
then
	FORMAT="BAM"
else
	FORMAT="BAMPE"
fi


echo  
echo
echo macs2 name is $MACS_o
echo 

term="_model.r"
modelterm=".model"
nomodelterm=".no_model"

PVALUE=$( echo "$6 > 0" | bc -l )
EXTSIZE=$( echo "$8 > 0" | bc -l )

if [[ "$3" == 'rn' ]];
then
	GENSIZE='1.87e9'
	# source: doi: 10.1038/srep33736
else
	GENSIZE=$3
fi

echo "PVALUE: $PVALUE"
echo "EXTSIZE: $EXTSIZE"
echo "GENSIZE: $GENSIZE"

if [ "$MACS_c" != "null" ]
then
	# BROAD PEAK CALL
	cd $1

	if [[ $EXTSIZE -eq 1 ]];
	then
		if [[ $PVALUE -eq 1 ]];
		then
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$nomodelterm" --broad -f $FORMAT \
				-g $GENSIZE --broad-cutoff $4 --bdg -p $6 $7 --nomodel --extsize=$8 --outdir $1 
			echo
		else
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$nomodelterm" --broad -f $FORMAT \
				-g $GENSIZE --broad-cutoff $4 --bdg -q $5 $7 --nomodel --extsize=$8 --outdir $1 
			echo
		fi
	else
		if [[ $PVALUE -eq 1 ]];
		then
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$modelterm" --broad -f $FORMAT \
				-g $GENSIZE --broad-cutoff $4 --bdg -p $6 $7 --outdir $1 
			echo
		else
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$modelterm" --broad -f $FORMAT \
				-g $GENSIZE --broad-cutoff $4 --bdg -q $5 $7 --outdir $1 
			echo
		fi
	fi	


	# NARROW PEAK CALL
	cd $2

	if [[ $EXTSIZE -eq 1 ]];
	then
		if [[ $PVALUE -eq 1 ]];
		then
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$nomodelterm" -f $FORMAT -g $GENSIZE \
				--call-summits --bdg -p $6 $7 --nomodel --extsize=$8 --outdir $2 
			echo
			
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$modelterm" -f $FORMAT -g $GENSIZE \
				--call-summits --bdg -p $6 --outdir $2 
			echo
		
			# Testing the --cutoff-analysis setting
	# 		srun --mem=$MEM -p gc128,gc256,gc512 -o "$LOG$MACS_o$nomodelterm"_N%j.out -e "$LOG$MACS_o$nomodelterm"_N%j.out \
	# 		macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$nomodelterm".cutoff -f $FORMAT -g $GENSIZE \
	# 			--bdg --cutoff-analysis $7 --nomodel --extsize=$8 --outdir $2
	# 			
	# 		srun --mem=$MEM -p gc128,gc256,gc512 -o "$LOG$MACS_o$modelterm"_N%j.out -e "$LOG$MACS_o$modelterm"_N%j.out \
	# 		macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$modelterm".cutoff -f $FORMAT -g $GENSIZE \
	# 			--bdg --cutoff-analysis --outdir $2
		
		else
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$nomodelterm" -f $FORMAT -g $GENSIZE \
				--call-summits --bdg -q $5 $7 --nomodel --extsize=$8 --outdir $2 
			echo
			
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$modelterm" -f $FORMAT -g $GENSIZE \
				--call-summits --bdg -q $5 --outdir $2 
			echo
		
			# Testing the --cutoff-analysis setting
	# 		srun --mem=$MEM -p gc128,gc256,gc512 -o "$LOG$MACS_o$nomodelterm"_N%j.out -e "$LOG$MACS_o$nomodelterm"_N%j.out \
	# 		macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$nomodelterm".cutoff -f $FORMAT -g $GENSIZE \
	# 			--bdg --cutoff-analysis $7 --nomodel --extsize=$8 --outdir $2
	# 			
	# 		srun --mem=$MEM -p gc128,gc256,gc512 -o "$LOG$MACS_o$modelterm"_N%j.out -e "$LOG$MACS_o$modelterm"_N%j.out \
	# 		macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$modelterm".cutoff -f $FORMAT -g $GENSIZE \
	# 			--bdg --cutoff-analysis --outdir $2
		
		fi
	else
		if [[ $PVALUE -eq 1 ]];
		then
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$modelterm" -f $FORMAT -g $GENSIZE \
				--call-summits --bdg -p $6 $7 --outdir $2 
			echo
		else
			macs2 callpeak -t $MACS_t -c $MACS_c -n "$MACS_o$modelterm" -f $FORMAT -g $GENSIZE \
				--call-summits --bdg -q $5 $7 --outdir $2 
			echo
		fi
	fi	

# 	wait 

	echo  

	echo  
	echo  ... completed broad peak call on $(date)
	echo  
	echo  

	echo  
	echo  ... completed narrow peak call on $(date)


	if [ -f $2$MACS_o$nomodelterm$term ];
	then
		Rscript $2$MACS_o$nomodelterm$term
		echo  
		echo 
		echo Completed pdf generation of narrow peak files on $(date)
		echo 

	fi

	if [ -f $2$MACS_o$modelterm$term ];
	then
		Rscript $2$MACS_o$modelterm$term
		echo  
		echo 
		echo Completed pdf generation of narrow peak files on $(date)
		echo 

	fi

	cd $1

	if [ -f $1$MACS_o$nomodelterm$term ];
	then
		Rscript $1$MACS_o$nomodelterm$term
		echo  
		echo 
		echo Completed pdf generation of broad peak files on $(date)
		echo 

	fi

	if [ -f $1$MACS_o$modelterm$term ];
	then
		Rscript $1$MACS_o$modelterm$term
		echo  
		echo 
		echo Completed pdf generation of broad peak files on $(date)
		echo 

	fi

	wait
fi
echo  
echo
echo 
echo Completed execution on $(date)
echo 
echo 
echo 

# changelog
# 1.4.0	introduced peak calling on pair-ended reads
# 1.3.0	introduced additional parameter settings
# 1.2.0	introduced soft coding for peak calling parameters
# 1.1.0	included single ChIP, inputDNA and NegCtrl peak calling
# 1.0.0	first working edition