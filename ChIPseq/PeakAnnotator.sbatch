#!/bin/bash
# PeakAnnotator.sh by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=PeakAnnotator

# Load Modules, and Run Programs

echo 
echo PeakAnnotator.sbatch version 1.1.0
echo ... started running on $(date)
echo
echo 
module load java
echo

echo
echo This job is running on
/bin/hostname
echo
echo


t=$((10 + 5 * $SLURM_ARRAY_TASK_ID));
c=$((11 + 5 * $SLURM_ARRAY_TASK_ID));
o=$((12 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((13 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((14 + 5 * $SLURM_ARRAY_TASK_ID));

MACS_o=${!o}     
MACS_sple=${!n}		
MACS_ctrl=${!m}		

term=".no_model_peaks.narrowPeak"
termG=".no_model_peaks.gappedPeak"

CHIP="$3$MACS_sple$term"
CHIP_B="$6$MACS_sple$term2"
CONTROL="$3$MACS_ctrl$term"
CTRL_B="$6$MACS_ctrl$term2"

OUTDIR="$2"
#OUTDIR=`echo "$OUTDIR" | rev | cut -c 2- | rev`

GENOME=$4
GENOMEFULL=$9

echo
cd $1

if [[ "$8" == "MACS2" ]];
then
    term2=".no_model_peaks.broadPeak"
    BROAD="$6$MACS_o$term2"
    NARROW="$3$MACS_o$term"
    GAPPED="$6$MACS_o$termG"
    
#     NEWNARROW="${NARROW}.tmp"
#     cp $NARROW $NEWNARROW
#     sed -i "s/$/\t0\t0/" $NEWNARROW
#     NEWNARROW=`cat $NEWNARROW | awk '{ $7 = "0"; $8 = "0"; $9 = "0"; print }'`
#     
#     NEWBROAD="${BROAD}.tmp"
#     cp $BROAD $NEWBROAD
#     sed -i "s/$/\t0\t0\t0/" $NEWBROAD
#     NEWBROAD=`cat $NEWBROAD | awk '{ $7 = "0"; $8 = "0"; $9 = "0"; print }'`
#     
#     java -jar -Xmx512m PeakAnnotator.jar -u NDG -p $NEWBROAD -a $GENOME -o $OUTDIR -g all -r $5
#     java -jar -Xmx512m PeakAnnotator.jar -u TSS -p $NEWBROAD -a $GENOME -o $OUTDIR -g all -r $5
#     
#     java -jar -Xmx512m PeakAnnotator.jar -u NDG -p $NEWNARROW -a $GENOME -o $OUTDIR -g all -r $5
#     java -jar -Xmx512m PeakAnnotator.jar -u TSS -p $NEWNARROW -a $GENOME -o $OUTDIR -g all -r $5
    
    java -jar -Xmx512m PeakAnnotator.jar -u NDG -p $GAPPED -a $GENOME -o $OUTDIR -g all -r $5
    java -jar -Xmx512m PeakAnnotator.jar -u TSS -p $GAPPED -a $GENOME -o $OUTDIR -g all -r $5
    
#     rm $NEWNARROW; rm $NEWBROAD;
else
    if [[ "$8" == "BayesPeak" ]];
    then
        term2=".bayesPeaks.bed"
        BROAD="$7$MACS_o$term2"
        java -jar -Xmx512m PeakAnnotator.jar -u NDG -p $BROAD -a $GENOME -o $OUTDIR -g all -r $5
        java -jar -Xmx512m PeakAnnotator.jar -u TSS -p $BROAD -a $GENOME -o $OUTDIR -g all -r $5
    else
        if [[ "$8" == "BOTH" ]];
        then
            term2=".no_model_peaks.broadPeak"
            BROAD="$6$MACS_o$term2"
    		NARROW="$3$MACS_o$term"
			GAPPED="$6$MACS_o$termG"
            
			java -jar -Xmx512m PeakAnnotator.jar -u NDG -p $GAPPED -a $GENOME -o $OUTDIR -g all -r $5
			java -jar -Xmx512m PeakAnnotator.jar -u TSS -p $GAPPED -a $GENOME -o $OUTDIR -g all -r $5
    
            term2=".bayesPeaks.bed"
            BROAD="$7$MACS_o$term2"
            java -jar -Xmx512m PeakAnnotator.jar -u NDG -p $BROAD -a $GENOME -o $OUTDIR -g all -r $5
            java -jar -Xmx512m PeakAnnotator.jar -u TSS -p $BROAD -a $GENOME -o $OUTDIR -g all -r $5
        fi
    fi
fi

echo "... starting gene symbol annotation"

perl /share/nordlab/users/rinaldo/chipseq/code/Gene_Annotator.pl $OUTDIR "${MACS_o}.no_model" $GENOMEFULL


echo
echo

CLOSEST="_closest_genes.bed"
bgraph="_coverage.bdg"
bwig="_coverage.bw"
noovlp="_gene_NO_overlap.bed"
ovlp="_gene_overlap.bed"
annot="_gene_overlap.bed_annotated"

term=".no_model_peaks.gappedPeak"
mv "$6$MACS_o$term$CLOSEST" $OUTDIR
mv "$6$MACS_o$term$bgraph" $OUTDIR
mv "$6$MACS_o$term$bwig" $OUTDIR
mv "$6$MACS_o$term$noovlp" $OUTDIR
mv "$6$MACS_o$term$ovlp" $OUTDIR
mv "$6$MACS_o$term$annot" $OUTDIR

# term="_peaks.broadPeak"
# mv "$6$MACS_o$term$CLOSEST" $OUTDIR
# mv "$6$MACS_o$term$bgraph" $OUTDIR
# mv "$6$MACS_o$term$bwig" $OUTDIR
# mv "$6$MACS_o$term$noovlp" $OUTDIR
# mv "$6$MACS_o$term$ovlp" $OUTDIR
# mv "$6$MACS_o$term$annot" $OUTDIR
# 
# term="_peaks.narrowPeak"
# mv "$3$MACS_o$term$CLOSEST" $OUTDIR
# mv "$3$MACS_o$term$bgraph" $OUTDIR
# mv "$3$MACS_o$term$bwig" $OUTDIR
# mv "$3$MACS_o$term$noovlp" $OUTDIR
# mv "$3$MACS_o$term$ovlp" $OUTDIR
# mv "$3$MACS_o$term$annot" $OUTDIR

echo
echo

echo
echo Completed execution on $(date)
echo

