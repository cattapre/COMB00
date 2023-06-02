#!/bin/bash
# BedTools.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=BedTools		   			# A single job name for the array

# Load Modules, and Run Programs

echo 
echo BedTools.sbatch version 1.1.0
echo 

echo ... started running on $(date)
echo
echo 
module load bedtools2/2.21.0
module load samtools/1.3
echo

echo
echo This job is running on
/bin/hostname
echo
echo


t=$(( 7 + 5 * $SLURM_ARRAY_TASK_ID));
c=$(( 8 + 5 * $SLURM_ARRAY_TASK_ID));
o=$(( 9 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((10 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((11 + 5 * $SLURM_ARRAY_TASK_ID));

MACS_o=${!o}
echo  

Bed_NAME=${!o}
TRT_file=${!t}
CTL_file=${!c}

TRT_name=${Bed_NAME/\_trim\-vs\.Input/\_ChIP}

Bed_suffix="_gene_overlap.bed"
ANNO="_annotated"
No_overlap="_gene_NO_overlap.bed"
Closest="_closest_genes.bed"
BedGraph="_coverage.bdg"
BigWig="_coverage.bw"

echo 
echo 
echo bedfile_name is $Bed_NAME
echo  


# Broad Peaks
cd $1

term=".no_model_peaks.gappedPeak"

INPUT_Bed=$1$Bed_NAME$term

if [[ -f "$INPUT_Bed" ]];
then
	echo  
	echo ... running on broad peak $INPUT_Bed
	echo 

	sort -k1,1 -k2,2n $INPUT_Bed > "${INPUT_Bed}.sorted"
	echo ... finished sorting on $(date)
	echo 
	echo 

	bedtools intersect -wb -a "${INPUT_Bed}.sorted" -b $3 > "${INPUT_Bed}$Bed_suffix"
	echo ... finished intersect on $(date)
	echo 

# 	bedtools intersect -wb -a "${INPUT_Bed}.sorted" -b $3 -v > "${INPUT_Bed}$No_overlap"
# 	echo ... finished intersect -v  on $(date)
# 	echo 

	bedtools closest -a "${INPUT_Bed}.sorted" -b $3 > "${INPUT_Bed}$Closest"
	echo ... finished closest on $(date)
	echo 

	samtools sort -o "${TRT_file}.sorted" -O bam -@ $6 "${TRT_file}"
	
	
	bedtools genomecov -bg -ibam "${TRT_file}.sorted" -g $4 > "$1${TRT_name}$BedGraph"
	echo ... finished genomecov on $(date)
	echo 

	bedtools annotate -both -i "${INPUT_Bed}.sorted" -files $3 > "${INPUT_Bed}${Bed_suffix}$ANNO"
	echo ... finished annotate on $(date)
	echo 

	/share/nordlab/codes/chipseq/bedGraphToBigWig "$1${TRT_name}$BedGraph" $4 "$1$TRT_name$BigWig"
	echo ... $TRT_name$BedGraph converted to $TRT_name$BigWig on $(date)
fi

# Narrow Peaks 
cd $2

term2=".no_model_peaks.narrowPeak"

INPUT_Bed2="$2${Bed_NAME}$term2"

if [[ -f "$INPUT_Bed2" ]];
then
	echo  
	echo
	echo  
	echo ... running on narrow peak $INPUT_Bed2
	echo 

	sort -k1,1 -k2,2n $INPUT_Bed2 > "${INPUT_Bed2}.sorted"
	echo ... finished sorting on $(date)
	echo 
	echo 

	bedtools intersect -wb -a "${INPUT_Bed2}.sorted" -b $3 > "${INPUT_Bed2}$Bed_suffix"
	echo ... finished intersect on $(date)
	echo 

# 	bedtools intersect -wb -a "${INPUT_Bed2}.sorted" -b $3 -v > "${INPUT_Bed2}$No_overlap"
# 	echo ... finished intersect -v  on $(date)
# 	echo 

	bedtools closest -a "${INPUT_Bed2}.sorted" -b $3 > "${INPUT_Bed2}$Closest"
	echo ... finished closest on $(date)
	echo  

	bedtools genomecov -bg -i "${INPUT_Bed2}.sorted" -g $4 > "${INPUT_Bed2}$BedGraph"
	echo ... finished genomecov on $(date)
	echo 

	bedtools annotate -both -i "${INPUT_Bed2}.sorted" -files $3 > "${INPUT_Bed2}${Bed_suffix}$ANNO"
	echo ... finished annotate on $(date)
	echo 

	/share/nordlab/codes/chipseq/bedGraphToBigWig "${INPUT_Bed2}$BedGraph" $4 "${INPUT_Bed2}$BigWig"
	echo ... $INPUT_Bed2$BedGraph converted to $INPUT_Bed2$BigWig on $(date)
	echo 

fi

# BayesPeaks
cd $5

term2=".bayesPeaks.bed"

INPUT_Bed2="$5${Bed_NAME}$term2"

if [[ -f "$INPUT_Bed2" ]];
then
	echo  
	echo
	echo  
	echo ... running on BayesPeaks $INPUT_Bed2
	echo 

	sort -k1,1 -k2,2n $INPUT_Bed2 > "${INPUT_Bed2}.sorted"
	echo ... finished sorting on $(date)
	echo 
	echo 

	bedtools intersect -wb -a "${INPUT_Bed2}.sorted" -b $3 > "${INPUT_Bed2}$Bed_suffix"
	echo ... finished intersect on $(date)
	echo 

# 	bedtools intersect -wb -a "${INPUT_Bed2}.sorted" -b $3 -v > "${INPUT_Bed2}$No_overlap"
# 	echo ... finished intersect -v  on $(date)
# 	echo 

	bedtools closest -a "${INPUT_Bed2}.sorted" -b $3 > "${INPUT_Bed2}$Closest"
	echo ... finished closest on $(date)
	echo  

	bedtools genomecov -bg -i "${INPUT_Bed2}.sorted" -g $4 > "${INPUT_Bed2}$BedGraph"
	echo ... finished genomecov on $(date)
	echo 

	bedtools annotate -both -i "${INPUT_Bed2}.sorted" -files $3 > "${INPUT_Bed2}${Bed_suffix}$ANNO"
	echo ... finished annotate on $(date)
	echo 

# 	/share/nordlab/codes/chipseq/bedGraphToBigWig "${INPUT_Bed2}$BedGraph" $4 "${INPUT_Bed2}$BigWig"
# 	echo ... $INPUT_Bed2$BedGraph converted to $INPUT_Bed2$BigWig on $(date)
# 	echo 

fi

echo  
echo
echo 
echo Completed execution on $(date)
echo 



# changelog
# 1.1.0 Implementation of BayesPeaks annotation
# 1.0.0 first production release
# 0.2.0	included annotate call
# 0.1.0	first working edition

