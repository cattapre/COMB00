#!/bin/bash -l

#SBATCH --job-name=deepTools
#SBATCH --time=3-00:00
#SBATCH --mem=32G
#SBATCH --partition=high
#SBATCH --cpus-per-task=16

echo
echo This job is running on
/bin/hostname
echo "crick"
echo
echo

module load bio3

echo
plotHeatmap --version
echo

SENSIT=$1       # pvalue
REGION=$2       # BG, CX or combined
BGSCOPE=$3      # 1 = all, 2 = excluding drivers, 3 = drivers
positionScope=$4  # 1 = all, 2 = distal, 3 = proximal 4 = directed
clusterID=$5

PROJECTID="COMB00"

PLATFORM="group"
THREADS=16
BLKLST="/${PLATFORM}/nordlab/libraries/externalDatasets/mm10/mm10.blacklist.bed"

TEMPDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/temp.folder"
BWADIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/bwa/"
BEDDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/data/bedfiles/"
DEEPTOOLS="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/deeptools/"

myWindow=1000
echo
echo $#
echo

if [ $positionScope == 1 ]
then
    if [[ $# -gt 4 ]]
    then
        myFile="${BEDDIR}FullMerge.p.${SENSIT}.${REGION}.noReps.bed.purged"
        DEEPDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/plots/DEEP/${REGION}.1.all.purged/"
    else
        myFile="${BEDDIR}FullMerge.p.${SENSIT}.${REGION}.noReps.bed"
        DEEPDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/plots/DEEP/${REGION}.1.all/"
    fi
    echo
    echo "$myFile"
    echo
    mkdir -p $DEEPDIR
elif [ $positionScope == 2 ]
then
    if [[ $# -gt 4 ]]
    then
        myFile="${BEDDIR}Distal.p.${SENSIT}.${REGION}.noReps.bed.purged"
        DEEPDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/plots/DEEP/${REGION}.1.distal.purged/"
    else
        myFile="${BEDDIR}Distal.p.${SENSIT}.${REGION}.noReps.bed"
        DEEPDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/plots/DEEP/${REGION}.1.distal/"
    fi
    echo
    echo "$myFile"
    echo
    mkdir -p $DEEPDIR
elif [ $positionScope == 3 ]
then
    if [[ $# -gt 4 ]]
    then
        myFile="${BEDDIR}Proximal.p.${SENSIT}.${REGION}.noReps.bed.purged"
        DEEPDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/plots/DEEP/${REGION}.1.proximal.purged/"
    else
        myFile="${BEDDIR}Proximal.p.${SENSIT}.${REGION}.noReps.bed"
        DEEPDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/plots/DEEP/${REGION}.1.proximal/"
    fi
    echo
    echo "$myFile"
    echo
    mkdir -p $DEEPDIR
elif [ $positionScope == 4 ]
then
    myFile="${BEDDIR}${clusterID}"
    DEEPDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/plots/DEEP/${clusterID}/"

    echo
    echo "$myFile"
    echo
    mkdir -p $DEEPDIR
fi


MYBWAallBG=('Arx_BG_e13-5_ChIP_rep3_R1_trimmed.fq.gz.srt.bam' 'Arx_BG_e13-5_Input_rep3_R1_trimmed.fq.gz.srt.bam' 'ChIP_ACSL1_BD_DSG_e13_S20_L002_R1_001_trimmed.fq.gz.srt.bam' 'input_DSG_e13_S19_L002_R1_001_trimmed.fq.gz.srt.bam' 'Dlx1-ChIP_BG_e13_Rep4_S124_L008_R1_001_trimmed.fq.gz.srt.bam' 'Dlx1-Input_BG_e13_Rep4_S123_L008_R1_001_trimmed.fq.gz.srt.bam' 'Dlx2-BG-e13_S78_L008_R1_001_trimmed.fq.gz.srt.bam' 'Input-Dlx2-BG-e13_S77_L008_R1_001_trimmed.fq.gz.srt.bam'  'ChIP-Dlx5_CCTCGG_L002_R1_trimmed.fq.gz.srt.bam' 'Dlx1_Dlx5_BG_e13-5_Input_rep1_trimmed.fq.gz.srt.bam' 'SLL11_bc05_AAGACG_L005_R1_001_trimmed.fq.gz.srt.bam' 'SLL11_bc04_TTCAGC_L005_R1_001_trimmed.fq.gz.srt.bam' 'SRR4035329_trimmed.fq.gz.srt.bam' 'SRR4035332_trimmed.fq.gz.srt.bam' 'SRR4035319_trimmed.fq.gz.srt.bam' 'SRR4035322_trimmed.fq.gz.srt.bam' 'SLL8_bc02_CCTTCA_L001_R1_001_trimmed.fq.gz.srt.bam' 'SLL8_bc01_AAGGGA_L001_R1_001_trimmed.fq.gz.srt.bam' 'OTX2_BG_E13.5_ChIP_rep3_trimmed.fq.gz.srt.bam' 'OTX2_BG_E13.5_Input_rep3_trimmed.fq.gz.srt.bam' 'SLL10_bc10_CACACA_L006_R1_001_trimmed.fq.gz.srt.bam' 'SLL10_bc09_ACACGA_L006_R1_001_trimmed.fq.gz.srt.bam' 'SP9-BG-e13_S87_L008_R1_001_trimmed.fq.gz.srt.bam' 'Input-SP9-BG-e13_S86_L008_R1_001_trimmed.fq.gz.srt.bam')

LABELSall0BG=('ARX_E13.5_BG_3' 'ASCL1_E13.5_BG_2' 'DLX1_E13.5_BG_N4' 'DLX2_E13.5_BG_N1' 'DLX5_E13.5_BG_O2' 'GSX2_E13.5_BG_1' 'LHX6_E13.5_BG_3' 'NKX2.1_E13.5_BG_2' 'NR2F1_E13.5_BG_1' 'OTX2_E13.5_BG_1' 'PBX1_E13.5_BG_1' 'SP9_E13.5_BG_N1')

MYBWAallCX=()

LABELSall0CX=()

if [ "$REGION" == "BG" ]
then
    LABELSall0=("${LABELSall0BG[@]}")
    MYBWAall=("${MYBWAallBG[@]}")
elif [ "$REGION" == "CX" ]
then
    LABELSall0=("${LABELSall0CX[@]}")
    MYBWAall=("${MYBWAallCX[@]}")
elif [ "$REGION" == "combined" ]
then
    LABELSall0=("${LABELSall0BG[@]}" "${LABELSall0CX[@]}")
    MYBWAall=("${MYBWAallBG[@]}" "${MYBWAallCX[@]}")
fi

for i in `echo ${LABELSall0[*]}`
do
    LABELSall="${LABELSall} ${i}_ChIP ${i}_input"
done

MYCOV=''
for i in `echo ${MYBWAall[*]}`
do
    MYCOV="${MYCOV} ${DEEPTOOLS}${i}.normalizedRPKM.coverage.bw"
done

MYBWA=''
for i in `echo ${MYBWAall[*]}`
do
    MYBWA="${MYBWA} ${BWADIR}${i}"
#    echo "... indexing $bwa"
#    samtools index "${BWADIR}$bwa"
done

echo
echo


MATRIX="${DEEPDIR}compute.matrix.${SENSIT}.${REGION}.gz"

# Only ChIP bams
MYBWAchipBG1=('Arx_BG_e13-5_ChIP_rep3_R1_trimmed.fq.gz.srt.bam'  'ChIP_ACSL1_BD_DSG_e13_S20_L002_R1_001_trimmed.fq.gz.srt.bam'  'Dlx1-ChIP_BG_e13_Rep4_S124_L008_R1_001_trimmed.fq.gz.srt.bam'  'Dlx2-BG-e13_S78_L008_R1_001_trimmed.fq.gz.srt.bam'  'ChIP-Dlx5_CCTCGG_L002_R1_trimmed.fq.gz.srt.bam'  'SLL11_bc05_AAGACG_L005_R1_001_trimmed.fq.gz.srt.bam'  'SRR4035329_trimmed.fq.gz.srt.bam' 'SRR4035319_trimmed.fq.gz.srt.bam'  'SLL8_bc02_CCTTCA_L001_R1_001_trimmed.fq.gz.srt.bam'  'OTX2_BG_E13.5_ChIP_rep3_trimmed.fq.gz.srt.bam'  'SLL10_bc10_CACACA_L006_R1_001_trimmed.fq.gz.srt.bam'  'SP9-BG-e13_S87_L008_R1_001_trimmed.fq.gz.srt.bam')

LABELSchipBG1=('ARX' 'ASCL1' 'DLX1' 'DLX2' 'DLX5' 'GSX2' 'LHX6' 'NKX2.1' 'NR2F1' 'OTX2' 'PBX1' 'SP9')

MYBWAchipCX1=()

LABELSchipCX1=()



if [ "$REGION" == "BG" ]
then
    if [ "$BGSCOPE" == "1" ]
    then
        LABELSchip=("${LABELSchipBG1[@]}")
        MYBWAchip=("${MYBWAchipBG1[@]}")
    elif [ "$BGSCOPE" == "2" ]
    then
        LABELSchip=("${LABELSchipBG2[@]}")
        MYBWAchip=("${MYBWAchipBG2[@]}")
    elif [ "$BGSCOPE" == "3" ]
    then
        LABELSchip=("${LABELSchipBG3[@]}")
        MYBWAchip=("${MYBWAchipBG3[@]}")
    fi
elif [ "$REGION" == "CX" ]
then
    if [ "$BGSCOPE" == "1" ]
    then
        LABELSchip=("${LABELSchipCX1[@]}")
        MYBWAchip=("${MYBWAchipCX1[@]}")
    elif [ "$BGSCOPE" == "2" ]
    then
        LABELSchip=("${LABELSchipCX2[@]}")
        MYBWAchip=("${MYBWAchipCX2[@]}")
    elif [ "$BGSCOPE" == "3" ]
    then
        LABELSchip=("${LABELSchipCX3[@]}")
        MYBWAchip=("${MYBWAchipCX3[@]}")
    fi
elif [ "$REGION" == "combined" ]
then
    LABELSchip=("${LABELSchipBG1[@]}" "${LABELSchipCX1[@]}")
    MYBWAchip=("${MYBWAchipBG1[@]}" "${MYBWAchipCX1[@]}")
fi

MYCOV=''
for i in `echo ${MYBWAchip[*]}`
do
    MYCOV="${MYCOV} ${DEEPTOOLS}${i}.normalizedRPKM.coverage.bw"
done

MYBWA=''
for i in `echo ${MYBWAchip[*]}`
do
    MYBWA="${MYBWA} ${BWADIR}${i}"
done

MATRIX="${DEEPDIR}compute.matrix.${SENSIT}.${REGION}.gz"

echo
echo "... starting multiBamSummary"

if [ "$REGION" == "BG" ] || [ "$REGION" == "CX" ]
then
    multiBamSummary bins \
                    --bamfiles $MYBWA \
                    --outFileName "${DEEPDIR}chip.matrix.${BGSCOPE}.npz" \
                    --blackListFileName $BLKLST \
                    --numberOfProcessors $THREADS \
                    --ignoreDuplicates

    echo "... finished multiBamSummary"

    echo
    echo "... starting plotCorrelation Pearson scatter"

    plotCorrelation --corData "${DEEPDIR}chip.matrix.${BGSCOPE}.npz" \
                    --corMethod pearson \
                    --skipZeros \
                    --plotTitle "Pearson Correlation of Average Scores" \
                    --whatToPlot scatterplot \
                    --plotFile "${DEEPDIR}scatterplot_PearsonCorr_bigwigScores.${BGSCOPE}.pdf"   \
                    --outFileCorMatrix "${DEEPDIR}PearsonCorr_bigwigScores.${BGSCOPE}.tab" \
                    --labels `echo ${LABELSchip[*]}` \
                    --zMin 0.5 \
                    --colorMap coolwarm &

    echo "... finished plotCorrelation Pearson scatter"

    echo
    echo "... starting plotCorrelation Pearson heat"

    plotCorrelation --corData "${DEEPDIR}chip.matrix.${BGSCOPE}.npz" \
                    --corMethod pearson \
                    --skipZeros \
                    --plotTitle "Pearson Correlation of Read Counts" \
                    --whatToPlot heatmap \
                    --colorMap BuGn \
                    --plotFile "${DEEPDIR}heatmap_PearsonCorr_readCounts.${BGSCOPE}.pdf"   \
                    --outFileCorMatrix "${DEEPDIR}SpearmanCorr_readCounts.${BGSCOPE}.tab" \
                    --labels `echo ${LABELSchip[*]}` \
                    --zMin 0.5 \
                    --colorMap coolwarm &

    echo "... finished plotCorrelation Pearson heat"

    echo
    echo "... starting plotCorrelation Spearman heat"

    plotCorrelation --corData "${DEEPDIR}chip.matrix.${BGSCOPE}.npz" \
                    --corMethod spearman \
                    --skipZeros \
                    --plotTitle "Spearman Correlation of Read Counts" \
                    --whatToPlot heatmap \
                    --colorMap BuGn \
                    --plotFile "${DEEPDIR}heatmap_SpearmanCorr_readCounts.${BGSCOPE}.pdf"   \
                    --outFileCorMatrix "${DEEPDIR}SpearmanCorr_readCounts.${BGSCOPE}.tab" \
                    --labels `echo ${LABELSchip[*]}` \
                    --zMin 0.5 --colorMap coolwarm &

    echo "... finished plotCorrelation Spearman heat"

    echo
    echo "... starting plotCorrelation Spearman heat"

    plotCorrelation --corData "${DEEPDIR}chip.matrix.${BGSCOPE}.npz" \
                    --corMethod spearman \
                    --skipZeros \
                    --plotTitle "Spearman Correlation of Average Scores" \
                    --whatToPlot scatterplot \
                    --plotFile "${DEEPDIR}scatterplot_SpearmanCorr_bigwigScores.${BGSCOPE}.pdf"   \
                    --outFileCorMatrix "${DEEPDIR}SpearmanCorr_bigwigScores.${BGSCOPE}.tab" \
                    --labels `echo ${LABELSchip[*]}` \
                    --zMin 0.5 \
                    --colorMap coolwarm &

    echo "... finished plotCorrelation Spearman scatter"

    echo
    echo "... starting plotPCA"

    plotPCA --corData "${DEEPDIR}chip.matrix.${BGSCOPE}.npz" \
            --plotFile "${DEEPDIR}PCA_readCounts.${BGSCOPE}.pdf" \
            --plotTitle "PCA of read counts"\
            --labels `echo ${LABELSchip[*]}`

    echo "... finished plotPCA"


    computeMatrix reference-point -S $MYCOV -R $myFile \
                    --downstream $myWindow \
                    --upstream $myWindow \
                    --numberOfProcessors $THREADS \
                    --referencePoint center \
                    --blackListFileName "$BLKLST" \
                    --outFileName "${MATRIX}.${REGION}.${BGSCOPE}.gz" \
                    --outFileSortedRegions "${DEEPDIR}COMB00.matrixSorted.${BGSCOPE}.bed" \
                    --samplesLabel `echo ${LABELSchip[*]}`

    echo "... finished computeMatrix"


    plotHeatmap --matrixFile "${MATRIX}.${REGION}.${BGSCOPE}.gz" \
                --outFileName "${DEEPDIR}Heatmap${BGSCOPE}.noclusters.${SENSIT}.pdf" \
                --outFileSortedRegions "${DEEPDIR}Heatmap${BGSCOPE}.noclusters.matrixSorted.${SENSIT}.bed" \
                --colorMap coolwarm \
                --xAxisLabel "distance (bp)" \
                --refPointLabel "0" \
                --samplesLabel `echo ${LABELSchip[*]}` &

    myMin=30

    for ncluster in {3..25}
        do
            plotHeatmap --matrixFile "${MATRIX}.${REGION}.${BGSCOPE}.gz" \
                        --outFileName "${DEEPDIR}Heatmap${BGSCOPE}.${ncluster}clusters.${SENSIT}.pdf" \
                        --outFileSortedRegions "${DEEPDIR}Heatmap${BGSCOPE}.${ncluster}clusters.matrixSorted.${SENSIT}.bed" \
                        --colorList black,yellow,darkred \
                        --alpha 1 \
                        --kmeans "$ncluster" \
                        --zMin $myMin \
                        --heatmapWidth 3 \
                        --interpolationMethod "gaussian" \
                        --xAxisLabel "distance (bp)" \
                        --refPointLabel "0" \
                        --whatToShow 'heatmap and colorbar' \
                        --samplesLabel `echo ${LABELSchip[*]}` &

            plotHeatmap --matrixFile "${MATRIX}.${REGION}.${BGSCOPE}.gz" \
                        --outFileName "${DEEPDIR}Heatmap${BGSCOPE}.${ncluster}clusters.${SENSIT}.grp.pdf" \
                        --outFileSortedRegions "${DEEPDIR}Heatmap${BGSCOPE}.${ncluster}clusters.matrixSorted.${SENSIT}.grp.bed" \
                        --outFileNameMatrix "${DEEPDIR}Heatmap${BGSCOPE}.${ncluster}clusters.Matrix.${SENSIT}.grp.txt" \
                        --colorList black,yellow,darkred \
                        --zMin $myMin \
                        --kmeans "$ncluster" \
                        --interpolationMethod "gaussian" \
                        --xAxisLabel "distance (bp)" \
                        --perGroup \
                        --refPointLabel "0" \
                        --whatToShow 'heatmap and colorbar' \
                        --samplesLabel `echo ${LABELSchip[*]}`


            echo "  "
            echo "     ... finished with ${ncluster} clusters"
        done
else
    multiBamSummary bins --bamfiles $MYBWA \
                    --outFileName "${DEEPDIR}chip.matrix.npz" \
                    --blackListFileName $BLKLST \
                    --numberOfProcessors $THREADS \
                    --ignoreDuplicates

    echo "... finished multiBamSummary"

    echo
    echo "... starting plotCorrelation Pearson scatter"

    plotCorrelation --corData "${DEEPDIR}chip.matrix.npz" \
                    --corMethod pearson \
                    --skipZeros \
                    --plotTitle "Pearson Correlation of Average Scores" \
                    --whatToPlot scatterplot \
                    --plotFile "${DEEPDIR}scatterplot_PearsonCorr_bigwigScores.pdf"   \
                    --outFileCorMatrix "${DEEPDIR}PearsonCorr_bigwigScores.tab" \
                    --labels `echo ${LABELSchip[*]}` \
                    --zMin 0.5 \
                    --colorMap coolwarm

    echo "... finished plotCorrelation Pearson scatter"

    echo
    echo "... starting plotCorrelation Pearson heat"

    plotCorrelation --corData "${DEEPDIR}chip.matrix.npz" \
                    --corMethod pearson --skipZeros \
                    --plotTitle "Pearson Correlation of Read Counts" \
                    --whatToPlot heatmap --colorMap BuGn \
                    --plotFile "${DEEPDIR}heatmap_PearsonCorr_readCounts.pdf"   \
                    --outFileCorMatrix "${DEEPDIR}SpearmanCorr_readCounts.tab" \
                    --labels `echo ${LABELSchip[*]}` \
                    --zMin 0.5 \
                    --colorMap coolwarm

    echo "... finished plotCorrelation Pearson heat"

    echo
    echo "... starting plotCorrelation Spearman heat"

    plotCorrelation --corData "${DEEPDIR}chip.matrix.npz" \
                    --corMethod spearman \
                    --skipZeros \
                    --plotTitle "Spearman Correlation of Read Counts" \
                    --whatToPlot heatmap --colorMap BuGn \
                    --plotFile "${DEEPDIR}heatmap_SpearmanCorr_readCounts.pdf"   \
                    --outFileCorMatrix "${DEEPDIR}SpearmanCorr_readCounts.tab" \
                    --labels `echo ${LABELSchip[*]}` \
                    --zMin 0.5 \
                    --colorMap coolwarm

    echo "... finished plotCorrelation Spearman heat"

    echo
    echo "... starting plotCorrelation Spearman heat"

    plotCorrelation --corData "${DEEPDIR}chip.matrix.npz" \
                    --corMethod spearman \
                    --skipZeros \
                    --plotTitle "Spearman Correlation of Average Scores" \
                    --whatToPlot scatterplot \
                    --plotFile "${DEEPDIR}scatterplot_SpearmanCorr_bigwigScores.pdf"   \
                    --outFileCorMatrix "${DEEPDIR}SpearmanCorr_bigwigScores.tab" \
                    --labels `echo ${LABELSchip[*]}` \
                    --zMin 0.5 \
                    --colorMap coolwarm

    echo "... finished plotCorrelation Spearman scatter"

    echo
    echo "... starting plotPCA"

    plotPCA --corData "${DEEPDIR}chip.matrix.npz" \
            --plotFile "${DEEPDIR}PCA_readCounts.pdf" \
            --plotTitle "PCA of read counts"\
            --labels `echo ${LABELSchip[*]}`

    echo "... finished plotPCA"


    computeMatrix reference-point -S $MYCOV -R $myFile \
                    --downstream $myWindow \
                    --upstream $myWindow \
                    --numberOfProcessors $THREADS \
                    --referencePoint center \
                    --blackListFileName "$BLKLST" \
                    --outFileName "$MATRIX" \
                    --outFileSortedRegions "${DEEPDIR}COMB00.matrixSorted.bed" \
                    --samplesLabel `echo ${LABELSchip[*]}`

    echo "... finished computeMatrix"


    plotHeatmap --matrixFile "${MATRIX}" \
                --outFileName "${DEEPDIR}Heatmap1.noclusters.${SENSIT}.pdf" \
                --outFileSortedRegions "${DEEPDIR}Heatmap1.noclusters.matrixSorted.${SENSIT}.bed" \
                --colorMap coolwarm \
                --xAxisLabel "distance (bp)" \
                --legendLocation "best" \
                --refPointLabel "0" \
                --samplesLabel `echo ${LABELSchip[*]}`

    for ncluster in {3..15}
    do
        plotHeatmap --matrixFile "${MATRIX}" \
                    --outFileName "${DEEPDIR}Heatmap1.${ncluster}clusters.${SENSIT}.pdf" \
                    --outFileSortedRegions "${DEEPDIR}Heatmap1.${ncluster}clusters.matrixSorted.${SENSIT}.bed" \
                    --colorMap coolwarm \
                    --kmeans "$ncluster" \
                    --xAxisLabel "distance (bp)" \
                    --legendLocation "best" \
                    --refPointLabel "0" \
                    --samplesLabel `echo ${LABELSchip[*]}`

        plotHeatmap --matrixFile "${MATRIX}" \
                    --outFileName "${DEEPDIR}Heatmap1.${ncluster}clusters.${SENSIT}.grp.pdf" \
                    --outFileSortedRegions "${DEEPDIR}Heatmap1.${ncluster}clusters.matrixSorted.${SENSIT}.grp.bed" \
                    --colorMap coolwarm \
                    --kmeans "$ncluster" \
                    --xAxisLabel "distance (bp)" \
                    --perGroup \
                    --legendLocation "best" \
                    --refPointLabel "0" \
                    --samplesLabel `echo ${LABELSchip[*]}`

        echo "  "
        echo "     ... finished with ${ncluster} clusters"
    done
fi


echo "... finished plotHeatmap"

