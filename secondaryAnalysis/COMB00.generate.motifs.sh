#!/bin/bash -l

#SBATCH --job-name=motifs
#SBATCH --time=1-00:00
#SBATCH --mem=40G
#SBATCH --partition=high
#SBATCH --cpus-per-task=12

echo
echo This job is running on
/bin/hostname
echo "crick"
echo
echo

module load homer

PROJECTID="COMB00"
PLATFORM="group"
GENOME="mm10"
THREADS=12

SOURCEBED_DIR=$1
FILEID=$2
PVALUE=$3


PREPARSEDIR="/${PLATFORM}/nordlab/libraries/externalDatasets/${GENOME}"
TEMPDIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/temp.folder"
HOMER_DIR="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/motifs"

mkdir -p "${HOMER_DIR}/p.${PVALUE}"
mkdir -p "${HOMER_DIR}/p.${PVALUE}/${FILEID}"

if [ -s "${SOURCEBED_DIR}/${FILEID}" ]
then
    mkdir -p "${HOMER_DIR}/p.${PVALUE}/${FILEID}/selected"
    mkdir -p $PREPARSEDIR
    mkdir -p "${TEMPDIR}/p.${PVALUE}/${FILEID}.tmp"
    cd "${TEMPDIR}/p.${PVALUE}/${FILEID}.tmp"

## De novo motif discovery

    findMotifsGenome.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" "${HOMER_DIR}/p.${PVALUE}/$FILEID" -size -300,300 \
                            -mask -gc -p $THREADS -preparsedDir $PREPARSEDIR -seqlogo

    findMotifsGenome.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" "${HOMER_DIR}/p.${PVALUE}/$FILEID" -size -300,300 \
                            -mask -gc -p $THREADS -preparsedDir $PREPARSEDIR

    allmotiffs="${HOMER_DIR}/p.${PVALUE}/${FILEID}/homerMotifs.all.motifs"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $allmotiffs -nmotifs \
                    -cpu $THREADS > "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.motif.counts"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $allmotiffs \
                     -matrix "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.motif.counts.cooccur" \
                     -cpu $THREADS > "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.motif.counts.cooccur.out"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $allmotiffs -size -300,300 -hist 10 \
                     -cpu $THREADS > "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.hist.txt"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $allmotiffs -size -300,300 -hist 10 -ghist \
                     -cpu $THREADS > "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.heatmap.hist.txt"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $allmotiffs -size -300,300 -hist 10 -CpG \
                     -cpu $THREADS > "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.CpG.hist.txt"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $allmotiffs \
                     -mbed "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.motif.instances.bed" \
                     -cpu $THREADS > "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.motif.instances.annotated"

    tail -n +2 "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.motif.instances.bed" | sort -k1,1 -k2,2n > \
                "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.motif.instances.sorted.bed"

    annotatePeaks.pl tss "${GENOME}" -m $allmotiffs -size -300,300 -hist 10 -cpu $THREADS > "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.tss.hist.txt"

    GO_dir="${HOMER_DIR}/p.${PVALUE}/${FILEID}"

#    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $allmotiffs -nmotifs -go $GO_dir -genomeOntology $GO_dir -cpu $THREADS > \
#                     "${HOMER_DIR}/p.${PVALUE}/${FILEID}/${FILEID}.motif.counts.GO"



## calculate enrichment over the selected motifs

    # chosen select motifs
    SEL_MOTIF_NAME='selected.motifs'
    TGT_DIR="${HOMER_DIR}/p.${PVALUE}/${FILEID}/selected"
    mkdir -p "$TGT_DIR"

    # create combined motifs
    SELECTEDMOTIF="/${PLATFORM}/nordlab/users/rinaldo/${PROJECTID}/ChIPseq/data/motifs/${SEL_MOTIF_NAME}"
    MYMOTIFS=`ls ${HOMER_DIR}/p.${PVALUE}/${FILEID}/knownResults/*.motif`
    cat $MYMOTIFS > $SELECTEDMOTIF

    findMotifsGenome.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" "$TGT_DIR" -size -300,300 \
                        -mask -gc -p $THREADS -preparsedDir $PREPARSEDIR -nomotif -mknown $SELECTEDMOTIF -seqlogo

    findMotifsGenome.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" "$TGT_DIR" -size -300,300 \
                        -mask -gc -p $THREADS -preparsedDir $PREPARSEDIR -nomotif -mknown $SELECTEDMOTIF

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $SELECTEDMOTIF -nmotifs -cpu $THREADS > \
                     "${TGT_DIR}/${FILEID}.motif.counts"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $SELECTEDMOTIF \
                     -matrix "${TGT_DIR}/${FILEID}.motif.counts.cooccur" \
                     -cpu $THREADS > "${TGT_DIR}/${FILEID}.motif.counts.cooccur.out"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $SELECTEDMOTIF -size -300,300 -hist 10 \
                     -cpu $THREADS > "${TGT_DIR}/${FILEID}.hist.txt"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $SELECTEDMOTIF -size -300,300 -hist 10 -ghist \
                     -cpu $THREADS > "${TGT_DIR}/${FILEID}.heatmap.hist.txt"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $SELECTEDMOTIF -size -300,300 -hist 10 -CpG \
                     -cpu $THREADS > "${TGT_DIR}/${FILEID}.CpG.hist.txt"

    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $SELECTEDMOTIF \
                     -mbed "${TGT_DIR}/${FILEID}.motif.instances.bed" \
                     -cpu $THREADS > "${TGT_DIR}/${FILEID}.motif.instances.annotated"

    tail -n +2 "${TGT_DIR}/${FILEID}.motif.instances.bed" | sort -k1,1 -k2,2n > \
               "${TGT_DIR}/${FILEID}.motif.instances.sorted.bed"

    annotatePeaks.pl tss "${GENOME}" -m $SELECTEDMOTIF -size -300,300 -hist 10 -cpu $THREADS > \
                     "${TGT_DIR}/${FILEID}.tss.hist.txt"

    GO_dir="${TGT_DIR}"

#    annotatePeaks.pl "${SOURCEBED_DIR}/${FILEID}" "${GENOME}" -m $SELECTEDMOTIF -nmotifs -go $GO_dir \
#                     -genomeOntology $GO_dir -cpu $THREADS > "${TGT_DIR}/${FILEID}.motif.counts.GO"

else
    echo "File ${SOURCEBED_DIR}/${FILEID} does not exist, or is empty"
fi

rm -rf "${TEMPDIR}/p.${PVALUE}/${FILEID}.tmp"

