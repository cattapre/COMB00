# COMB00 master table
# R. Catta-Preta [2019-05-05 version]
starttime = Sys.time()

# Datasets
{
  source("COMB00.datasets.R")  
}

# Global parameters
{
  region = "BG"   # (options are BG, CX or combined)
  source("COMB00.global.parameters.R")
}

# Directories
{
  source("COMB00.directories.R")
}

# Functions
{
  source(paste0(base.dir, "make.master.table.R"))
  source(paste0(base.dir, "format.table.R"))
  source(paste0(base.dir, "export.to.cluster.R"))
  source(paste0(base.dir, "generate.bed.files.R"))
}

# Initialize objects
{
  library(GenomicFeatures)
  base.peak.set = GRanges()
  bedfiles.to.transfer = c()
}

# Run codes
{
  base.peak.set = make.master.table(
    dataset, dataset.histones, dataset.atacseq, dataset.placseq, 
    filter.peaks.by.FC, filter.peaks.by.PValue, 
    use.histones,
    diff.histones,
    valid.enhancers,
    DEX,
    plac,
    atac,
    chromatin
  )
  
  base.peak.set = format.table(base.peak.set)
  
  # Filter for significant peaks
  base.peak.set = base.peak.set[mcols(base.peak.set)$Pass.Filter == 1]
  
  bedfiles.to.transfer = generate.bed.file(base.peak.set)
  tempFile = data.frame(base.peak.set, stringsAsFactors = FALSE)
  if (filter.peaks.by.peakHeight) {
    file.name = paste0(output.dir, "FullMerge.h.", pkHeight.thrsh, ".", region, ".txt")
  } else if (filter.peaks.by.PValue) {
    file.name = paste0(output.dir, "FullMerge.p.", pkHeight.thrsh, ".", region, ".txt")
  }
  write.table(tempFile, file.name, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  
  if (!allReps) {
    peaks = data.frame(base.peak.set, stringsAsFactors = FALSE)
    TF.Columns = colnames(peaks)[grep("normHeight", colnames(peaks))]
    TFs = unique(sub("_.*", "", TF.Columns))
    myFilter = get(paste0(region, "Reps"))
    for (myTF in TFs) {
      myColumns = TF.Columns[grep(myTF, TF.Columns)]
      for (myColumn in myColumns) {
        myRep = sub("rep", "", sub("\\..*", "", sub(".*_", "", myColumn)))
        if (myRep != myFilter$rep[myFilter$TF == myTF]) {
          columns.to.exclude = sub(".normHeight", "", myColumn)
          columns.to.exclude = colnames(peaks)[grep(columns.to.exclude, colnames(peaks))]
          peaks[, columns.to.exclude] = NULL
        }
      }
    }
    TF.Columns = colnames(peaks)[grep("normHeight", colnames(peaks))]
    peaks = peaks[rowSums(peaks[, TF.Columns]) > 0, ]
    if (filter.peaks.by.PValue) {
      qc.test = colnames(peaks)[grep("negLog10p_value", colnames(peaks))]
    } else if (filter.peaks.by.peakHeight) {
      qc.test = colnames(peaks)[grep("peakHeight", colnames(peaks))]
    }
    qc.test = qc.test[grep(paste(TFs, collapse = "|"), qc.test)]
    qc.test = peaks[, qc.test]
    qc.test = qc.test > ifelse(filter.peaks.by.PValue, PV.thrsh, pkHeight.thrsh)
    qc.test = apply(qc.test, 2, as.integer)
    qc.test = as.data.frame(qc.test)
    qc.test$test = rowSums(qc.test) > 0
    qc.test = qc.test$test
    peaks = peaks[qc.test, ]
    remove(qc.test)

    if (filter.peaks.by.peakHeight) {
      write.table(peaks, paste0(output.dir, "FullMerge.h.", pkHeight.thrsh, ".", region, ".noReps.txt"),
                  sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
      write.table(peaks[, c("seqnames", "start", "end")], 
                  paste0(output.dir, "FullMerge.h.", pkHeight.thrsh, ".", region, ".noReps.bed"),
                  sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
    } else if (filter.peaks.by.PValue) {
      write.table(peaks, paste0(output.dir, "FullMerge.p.", PV.thrsh, ".", region, ".noReps.txt"),
                  sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
      write.table(peaks[, c("seqnames", "start", "end")], 
                  paste0(output.dir, "FullMerge.p.", PV.thrsh, ".", region, ".noReps.bed"),
                  sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
    }
  }
  
  if (use.export.to.cluster) {
    temp = export.to.cluster(bedfiles.to.transfer)
    remove(temp)
  }
}

print(noquote(Sys.time() - starttime))
