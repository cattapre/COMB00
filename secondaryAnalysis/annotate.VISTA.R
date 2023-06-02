annotate.VISTA.enhancers = function(base.peak.set) {
  print(noquote(" "))
  print(noquote("   ... VISTA enhancers"))
  file.name = paste0(enhancer.dir, enhancer.file)
  if (file.exists(file.name)) {
    InVivo = GRanges(read.table(file.name, sep='\t', header=TRUE,
                                comment.char="#", as.is=TRUE, stringsAsFactors=FALSE))
  } else {
    library(rtracklayer)
    InVivo = GRanges(read.table("/Volumes/GoogleDrive/Team Drives/Nord Lab - Personal Folders/Rinaldo/VISTA_enhancers/VISTA.enhancers.20180710.mm9.txt", sep='\t', header=TRUE,
                                comment.char="#", as.is=TRUE, stringsAsFactors=FALSE))
    chain = import.chain(paste0(genome.dir, "mm9ToMm10.over.chain"))
    seqlevelsStyle(InVivo) = "UCSC"
    temp.feature = unlist(liftOver(InVivo, chain))
    write.table(data.frame(InVivo, stringsAsFactors = FALSE), file.name, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  }

  ovlp = as.data.frame(findOverlaps(base.peak.set, InVivo))
  mcols(base.peak.set)[ovlp$queryHits, "InVivo.Name"] = InVivo$name[ovlp$subjectHits]
  mcols(base.peak.set)[ovlp$queryHits, "InVivo.Exp.Validated"] = InVivo$exp.valid[ovlp$subjectHits]
  mcols(base.peak.set)[ovlp$queryHits, "InVivo.Tissue"] = InVivo$tissue[ovlp$subjectHits]
  mcols(base.peak.set)[ovlp$queryHits, "InVivo.Same"] = 0
  same.condition = mcols(base.peak.set)[ovlp$queryHits, "InVivo.Tissue"] %in% c("forebrain", "basal_ganglia") |
    mcols(base.peak.set)[ovlp$queryHits, "InVivo.Exp.Validated"] == "positive"
  mcols(base.peak.set)[ovlp$queryHits, "InVivo.Same"][same.condition] = 1
  mcols(base.peak.set)[ovlp$queryHits, "InVivo.Same"] = 0
  same.condition = grepl("forebrain", mcols(base.peak.set)[ovlp$queryHits, "InVivo.Tissue"])
  mcols(base.peak.set)[ovlp$queryHits, "InVivo.Same"][same.condition] = 1
  mcols(base.peak.set)[ovlp$queryHits, "InVivo.NotSame"] = 0
  same.not.condition = ! grepl("forebrain|basal\\ ganglia", mcols(base.peak.set)[ovlp$queryHits, "InVivo.Tissue"])
  mcols(base.peak.set)[ovlp$queryHits, "InVivo.NotSame"][same.not.condition] = 1

  return(base.peak.set)
}