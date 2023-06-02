apply.genomic.filter = function(base.peak.set, peak.names) {
  # Annotate gaps
  file.name = paste0(genome.dir, "mm10.gaps.txt")
  gaps = read.table(file.name, sep="\t", header=FALSE, comment.char = "#", col.names = c("chr", "start", "end", 
                                                                                         "un.known.size", "size", 
                                                                                         "GapName", "bridge"))
  gaps$GapName[gaps$GapName == "other"] = "fragment"
  gaps$GapName = paste0(gaps$GapName, "_", seq(1, nrow(gaps)))
  gaps = GRanges(gaps)
  
  ovlp = as.data.frame(findOverlaps(base.peak.set, gaps))
  mcols(base.peak.set)[ovlp$queryHits, "gaps"] = as.character(mcols(gaps)[ovlp$subjectHits, "GapName"])
  overlaps = pintersect(base.peak.set[ovlp$queryHits], gaps[ovlp$subjectHits])
  percentOverlap = width(overlaps) / width(base.peak.set[ovlp$queryHits])
  mcols(base.peak.set)$Pass.Filter[ovlp$queryHits] = mcols(base.peak.set)$Pass.Filter[ovlp$queryHits] - percentOverlap
  mcols(base.peak.set)$Pass.Filter[!is.na(mcols(base.peak.set)$gaps) &
                                     mcols(base.peak.set)$Pass.Filter < filter.proportion] = 0
  mcols(base.peak.set)$Pass.Filter = ceiling(mcols(base.peak.set)$Pass.Filter)
  
  # Annotate blacklisted regions
  file.name = paste0(genome.dir, "mm10.blacklist.bed")
  blacklist = GRanges(read.table(file.name, sep="\t", header=FALSE, stringsAsFactors = FALSE, col.names = c("chr", "start", "end")))
  blacklist = sort(blacklist)
  mcols(blacklist)$BlacklistName = paste0("Blacklist_", seq(1, length(blacklist)))
  ovlp = as.data.frame(findOverlaps(base.peak.set, blacklist))
  mcols(base.peak.set)[ovlp$queryHits, "Blacklist"] = as.character(mcols(blacklist)$BlacklistName[ovlp$subjectHits])
  overlaps = pintersect(base.peak.set[ovlp$queryHits], blacklist[ovlp$subjectHits])
  percentOverlap = width(overlaps) / width(base.peak.set[ovlp$queryHits])
  mcols(base.peak.set)$Pass.Filter[ovlp$queryHits] = mcols(base.peak.set)$Pass.Filter[ovlp$queryHits] - percentOverlap
  mcols(base.peak.set)$Pass.Filter[!is.na(mcols(base.peak.set)$Blacklist)] = 0
  mcols(base.peak.set)$Pass.Filter = ceiling(mcols(base.peak.set)$Pass.Filter)
  
  # Annotate Satelite/Repeat regions
  Repeats = read.table(paste0(genome.dir, "mm10.rmsk.txt"), sep='\t', header=FALSE,
                       comment.char="#", as.is=TRUE, stringsAsFactors=FALSE)
  Repeats = Repeats[, -c(1:5, 9:10, 13:16)]
  colnames(Repeats) = c("chr", "start", "end", "description", "type", "RepeatName")
  
  Repeats$RepeatName = paste0("Satellite_", seq(1, nrow(Repeats)))
  Repeats = Repeats[Repeats$type == "Satellite", ]
  Repeats = GRanges(Repeats)
  
  ovlp = as.data.frame(findOverlaps(base.peak.set, Repeats))
  mcols(base.peak.set)[ovlp$queryHits, "Repeats"] = as.character(Repeats$RepeatName[ovlp$subjectHits])
  mcols(base.peak.set)[ovlp$queryHits, "Repeat.Description"] = as.character(Repeats$description[ovlp$subjectHits])
  overlaps = pintersect(base.peak.set[ovlp$queryHits], Repeats[ovlp$subjectHits])
  percentOverlap = width(overlaps) / width(base.peak.set[ovlp$queryHits])
  mcols(base.peak.set)$Pass.Filter[ovlp$queryHits] = mcols(base.peak.set)$Pass.Filter[ovlp$queryHits] - percentOverlap
  mcols(base.peak.set)$Pass.Filter[!is.na(mcols(base.peak.set)$Repeats) &
                                     mcols(base.peak.set)$Pass.Filter < filter.proportion] = 0
  mcols(base.peak.set)$Pass.Filter = ceiling(mcols(base.peak.set)$Pass.Filter)
  columns.to.exclude = grep(paste(peak.names, collapse = "|"), names(mcols(base.peak.set)))
  mcols(base.peak.set)[, columns.to.exclude] = NULL

  return(base.peak.set)
}
