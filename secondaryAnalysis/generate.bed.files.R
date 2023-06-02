generate.bed.file = function(base.peak.set) {
  peaks = data.frame(base.peak.set, stringsAsFactors = FALSE)
  peaks = peaks[, c("seqnames", "start", "end")]
  file.name = paste0(output.dir, "FullMerge.", 
                     ifelse(filter.peaks.by.PValue, paste0("p.", PV.thrsh), paste0("h.", pkHeight.thrsh)), 
                     ".", region, ".bed")
  write.table(peaks, file.name, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  bedfiles.to.transfer = c(bedfiles.to.transfer, file.name)
  return(bedfiles.to.transfer)
}