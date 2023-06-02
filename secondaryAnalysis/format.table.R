format.table = function(base.peak.set) {
  temp.frame = data.frame(base.peak.set, stringsAsFactors = FALSE)
  first.column = which(colnames(temp.frame) == "PeakID")
  second.column = which(colnames(temp.frame) == "Gene")
  temp.frame = cbind(PeakID=temp.frame[, first.column], Gene=temp.frame[, second.column], temp.frame[, -c(first.column, second.column)])

  return(base.peak.set)
}