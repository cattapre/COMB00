include.peaks = function(peak.files, peak.names, raw.data.dir, dset, base.peak.set) {
  norm.factors = c()
  for (peak.file in peak.files) {
    file.index = which(peak.files == peak.file)
    if (file.index == 1) {
      print(noquote(" "))
      print(noquote("... including peak sets into base data set"))
    }
    print(noquote(peak.file))
    file.name = paste0(raw.data.dir, peak.file)
    
    if (grepl("histone.broad", dset)) {
      peak.data = read.table(file.name, sep='\t', header=TRUE, comment.char="#", as.is=TRUE, stringsAsFactors=FALSE)
      colnames(peak.data)[c(5,6,8,9)] = c("peak.height", "negLog10p_value", "negLog10q_value", "PeakID")
    } else if (grepl("TF.narrow", dset)) {
      peak.data = read.table(file.name, sep='\t', header=TRUE, comment.char="#", as.is=TRUE, stringsAsFactors=FALSE)
      colnames(peak.data)[c(6,7,9,10)] = c("peak.height", "negLog10p_value", "negLog10q_value", "PeakID")
      norm.factor = readLines(file.name)
      norm.factor = norm.factor[grep("tags after filtering in treatment: ", norm.factor)]
      norm.factor = unlist(strsplit(norm.factor, ": "))
      norm.factor = as.integer(norm.factor[2])
      print(norm.factor)
      norm.factors = c(norm.factors, norm.factor)
    } else if (grepl("histone.diff", dset)) {
      peak.data = read.table(file.name, sep='\t', header=FALSE, comment.char="#", as.is=TRUE, skip = 1,
                             stringsAsFactors=FALSE, col.names = c("chr", "start", "end", "PeakID", "log10.likelihood.ratio"))
      peak.data$length = peak.data$end - peak.data$start + 1
    }
    
    the.data = peak.data[peak.data$chr != "chrM", ]
    the.data = the.data[grep("RANDOM", the.data$chr, invert = TRUE),]
    the.data = the.data[grep("random", the.data$chr, invert = TRUE),]
    
    # filter for highest peak within complex region
    if (! grepl("diff|FC", peak.file)) {
      library(dplyr)
      the.data.tbl = as_tibble(the.data)
      the.data.tbl = the.data.tbl %>% group_by(chr, start, end) %>% summarise(peak.height = max(peak.height))
      the.data.tbl = as.data.frame(the.data.tbl)
      the.data = merge(the.data.tbl, the.data, by.x = c("chr", "start", "end", "peak.height"),
                       by.y = c("chr", "start", "end", "peak.height"))
      the.data = the.data[!duplicated(the.data[, c("chr", "start", "end")]),]
    }
    
    data.gr = GRanges(the.data)
    start(data.gr) = start(data.gr) + 1
    diff.prefix = sub(".bed", "", peak.file)
    diff.prefix = sub(paste0("_c", "[[:digit:]].[[:digit:]]"), "", diff.prefix)
    mcols(data.gr)$PeakID = sub(diff.prefix, peak.names[file.index], mcols(data.gr)$PeakID)
    mcols(data.gr)$PeakID = sub(sub("_peaks.xls", "", peak.file), peak.names[file.index], mcols(data.gr)$PeakID)
    mcols(data.gr)[, paste0(peak.names[file.index], ".PeakID")] = mcols(data.gr)$PeakID
    
    if (grepl("diff", dset)) {
      mcols(data.gr)[, paste0(peak.names[file.index], ".logLR")] = as.numeric(mcols(data.gr)$log10.likelihood.ratio)
      mcols(data.gr)[, c("length", "log10.likelihood.ratio")] = NULL
    } else if (grepl("FC", peak.file)) {
      mcols(data.gr)[, "length"] = NULL
    } else {
      mcols(data.gr)[, paste0(peak.names[file.index], ".peakHeight")] = as.numeric(mcols(data.gr)$peak.height)
      mcols(data.gr)[, paste0(peak.names[file.index], ".FC")] = as.numeric(mcols(data.gr)$fold_enrichment)
      mcols(data.gr)[, paste0(peak.names[file.index], ".negLog10p_value")] = as.numeric(mcols(data.gr)$negLog10p_value)
      mcols(data.gr)[, c("length", "peak.height", "negLog10p_value", "fold_enrichment", "negLog10q_value", "PeakID")] = NULL
      if (! grepl("broad", dset)) mcols(data.gr)[, "abs_summit"] = NULL
    }
    
    if (grepl("diff|FC", peak.file)) {
      ovlp = as.data.frame(findOverlaps((base.peak.set + peak.buffer), data.gr))
      dups.peak.tbl = as_tibble(ovlp)
      dups.peak.tbl = dups.peak.tbl %>% mutate(strength = mcols(data.gr)[subjectHits, 2])
      ovlp = dups.peak.tbl %>% group_by(queryHits) %>% summarise(strength = max(strength))
      ovlp = dplyr::semi_join(dups.peak.tbl, ovlp, by = c("queryHits", "strength"))
      ovlp = as.data.frame(ovlp[, c("queryHits", "subjectHits")])
    } else {
      ovlp = as.data.frame(findOverlaps(base.peak.set, data.gr))
    }
    
    mcols(base.peak.set)[paste0(peak.names[file.index], ".PeakID")] = ""
    mcols(base.peak.set)[ovlp$queryHits, paste0(peak.names[file.index], ".PeakID")] = mcols(data.gr)[ovlp$subjectHits, paste0(peak.names[file.index], ".PeakID")]
    
    if (grepl("diff", peak.file)) {
      mcols(base.peak.set)[, paste0(peak.names[file.index], ".logLR")] = 0
      mcols(base.peak.set)[ovlp$queryHits, paste0(peak.names[file.index], ".logLR")] = mcols(data.gr)[ovlp$subjectHits, paste0(peak.names[file.index], ".logLR")]
    } else if (grepl("FC", peak.file)) {
      mcols(base.peak.set)[, paste0(peak.names[file.index], "fold_change")] = 0
      mcols(base.peak.set)[ovlp$queryHits, paste0(peak.names[file.index], "fold_change")] = mcols(data.gr)[ovlp$subjectHits, "fold_change"]
    } else if (grepl("TF.narrow", dset)) {
      mcols(base.peak.set)[, paste0(peak.names[file.index], ".peakHeight")] = 0
      mcols(base.peak.set)[ovlp$queryHits, paste0(peak.names[file.index], ".peakHeight")] = mcols(data.gr)[ovlp$subjectHits, paste0(peak.names[file.index], ".peakHeight")]
      mcols(base.peak.set)[, paste0(peak.names[file.index], ".normHeight")] = mcols(base.peak.set)[, paste0(peak.names[file.index], ".peakHeight")] / norm.factor
      #print(noquote(paste(peak.file, "normalized")))
      #print(mcols(base.peak.set)[, paste0(peak.names[file.index], ".normHeight")])
      mcols(base.peak.set)[, paste0(peak.names[file.index], ".FC")] = 0
      mcols(base.peak.set)[ovlp$queryHits, paste0(peak.names[file.index], ".FC")] = mcols(data.gr)[ovlp$subjectHits, paste0(peak.names[file.index], ".FC")]
      mcols(base.peak.set)[, paste0(peak.names[file.index], ".negLog10p_value")] = 0
      mcols(base.peak.set)[ovlp$queryHits, paste0(peak.names[file.index], ".negLog10p_value")] = mcols(data.gr)[ovlp$subjectHits, paste0(peak.names[file.index], ".negLog10p_value")]
    } else {
      # mcols(data.gr)[, paste0(peak.names[file.index], ".peakHeight")] = as.numeric(mcols(data.gr)$peak.height)
      mcols(base.peak.set)[, paste0(peak.names[file.index], ".FC")] = 0
      mcols(base.peak.set)[ovlp$queryHits, paste0(peak.names[file.index], ".FC")] = mcols(data.gr)[ovlp$subjectHits, paste0(peak.names[file.index], ".FC")]
      mcols(base.peak.set)[, paste0(peak.names[file.index], ".negLog10p_value")] = 0
      mcols(base.peak.set)[ovlp$queryHits, paste0(peak.names[file.index], ".negLog10p_value")] = mcols(data.gr)[ovlp$subjectHits, paste0(peak.names[file.index], ".negLog10p_value")]
    }
    
    
    
  }
  
  # >>>>>>>>>
  if (dset == "TF.narrow") {
    for (peak.name in peak.names) {
      file.index = which(peak.names == peak.name)
      #print(mcols(base.peak.set)[, paste0(peak.names[file.index], ".normHeight")])
      mcols(base.peak.set)[, paste0(peak.names[file.index], ".normHeight")] = mcols(base.peak.set)[, paste0(peak.names[file.index], ".normHeight")] * min(norm.factors)
    }
  }

  diff.columns = sort(names(mcols(base.peak.set))[grep("_vs_.*cond..logLR", names(mcols(base.peak.set)))])
  diff.names = unique(sub("_c[[:digit:]].*", "", diff.columns))
  
  for (diff in diff.names) {
    chosen.diff = diff.columns[grep(diff, diff.columns)]
    mcols(base.peak.set)[, paste0(diff, ".strength")] = mcols(base.peak.set)[, chosen.diff[1]] - mcols(base.peak.set)[, chosen.diff[2]]
  }
  
  strength.columns = names(mcols(base.peak.set))[grep("strength", names(mcols(base.peak.set)))]
  for (str.column in strength.columns) {
    mcols(base.peak.set)[, str.column][10 ^ (mcols(base.peak.set)[, str.column]) == 0] = 0
  }
  
  print(noquote("Finished including peaks"))
  return(base.peak.set)
}
