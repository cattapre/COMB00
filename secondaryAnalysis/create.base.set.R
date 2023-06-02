create.base.set = function(peak.files, peak.names, filter.peaks, raw.data.dir, dset) {
  library(tidyverse)
  library(paste0("TxDb.Mmusculus.UCSC.", genome.used, ".knownGene"), character.only = TRUE)
  
  for (peak.file in peak.files) {
    file.index <<- which(peak.files == peak.file)
    if (file.index == 1) print(noquote("... creating base data set"))
    print(noquote(peak.file))
    file.name = paste0(raw.data.dir, peak.file)
    if (grepl("histone.broad", dset)) {
      peak.data = read.table(file.name, sep='\t', header=TRUE, comment.char="#", as.is=TRUE, stringsAsFactors=FALSE)
      colnames(peak.data) = c("chr", "start", "end", "length", "peak.height", "negLog10p_value", "fold_enrichment", 
                              "negLog10q_value", "PeakID")
    } else {
      peak.data = read.table(file.name, sep='\t', header=TRUE, comment.char="#", as.is=TRUE, stringsAsFactors=FALSE)
      colnames(peak.data) = c("chr", "start", "end", "length", "abs_summit", "peak.height", "negLog10p_value", 
                              "fold_enrichment", "negLog10q_value", "PeakID")
    }
    the.data = peak.data[peak.data$chr != "chrM", ]
    the.data = the.data[grep("RANDOM", the.data$chr, invert = TRUE),]
    the.data = the.data[grep("random", the.data$chr, invert = TRUE),]
    
    # filter for largest peak within complex region
    the.data.tbl = as_tibble(the.data)
    the.data.tbl = the.data.tbl %>% 
      group_by(chr, start, end) %>% 
      summarise(peak.height = max(peak.height))
    the.data.tbl = as.data.frame(the.data.tbl)
    the.data = merge(the.data.tbl, the.data, by.x = c("chr", "start", "end", "peak.height"),
                     by.y = c("chr", "start", "end", "peak.height"))
    the.data = the.data[!duplicated(the.data[, c("chr", "start", "end")]),]
    
    # filter for significant peaks
    if (filter.peaks) {
      if (filter.peaks.by.peakHeight) {
        peak.condition = (the.data$length < peak.width.thrsh & the.data$peak.height > pkHeight.thrsh) |
                         (the.data$length >= peak.width.thrsh & the.data$peak.height >= pkHeight.thrsh)
        the.data = the.data[peak.condition, ]
      } else if (filter.peaks.by.PValue) {
        peak.condition = (the.data$length < peak.width.thrsh & the.data$negLog10p_value > PV.thrsh) |
                         (the.data$length >= peak.width.thrsh & the.data$negLog10p_value >= PV.thrsh)
        the.data = the.data[peak.condition, ]
      }
    }
    
    init = ncol(the.data) + 1
    the.data[, (init):(init + 4*length(peak.names) - 1)] = NA
    columns.yes = unlist(lapply(peak.names, function(x){c(paste0(x, ".PeakID"), paste0(x, ".peakHeight"), paste0(x, ".FC"), paste0(x, ".negLog10p_value"))}))
    colnames(the.data)[(init):(init+4*length(peak.names)-1)] = columns.yes
    
    data.gr = GRanges(the.data)
    
    # Convert macs2 BED 0-base to 1-base
    start(data.gr) = start(data.gr) + 1
    
    mcols(data.gr)[, paste0(peak.names[file.index], ".PeakID")] = sub(sub("_peaks.xls", "", peak.file), peak.names[file.index], mcols(data.gr)$PeakID)
    mcols(data.gr)[, paste0(peak.names[file.index], ".peakHeight")] = as.numeric(mcols(data.gr)$peak.height)
    mcols(data.gr)[, paste0(peak.names[file.index], ".FC")] = as.numeric(mcols(data.gr)$fold_enrichment)
    mcols(data.gr)[, paste0(peak.names[file.index], ".negLog10p_value")] = as.numeric(mcols(data.gr)$negLog10p_value)
    mcols(data.gr)[, c("peak.height", "abs_summit", "length", "PeakID", "fold_enrichment", "negLog10p_value", "negLog10q_value")] = NULL
    
    columns.to.merge = c()
    for (mycolumn in seq_along(mcols(data.gr))) {
      if (sum(!is.na(mcols(data.gr))[, mycolumn]) > 0) {
        columns.to.merge = c(columns.to.merge, mycolumn)
      }
    }
    
    if (file.index == 1) {
      data.gr.reduced = GenomicRanges::reduce(data.gr, with.revmap = FALSE, min.gapwidth = 10)
      base.peak.set = sort(data.gr[data.gr %over% data.gr.reduced])
      
    } else {
      data.gr.reduced = GenomicRanges::reduce(data.gr, with.revmap = FALSE, min.gapwidth = 10)
      data.gr.annot = data.gr[data.gr %over% data.gr.reduced]
      
      new.base.peak = c(base.peak.set, data.gr.annot)
      base.set.reduced = GenomicRanges::reduce(new.base.peak, with.revmap = TRUE, min.gapwidth = 10)
      revmap = as.data.frame(mcols(base.set.reduced)$revmap)
      prior.lines = revmap[revmap$value <= length(base.peak.set), ]
      prior.lines$calcsum = rowSums(as.data.frame(mcols(new.base.peak)[prior.lines$value, c(FALSE, TRUE, FALSE, FALSE)], row.names = NULL), na.rm = TRUE)
      prior.prgd = prior.lines %>%
        group_by(group, value) %>%
        summarise(max.FC = max(calcsum))
      prior.prgd = prior.prgd[!duplicated(prior.prgd$group), ]
      
      postr.lines = revmap[revmap$value > length(base.peak.set), ]
      postr.lines = as_tibble(postr.lines)
      postr.lines$calcsum = rowSums(as.data.frame(mcols(new.base.peak)[postr.lines$value, c(FALSE, TRUE, FALSE, FALSE)], row.names = NULL), na.rm = TRUE)
      postr.prgd = postr.lines %>%
        group_by(group, value) %>%
        summarise(max.FC = max(calcsum))
      postr.prgd = postr.prgd[!duplicated(postr.prgd$group), ]
      
      # metadata columns from previous interation
      base.set.reduced = base.set.reduced[, -1]
      mcols(base.set.reduced)[, columns.yes] = NA
      mcols(base.set.reduced)[, c(TRUE, FALSE, FALSE, FALSE)] = ""
      mcols(base.set.reduced)[, c(FALSE, TRUE, TRUE, TRUE)] = 0
      
      mcols(base.set.reduced)[prior.prgd$group, columns.yes[1:(columns.to.merge[1]-1)]] =
        mcols(new.base.peak)[prior.prgd$value, columns.yes[1:(columns.to.merge[1]-1)]]
      mcols(base.set.reduced)[postr.prgd$group, columns.yes[columns.to.merge]] =
        mcols(new.base.peak)[postr.prgd$value, columns.yes[columns.to.merge]]
      base.peak.set = base.set.reduced
    }
  }
  
  # Add column for filtering
  mcols(base.peak.set)$Pass.Filter = 1

  # Create filters
  source(paste0(base.dir, "annotate.feature.R"))
  feature.types = c("genomic_features", "conservancy")
  base.peak.set = annotate.feature(base.peak.set, feature.types)
  source(paste0(base.dir, "filter.blklst.gaps.R"))
  base.peak.set = apply.genomic.filter(base.peak.set, peak.names)
  source(paste0(base.dir, "annotate.VISTA.R"))
  base.peak.set = annotate.VISTA.enhancers(base.peak.set)

  arguments = list(peak.files, peak.names, raw.data.dir, dset, base.peak.set)
  base.peak.set = do.call(include.peaks, arguments)
  
  print(noquote("Finished creating, annotating primary table"))
  
  return(base.peak.set)
}
