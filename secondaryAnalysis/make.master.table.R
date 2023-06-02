make.master.table = function(dataset = dataset, filter.peaks.by.peakHeight, filter.peaks.by.PValue, ...) {
  source(paste0(base.dir, "create.base.set.R"))
  source(paste0(base.dir, "include.peaks.R"))
  
  if (diff.histones) source(paste0(base.dir, "include.differential.peaks.R"))
  if (valid.enhancers) source(paste0(base.dir, "annotate.validated.enhancers.R"))
  if (DEX) source(paste0(base.dir, "annotate.DEX.R"))
  if (atac) source(paste0(base.dir, "annotate.atac.seq.R"))
  if (plac) source(paste0(base.dir, "annotate.plac.seq.R"))
  if (chromatin) source(paste0(base.dir, "associate.chromatin.states.R"))
  
  for (dset in scopes.dataset) {
    if (dset == "TF.narrow") {
      print(noquote(" "))
      print(noquote("Starting generation of Master Table"))
      raw.data.dir = narrow.peak.dir
      peak.files = list.files(path = raw.data.dir, pattern = "^.*vs.Input.no_model_peaks.xls")
      peak.files = intersect(peak.files, paste0(dataset, "_trim-vs.Input.no_model_peaks.xls"))
      peak.names = gsub('_trim-vs.Input.no_model_peaks.xls', "", peak.files)
      base.peak.set = create.base.set(peak.files, peak.names, filter.peaks, raw.data.dir, dset)
      
    } else if (dset == "histone.broad" & use.histones) {
      print(noquote(" "))
      print(noquote("... Histone broadPeak files"))
      raw.data.dir = broad.peak.dir
      peak.files = list.files(path = raw.data.dir,
                              pattern = "^H3K.*vs.Input.no_model_peaks.xls")
      peak.files = peak.files[grep("KO", peak.files, invert = TRUE)]
      peak.files = peak.files[grep("Mut", peak.files, invert = TRUE)]
      peak.files = intersect(peak.files, paste0(dataset.histones, "_trim-vs.Input.no_model_peaks.xls"))
      peak.names = sub("_trim-vs.Input.no_model_peaks.xls", "", peak.files)
      arguments = list(peak.files, peak.names, raw.data.dir, dset, base.peak.set)
      base.peak.set = do.call(include.peaks, arguments)

    } else if (dset == "histone.diff" & diff.histones) {
      print(noquote(" "))
      print(noquote(paste0("... ", "Differential histone calls with bdgdiff")))
      raw.data.dir = diff.hist.dir
      peak.files = list.files(path = raw.data.dir, pattern = "trim_diff_BROAD_c[[:digit:]]\\.[[:digit:]]_co")
      peak.files = peak.files[grep("H3K", peak.files, invert = FALSE)]
      # temp.files = c()
      # for (histone in histones) {
      #   temp.file = peak.files[grep(paste0(histone, ".*WT.*", histone, ".*Mut.*"), peak.files)]
      #   temp.files = c(temp.files, temp.file)
      # }
      # peak.files2 = peak.files
      # peak.files = temp.files
      # peak.names = gsub(".bed", "", peak.files)
      # peak.names = gsub("_trim", "", peak.names)
      # peak.names = gsub("_Cx", "", peak.names)
      # peak.names = gsub("_diff_BROAD", "", peak.names)
      # # peak.names = gsub("vs_.*_KO", "vs_KO", peak.names)
      # # peak.names = gsub("vs_.*_WT", "vs_WT", peak.names)
      # # peak.names = gsub("e12.5_", "", peak.names)
      # # peak.names = gsub("rep1_", "", peak.names)
      
      # arguments = list(peak.files, peak.names, raw.data.dir, dset, base.peak.set)
      # base.peak.set = do.call(include.peaks, arguments)
      
      # base.peak.set = include.peaks(peak.files, peak.names, base.peak.set, raw.data.dir, dset)
      # 
      # temp.files = c()
      # for (histone in histones) {
      #   temp.file = peak.files2[grep(paste0(histone, ".*Mut.*", histones[histones != histone], ".*Mut.*"), peak.files2)]
      #   temp.files = c(temp.files, temp.file)
      #   temp.file = peak.files2[grep(paste0(histone, ".*WT.*", histones[histones != histone], ".*WT.*"), peak.files2)]
      #   temp.files = c(temp.files, temp.file)
      # }
      # peak.files2 = temp.files
    }
  }
  
  if (atac) base.peak.set = annotate.atac.seq(base.peak.set)
  if (plac) base.peak.set = suppressWarnings(annotate.plac.seq(base.peak.set))
  if (chromatin) base.peak.set = associate.chromatin.states(base.peak.set)

  return(base.peak.set)
}


