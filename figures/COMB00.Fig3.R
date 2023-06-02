# Load libraries
{
  initial.packages = c("GenomicFeatures", "tidyverse", "pheatmap", "rtracklayer", "ggcorrplot", 
                       "ComplexHeatmap", "randomcoloR", "org.Mm.eg.db", "clusterProfiler", "ggpubr", 
                       "STRINGdb", "simpIntLists", "ggrepel", "ggsci")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (myPackage in initial.packages) {
    if (!require(myPackage, character.only = TRUE)) BiocManager::install(myPackage)
    library(myPackage, character.only = TRUE)
  }
  remove(list = c("initial.packages", "myPackage"))
}

# Load functions
{
  source("annotateHMM.R")
  source("annotateClusters.R")
  source("annotateEnsembles.R")
  source("ovlpEnhancers.R")
}

# Load parameters
{
  allReps = FALSE
  region = "BG"
  nPerm = 10000
  samplingRatio = 0.75
  minTargetThrsh = 10
  strippedVersion = TRUE
  bgTFs = c("ARX", "ASCL1", "DLX1", "DLX2", "DLX5", "GSX2", "LHX6", "NKX2.1", "NR2F1", "OTX2", "PBX1", "SP9")
  cxTFs = c("EMX2", "LHX2", "NR2F1", "PAX6", "PBX1")
  histones = c("H3K27ac", "H3K27me3", "H3K4me3", "H3K4me1")
  source("COMB00.directories.R")
  hmmDir = paste0(working.dir, "ChromHMM/Output/WT/")
  dexDir = paste0(working.dir, "DEX/")
  interactDir = paste0(working.dir, "interactome/")
  motif.dir = paste0(motif.dir, "noRepsStripped/")
  numHMM = 9
  promoterLimits = 2000
  chainMM9_10 = paste0("/Volumes/GoogleDrive/Shared drives/NordLabRinaldo/mm10/", 
                       "mm9ToMm10.over.chain")
  scaleLimit = 3
  minDiff = 2
  fontSize = 12
}

# Load data
{
  # bgPeaks
  {
    bgPeaks = paste0(output.dir, "FullMerge.p.0.BG.noReps.txt")
    if (file.exists(bgPeaks)) {
      bgPeaks = read.table(bgPeaks, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
      bgPeaks = bgPeaks[bgPeaks$Gene != 'Cwc22', ]
      bgPeaks$bgCluster = NA
    } else {
      stop(" ... Peak file not found")
    }
    bgPeaks = do.call(annotateClusters, list(bgPeaks, "bg"))
    bgPeaks = do.call(annotateEnsembles, list(bgPeaks, "bg"))
    bgPeaks = do.call(annotateHMM, list(bgPeaks, numHMM))
  }
  
  # cxPeaks
  {
    cxPeaks = paste0(output.dir, "FullMerge.p.0.CX.noReps.txt")
    cxPeaks = read.table(cxPeaks, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    cxPeaks = do.call(annotateClusters, list(cxPeaks, "cx"))
    
  }
  
  # randomSample
  {
    randomSample = paste0(output.dir, "randomSampleAnnotated.txt")
    if (file.exists(randomSample)) {
      randomSample = read.table(randomSample, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    } else {
      randomSample = paste0(output.dir, "randomSample.txt")
      randomSample = read.table(randomSample, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
      randomSample$seqnames[randomSample$seqnames == 'chrM'] = 'chr1'
      
      source("annotate.feature.R")
      myArguments = list(GRanges(randomSample), c("genomic_features", "conservancy"))
      randomSample = data.frame(do.call(annotate.feature, myArguments), stringsAsFactors = FALSE)
      source("annotate.VISTA.R")
      randomSample = data.frame(annotate.VISTA.enhancers(GRanges(randomSample)), stringsAsFactors = FALSE)
      
      # overlap with peaks
      ovlpRandom = as.data.frame(findOverlaps(GRanges(randomSample), GRanges(bgPeaks)))
      randomSample[ovlpRandom$queryHits, colnames(bgPeaks)[35:ncol(bgPeaks)]] = bgPeaks[ovlpRandom$subjectHits, colnames(bgPeaks)[35:ncol(bgPeaks)]]
      randomSample$bgCluster = 'random'
      
      file.name = paste0(output.dir, "randomSampleAnnotated.txt")
      write.table(randomSample, file.name, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
    randomSample$PeakID = sub("Merged", "Random", randomSample$PeakID)
    colnames(randomSample)[which(colnames(randomSample) == "cluster")] = "bgCluster"
    randomSample = do.call(annotateEnsembles, list(randomSample, "bg"))
    randomSample = do.call(annotateHMM, list(randomSample, numHMM))
  }
  
  # bgEmission probabilities (HMM)
  {
    bgEmission = paste0(hmmDir, "DLX/emissions_", numHMM, ".txt")
    bgEmission = read.table(bgEmission, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    colnames(bgEmission)[1] = "bgState"
  }
  
  # expressed genes in BG E13.5 (from DLX study)
  {
    expressedGenes = paste0(output.dir, "expressedGenes.txt")
    if (file.exists(expressedGenes)) {
      expressedGenes = read.table(expressedGenes, sep = '\t', quote = "", header = TRUE, stringsAsFactors = FALSE)
    } else {
      expressedGenes = paste0(dexDir, "DLX_RubensteinProject_mm9_fC_UN.txt")
      expressedGenes = read.table(expressedGenes, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
      myExcels = grep("-Sep|-Mar", expressedGenes$Gene)
      expressedGenes$Gene[myExcels] = paste0("Sept", sub("-Sep", "", expressedGenes$Gene[myExcels]))
      expressedGenes$Gene[myExcels] = paste0("Mar", sub("-Mar", "", expressedGenes$Gene[myExcels]))
      expressedGenes$Gene[myExcels] = paste0("Oct", sub("-Oct", "", expressedGenes$Gene[myExcels]))
      expressedGenes$meanCounts = rowMeans(expressedGenes[, grep("WT", colnames(expressedGenes))])
      expressedGenes = expressedGenes[, c("Gene", "meanCounts")]
      
      # gene aliases
      geneAlias = paste0("/Volumes/GoogleDrive/Shared drives/NordLabRinaldo/mm10/MGI_MRK_List2.rpt.txt")
      geneAlias = read.table(geneAlias, sep = '\t', header = TRUE, stringsAsFactors = FALSE, quote = "", fill = TRUE)
      geneAlias = geneAlias[, c("MGI.Accession.ID", "Marker.Symbol", "Marker.Name", "Marker.Type", "Feature.Type", "Marker.Synonyms..pipe.separated.")]
      geneAlias0 = geneAlias[grep("\\|", geneAlias$Marker.Synonyms..pipe.separated.), ]
      geneAlias = geneAlias[grep("\\|", geneAlias$Marker.Synonyms..pipe.separated., invert = TRUE), ]
      geneAlias = as.data.frame(geneAlias[, c("Marker.Symbol", "Marker.Synonyms..pipe.separated.")])
      colnames(geneAlias) = c("Gene", "Alias")
      geneAlias = geneAlias[geneAlias$Alias != "", ]
      
      geneAlias0 = as.data.frame(geneAlias0[, c("Marker.Symbol", "Marker.Synonyms..pipe.separated.")])
      colnames(geneAlias0) = c("Gene", "Alias")
      for (myLine in seq_along(geneAlias0)) {
        tempLine = unlist(strsplit(geneAlias0$Alias[myLine], "\\|"))
        for (myItem in seq_along(tempLine)) {
          geneAlias = rbind(geneAlias, data.frame(Gene = geneAlias0$Gene[myLine], 
                                                  Alias = tempLine[myItem], 
                                                  stringsAsFactors = FALSE))
        }
      }
      remove(geneAlias0)
      
      geneAlias$meanCounts = expressedGenes$meanCounts[match(geneAlias$Gene, expressedGenes$Gene)]
      geneAlias$meanCounts[is.na(geneAlias$meanCounts)] = expressedGenes$meanCounts[match(geneAlias$Alias[is.na(geneAlias$meanCounts)], expressedGenes$Gene)]
      geneAlias = geneAlias[!is.na(geneAlias$meanCounts), ]
      expressedGenes = do.call("rbind", 
                               list(expressedGenes, 
                                    geneAlias[, c("Gene", "meanCounts")], 
                                    data.frame(Gene = geneAlias$Alias, meanCounts = geneAlias$meanCounts, stringsAsFactors = FALSE)))
      expressedGenes = expressedGenes[!duplicated(expressedGenes$Gene), ]
      expressedGenes = expressedGenes[grep("\\[", expressedGenes$Gene, invert = TRUE), ]
      expressedGenes = expressedGenes[grep("\\ ", expressedGenes$Gene, invert = TRUE), ]
      expressedGenes$Gene = gsub("\\ ", "", expressedGenes$Gene)
      expressedGenes$Gene = gsub("\\\t", "", expressedGenes$Gene)
      
      file.name = paste0(output.dir, "expressedGenes.txt")
      write.table(expressedGenes, file.name, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
    
    geneExpression = paste0(output.dir, "geneExpression.txt")
    if (file.exists(geneExpression)) {
      geneExpression = read.table(geneExpression, sep = '\t', quote = "", header = TRUE, stringsAsFactors = FALSE)
    } else {
      geneExpression = paste0(dexDir, "dlx.DE.txt")
      geneExpression = read.table(geneExpression, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
      geneExpression$Gene = rownames(geneExpression)
      myExcels = grep("-Sep|-Mar", expressedGenes$Gene)
      geneExpression$Gene[myExcels] = paste0("Sept", sub("-Sep", "", geneExpression$Gene[myExcels]))
      geneExpression$Gene[myExcels] = paste0("Mar", sub("-Mar", "", geneExpression$Gene[myExcels]))
      geneExpression$Gene[myExcels] = paste0("Oct", sub("-Oct", "", geneExpression$Gene[myExcels]))
      geneExpression$meanRPKM = rowMeans(geneExpression[, grep("RPKM.WT", colnames(geneExpression))])
      geneExpression = geneExpression[, c("Gene", "meanRPKM")]
      
      # gene aliases
      geneAlias = paste0("/Volumes/GoogleDrive/Shared drives/NordLabRinaldo/mm10/MGI_MRK_List2.rpt.txt")
      geneAlias = read.table(geneAlias, sep = '\t', header = TRUE, stringsAsFactors = FALSE, quote = "", fill = TRUE)
      geneAlias = geneAlias[, c("MGI.Accession.ID", "Marker.Symbol", "Marker.Name", "Marker.Type", "Feature.Type", "Marker.Synonyms..pipe.separated.")]
      geneAlias0 = geneAlias[grep("\\|", geneAlias$Marker.Synonyms..pipe.separated.), ]
      geneAlias = geneAlias[grep("\\|", geneAlias$Marker.Synonyms..pipe.separated., invert = TRUE), ]
      geneAlias = as.data.frame(geneAlias[, c("Marker.Symbol", "Marker.Synonyms..pipe.separated.")])
      colnames(geneAlias) = c("Gene", "Alias")
      geneAlias = geneAlias[geneAlias$Alias != "", ]
      
      geneAlias0 = as.data.frame(geneAlias0[, c("Marker.Symbol", "Marker.Synonyms..pipe.separated.")])
      colnames(geneAlias0) = c("Gene", "Alias")
      for (myLine in seq_along(geneAlias0)) {
        tempLine = unlist(strsplit(geneAlias0$Alias[myLine], "\\|"))
        for (myItem in seq_along(tempLine)) {
          geneAlias = rbind(geneAlias, data.frame(Gene = geneAlias0$Gene[myLine], 
                                                  Alias = tempLine[myItem], 
                                                  stringsAsFactors = FALSE))
        }
      }
      remove(geneAlias0)
      
      geneAlias$meanRPKM = geneExpression$meanRPKM[match(geneAlias$Gene, geneExpression$Gene)]
      geneAlias$meanRPKM[is.na(geneAlias$meanRPKM)] = geneExpression$meanRPKM[match(geneAlias$Alias[is.na(geneAlias$meanRPKM)], geneExpression$Gene)]
      geneAlias = geneAlias[!is.na(geneAlias$meanRPKM), ]
      
      geneExpression = do.call("rbind", 
                               list(geneExpression, 
                                    geneAlias[, c("Gene", "meanRPKM")], 
                                    data.frame(Gene = geneAlias$Alias, meanRPKM = geneAlias$meanRPKM, stringsAsFactors = FALSE)))
      geneExpression = geneExpression[!duplicated(geneExpression$Gene), ]
      
      file.name = paste0(output.dir, "geneExpression.txt")
      write.table(geneExpression, file.name, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
  }
  
  # ensembleGenes
  {
    ensembleGenes = paste0(output.dir, "ensembleGenesLong.txt")
    ensembleGenes = read.table(ensembleGenes, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  }
  
  # plac-seq data
  {
    placReduced = paste0(output.dir, "consolidatedEnsemblesSummary.txt")
    placReduced = read.table(placReduced, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    placRaw = paste0(plac.seq.dir, "IJ010.IJ057.IJ058.peaks/IJ010.IJ057.IJ058.10k.2.peaks.bedpe")
    placRaw = read.table(placRaw, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    placFiltered = placRaw %>% 
      dplyr::mutate(interactionID = row_number()) %>% 
      as.data.frame()
    placFiltered = rbind(
      data.frame(seqnames = placFiltered$chr1, start = placFiltered$start1, end = placFiltered$end1, negLogFDR = -log10(placFiltered$fdr), negLog10P = placFiltered$ClusterNegLog10P, interactionID = placFiltered$interactionID, side = 1), 
      data.frame(seqnames = placFiltered$chr2, start = placFiltered$start2, end = placFiltered$end2, negLogFDR = -log10(placFiltered$fdr), negLog10P = placFiltered$ClusterNegLog10P, interactionID = placFiltered$interactionID, side = 2)
    )
    ovlp = as.data.frame(findOverlaps(GRanges(placFiltered), GRanges(placReduced)))
    placFiltered$bgEnsemble = placReduced$ensemble[ovlp$subjectHits]
  }
  
  # DLX/NKX REs
  {
    ## DLX
    {
      dlx = paste0(gene.grad.dir, "DLX_table.txt")
      dlx = read.table(dlx, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
      dlx = GRanges(dlx[, c("seqnames", "start", "end", "a.RE.bdg", "r.RE.bdg")])
      myChain = import.chain(chainMM9_10)
      dlx = data.frame(sort(unlist(liftOver(dlx, myChain))), stringsAsFactors = FALSE)
      
      ovlp = as.data.frame(findOverlaps(GRanges(placReduced), GRanges(dlx)))
      placReduced$dlx_aRE = placReduced$dlx_rRE = placReduced$dlx_RE = NA
      placReduced$dlx_aRE[ovlp$queryHits] = dlx$a.RE.bdg[ovlp$subjectHits]
      placReduced$dlx_rRE[ovlp$queryHits] = dlx$r.RE.bdg[ovlp$subjectHits]
      placReduced$dlx_RE = rowSums(placReduced[, c("dlx_aRE", "dlx_rRE")])
      placReduced$dlx_RE[placReduced$dlx_RE == 0] = "noSignal"
      placReduced$dlx_RE[placReduced$dlx_RE == 2] = "noSignal"
      placReduced$dlx_RE[placReduced$dlx_RE == 1 & placReduced$dlx_aRE == 1] = "a.RE"
      placReduced$dlx_RE[placReduced$dlx_RE == 1 & placReduced$dlx_rRE == 1] = "r.RE"
      placReduced$dlx_aRE = placReduced$dlx_rRE = NULL
    }
    
    ## NKX
    {
      nkx = paste0(gene.grad.dir, "NKX_table.txt")
      nkx = read.table(nkx, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
      nkx = GRanges(nkx[, c("chr", "start", "end", "RE.class")])
      nkx = data.frame(sort(unlist(liftOver(nkx, myChain))), stringsAsFactors = FALSE)
      
      ovlp = as.data.frame(findOverlaps(GRanges(placReduced), GRanges(nkx)))
      placReduced$nkx_RE = NA
      placReduced$nkx_RE[ovlp$queryHits] = nkx$RE.class[ovlp$subjectHits]
      placReduced$nkx_RE[placReduced$nkx_RE == 'NoChange'] = 'noSignal'
      placReduced$nkx_RE[placReduced$nkx_RE == 'aRE'] = 'a.RE'
      placReduced$nkx_RE[placReduced$nkx_RE == 'rRE'] = 'r.RE'
    }
  }
}

# overlap tested enhancers with cxPeaks, bgPeaks
{
  cxOvlp = do.call(ovlpEnhancers, list(cxPeaks, "cx", "no"))
  bgOvlp = do.call(ovlpEnhancers, list(bgPeaks, "bg", "yes"))
}

# annotate DKX/NKX REs data with bgPeaks
{
  ovlp = as.data.frame(findOverlaps(GRanges(bgPeaks), GRanges(dlx)))
  bgPeaks$dlx_aRE = bgPeaks$dlx_rRE = bgPeaks$dlx_RE = NA
  bgPeaks$dlx_aRE[ovlp$queryHits] = dlx$a.RE.bdg[ovlp$subjectHits]
  bgPeaks$dlx_rRE[ovlp$queryHits] = dlx$r.RE.bdg[ovlp$subjectHits]
  bgPeaks$dlx_RE = rowSums(bgPeaks[, c("dlx_aRE", "dlx_rRE")])
  bgPeaks$dlx_RE[bgPeaks$dlx_RE == 0] = "noSignal"
  bgPeaks$dlx_RE[bgPeaks$dlx_RE == 2] = "noSignal"
  bgPeaks$dlx_RE[bgPeaks$dlx_RE == 1 & bgPeaks$dlx_aRE == 1] = "a.RE"
  bgPeaks$dlx_RE[bgPeaks$dlx_RE == 1 & bgPeaks$dlx_rRE == 1] = "r.RE"
  bgPeaks$dlx_aRE = bgPeaks$dlx_rRE = NULL
  
  ovlp = as.data.frame(findOverlaps(GRanges(bgPeaks), GRanges(nkx)))
  bgPeaks$nkx_RE = NA
  bgPeaks$nkx_RE[ovlp$queryHits] = nkx$RE.class[ovlp$subjectHits]
  bgPeaks$nkx_RE[bgPeaks$nkx_RE == "aRE"] = "a.RE"
  bgPeaks$nkx_RE[bgPeaks$nkx_RE == "rRE"] = "r.RE"
  bgPeaks$nkx_RE[bgPeaks$nkx_RE == "NoChange"] = "noSignal"
}

# Annotate dlx_DE / nkx_DE to ensembleGenes and mm10Promoters
{
  ensembleGenes$dlx_DE = dlxDE$logFC[match(ensembleGenes$gene0, dlxDE$gene)]
  ensembleGenes$nkx_DE = nkxDE$logFC[match(ensembleGenes$gene0, nkxDE$gene)]
  ensembleGenes$dlx_DE[is.na(ensembleGenes$dlx_DE)] = 0
  ensembleGenes$nkx_DE[is.na(ensembleGenes$nkx_DE)] = 0
  
  mm10Promoters$dlx_DE = dlxDE$logFC[match(mm10Promoters$gene_id, dlxDE$gene)]
  mm10Promoters$nkx_DE = nkxDE$logFC[match(mm10Promoters$gene_id, nkxDE$gene)]
  mm10Promoters$dlx_DE[is.na(mm10Promoters$dlx_DE)] = 0
  mm10Promoters$nkx_DE[is.na(mm10Promoters$nkx_DE)] = 0
  bgGenes = mm10Promoters[!is.na(mm10Promoters$fpkm), ]
}


### Main figure
# 3A - - Heatmap bgClusters versus HMM states
{
  fontSize = 10
  # distal
  {
    distal = bgPeaks %>% 
      dplyr::filter(Feature != 'Promoter') %>% 
      dplyr::select(bgCluster, bgState) %>% 
      dplyr::filter(!is.na(bgCluster)) %>% 
      dplyr::group_by(bgCluster, bgState) %>% 
      tally() %>% 
      dplyr::group_by(bgCluster) %>% 
      dplyr::mutate(pct = 100 * n / sum(n)) %>% 
      dplyr::select(-n) %>% 
      tidyr::spread(key = bgCluster, value = pct, fill = 0) %>% 
      as.data.frame()
    rownames(distal) = distal$bgState
    distal$bgState = NULL
    distal = as.matrix(distal)
    myColors = colorRamp2(c(0, max(distal) / 2, max(distal)), c("black", "yellow", "red"))
    distalPlot = Heatmap(distal
                         , col = myColors
                         , cluster_rows = FALSE
                         , rect_gp = gpar(col = "grey30")
                         , cluster_columns = TRUE
                         , column_names_rot = 45
                         , column_names_gp = gpar(fontsize = fontSize) 
                         , column_names_side = "bottom"
                         , row_names_gp = gpar(fontsize = fontSize + 1) 
                         , row_names_side = "left"
                         , name = " ")
  }
  
  # random
  {
    random = randomSample %>% 
      dplyr::select(bgCluster, bgState) %>% 
      dplyr::group_by(bgCluster, bgState) %>% 
      tally() %>% 
      dplyr::group_by(bgCluster) %>% 
      dplyr::mutate(pct = 100 * n / sum(n)) %>% 
      dplyr::select(-n) %>% 
      tidyr::spread(key = bgCluster, value = pct, fill = 0) %>% 
      as.data.frame()
    rownames(random) = random$bgState
    random$bgState = NULL
    random = as.matrix(random)
    randomPlot = Heatmap(random
                         , col = myColors
                         , cluster_rows = FALSE
                         , rect_gp = gpar(col = "grey30")
                         , cluster_columns = TRUE
                         , column_names_rot = 45
                         , column_names_gp = gpar(fontsize = fontSize) 
                         , column_names_side = "bottom"
                         , row_names_side = "left"
                         , name = " ")
  }
  
  # proximal
  {
    proximal = bgPeaks %>% 
      dplyr::filter(Feature == 'Promoter') %>% 
      dplyr::select(bgCluster, bgState) %>% 
      dplyr::filter(!is.na(bgCluster)) %>% 
      dplyr::group_by(bgCluster, bgState) %>% 
      tally() %>% 
      dplyr::group_by(bgCluster) %>% 
      dplyr::mutate(pct = 100 * n / sum(n)) %>% 
      dplyr::select(-n) %>% 
      tidyr::spread(key = bgCluster, value = pct, fill = 0) %>% 
      as.data.frame()
    rownames(proximal) = proximal$bgState
    proximal$bgState = NULL
    proximal = as.matrix(proximal)
    
    # histones
    {
      myColors2 = colorRamp2(c(0, 1), c("#FFFFFF", "#000099"))
      histonePlot = rowAnnotation(
        H3K27me3 = bgEmission$H3K27me3
        , H3K27ac = bgEmission$H3K27ac
        , H3K4me3 = bgEmission$H3K4me3
        , H3K4me1 = bgEmission$H3K4me1
        , gp = gpar(col = "grey50")
        , show_legend = FALSE
        , annotation_name_rot = 45
        , annotation_name_gp = gpar(fontsize = fontSize) 
        , simple_anno_size = unit(0.45, "cm")
        , col = list(
          H3K27me3 = myColors2, 
          H3K27ac = myColors2, 
          H3K4me3 = myColors2, 
          H3K4me1 = myColors2
        )
      )
    }
    
    proximalPlot = Heatmap(proximal
                           , col = myColors
                           , cluster_rows = FALSE
                           , rect_gp = gpar(col = "grey30")
                           , cluster_columns = TRUE
                           , column_names_rot = 45
                           , column_names_gp = gpar(fontsize = fontSize) 
                           , column_names_side = "bottom"
                           , row_names_side = "left"
                           , name = " "
                           , right_annotation = histonePlot)
  }
  
  # compose plot
  {
    file.name = paste0(output.dir, "bgClusters.vs.HMM.pdf")
    pdf(file.name, height = 3, width = 7)
    totalPlot = distalPlot + randomPlot + proximalPlot
    draw(totalPlot, padding = unit(c(0.25, 0.25, 0.1, 0.05), "in"))
    dev.off()
  }
}

# 3B - GO analysis by proximity (GREAT)
{
  goToInclude = paste0(gene.grad.dir, "PlotGO.txt")
  goToInclude = read.table(goToInclude, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  goToInclude = goToInclude$ID[!is.na(goToInclude$Use)]
  
  greatGO = paste0(output.dir, "GREAT.GOterm.clusters.txt")
  greatGO = read.table(file.name, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
  
  # pre-process data
  {
    ceiling.fold = 4
    floor.fold = 1
    
    shortGREAT = greatGO %>% 
      dplyr::filter(ID %in% goToInclude, 
                    Cluster != "6_D") %>% 
      dplyr::mutate(GOterm = paste(ID, name)) %>% 
      dplyr::select(Cluster, name, ID, Binom_Fold_Enrichment, Binom_Adjp_BH) %>% 
      dplyr::mutate(Binom_Fold_Enrichment = pmin(Binom_Fold_Enrichment, ceiling.fold)) %>% 
      dplyr::mutate(Binom_Fold_Enrichment = pmax(Binom_Fold_Enrichment, floor.fold)) %>% 
      dplyr::mutate(log2Fold = log2(Binom_Fold_Enrichment)) %>% 
      as.data.frame()
    shortGREAT$name = sub("RNA polymerase", "RNApol", shortGREAT$name)
    shortGREAT$name = sub("transcription factor", "TF", shortGREAT$name)
    shortGREAT$name = sub("sequence", "seq", shortGREAT$name)
    shortGREAT$name = sub("specific", "specif", shortGREAT$name)
    shortGREAT$name = sub("binding", "bind", shortGREAT$name)
    shortGREAT$name = sub("distal enhancer", "dist enh", shortGREAT$name)
    shortGREAT$name = sub("proliferation", "prolif", shortGREAT$name)
    shortGREAT$name = sub("differentiation", "diff.", shortGREAT$name)
    shortGREAT$name = sub("activity", "activ", shortGREAT$name)
  }
  
  file.name = paste0(output.dir, "asd.GREAT.GO_heatmap.pdf")
  {
    heatData = shortGREAT %>% 
      dplyr::select(Cluster, name, Binom_Fold_Enrichment) %>% 
      spread(key = Cluster, value = Binom_Fold_Enrichment) %>% 
      as.data.frame()
    rownames(heatData) = heatData$name
    heatData$name = NULL
    pheatmap::pheatmap(heatData
                       , cluster_rows = TRUE
                       , cluster_cols = TRUE
                       , show_rownames = TRUE
                       , show_colnames = TRUE
                       , color = colorRampPalette(c("black", "yellow", "firebrick3"))(50)
                       , border_color = NA
                       , na_col = "black"
                       , angle_col = 45
                       , fontsize_row = 10
                       , fontsize_col = 9
                       , height = 5.5
                       , width = 8.25
                       , filename = file.name)  
  }
}

# 3C - Cluster cell type inference from Preissl et al. (2018) (https://www.nature.com/articles/s41593-018-0079-3)
{
  # by Alex
  # data input
  {
    total.bg = bgPeaks %>% 
      dplyr::filter(!is.na(bgCluster)) %>% 
      dplyr::group_by(bgCluster) %>% 
      tally() %>% 
      as.data.frame()
    bgNames = total.bg$bgCluster
    total.bg = as.list(total.bg$n)
    names(total.bg) = bgNames
    
    my.data = read.table(paste0(output.dir, "fileForAlex.all.txt"), sep = "\t", header = TRUE)
    total.snATAC = list(
      "K1" = 880
      , "K10" = 764
      , "K11" = 1238
      , "K12" = 623
      , "K13" = 1438
      , "K14" = 1112
      , "K15" = 698
      , "K2" = 1660
      , "K3" = 1838
      , "K4" = 1015
      , "K5" = 1276
      , "K6" = 444
      , "K7" = 1073
      , "K8" = 1263
      , "K9" = 1042
    )
    my.int.temp = table(my.data$bgCluster, my.data$atacCluster)
    my.int.bg = cbind(my.int.temp/unlist(total.bg), All = rowSums(my.int.temp) / unlist(total.bg))
    my.int.temp = table(my.data$bgCluster, my.data$atacCluster)
    my.int.bg = my.int.temp/unlist(total.bg)
    my.int.dev = table(my.data$bgCluster) / unlist(total.bg)
    
    my.int.temp = table(my.data$atacCluster, my.data$bgCluster)
    my.int.snATAC = cbind(my.int.temp / unlist(total.snATAC), All = rowSums(my.int.temp) / unlist(total.snATAC))
    my.int.dev.atac = table(my.data$atacCluster) / unlist(total.snATAC)
    
    my.data = read.table(paste0(output.dir, "fileForAlex.all.P56.txt"), sep = "\t", header = TRUE)
    total.snATAC = list(
      "K1" = 529
      , "K10" = 434
      , "K2" = 586
      , "K3" = 737
      , "K4" = 270
      , "K5" = 601
      , "K6" = 513
      , "K7" = 538
      , "K8" = 490
      , "K9" = 282
    )
    my.int.temp = table(my.data$bgCluster, my.data$atacCluster)
    my.int.bg = cbind(my.int.temp / unlist(total.bg), All = rowSums(my.int.temp) / unlist(total.bg))
    my.int.temp = table(my.data$bgCluster, my.data$atacCluster)
    my.int.bg = my.int.temp / unlist(total.bg)
    my.int.adult = table(my.data$bgCluster) / unlist(total.bg)
    
    my.int.temp = table(my.data$atacCluster, my.data$bgCluster)
    my.int.snATAC = cbind(my.int.temp / unlist(total.snATAC), All = rowSums(my.int.temp) / unlist(total.snATAC))
    
    my.data = read.table(paste0(output.dir, "fileForAlex.all.AllPeaks.txt"), sep = "\t", header = TRUE)
    my.int.temp = table(my.data$bgCluster)
    my.int.bg = my.int.temp / unlist(total.bg)
    my.int.all = my.int.temp / unlist(total.bg)
    sum.rep = rbind(
      All_ATAC = c(my.int.all, All = 0.6307127),
      CellType_Dev = c(my.int.dev, All = 0.2604835),
      CellType_P56 = c(my.int.adult, All = 0.02955555)
    )
    temp.nondev = sum.rep[1, ] - sum.rep[2, ]
    prom.cols = c(2, 13, 15, 17, 19, 21, 23)

    my.data = read.table(paste0(output.dir, "fileForAlex.all.txt"), sep = "\t", header = TRUE)
    my.int.temp = table(my.data$bgCluster, my.data$atacCluster)
    my.int.bg = my.int.temp / unlist(total.bg)
    select.rows = c("1_D", "2_D", "3_D", "4_D", "5_D", "7_D", "8_D", "9_D", "10_D", "12_D", "14_D", "15_D")
    
    use.rows = c(1, 12, 14, 16, 18, 20, 22, 24, 25, 3:11)
    my.int.final = cbind(my.int.bg, NonDev = unlist(sum.rep[1, ] - sum.rep[2, ])[-ncol(sum.rep)])
    my.int.tibble = as_tibble(cbind(BG_Cluster=rownames(my.int.final)[use.rows], my.int.final[use.rows, ]))
    my.int.tibble = my.int.tibble %>%
      mutate(across(2:ncol(my.int.tibble), as.numeric)) %>%
      mutate(BG_Cluster, fct_relevel(BG_Cluster, rownames(my.int.bg)[use.rows]))
    bg.order = c("1_D", "2_D", "12_D", "4_D", "10_D", "11_D", "13_D", "17_D", "9_D", "3_D", "7_D", "8_D", "14_D", "15_D",  "16_D", "5_D", "18_D", "6_D")
    sn.order = c("K1", "K2", "K3", "K4", "K5", "K6", "K7", "K8", "K9", "K10", "K11", "K12", "K13", "K14", "K15", "NonDev")
    my.int.tibble.long = my.int.tibble %>%
      pivot_longer(2:17, names_to = "snATAC_Cluster", values_to = "Proportion") %>%
      filter(BG_Cluster %in% bg.order) 
    
    # fix cluster names
    my.int.tibble.final = my.int.tibble.long %>%
      mutate(BG_Cluster=fct_relevel(BG_Cluster, rev(bg.order))) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K1", "Other", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K2", "Other", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K3", "Pal PM to Ex/AC transition (K3)", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K4", "Pal emerging Ex neuron (K4)", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K5", "SP+P neural progenitor (K5)", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K6", "Other", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K7", "Other", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K8", "SP early MSN (K8)", snATAC_Cluster)) %>%
      #  mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K8", "Other", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K9", "SP early IN neuron (K9)", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K10", "SP PM to neuron transition (K10)", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K11", "SP diff. IN neuron (K11)", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K12", "Pal early Ex neuron (K12)", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K13", "Pal diff. Ex neuron (K13)", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K14", "Other", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "K15", "Other", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster == "Other", "Other developmental cell type specific", snATAC_Cluster)) %>%
      mutate(snATAC_Cluster=ifelse(snATAC_Cluster =="NonDev", "snATAC, not developmental cell type", snATAC_Cluster)) %>%
      #  mutate(snATAC_Cluster = fct_relevel(snATAC_Cluster, c("Other", rev(sn.order))))
      mutate(snATAC_Cluster=fct_relevel(snATAC_Cluster, 
                                        rev(
                                          c("SP+P neural progenitor (K5)"
                                            ,"SP PM to neuron transition (K10)"
                                            ,"SP early IN neuron (K9)"
                                            ,"SP diff. IN neuron (K11)"
                                            ,"SP early MSN (K8)"
                                            ,"Pal PM to Ex/AC transition (K3)"
                                            ,"Pal emerging Ex neuron (K4)"
                                            ,"Pal early Ex neuron (K12)"
                                            ,"Pal diff. Ex neuron (K13)"
                                            ,"Other developmental cell type specific"
                                            ,"snATAC, not developmental cell type"
                                          ))))
  }
  
  # plot snATAC-seq against clusters
  {
    use.colors = c(
      # generic neural progenitor (K5)
      "#8A2BE2", #blue violet (mix red+blue)
      # SP palette blues  
      "#0000FF", #  blue PM to neuron (K10)
      "#4682B4", # steel blue early IN neuron (K9)
      "#1E90FF", # dodger blue early IN neuron (K9)
      "#87CEFA", # light sky blue early MSN  (K8)
      # Pallial pallette reds
      "#800000", # maroon, Pal progen K3
      "#DC143C", # crimson, Pal emerging neuron K4
      "#FF7F50", # coral, Pal early neuron K12
      "#FFA07A", # light salmon Pal diff neuron K13
      # other  
      "#000000", # black (Dev)
      "#D3D3D3" # light grey (non-Dev) 
    )
    
    ggplot(my.int.tibble.final, aes(x = BG_Cluster,  y = Proportion, fill = snATAC_Cluster)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      xlab("cluster") +
      ylab("Proportion overlapping snATAC-seq") + 
      scale_fill_manual(values = rev(use.colors)) + 
      guides(fill = guide_legend(reverse = TRUE)) + 
      theme_pubr() +
      theme(text = element_text(size = 14), 
            panel.grid.major.x = element_line(linetype = "dashed", color = "black"), 
            legend.position = "right")
    ggsave(filename = paste0(output.dir, "Curated_BG_to_snATAC.pdf"), width = 8, height = 5)
  }
}

# 3D - Percent peaks in PSC vs. ensemble size per cluster
{
  tempObj0 = bgPeaks %>% 
    dplyr::filter(!is.na(bgCluster)) %>% 
    dplyr::select(bgCluster, location, bgEnsembleSize) %>% 
    dplyr::group_by(bgCluster, location) %>% 
    tally() %>% 
    dplyr::group_by(bgCluster) %>% 
    dplyr::mutate(pctInPSC = 100 * n / sum(n)) %>% 
    dplyr::filter(location == "PSC") %>% 
    ungroup() %>% 
    dplyr::select(bgCluster, pctInPSC) %>% 
    as.data.frame()
  tempObj1 = bgPeaks %>% 
    dplyr::filter(!is.na(bgCluster)) %>% 
    dplyr::select(bgCluster, location, bgEnsembleSize) %>% 
    dplyr::filter(location == "PSC") %>% 
    dplyr::group_by(bgCluster) %>% 
    dplyr::summarise(meanSize = mean(bgEnsembleSize)) %>% 
    as.data.frame()
  full_join(tempObj0, tempObj1, by = "bgCluster") %>% 
    mutate(loc = ifelse(grepl("_D", bgCluster), "distal", "proximal"), 
           bgCluster = factor(bgCluster, levels = gtools::mixedsort(unique(bgCluster)))) %>% 
    ggplot(aes(x = pctInPSC, y = meanSize, color = bgCluster, label = bgCluster)) +
    geom_point(aes(shape = loc), size = 4) +
    geom_text_repel(box.padding = 1, size = 5) +
    theme_classic() + 
    xlab("percent peaks in PSC") +
    ylab("mean ensemble size") +
    theme(legend.position = "none", 
          legend.title = element_blank(), 
          axis.text = element_text(size = 20), 
          axis.title = element_text(size = 20))
  file.name = paste0(output.dir, "ensembleSize.vs.pctPeaksInPSC.pdf")
  ggsave(file.name, width = 6, height = 6.5)
  remove(list = c("tempObj0", "tempObj1"))
}

# 3E - BG vertebrate evolutionary conservation
{
  bgPhastCons = data.frame()
  for (myCluster in unique(bgPeaks$bgCluster[!is.na(bgPeaks$bgCluster)])) {
    tempObj = bgPeaks[bgPeaks$bgCluster == myCluster, c("VertebratePhastcons", "Feature")] %>% 
      dplyr::mutate(Feature = ifelse(Feature == "Promoter", "Proximal", "Distal")) %>% 
      as.data.frame()
    mySample = replicate(nPerm, mean(sample(tempObj$VertebratePhastcons, 
                                            size = samplingRatio * nrow(tempObj), 
                                            replace = FALSE), na.rm = TRUE))
    tempObj = data.frame(Feature = tempObj$Feature[1], cluster = myCluster, conservMeans = mySample, stringsAsFactors = FALSE)
    if (nrow(bgPhastCons) == 0) {
      bgPhastCons = tempObj
    } else {
      bgPhastCons = rbind(bgPhastCons, tempObj)
    }
    print(noquote(myCluster))
    remove(list = c("tempObj", "mySample"))
  }
  
  # random sample
  for (myCluster in unique(randomSample$bgCluster)) {
    tempObj = randomSample[randomSample$bgCluster == myCluster, c("VertebratePhastcons", "Feature")] %>% 
      dplyr::mutate(Feature = ifelse(Feature == "Promoter", "Proximal", "Distal")) %>% 
      as.data.frame()
    for (myLoc in c("Proximal", "Distal")) {
      mySize = bgPeaks %>% mutate(Feature = ifelse(Feature == "Promoter", "Proximal", "Distal")) %>% 
        dplyr::filter(!is.na(bgCluster), Feature == myLoc) %>% dplyr::group_by(bgCluster) %>% tally() %>% arrange(desc(n)) %>% dplyr::select(n) %>% unlist() %>% mean() %>% round(0)
      mySample = replicate(nPerm, mean(sample(tempObj$VertebratePhastcons[tempObj$Feature == myLoc], 
                                              size = samplingRatio * (ifelse(mySize < length(tempObj$VertebratePhastcons[tempObj$Feature == myLoc]), mySize, length(tempObj$VertebratePhastcons[tempObj$Feature == myLoc]))), 
                                              replace = FALSE), na.rm = TRUE))
      tempObj0 = data.frame(Feature = myLoc, cluster = paste(myLoc, myCluster), conservMeans = mySample, stringsAsFactors = FALSE)
      if (nrow(bgPhastCons) == 0) {
        bgPhastCons = tempObj0
      } else {
        bgPhastCons = rbind(bgPhastCons, tempObj0)
      }
      remove(tempObj0)
    }
    
    print(noquote(myCluster))
    remove(list = c("tempObj", "mySample"))
  }
  
  myLevels = bgPhastCons %>% dplyr::group_by(cluster, Feature) %>% summarise(mean = mean(conservMeans)) %>% arrange(mean, cluster) %>% dplyr::select(cluster) %>% unlist()
  bgPhastCons$cluster = factor(bgPhastCons$cluster, levels = myLevels)
  
  ggplot(data = bgPhastCons, aes(x = cluster, y = conservMeans, fill = cluster)) +
    geom_boxplot(outlier.size = 0.25, size = 0.25) +
    theme_classic() +
    ylab("means of PhastCons scores") +
    xlab("cluster") +
    ylim(c(350, 700)) +
    coord_flip() +
    theme(axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 10), 
          axis.text.x = element_text(size = 10), 
          axis.line = element_line(size = 0.5), 
          legend.position = "none") +
    facet_wrap(~ Feature, ncol = 1, scales = "free_y") +
    theme(strip.background = element_rect(color = "white", fill = "navyblue", size = 2, linetype = "solid"), 
          strip.placement = "outside", 
          strip.text = element_text(face = "bold", colour = "white"))
  file.name = paste0(output.dir, "bgVertebratePhastCons.pdf")
  ggsave(file.name, height = 6, width = 4)
}


### Supplementary
# S3a - GO term analysis in PSCs and inside loops (curated)
{
  # PSCs and inside loops
  {
    # Peak data
    {
      # from bgPeaks
      tempObjP = bgPeaks %>% 
        dplyr::select(bgCluster, inEnsemble, location) %>% 
        dplyr::filter(location != "outside", 
                      !is.na(bgCluster)) %>% 
        as.data.frame()
      
      # from randomSample
      tempObjR = randomSample %>% 
        dplyr::select(bgCluster, inEnsemble, location) %>% 
        dplyr::filter(location != "outside", 
                      !is.na(bgCluster)) %>% 
        as.data.frame()
      
      # combine sets
      tempObj = rbind(tempObjP, tempObjR)
      remove(list = c("tempObjP", "tempObjR"))
    }
    
    # GO analysis
    totalGO = data.frame()
    for (myLoc in unique(tempObj$location)) {
      for (myCluster in gtools::mixedsort(unique(tempObj$bgCluster))) {
        print(noquote(paste(myLoc, myCluster)))
        if (myCluster == 'random') {
          tempGO = tempObj[tempObj$bgCluster == myCluster, ]
        } else {
          tempGO = tempObj[tempObj$bgCluster == myCluster & tempObj$location == myLoc, ]
        }
        tempGO = tempGO[!is.na(tempGO[,1]), ]
        if (nrow(tempGO) == 0) next
        GOgenes = mm10Promoters$gene_id[mm10Promoters$inEnsemble %in% unique(tempGO$inEnsemble)]
        myEnrich = enrichGO(gene = GOgenes
                            , universe = unique(mm10Promoters$gene_id)
                            , OrgDb = org.Mm.eg.db
                            , keyType = "SYMBOL"
                            , ont = "BP"
                            , pAdjustMethod = "BH"
                            , pvalueCutoff = 0.05
                            , qvalueCutoff = 0.05
                            , minGSSize = 10
                            , maxGSSize = 500)
        
        if (!is.null(myEnrich)) enrichSimple = clusterProfiler::simplify(myEnrich, 
                                                                         cutoff     = 0.7,
                                                                         by         = "p.adjust",
                                                                         select_fun = min, 
                                                                         measure = "Wang", 
                                                                         semData = NULL)
        if (nrow(enrichSimple@result) > 0) {
          tempGO = data.frame(cluster = myCluster, loc = myLoc, enrichSimple@result, stringsAsFactors = FALSE)
          if (nrow(totalGO) == 0) {
            totalGO = tempGO
          } else {
            totalGO = rbind(totalGO, tempGO)
          }
        }
      }
    }
    
    totalGO.tbl = totalGO %>% 
      dplyr::filter(p.adjust <= 0.05) %>% 
      rowwise() %>% 
      dplyr::mutate(GeneRatio = eval(parse(text = GeneRatio)), 
                    BgRatio = eval(parse(text = BgRatio))) %>% 
      dplyr::mutate(Enrichment = GeneRatio / BgRatio, 
                    GOterm = paste(ID, Description), 
                    totCluster = paste(cluster, loc)) %>% 
      dplyr::select(totCluster, GOterm, Enrichment) %>% 
      tidyr::spread(key = totCluster, value = Enrichment, fill = 0) %>% 
      dplyr::mutate(include = 0) %>% 
      as.data.frame()
    
    for (myLoc in c("PSC", "inside")) {
      goToInclude = paste0(output.dir, "goTermsToChoosePLAC.curated.txt")
      if (FALSE) write.table(totalGO.tbl, goToInclude, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
      
      if (file.exists(goToInclude)) {
        goToInclude = read.table(goToInclude, sep = '\t', header = TRUE, quote = "", stringsAsFactors = FALSE)
        goToInclude = goToInclude$GOterm[goToInclude$include == 1]
      } else {
        stop("include GO terms to include")
      }
      ## common GO terms
      {
        totalGO.tbl = totalGO %>% 
          dplyr::filter(p.adjust <= 0.05, 
                        cluster != "6_D", 
                        loc == myLoc | cluster == 'random') %>% 
          rowwise() %>% 
          dplyr::mutate(GeneRatio = eval(parse(text = GeneRatio)), 
                        BgRatio = eval(parse(text = BgRatio))) %>% 
          dplyr::mutate(Enrichment = GeneRatio / BgRatio, 
                        GOterm = paste(ID, Description)) %>% 
          dplyr::filter(GOterm %in% goToInclude) %>% 
          dplyr::select(cluster, GOterm, Enrichment) %>% 
          dplyr::group_by(cluster, GOterm) %>% 
          summarise(Enrichment = mean(Enrichment)) %>% 
          tidyr::spread(key = cluster, value = Enrichment, fill = 0) %>% 
          as.data.frame()
        totalGO.tbl$chk = rowSums(totalGO.tbl[, -1]) > (ncol(totalGO.tbl) - 0) / 2
        totalGO.tbl = totalGO.tbl[totalGO.tbl$chk, grep("chk", colnames(totalGO.tbl), invert = TRUE)]
        rownames(totalGO.tbl) = totalGO.tbl$GOterm
        totalGO.tbl$GOterm = NULL
        
        file.name = paste0(output.dir, "placSeq.GOheatmap.", myLoc, ".curated.common.pdf")
        pheatmap::pheatmap(totalGO.tbl
                           , cluster_rows = TRUE
                           , cluster_cols = TRUE
                           , show_rownames = TRUE
                           , show_colnames = TRUE
                           , color = colorRampPalette(c("black", "yellow", "red"), bias = 1)(100)
                           , border_color = "grey10"
                           , angle_col = 45
                           , fontsize_row = 10
                           , fontsize_col = 10
                           , height = ifelse(myLoc == 'PSC', 2.75, 2)
                           , width = ifelse(myLoc == 'PSC', 13, 12)
                           , filename = file.name)  
      }
      
      ## distinguishing GO terms
      {
        totalGO.tbl = totalGO %>% 
          dplyr::filter(p.adjust <= 0.05, 
                        cluster != "6_D", 
                        loc == myLoc | cluster == 'random') %>% 
          rowwise() %>% 
          dplyr::mutate(GeneRatio = eval(parse(text = GeneRatio)), 
                        BgRatio = eval(parse(text = BgRatio))) %>% 
          dplyr::mutate(Enrichment = GeneRatio / BgRatio, 
                        GOterm = paste(ID, Description)) %>% 
          dplyr::filter(GOterm %in% goToInclude) %>% 
          dplyr::select(cluster, GOterm, Enrichment) %>% 
          dplyr::group_by(cluster, GOterm) %>% 
          summarise(Enrichment = mean(Enrichment)) %>% 
          tidyr::spread(key = cluster, value = Enrichment, fill = 0) %>% 
          as.data.frame()
        totalGO.tbl$chk = rowSums(totalGO.tbl[, -1]) < (ncol(totalGO.tbl) - 0) / 2
        totalGO.tbl = totalGO.tbl[totalGO.tbl$chk, grep("chk", colnames(totalGO.tbl), invert = TRUE)]
        rownames(totalGO.tbl) = totalGO.tbl$GOterm
        totalGO.tbl$GOterm = NULL
        totalGO.tbl = totalGO.tbl[, colSums(totalGO.tbl[, ]) > 0]
        
        file.name = paste0(output.dir, "placSeq.GOheatmap.", myLoc, ".curated.unique.pdf")
        pheatmap::pheatmap(totalGO.tbl
                           , cluster_rows = TRUE
                           , cluster_cols = TRUE
                           , show_rownames = TRUE
                           , show_colnames = TRUE
                           , color = colorRampPalette(c("black", "yellow", "red"), bias = 0.875)(100)
                           , border_color = "grey10"
                           , angle_col = 45
                           , fontsize_row = 10
                           , fontsize_col = 10
                           , height = ifelse(myLoc == 'PSC', 4, 4)
                           , width = ifelse(myLoc == 'PSC', 8.5, 7.5)
                           , filename = file.name)  
      }
      
      ## reduced overall GO terms (curated by Alex)
      {
        goToInclude = paste0(output.dir, "goTermsToChoosePLAC.curatedAlex.txt")
        goToInclude = read.table(goToInclude, sep = '\t', header = TRUE, quote = "", stringsAsFactors = FALSE)
        goToInclude = goToInclude$GOterm[goToInclude$include == 1]
        
        totalGO.tbl = totalGO %>% 
          dplyr::filter(p.adjust <= 0.05, 
                        cluster != "6_D", 
                        loc == myLoc | cluster == 'random') %>% 
          rowwise() %>% 
          dplyr::mutate(GeneRatio = eval(parse(text = GeneRatio)), 
                        BgRatio = eval(parse(text = BgRatio))) %>% 
          dplyr::mutate(Enrichment = GeneRatio / BgRatio, 
                        GOterm = paste(ID, Description)) %>% 
          dplyr::filter(GOterm %in% goToInclude) %>% 
          dplyr::select(cluster, GOterm, Enrichment) %>% 
          dplyr::group_by(cluster, GOterm) %>% 
          summarise(Enrichment = mean(Enrichment)) %>% 
          tidyr::spread(key = cluster, value = Enrichment, fill = 0) %>% 
          as.data.frame()
        rownames(totalGO.tbl) = totalGO.tbl$GOterm
        totalGO.tbl$GOterm = NULL
        
        file.name = paste0(output.dir, "placSeq.GOheatmap.", myLoc, ".curated.Alex.pdf")
        pheatmap::pheatmap(totalGO.tbl
                           , cluster_rows = TRUE
                           , cluster_cols = TRUE
                           , show_rownames = TRUE
                           , show_colnames = TRUE
                           , color = colorRampPalette(c("black", "yellow", "red"), bias = 1)(100)
                           , border_color = "grey10"
                           , angle_col = 45
                           , fontsize_row = 10
                           , fontsize_col = 10
                           , height = ifelse(myLoc == 'PSC', 3.75, 3.75)
                           , width = ifelse(myLoc == 'PSC', 9.25, 8)
                           , filename = file.name)  
      }
    }
  }
}

# S3b - Representation of Oct4, Sox2, Nanog, Klf4 motif relative enrichment
{
  # Alex's plot
}

# S3c - REs on clusters from DLX, NKX studies
{
  file.name = paste0(output.dir, "bothREsInClusters.pdf")
  ovlp = as.data.frame(findOverlaps(GRanges(bgPeaks), GRanges(dlx)))
  bgPeaks$dlx_aRE = bgPeaks$dlx_rRE = bgPeaks$dlx_RE = NA
  bgPeaks$dlx_aRE[ovlp$queryHits] = dlx$a.RE.bdg[ovlp$subjectHits]
  bgPeaks$dlx_rRE[ovlp$queryHits] = dlx$r.RE.bdg[ovlp$subjectHits]
  bgPeaks$dlx_RE = rowSums(bgPeaks[, c("dlx_aRE", "dlx_rRE")])
  bgPeaks$dlx_RE[bgPeaks$dlx_RE == 0] = "noSignal"
  bgPeaks$dlx_RE[bgPeaks$dlx_RE == 2] = "noSignal"
  bgPeaks$dlx_RE[bgPeaks$dlx_RE == 1 & bgPeaks$dlx_aRE == 1] = "a.RE"
  bgPeaks$dlx_RE[bgPeaks$dlx_RE == 1 & bgPeaks$dlx_rRE == 1] = "r.RE"
  bgPeaks$dlx_aRE = bgPeaks$dlx_rRE = NULL
  
  ovlp = as.data.frame(findOverlaps(GRanges(bgPeaks), GRanges(nkx)))
  bgPeaks$nkx_RE = NA
  bgPeaks$nkx_RE[ovlp$queryHits] = nkx$RE.class[ovlp$subjectHits]
  bgPeaks$nkx_RE[bgPeaks$nkx_RE == 'NoChange'] = 'noSignal'
  bgPeaks$nkx_RE[bgPeaks$nkx_RE == 'aRE'] = 'a.RE'
  bgPeaks$nkx_RE[bgPeaks$nkx_RE == 'rRE'] = 'r.RE'
  
  tempObj = bgPeaks %>% 
    dplyr::filter(!is.na(bgCluster)) %>% 
    as.data.frame()
  dlxObj = tempObj %>% 
    dplyr::filter(!is.na(dlx_RE)) %>% 
    dplyr::select(bgCluster, dlx_RE) %>% 
    group_by(bgCluster, dlx_RE) %>% 
    tally() %>% 
    as.data.frame()
  dlxObj$dlx_RE = paste0("dlx_", dlxObj$dlx_RE)
  dlxObj = dlxObj %>% 
    group_by(bgCluster, dlx_RE) %>% 
    dplyr::summarise(n = sum(n)) %>%
    tidyr::spread(key = dlx_RE, value = n, fill = 0) %>% 
    rowwise() %>% 
    dplyr::mutate(totalDLX = dlx_a.RE + dlx_noSignal + dlx_r.RE) %>% 
    as.data.frame()
  colnames(dlxObj)[c(2, 4)] = c("obs_dlxA", "obs_dlxR")
  dlxObj$exp_dlxA = round(sum(dlxObj$obs_dlxA) / sum(dlxObj$totalDLX) * dlxObj$totalDLX, 0)
  dlxObj$exp_dlxR = round(sum(dlxObj$obs_dlxR) / sum(dlxObj$totalDLX) * dlxObj$totalDLX, 0)
  dlxObj$dlx_a.RE = ifelse(abs(dlxObj$obs_dlxA - dlxObj$exp_dlxA) > minDiff, 
                           log2((dlxObj$obs_dlxA + 0.001) / (dlxObj$exp_dlxA + 0.001)), NA)
  dlxObj$dlx_r.RE = ifelse(abs(dlxObj$obs_dlxR - dlxObj$exp_dlxR) > minDiff, 
                           log2((dlxObj$obs_dlxR + 0.001) / (dlxObj$exp_dlxR + 0.001)), NA)
  dlxObj$obs_dlxA = dlxObj$dlx_noSignal = dlxObj$obs_dlxR = dlxObj$totalDLX = dlxObj$exp_dlxA = dlxObj$exp_dlxR = NULL
  
  nkxObj = tempObj %>% 
    dplyr::filter(!is.na(nkx_RE)) %>% 
    dplyr::select(bgCluster, nkx_RE) %>% 
    group_by(bgCluster, nkx_RE) %>% 
    tally() %>% 
    as.data.frame()
  nkxObj$nkx_RE = paste0("nkx_", nkxObj$nkx_RE)
  nkxObj = nkxObj %>% 
    group_by(bgCluster, nkx_RE) %>% 
    dplyr::summarise(n = sum(n)) %>%
    tidyr::spread(key = nkx_RE, value = n, fill = 0) %>% 
    mutate(totalNKX = nkx_a.RE + nkx_noSignal + nkx_r.RE) %>% 
    as.data.frame()
  colnames(nkxObj)[c(2, 4)] = c("obs_nkxA", "obs_nkxR")
  nkxObj$exp_nkxA = round(sum(nkxObj$obs_nkxA) / sum(nkxObj$totalNKX) * nkxObj$totalNKX, 0)
  nkxObj$exp_nkxR = round(sum(nkxObj$obs_nkxR) / sum(nkxObj$totalNKX) * nkxObj$totalNKX, 0)
  nkxObj$nkx_a.RE = ifelse(abs(nkxObj$obs_nkxA - nkxObj$exp_nkxA) > minDiff, 
                           log2((nkxObj$obs_nkxA + 0.001) / (nkxObj$exp_nkxA + 0.001)), NA)
  nkxObj$nkx_r.RE = ifelse(abs(nkxObj$obs_nkxR - nkxObj$exp_nkxR) > minDiff, 
                           log2((nkxObj$obs_nkxR + 0.001) / (nkxObj$exp_nkxR + 0.001)), NA)
  nkxObj$obs_nkxA = nkxObj$nkx_noSignal = nkxObj$obs_nkxR = nkxObj$totalNKX = nkxObj$exp_nkxA = nkxObj$exp_nkxR = NULL
  
  tempObj = full_join(dlxObj, nkxObj, by = "bgCluster")
  myLevels = c(
    gtools::mixedsort(tempObj$bgCluster[grep("_D", tempObj$bgCluster)]), 
    gtools::mixedsort(tempObj$bgCluster[grep("_P", tempObj$bgCluster)])
  )
  tempObj$bgCluster = factor(tempObj$bgCluster, levels = myLevels)
  tempObj = tempObj[order(tempObj$bgCluster), ]
  rownames(tempObj) = tempObj$bgCluster
  tempObj$bgCluster = NULL
  myColor = colorRampPalette(c("blue", "white", "red"))(100)
  myBreaks = c(seq(-scaleLimit, 0, length.out = 51), seq(scaleLimit/100, scaleLimit, length.out = 50))
  pheatmap::pheatmap(tempObj
                     , height = 6
                     , width = 3
                     , color = myColor
                     , breaks = myBreaks
                     , border_color = "grey30"
                     , na_col = "grey80"
                     , fontsize_row = 11
                     , fontsize_col = 11
                     , cluster_rows = FALSE
                     , cluster_cols = FALSE
                     , show_rownames = TRUE
                     , angle_col = 45
                     , filename = file.name)
}

# S3d - Likelihood of DLX/NKX REs in ensembles
{
  file.name = paste0(output.dir, "RE.effect.onEnsembles.pdf")
  rownames(bgPeaks) = NULL
  mySample = as.data.frame(replicate(1000, sample(rownames(bgPeaks), size = 0.75 * nrow(bgPeaks), replace = FALSE)))
  InOuts = data.frame()
  for (myColumn in colnames(mySample)) {
    # DLX
    {
      tempObj = bgPeaks[mySample[, myColumn], c("PeakID", "location", "dlx_RE")]
      tempObj = tempObj %>% 
        dplyr::group_by(location, dlx_RE) %>% 
        tally() %>% 
        spread(key = location, value = n) %>% 
        mutate(pctIn = inside / sum (inside), 
               pctOut = outside / sum (outside)) %>% 
        mutate(logIn_Out = log2(pctIn / pctOut)) %>% 
        dplyr::select(dlx_RE, logIn_Out) %>% 
        gather(key = "dlx_RE", value = "logIn_Out") %>% 
        dplyr::filter(!is.na(dlx_RE)) %>% 
        as.data.frame()
      tempObj$dlx_RE = paste0("DLX.", tempObj$dlx_RE)
      colnames(tempObj)[1] = "feature"
      if (nrow(InOuts) == 0) {
        InOuts = tempObj
      } else {
        InOuts = rbind(InOuts, tempObj)
      }
    }
    
    # NKX
    {
      tempObj = bgPeaks[mySample[, myColumn], c("PeakID", "location", "nkx_RE")]
      tempObj = tempObj %>% 
        dplyr::group_by(location, nkx_RE) %>% 
        tally() %>% 
        spread(key = location, value = n) %>% 
        mutate(pctIn = inside / sum (inside), 
               pctOut = outside / sum (outside)) %>% 
        mutate(logIn_Out = log2(pctIn / pctOut)) %>% 
        dplyr::select(nkx_RE, logIn_Out) %>% 
        gather(key = "nkx_RE", value = "logIn_Out") %>% 
        dplyr::filter(!is.na(nkx_RE)) %>% 
        as.data.frame()
      tempObj$nkx_RE = paste0("NKX.", tempObj$nkx_RE)
      colnames(tempObj)[1] = "feature"
      if (nrow(InOuts) == 0) {
        InOuts = tempObj
      } else {
        InOuts = rbind(InOuts, tempObj)
      }
    }
    remove(tempObj)
  }
  myLevels = c("NKX.r.RE", "NKX.noSignal", "NKX.a.RE", "DLX.r.RE", "DLX.noSignal", "DLX.a.RE")
  InOuts$feature = factor(InOuts$feature, levels = myLevels)
  InOuts = InOuts[order(InOuts$feature), ]
  
  ggplot(data = InOuts, aes(x = feature, y = logIn_Out, fill = feature)) +
    geom_boxplot(alpha = 0.85, color = 'black') +
    xlab(" ") +
    ylab("log2 likelihood in ensembles") +
    scale_fill_manual(values = rep(c("red", "blue", "green"), 2)) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title = element_text(size = 14), 
          axis.text = element_text(size = 14))
  ggsave(file.name, height = 3, width = 4)
  remove(InOuts)
}

# S3e - Percent CX peaks in BG clusters
{
  file.name = paste0(output.dir, "CX_in_BG.ensembles.pdf")
  ovlp = as.data.frame(findOverlaps(GRanges(bgPeaks), GRanges(cxPeaks)))
  cxPeaks$bgCluster[ovlp$subjectHits] = bgPeaks$bgCluster[ovlp$queryHits]
  cxPeaks$bgEnsemble[ovlp$subjectHits] = bgPeaks$bgEnsemble[ovlp$queryHits]
  
  ensembleRatios = cxPeaks %>% 
    dplyr::group_by(!is.na(bgEnsemble)) %>% 
    tally() %>% 
    mutate(pct = round(100 * n / sum(n), 1)) %>% 
    as.data.frame()
  colnames(ensembleRatios)[1] = 'peaks in BG contacts'
  
  fontSize = 11
  ensembleRatios %>% 
    dplyr::arrange(desc(`peaks in BG contacts`)) %>% 
    dplyr::mutate(prop = pct / sum(pct) * 100) %>%
    dplyr::mutate(ypos = cumsum(prop)- 0.5 * prop) %>% 
    ggplot(aes(x = "", y = pct, fill = `peaks in BG contacts`)) +
    geom_bar(stat = "identity", position = "stack", color = "grey50") +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    coord_polar("y", start = 0) +
    theme_void() + 
    theme(legend.position = "bottom", 
          legend.text = element_text(size = fontSize + 10), 
          legend.title = element_text(size = fontSize + 10)) +
    geom_text(aes(y = ypos, label = paste0(pct, "%")), size = fontSize) +
    guides(guide_legend(nrow = 2))
  ggsave(file.name, height = 4, width = 7)
}

# S3f - Percent cxPeaks in BG clusters
{
  ovlp = as.data.frame(findOverlaps(GRanges(bgPeaks), GRanges(cxPeaks)))
  cxPeaks$bgCluster[ovlp$subjectHits] = bgPeaks$bgCluster[ovlp$queryHits]
  cxPeaks$bgEnsemble[ovlp$subjectHits] = bgPeaks$bgEnsemble[ovlp$queryHits]
  
  ensemblePercent = cxPeaks %>% 
    dplyr::filter(!is.na(bgEnsemble), !is.na(bgCluster)) %>% 
    dplyr::group_by(bgCluster) %>% 
    tally() %>% 
    mutate(pct = round(100 * n / sum(n), 1)) %>% 
    as.data.frame()
  myLevels = ensemblePercent %>% arrange(desc(pct)) %>% dplyr::select(bgCluster) %>% unlist()
  ensemblePercent$bgCluster = factor(ensemblePercent$bgCluster, levels = myLevels)
  
  file.name = paste0(output.dir, "CX.BGclusters.onBGensembles.landscape.pdf")
  ggplot(data = ensemblePercent, aes(x = bgCluster, y = pct)) +
    geom_bar(stat = "identity", fill = "lightgreen", color = "grey20") +
    xlab("cluster") +
    ylab("percent of CX peaks in BG clusters") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
  ggsave(file.name, height = 3, width = 5)
  
  myLevels = ensemblePercent %>% arrange(pct) %>% dplyr::select(bgCluster) %>% unlist()
  ensemblePercent$bgCluster = factor(ensemblePercent$bgCluster, levels = myLevels)
  
  file.name = paste0(output.dir, "CX.BGclusters.onBGensembles.portrait.pdf")
  ggplot(data = ensemblePercent, aes(x = bgCluster, y = pct)) +
    geom_bar(stat = "identity", fill = "lightgreen", color = "grey20") +
    xlab("cluster") +
    ylab("percent of CX peaks in BG clusters") +
    coord_flip() +
    theme_classic() +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 14))
  ggsave(file.name, height = 5, width = 3)
}

# S3g - Percent of BG clusters sharing cxPeaks
{
  ovlp = as.data.frame(findOverlaps(GRanges(bgPeaks), GRanges(cxPeaks)))
  cxPeaks$bgCluster[ovlp$subjectHits] = bgPeaks$bgCluster[ovlp$queryHits]
  cxPeaks$bgEnsemble[ovlp$subjectHits] = bgPeaks$bgEnsemble[ovlp$queryHits]
  
  numBGclusters = bgPeaks %>% 
    dplyr::filter(!is.na(bgCluster)) %>% 
    dplyr::group_by(bgCluster) %>% 
    tally() %>% 
    as.data.frame()
  
  ensemblePercent = cxPeaks %>% 
    dplyr::filter(!is.na(bgEnsemble), !is.na(bgCluster)) %>% 
    dplyr::group_by(bgCluster) %>% 
    tally() %>% 
    mutate(numBG = numBGclusters$n[match(bgCluster, numBGclusters$bgCluster)]) %>% 
    mutate(pct = round(100 * n / numBG, 1)) %>% 
    as.data.frame()
  
  myLevels = ensemblePercent %>% arrange(desc(pct)) %>% dplyr::select(bgCluster) %>% unlist()
  ensemblePercent$bgCluster = factor(ensemblePercent$bgCluster, levels = myLevels)
  
  file.name = paste0(output.dir, "sharedPropCX.BGclusters.onBGensembles.landscape.pdf")
  ggplot(data = ensemblePercent, aes(x = bgCluster, y = pct)) +
    geom_bar(stat = "identity", fill = "lightgreen", color = "grey20") +
    xlab("cluster") +
    ylab("percent of BG peaks in shared BG/CX loci") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
  ggsave(file.name, height = 3, width = 5)
  
  myLevels = ensemblePercent %>% arrange(pct) %>% dplyr::select(bgCluster) %>% unlist()
  ensemblePercent$bgCluster = factor(ensemblePercent$bgCluster, levels = myLevels)
  
  file.name = paste0(output.dir, "sharedPropCX.BGclusters.onBGensembles.portrait.pdf")
  ggplot(data = ensemblePercent, aes(x = bgCluster, y = pct)) +
    geom_bar(stat = "identity", fill = "lightgreen", color = "grey20") +
    xlab("cluster") +
    ylab("percent of BG peaks in shared BG/CX loci") +
    coord_flip() +
    theme_classic() +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 14))
  ggsave(file.name, height = 5, width = 3)
}
