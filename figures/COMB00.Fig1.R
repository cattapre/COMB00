# Load packages
{
  initial.packages = c("GenomicFeatures", "tidyverse", "pheatmap", "ggcorrplot", "ggseqlogo", "smooth", 
                       "Mcomp", "XML", "rlist", "ggpubr", "ggimage", "webr","ComplexHeatmap")
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
}

# Load parameters
{
  allReps = FALSE
  region = "BG"
  nPerm = 10000
  samplingRatio = 0.75
  minTargetThrsh = 10
  strippedVersion = TRUE
  numHMM = 9
  source("COMB00.directories.R")
  
  bgTFs = c("ARX", "ASCL1", "DLX1", "DLX2", "DLX5", "GSX2", "LHX6", "NKX2.1", "NR2F1", "OTX2", "PBX1", "SP9")
}

# Load data
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
  
  rawPLAC = paste0(plac.seq.dir, "IJ010.IJ057.IJ058.peaks/IJ010.IJ057.IJ058.10k.2.peaks.bedpe")
  rawPLAC = read.table(rawPLAC, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  placReduced = paste0(output.dir, "consolidatedEnsemblesSummary.txt")
  placReduced = read.table(placReduced, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  placFiltered = rawPLAC %>% 
    dplyr::mutate(interactionID = row_number()) %>% 
    as.data.frame()
  placFiltered = rbind(
    data.frame(seqnames = placFiltered$chr1, start = placFiltered$start1, end = placFiltered$end1, negLogFDR = -log10(placFiltered$fdr), negLog10P = placFiltered$ClusterNegLog10P, interactionID = placFiltered$interactionID, side = 1), 
    data.frame(seqnames = placFiltered$chr2, start = placFiltered$start2, end = placFiltered$end2, negLogFDR = -log10(placFiltered$fdr), negLog10P = placFiltered$ClusterNegLog10P, interactionID = placFiltered$interactionID, side = 2)
  )
  ovlp = as.data.frame(findOverlaps(GRanges(placFiltered), GRanges(placReduced)))
  placFiltered$bgEnsemble = placReduced$ensemble[ovlp$subjectHits]
  
  ensembleGenes = paste0(output.dir, "ensembleGenesLong.txt")
  ensembleGenes = read.table(ensembleGenes, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  bgEmissions = paste0(working.dir, "ChromHMM/Output/WT/DLX/emissions_", numHMM, ".edited.txt")
  bgEmissions = read.table(bgEmissions, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  bgHMM = paste0(working.dir, "ChromHMM/Output/WT/DLX/BG_WT_DLX_", numHMM, "_overlap.edited.txt")
  bgHMM = read.table(bgHMM, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  colnames(bgHMM) = c("state", 'Genome %', 'CpG Island', 'Exon', 'Gene', 'TES', 'TSS', 'TSS 2kb')
}


### Main figure
# 1A - summary diagram
{
  # composed from other sources, e.g. Biorender
}

# 1B - Positional summary of TF binding
{
  # Distribution of TFs distal/proximal
  {
    file.name = paste0(output.dir, "numTFs.vs.DistalProx.pdf")
    myColumns = colnames(bgPeaks)[grep("normHeight", colnames(bgPeaks))]
    tempObj = bgPeaks[, c("Feature", myColumns)]
    tempObj$numTFs = apply(tempObj[, myColumns], 1, function(X){sum(X > 0)})
    tempObj$Feature = ifelse(tempObj$Feature == "Promoter", "Proximal", "Distal")
    tempObj$numTFs = factor(tempObj$numTFs, levels = sort(unique(tempObj$numTFs), decreasing = TRUE))
    tempObj %>% 
      dplyr::select(Feature, numTFs) %>% 
      dplyr::group_by(numTFs, Feature) %>% 
      tally() %>% 
      dplyr::group_by(numTFs) %>% 
      dplyr::mutate(pct = 100 * n / sum(n), 
                    label = paste0("(", sum(n), ")")) %>% 
      ggplot(aes(x = numTFs, y = pct, fill = Feature)) +
      geom_bar(position = "stack", stat = "identity", color = "black", alpha = 0.75) +
      geom_text(aes(y = 99, label = label), hjust = 1, color = "grey30", size = 5) +
      theme_classic() +
      xlab("no. of TFs binding to locus") +
      ylab("percent peaks") +
      theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14), 
            legend.title = element_blank(),
            legend.text = element_text(size = 12), 
            legend.position = c(0.45, 0.5)) +
      coord_flip()
    ggsave(file.name, height = 4, width = 3)
    
    remove(list = c("tempObj", "myColumns"))
  }
  
  # Percentage of peaks in distance segments
  {
    file.name = paste0(output.dir, "distanceToTSS.bar.pdf")
    tempObj = bgPeaks[, c("Feature", "PeakID", "dist.nearest.promoter")]
    tempObj$dist.nearest.promoter = tempObj$dist.nearest.promoter / 1000
    tempObj$Proximal = as.integer(tempObj$Feature == "Promoter")
    tempObj$`> 50 kb` = as.integer(tempObj$dist.nearest.promoter > 50)
    tempObj$`10 - 50 kb` = as.integer(tempObj$dist.nearest.promoter > 10)
    tempObj$`10 - 50 kb`[tempObj$`> 50 kb` == 1] = 0
    tempObj$`2 - 10 kb` = as.integer(rowSums(tempObj[, 4:6]) == 0)
    tempObj %>% 
      dplyr::select(-c(Feature, PeakID, dist.nearest.promoter)) %>% 
      gather(key = "cat", value = "ct") %>% 
      dplyr::filter(ct > 0) %>% 
      dplyr::group_by(cat) %>% 
      tally() %>% 
      dplyr::mutate(pct = round(100 * n / sum(n), 0)) %>% 
      dplyr::mutate(cat = factor(cat, levels = c('> 50 kb', '10 - 50 kb', '2 - 10 kb', "Proximal"))) %>% 
      ggplot(aes(x = 1, y = pct, fill = cat)) +
      geom_bar(stat = "identity", position = "stack", color = NA, width = 0.85, alpha = 0.95) +
      geom_text(aes(y = 100 - cumsum(pct) + 0.5 * pct, label = pct), vjust = 0.5, 
                color = "black", size = 6) +
      ylab("percent of peaks") +
      theme_classic() + 
      scale_fill_manual(values = rev(c("darkblue", "darkgreen", "green", "lightgreen"))) +
      coord_flip() + 
      theme(axis.title.y = element_blank(), 
            axis.title.x = element_text(size = 18), 
            axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            axis.line = element_blank(), 
            legend.position = "bottom", 
            legend.title = element_blank(), 
            legend.text = element_text(size = 16)) +
      guides(fill = guide_legend(ncol = 2, reverse = TRUE, byrow = TRUE))
    ggsave(file.name, width = 4, height = 1.7)
    remove(list = c("tempObj"))
  }
}

# 1C - Percentage pie loci in contacts
{
  file.name = paste0(output.dir, "lociInPSC.pdf")
  bgPeaks[, c("PeakID", "location")] %>% 
    dplyr::group_by(location) %>% 
    tally() %>% 
    dplyr::mutate(pct = round(100 * n / sum(n), 0)) %>% 
    dplyr::mutate(cumsum = 100 - cumsum(pct)) %>% 
    ggplot(aes(x = "", y = pct, fill = location)) +
    geom_bar(width = 1, stat = "identity", color = "grey50") +
    scale_fill_manual(values=c("#E69F00", "#999999", "#56B4E9")) +
    coord_polar("y", start = 0) +
    theme_void() + 
    theme(legend.text = element_text(size = 16), 
          legend.title = element_text(size = 16), 
          legend.position = "bottom") +
    geom_text(aes(y = pct/2 + cumsum, label = paste0(pct, "%")), size = 12) +
    guides(fill = guide_legend(title=" "))
  
  ggsave(file.name, height = 4, width = 7)
}

# 1D - UCSC browser track
{
  # manual editing from UCSC browser printout
  
  # contact map for plac-seq display
  {
    myTF = 'Sp9'
    myChrom = 'chr2'
    windowStart = 72008513
    windowEnd = 73432142
    
    file.name = paste0(output.dir, myTF, ".placSeq.Contacts.pdf")
    tempPlac = placFiltered %>% 
      dplyr::filter(seqnames == myChrom, 
                    start >= windowStart, 
                    end <= windowEnd) %>% 
      as.data.frame()
    
    myRange = GRanges(paste0(myChrom, ':', windowStart, '-', windowEnd))
    myStep = as.integer((min(c(tempPlac$start, tempPlac$end)) - windowStart) / 10000)
    start0 = min(c(tempPlac$start, tempPlac$end)) - 10000 * myStep
    myStep = as.integer((windowEnd - max(c(tempPlac$start, tempPlac$end))) / 10000)
    end0 = max(c(tempPlac$start, tempPlac$end)) + 10000 * myStep
    end(myRange) = end0
    start(myRange) = start0
    myRange = tile(myRange, width = 10000) %>% 
      as.data.frame()
    myRange$center = as.integer(rowMeans(myRange[, c("start", "end")]))
    
    ovlp = as.data.frame(findOverlaps(GRanges(myRange), GRanges(tempPlac)))
    ovlp$center = myRange$center[ovlp$queryHits]
    ovlp[, c("negLogFDR", "negLog10P", "interactionID", "side")] = tempPlac[ovlp$subjectHits, c("negLogFDR", "negLog10P", "interactionID", "side")]
    
    myStep = nrow(myRange)
    tempObj = data.frame(matrix(ncol = myStep, nrow = myStep))
    tempObj = tempObj %>% 
      mutate(row = row_number()) %>% 
      gather(key = 'column', value = 'value', -row) %>% 
      mutate(column = sub('X', '', column)) %>% 
      mutate(label = paste(row, column), 
             invLabel = paste(column, row)) %>% 
      as.data.frame()
    for (myInter in unique(ovlp$interactionID)) {
      tempElem1 = ovlp$queryHits[ovlp$interactionID == myInter & ovlp$side == 1]
      tempElem2 = ovlp$queryHits[ovlp$interactionID == myInter & ovlp$side == 2]
      tempTable = data.frame(t(combn(c(tempElem1, tempElem2), 2)))
      tempTable$value = unique(ovlp$negLogFDR[ovlp$interactionID == myInter & ovlp$side == 1])
      colnames(tempTable)[1:2] = c('row', 'column')
      tempTable$label = paste(tempTable$row, tempTable$column)
      myTemp = cbind(tempTable[, c("label", "value")], tempObj[match(tempTable$label, tempObj$label), c("label", "value")])
      myTemp$avg = rowMeans(myTemp[, c(2,4)], na.rm = TRUE)
      tempObj$value[tempObj$label %in% myTemp$label] = myTemp$avg
      tempObj$value[tempObj$invLabel %in% myTemp$label] = myTemp$avg
    }
    tempObj = tempObj %>% 
      dplyr::select(-c(label, invLabel)) %>% 
      tidyr::spread(key = column, value = value, fill = 0) %>% 
      as.data.frame()
    tempObj$row = NULL
    tempObj = tempObj[, gtools::mixedsort(colnames(tempObj))]
    
    pheatmap::pheatmap(tempObj
                       , cluster_rows = FALSE
                       , cluster_cols = FALSE
                       , show_rownames = FALSE
                       , show_colnames = FALSE
                       , color = colorRampPalette(c("lightyellow", "red", "darkred"), bias = 1.5)(100)
                       , border_color = "yellow"
                       , fontsize_row = 10
                       , fontsize_col = 10
                       , height = 5
                       , width = 5
                       , legend = FALSE
                       , filename = file.name)  
    # print ensemble size
    print(noquote(paste("n =", unique(placReduced$ensembleSize[placReduced$ensemble == ensembleGenes$ensemble[ensembleGenes$gene0 == myTF]]))))
  }
}

# 1E - Binding summary plot
{
  # daisy wheel of binding
  {
    myColumns = colnames(bgPeaks)[grep("peakHeight", colnames(bgPeaks))]
    myTitle = paste(nrow(bgPeaks), "loci")
    tempObj0 = bgPeaks[, c("Feature", myColumns)] %>% 
      dplyr::mutate(Feature = ifelse(Feature == "Promoter", "proximal", "distal")) %>% 
      dplyr::rename_with(~ sub("_.*", "", .x)) %>% 
      tidyr::gather(key = "TF", value = "height", -Feature) %>% 
      dplyr::filter(height > 0) %>% 
      dplyr::select(-height) %>% 
      dplyr::group_by(TF, Feature) %>% 
      tally() %>% 
      dplyr::group_by(TF) %>% 
      dplyr::mutate(pct = round(n / sum(n) * 100, 0)) %>% 
      dplyr::select(-n) %>% 
      as.data.frame()
    tempObj = data.frame()
    for (myLine in seq_along(tempObj0[, 1])) {
      nLines = tempObj0$pct[myLine]
      myTemp = data.frame(TF = tempObj0$TF[myLine], Feature = tempObj0$Feature[myLine], seq = seq(1, nLines))
      if (nrow(tempObj) == 0) {
        tempObj = myTemp
      } else {
        tempObj = rbind(tempObj, myTemp)
      }
      remove(list = c("nLines", "myTemp"))
    }
    
    file.name = paste0(output.dir, "daisyWheel.pdf")
    pdf(file.name)
    PieDonut(tempObj
             , aes(pies = TF, donuts = Feature)
             , selected = 1
             , labelposition = 0
             , explode = seq(1, length(myColumns))
             , showPieName = F)
    dev.off()
    
    remove(list = c("tempObj"))
  }
  
  # accessory motif plots
  {
    # setting parameters
    {
      deNovo = data.frame(TF = c("ARX", "ASCL1", "DLX1", "DLX2", "DLX5", "GSX2", "LHX6", "NKX2.1", "NR2F1", 
                                 "OTX2", "PBX1", "SP9"), 
                          motifNumber = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                          1, 1, 16), 
                          dataset = c("", "strict", "", "strict", "strict", "strict", "strict", "strict", "strict", 
                                      "strict", "", "strict"), 
                          direction = c("", "", "", "", "", "", "", "", "RV", 
                                        "", "", ""), 
                          stringsAsFactors = FALSE)
      # Load TF colors
      {
        tf.colors = data.frame(TF = c("ARX", "ASCL1", "DLX1", "DLX2", "DLX5", "GSX2", "LHX6", "NKX2.1", "NR2F1", 
                                      "OTX2", "PBX1", "RBBP4", "SP9"), 
                               tfcolor = c("153,76,0", "127,0,255", "51,153,255", "0,204,0", "204,102,0", 
                                           "102,102,0", "127,0,255", "102,0,0", "51,255,255", "255,51,51", 
                                           "0,102,102", "96,96,96", "102,0,102"), 
                               stringsAsFactors = FALSE)
        tf.colors$r = sapply(tf.colors$tfcolor, function(X) { unlist(strsplit(X, "\\,"))[1] })
        tf.colors$g = sapply(tf.colors$tfcolor, function(X) { unlist(strsplit(X, "\\,"))[2] })
        tf.colors$b = sapply(tf.colors$tfcolor, function(X) { unlist(strsplit(X, "\\,"))[3] })
        tf.colors$colorcode = rgb(tf.colors$r, tf.colors$g, tf.colors$b, max = 255)
      }
    }
    
    # TF order
    {
      TF.order = data.frame(TF = bgTFs,
                            type = c("homeo", "bHLH", "homeo", "homeo", "homeo", "homeo", "homeo", "HTH", 
                                     "NRC4ZF", "homeo2", "HTH", "C2H2ZF"), 
                            order = c(1, 3, 1, 1, 1, 1, 1, 4, 5, 2, 4, 6), 
                            stringsAsFactors = FALSE)
    }
    
    # generate motif logos
    {
      # SP9
      {
        sp9 = paste0(gene.grad.dir, "motifLogos/SP9_MA1564.1.jaspar")
        sp9 = read.table(sp9, sep = '\t', header = FALSE, skip = 1, stringsAsFactors = FALSE)
        tempSp9 = data.frame()
        for (myLine in seq_along(sp9[,1])) {
          tempLine = unlist(strsplit(sp9[myLine, 1], " "))
          tempLine = tempLine[tempLine != ""]
          tempLine = tempLine[! tempLine %in% c("[", "]")]
          if (nrow(tempSp9) == 0) {
            tempSp9 = data.frame(as.numeric(tempLine[-1]), stringsAsFactors = FALSE)
            colnames(tempSp9) = tempLine[1]
          } else {
            tempSp9 = cbind(tempSp9, as.numeric(tempLine[-1]))
            colnames(tempSp9)[ncol(tempSp9)] = tempLine[1]
          }
        }
        myLength = nrow(tempSp9)
        sp9 = t(tempSp9)
        sp9 = ggseqlogo(sp9, method = 'prob') +
          theme(axis.title = element_blank(), 
                axis.text.y = element_blank(), 
                axis.text.x = element_blank())
        file.name = paste0(gene.grad.dir, "motifLogos/", "SP9.pdf")
        ggsave(file.name, height = 0.5, width = 1.2/10*myLength)
        file.name = paste0(gene.grad.dir, "motifLogos/", "SP9.png")
        ggsave(file.name, height = 0.5, width = 1.2/10*myLength)
        remove(list = c("sp9", "tempSp9"))
      }
      
      # all other
      for (myTF in bgTFs) {
        if (!myTF %in% c(" ")) {
          myFile = list.dirs(path = homer.dir, full.names = FALSE, recursive = FALSE)
          myFile = myFile[grep(myTF, myFile)]
          myFile = myFile[grep("complex", myFile, invert = TRUE)]
          myFile = ifelse(deNovo$dataset[deNovo$TF == myTF] == "strict", myFile[grep("strict", myFile)], myFile[grep("strict", myFile, invert = TRUE)])
          file.name = paste0(homer.dir, myFile, "/homerResults/motif", 
                             deNovo$motifNumber[deNovo$TF == myTF], 
                             deNovo$direction[deNovo$TF == myTF], 
                             ".motif")
          myMotif = read.table(file.name, sep = '\t', header = FALSE, skip = 1, stringsAsFactors = FALSE, 
                               col.names = c("A", "C", "G", "T"))
          myLength = nrow(myMotif)
          myMotif = t(myMotif)
          myMotif = ggseqlogo(myMotif, method = 'prob') +
            theme(axis.title = element_blank(), 
                  axis.text.y = element_blank(), 
                  axis.text.x = element_blank())
          file.name = paste0(gene.grad.dir, "motifLogos/", myTF, ".pdf")
          ggsave(file.name, height = 0.5, width = 1.2/10*myLength)
          file.name = paste0(gene.grad.dir, "motifLogos/", myTF, ".png")
          ggsave(file.name, height = 0.75, width = 1.2/10*myLength)
          remove(myMotif)
        }
      }
    }
    
    # facet_wrap of TFname/enrich/targets/motifCoverages, distal/prox, distance, widths
    {
      myColumns = colnames(bgPeaks)[grep("peakHeight", colnames(bgPeaks))]
      TFsizeLabels = bgPeaks[, myColumns]
      TFsizeLabels[TFsizeLabels > 0] = 1
      colnames(TFsizeLabels) = sub("_.*", "", colnames(TFsizeLabels))
      TFsizeLabels = TFsizeLabels %>% 
        gather(key = "TF", value = "hits") %>% 
        group_by(TF) %>% 
        summarise(numPeaks = sum(hits)) %>% 
        as.data.frame()
      
      myFeatures = data.frame()
      peakWidths = data.frame()
      TFlabel = data.frame()
      Mcoverage = data.frame()
      
      for (myTF in bgTFs) {
        myColumn = myColumns[grep(myTF, myColumns)]
        tempObj = bgPeaks[bgPeaks[, myColumn] > 0, colnames(bgPeaks)[1:34]]
        
        # b(top). De novo motif enrichment/targets (label)
        {
          myMotif = list.dirs(path = homer.dir, full.names = FALSE, recursive = FALSE)
          myMotif = myMotif[grep(myTF, myMotif)]
          myMotif = myMotif[grep("complex", myMotif, invert = TRUE)]
          myMotif = ifelse(deNovo$dataset[deNovo$TF == myTF] == "strict", myMotif[grep("strict", myMotif)], myMotif[grep("strict", myMotif, invert = TRUE)])
          myMotif = paste0(homer.dir, myMotif, "/", "homerResults.html")
          myMotif = readHTMLTable(myMotif)
          myMotif = list.clean(myMotif, fun = is.null, recursive = FALSE)
          myMotif = myMotif$`NULL`
          myMotif = myMotif[grep("\\*", myMotif[, 1], invert = TRUE), ]
          myMotif$`% of Targets` = as.numeric(sub("\\%", "", as.character(myMotif$`% of Targets`)))
          myMotif$`% of Background` = as.numeric(sub("\\%", "", as.character(myMotif$`% of Background`)))
          myMotif = myMotif[myMotif$`% of Targets` >= minTargetThrsh, ]
          myMotif$Enrichment = myMotif$`% of Targets` / myMotif$`% of Background`
          myMotif = myMotif[, grep("Motif|STD", colnames(myMotif), invert = TRUE)]
          myMotif$`log P-pvalue` = as.numeric(as.character(myMotif$`log P-pvalue`))
          myMotif$`P-value` = NULL
          myMotif$motif = as.character(myMotif$`Best Match/Details`)
          myMotif$`Best Match/Details` = NULL
          myMotif$motif = sub("\\(.*", "", sub("/.*", "", myMotif$motif))
          myMotif$motif = sub("^PB[[:digit:]]{4}.._", "", myMotif$motif)
          myMotif$motif = sub("^PH[[:digit:]]{4}.._", "", myMotif$motif)
          if (myTF != "NR2F1") {
            myMotif = myMotif[myMotif$Rank == deNovo$motifNumber[deNovo$TF == myTF], ]
          } else {
            mm1 = myMotif[myMotif$Rank == deNovo$motifNumber[deNovo$TF == myTF], ]
            mm2 = myMotif[myMotif$Rank == deNovo$motifNumber[deNovo$TF == myTF]+1, ]
            myMotif$`% of Targets` = mm1$`% of Targets` + mm2$`% of Targets`
            myMotif$`% of Background` = mm1$`% of Background` + mm2$`% of Background`
            myMotif$Enrichment = myMotif$`% of Targets` / myMotif$`% of Background`
            myMotif = myMotif[myMotif$Rank == deNovo$motifNumber[deNovo$TF == myTF], ]
          }
          myLabel = paste0(" (", round(myMotif$Enrichment, 2), ", ", round(myMotif$`% of Targets`, 0), "%)")
          
          if (nrow(TFlabel) == 0) {
            TFlabel = data.frame(TF = myTF, label = myLabel, stringsAsFactors = FALSE)
          } else {
            TFlabel = rbind(TFlabel, data.frame(TF = myTF, label = myLabel, stringsAsFactors = FALSE))
          }
          remove(list = c("myMotif", "myLabel"))
        }
        
        # c. motif coverage average plot
        {
          myMotif = list.dirs(path = homer.dir, full.names = FALSE, recursive = FALSE)
          myMotif = myMotif[grep(myTF, myMotif)]
          myMotif = myMotif[grep("complex", myMotif, invert = TRUE)]
          myMotif = ifelse(deNovo$dataset[deNovo$TF == myTF] == "strict", 
                           myMotif[grep("strict", myMotif)], 
                           myMotif[grep("strict", myMotif, invert = TRUE)])
          myMotif = paste0(homer.dir, myMotif, "/", myMotif, ".hist.txt")
          myMotif = read.table(myMotif, sep ='\t', header = TRUE, stringsAsFactors = FALSE)
          colnames(myMotif)[1] = "distance"
          myMotif = myMotif[, c("distance", colnames(myMotif)[grep("total.sites", colnames(myMotif))])]
          
          col1 = grep(paste0("^X", deNovo$motifNumber[deNovo$TF == myTF], "\\."), colnames(myMotif))[1]
          if (myTF != "NR2F1") {
            myMotif = myMotif[, c(1, col1)]
          } else {
            col2 = col1 + 1
            myMotif = myMotif[, c(1, col1, col2)]
            remove(list = c("col2"))
          }
          remove(list = c("col1"))
          
          myMotif$mean = apply(myMotif, 1, function(X){ ifelse(length(X) == 2, X[2], mean(X[-1])) })
          myMotif = myMotif[, c("distance", "mean")]
          myMotif$smoothed = sma(myMotif$mean, order = 3, h = 2)$fitted
          if (nrow(Mcoverage) == 0) {
            Mcoverage = data.frame(TF = myTF, myMotif, stringsAsFactors = FALSE)
          } else {
            Mcoverage = rbind(Mcoverage, data.frame(TF = myTF, myMotif, stringsAsFactors = FALSE))
          }
          remove(myMotif); 
        }
        
        # d. Genomic feature distribution
        {
          myFeature = tempObj %>% 
            dplyr::mutate(Feature = ifelse(Feature == "Promoter", "Proximal", "Distal")) %>% 
            dplyr::group_by(Feature) %>% 
            tally() %>% 
            dplyr::mutate(pct = 100 * n / sum(n)) %>% 
            dplyr::select(-n) %>% 
            dplyr::mutate(TF = myTF) %>% 
            as.data.frame()
          if (nrow(myFeatures) == 0) {
            myFeatures = myFeature
          } else {
            myFeatures = rbind(myFeatures, myFeature)
          }
          remove(myFeature)
        }
        
        # e. peak widths
        {
          myWidths = data.frame(TF = myTF, myWidths = tempObj$width / 1000, stringsAsFactors = FALSE)
          if (nrow(peakWidths) == 0) {
            peakWidths = myWidths
          } else {
            peakWidths = rbind(peakWidths, myWidths)
          }
          remove(myWidths)
        }
      }
      
      fontSize = 8
      newTFcolors = tf.colors$colorcode
      names(newTFcolors) = tf.colors$TF
      newTFcolors = c(newTFcolors, Distal = "yellow", Proximal = "blue")
      
      # Draw plot
      {
        # Zeroth column: TF caption, label and logo  <<<--- TODO
        {
          myZero = myFeatures[myFeatures$Feature == "Distal", ]
          myZero$label = TFlabel$label[match(myZero$TF, TFlabel$TF)]
          myZero$logo = paste0(gene.grad.dir, "motifLogos/", myZero$TF, ".png")
          
          zerothColumn = ggbarplot(myZero, x = "Feature", y = 1, fill = "white", color = "white") +
            theme(legend.position = "none", 
                  legend.title = element_blank(), 
                  axis.line = element_blank(),
                  axis.text = element_text(color = "white"), 
                  axis.title = element_text(color = "white"), 
                  axis.ticks = element_blank()) +
            ylim(c(0, 2)) +
            geom_image(x = 0.25, y = 0.65, aes(image = logo), size = 1.25, hjust = 0) +
            geom_text(data = myZero, x = 0.5, y = 1.75, aes(label = label), hjust = 0, size = fontSize*0.5) +
            ylab(" ") +
            xlab(" \n ")
          
          zerothColumn = facet(zerothColumn, facet.by = "TF", ncol = 1, scales = "fixed", strip.position = "left") +
            theme(panel.border = element_blank(), 
                  strip.text = element_text(size = 12), 
                  legend.title = element_blank())
        }
        
        # First column: motif coverage
        {
          Mcoverage$TF = factor(Mcoverage$TF, levels = TF.order$TF)
          Mcoverage$smoothed = as.numeric(Mcoverage$smoothed)
          
          firstColumn = ggline(Mcoverage, x = "distance", y = "smoothed", color = "TF", palette = newTFcolors, 
                               size = 0.75, numeric.x.axis = TRUE, plot_type = "l") +
            xlab("distance (bp)") +
            rremove("y.title") + 
            rremove("y.ticks") + 
            rremove("x.axis") +
            rremove("y.text") +
            scale_x_continuous(breaks = c(-200,0,200)) +
            theme(legend.position = "none", 
                  legend.title = element_blank(), 
                  axis.text.x = element_text(size = 4+fontSize), 
                  axis.title.x = element_text(size = 4+fontSize)) +
            geom_hline(yintercept = 0)
          
          firstColumn = facet(firstColumn, facet.by = "TF", ncol = 1, scales = "free_y") +
            theme(panel.border = element_blank(), 
                  strip.background = element_blank(),
                  strip.text = element_blank(), 
                  legend.title = element_blank())
        }
        
        # Second column: distal/proximal
        {
          myFeatures$size = TFsizeLabels$numPeaks[match(myFeatures$TF, TFsizeLabels$TF)]
          myFeatures$x = 1
          
          secondColumn = ggbarplot(myFeatures, x = "x", y = "pct", 
                                   fill = "Feature", 
                                   palette = newTFcolors, 
                                   ylab = "percent of peaks",
                                   xlab = "  ",
                                   rotate = TRUE, 
                                   position = position_stack()) +
            scale_y_continuous(limits = c(0,150), breaks = c(0, 50, 100)) +
            scale_x_discrete(breaks = c(0, 1, 2)) +
            theme(legend.position = "none", 
                  legend.title = element_blank(), 
                  axis.line.x = element_blank(),
                  axis.title.x = element_text(size = fontSize+5, angle = 0),
                  axis.title.y = element_text(colour = "white", size = 3),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.line.y = element_blank(),
                  axis.text.x = element_text(size = fontSize+4, angle = 0))
          secondColumn = secondColumn +
            geom_text(data = myFeatures, size = fontSize*0.5, aes(label = paste0("n = ", size)), x = 1.1, y = 98, hjust = 1)
          # geom_text(data = myFeatures, size = fontSize*0.7, aes(label = paste0(round(pct, 0), "%"), x = 1, y = ifelse(pct > 50, 85, 12 + pct)))
          
          secondColumn = facet(secondColumn, facet.by = "TF", ncol = 1, scales = "fixed") +
            theme(panel.border = element_blank(), 
                  strip.text = element_blank(), 
                  legend.title = element_blank())
          
        }
        
        # Third column: peak widths
        {
          peakWidths$TF = factor(peakWidths$TF, levels = TF.order$TF)
          fourthColumn = ggdensity(peakWidths, x = "myWidths", y = "..density..", fill = "TF", palette = newTFcolors, 
                                   size = 0.5, numeric.x.axis = TRUE, adjust = 1) +
            geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
            coord_cartesian(xlim = c(0, 3)) +
            xlab("loci width (kb)") +
            rremove("y.title") + 
            rremove("y.ticks") + 
            rremove("x.axis") +
            rremove("y.text") +
            scale_x_continuous(breaks = seq(0, 3, 1)) +
            theme(legend.position = "none", 
                  axis.title.x = element_text(size = 5+fontSize), 
                  legend.title = element_blank()) +
            geom_hline(yintercept = 0)
          fourthColumn = facet(fourthColumn, facet.by = "TF", ncol = 1, scales = "fixed", 
                               panel.labs.background = ) +
            theme(panel.border = element_blank(), 
                  strip.text = element_blank(), 
                  legend.title = element_blank())
        }
        
        # Compose
        {
          myPlot = ggarrange(zerothColumn, firstColumn, secondColumn, fourthColumn, 
                             nrow = 1, widths = c(2, 0.75, 1.25, 0.85))
          file.name = paste0(output.dir, "indivTFSummary.pdf")
          ggsave(myPlot, filename = file.name, height = 10, width = 7)
        }
      }
    }
    
    remove(list = c("deNovo", "TF.order", "myColumns", "firstColumn", "secondColumn", "fourthColumn", "Mcoverage", 
                    "mm1", "mm2", "myFeatures", "myPlot", "myZero", "peakWidths", "tf.colors", "TF.order", "TFlabel", 
                    "TFsizeLabels", "zerothColumn"))
  }
}


### Supplementary
# S1a - Analysis diagram
{
  # composed in MS Powerpoint
}

# S1b - In situ TF expression
{
  # composed by Athena from https://developingmouse.brain-map.org/
}

# S1c - scRNA-seq expression trajectory in pseudotime
{
  # composed by Linda from diffusion data from https://doi.org/10.1073/pnas.2108760119
}

# S1d - Means of widths by numTFs per location
{
  for (myLoc in c("distal", "proximal")) {
    file.name = paste0(output.dir, "distributionWidthsMeans.", myLoc, ".pdf")
    myColumns = colnames(bgPeaks)[grep("normHeight", colnames(bgPeaks))]
    tempObj = bgPeaks[, c("width", "Feature", myColumns)]
    tempObj$numTFs = rowSums(tempObj[, -c(1,2)] > 0)
    
    tempObj = tempObj %>% 
      select(numTFs, Feature, width) %>% 
      mutate(Feature = ifelse(Feature == "Promoter", "proximal", "distal")) %>% 
      filter(Feature == myLoc) %>% 
      as.data.frame()
    
    distrObj = data.frame()
    for (myNum in unique(tempObj$numTFs)) {
      tempObj0 = tempObj[tempObj$numTFs == myNum, ]
      tempObj0 = tempObj0[!is.na(tempObj0[, 1]), ]
      tempObj0 = replicate(nPerm, mean(sample(tempObj0$width, size = samplingRatio * length(tempObj0$width), replace = FALSE)))
      tempObj0 = data.frame(numTFs = myNum, width = tempObj0)
      if (nrow(distrObj) == 0) {
        distrObj = tempObj0
      } else {
        distrObj = rbind(distrObj, tempObj0)
      }
      remove(tempObj0)
    }
    
    distrObj %>% ggplot(aes(x = as.factor(numTFs), y = width, fill = as.factor(numTFs))) +
      geom_boxplot(outlier.size = 0.1) +
      scale_y_log10(breaks = c(400, 600, 1000, 1500, 2000, 3000)) +
      xlab("number of TFs") +
      ylab("width (bp)") +
      ggtitle(" ", subtitle = myLoc) +
      theme_classic() +
      theme(legend.position = "none",
            title = element_text(size = 14), 
            axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14), 
            panel.grid.minor.y = element_line(color = "blue", linetype = "dotted"), 
            panel.grid.major.y = element_line(color = "blue", linetype = "dotted"))
    ggsave(file.name, height = 4, width = 2.75)
  }
  remove(list = c("tempObj", "myColumns", "distrObj"))
}

# S1e / S1f - Individual TFs hitting PSC/loops segmented by distal/proximal
{
  for (myLoc in c("PSC", "inside")) {
    file.name = paste0(output.dir, "TFsIn", myLoc, ".pdf")
    myColumns = colnames(bgPeaks)[grep("normHeight", colnames(bgPeaks))]
    tempObj = bgPeaks %>% 
      dplyr::select(myColumns, Feature, location) %>% 
      dplyr::rename_at(.vars = myColumns, .funs = function(X){sub("_.*", "", X)}) %>% 
      dplyr::mutate_at(.vars = bgTFs, .funs = function(X){ as.integer(X > 0)}) %>% 
      dplyr::mutate(Feature = ifelse(Feature == "Promoter", "proximal", "distal")) %>% 
      dplyr::group_by(Feature, location) %>% 
      dplyr::summarise_at(.vars = bgTFs, sum) %>% 
      gather(key = "TF", value = "n", -c(Feature, location)) %>% 
      dplyr::group_by(TF) %>% 
      dplyr::mutate(pct = 100 * n / sum(n)) %>% 
      ungroup() %>% 
      dplyr::filter(location == myLoc) %>% 
      dplyr::select(-c(n, location)) %>% 
      as.data.frame()
    myLevels = tempObj %>% filter(Feature == "distal") %>% arrange(pct) %>% select(TF) %>% unlist()
    tempObj %>% 
      dplyr::mutate(TF = factor(TF, levels = myLevels)) %>% 
      ggplot(aes(x = TF, y = pct, fill = Feature)) +
      geom_bar(stat = "identity", color = "grey30", position = "dodge") +
      scale_fill_manual(values = c("yellow", "royalblue")) +
      ylab(paste0("percent peaks \nin ", ifelse(myLoc == "PSC", myLoc, "loops"))) +
      theme_classic() + 
      coord_flip() +
      theme(legend.position = "bottom", 
            legend.title = element_blank(), 
            legend.text = element_text(size = 12), 
            axis.text = element_text(size = 14), 
            axis.title = element_text(size = 14), 
            axis.title.y = element_blank()) +
      guides(fill = guide_legend(ncol = 1))
    ggsave(file.name, height = 5.5, width = 2)
    remove(tempObj)
  }
}

# S1g - Occurrence of ensemble sizes
{
  placReduced %>% 
    dplyr::group_by(ensemble) %>% 
    dplyr::filter(row_number() == 1) %>% 
    group_by(ensembleSize) %>% 
    tally() %>% 
    ggplot(aes(x = ensembleSize, y = n)) +
    geom_bar(stat = "identity", fill = "orange3") +
    xlab("interaction ensemble size") +
    ylab("number of occurrences") +
    theme_classic() +
    theme(text = element_text(size = 14))
  file.name = paste0(output.dir, "totalEnsembleSize Distribution.pdf")
  ggsave(file.name, width = 2.5, height = 2.5)
}

# S1h - Cumulative distribution of TF binding in plac-seq ensembles
{
  file.name = paste0(output.dir, "boundEnsembleSize Percent Distribution.pdf")
  placReduced %>% 
    dplyr::group_by(ensemble) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::mutate(inTFset = as.integer(ensemble %in% bgPeaks$bgEnsemble)) %>% 
    dplyr::mutate(inTFset = ifelse(inTFset == 1, "TFs bound", "other")) %>% 
    group_by(ensembleSize, inTFset) %>% 
    tally() %>% 
    dplyr::group_by(ensembleSize) %>% 
    dplyr::mutate(pct = n / sum(n)) %>% 
    dplyr::filter(inTFset == "TFs bound") %>% 
    ggplot(aes(x = ensembleSize, y = pct, color = as.factor(inTFset))) +
    geom_smooth(se = 0, method = "glm", method.args = list(family = "binomial"), 
                linetype = "dashed", color = 'red2') +
    geom_point(size = 1) +
    scale_color_manual(values = c("blue", "yellow3")) +
    coord_cartesian(xlim = c(0, 25)) +
    xlab("interaction ensemble \nsize") +
    ylab("occurrences with TFs bound") +
    scale_y_continuous(labels = scales::percent, breaks = c(0.4, 0.5, 0.60, 0.70, 0.80, 0.90, 1.00)) +
    theme_classic() +
    theme(axis.text = element_text(size = 13), 
          axis.title = element_text(size = 12), 
          legend.title = element_blank(), 
          legend.position = "none")
  ggsave(file.name, width = 2.2, height = 2.3)
}

# S1i - HMM genome-wide summary
{
  file.name = paste0(output.dir, "genomeHMM.pdf")
  columnPos = "bottom"
  
  # Left
  tempObj = bgHMM[-nrow(bgHMM), ]
  if (exists("ht")) remove(ht)
  for (myColumn in 2:ncol(tempObj)) {
    tempObj[, myColumn] = tempObj[, myColumn] / bgHMM[nrow(bgHMM), myColumn]
    myTemp = Heatmap(as.matrix(tempObj[, myColumn])
                     , col = circlize::colorRamp2(c(min(tempObj[, myColumn]), max(tempObj[, myColumn])), c("white", "royalblue"))
                     , cluster_rows = FALSE
                     , row_labels = tempObj$state
                     , column_labels = colnames(tempObj)[myColumn]
                     , cluster_columns = FALSE
                     , show_column_names = TRUE
                     , show_row_names = TRUE
                     , row_names_side = "left"
                     , column_names_side = columnPos
                     , column_names_rot = 45
                     , show_heatmap_legend = FALSE
    )
    if (exists("ht")) {
      ht = add_heatmap(ht, myTemp, direction = "horizontal")
    } else {
      ht = myTemp
    }
    remove(list = c("myTemp"))
  }
  leftSide = length(2:ncol(tempObj)) - 1
  
  # Right
  tempObj = bgEmissions
  if (exists("ht2")) remove(ht2)
  for (myColumn in 2:ncol(tempObj)) {
    myTemp = Heatmap(as.matrix(tempObj[, myColumn])
                     , col = circlize::colorRamp2(c(min(tempObj[, myColumn]), max(tempObj[, myColumn])), c("white", "royalblue"))
                     , cluster_rows = FALSE
                     , column_labels = colnames(tempObj)[myColumn]
                     , cluster_columns = FALSE
                     , show_column_names = TRUE
                     , show_row_names = FALSE
                     , row_names_side = "left"
                     , column_names_side = columnPos
                     , column_names_rot = 45
                     , show_heatmap_legend = FALSE
    )
    if (exists("ht2")) {
      ht2 = add_heatmap(ht2, myTemp, direction = "horizontal")
    } else {
      ht2 = myTemp
    }
    remove(list = c("myTemp"))
  }
  rightSide = length(2:ncol(tempObj)) - 1
  
  pdf(file.name, width = 4, height = 3)
  draw(ht + ht2
       , ht_gap = unit(c(rep(0, leftSide), 0.2, rep(0, rightSide)), "in")
       , padding = unit(c(0.1, 0.5, 0.1, 0.5), "in"))
  dev.off()
}

# S1j - Heatmap HMM vs PLAC-seq summary
{
  file.name = paste0(output.dir, "placContacts.HMM_Distribution.overall.Location.pdf")
  hmmDir = paste0(working.dir, "ChromHMM/Output/WT/")
  hmm = paste0(hmmDir, "DLX/BG_WT_DLX_", numHMM, "_dense.bed")
  hmm = read.table(hmm, sep = '\t', skip = 1, header = FALSE, stringsAsFactors = FALSE, 
                   col.names = c("seqnames", "start", "end", "bgState", "rank", "strand", "blockStart", 
                                 "blockEnd", "color"))
  hmm$hmmPeakID = paste0("hmm.", seq(1, nrow(hmm)))
  
  newHMM = disjoin(c(GRanges(bgPeaks), GRanges(hmm)), with.revmap = FALSE, ignore.strand = TRUE)
  newHMM = newHMM[width(newHMM) > 1]
  ovlp = as.data.frame(findOverlaps(newHMM, GRanges(hmm)))
  ovlp$seqnames = data.frame(newHMM, stringsAsFactors = FALSE)[ovlp$queryHits, "seqnames"]
  ovlp$start = data.frame(newHMM, stringsAsFactors = FALSE)[ovlp$queryHits, "start"]
  ovlp$end = data.frame(newHMM, stringsAsFactors = FALSE)[ovlp$queryHits, "end"]
  ovlp$hmmWidth = data.frame(newHMM, stringsAsFactors = FALSE)[ovlp$queryHits, "width"]
  ovlp[, colnames(hmm)[-c(1:3)]] = hmm[ovlp$subjectHits, colnames(hmm)[-c(1:3)]]
  newHMM = ovlp[, -c(1,2)]
  ovlp = as.data.frame(findOverlaps(GRanges(bgPeaks), GRanges(newHMM)))
  ovlp$PeakID = bgPeaks$PeakID[ovlp$queryHits]
  ovlp$width = bgPeaks$width[ovlp$queryHits]
  ovlp$bgCluster = bgPeaks$bgCluster[ovlp$queryHits]
  ovlp$bgEnsemble = bgPeaks$bgEnsemble[ovlp$queryHits]
  ovlp$bgLocation = bgPeaks$location[ovlp$queryHits]
  ovlp[, colnames(newHMM)[-c(1:3)]] = newHMM[ovlp$subjectHits, colnames(newHMM)[-c(1:3)]]
  ovlp = ovlp %>% 
    group_by(bgLocation, bgState) %>% 
    tally() %>% 
    group_by(bgLocation) %>% 
    mutate(pct = 100 * n / sum(n)) %>% 
    dplyr::select(-n) %>% 
    tidyr::spread(key = bgLocation, value = pct, fill = 0) %>% 
    as.data.frame()
  rownames(ovlp) = ovlp$bgState
  ovlp$bgState = NULL
  pheatmap::pheatmap(ovlp
                     , cluster_rows = TRUE
                     , cluster_cols = TRUE
                     , show_rownames = TRUE
                     , show_colnames = TRUE
                     , color = colorRampPalette(c("black", "yellow", "red"), bias = 1)(100)
                     , border_color = "grey10"
                     , angle_col = 45
                     , fontsize_row = 10
                     , fontsize_col = 10
                     , height = 2.5
                     , width = 2
                     , filename = file.name)  
}
