# Load packages
{
  initial.packages = c("GenomicFeatures", "tidyverse", "pheatmap", "ComplexHeatmap", "circlize", 
                       "smooth", "rtracklayer", "clusterProfiler", "org.Mm.eg.db", "biomaRt", "ggpubr", "ggsci")
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
  bgTFs = c("ARX", "ASCL1", "DLX1", "DLX2", "DLX5", "GSX2", "LHX6", "NKX2.1", "NR2F1", "OTX2", "PBX1", "SP9")
  source("COMB00.directories.R")
  hmmDir = paste0(working.dir, "ChromHMM/Output/WT/")
  dexDir = paste0(working.dir, "DEX/")
  numHMM = 9
  promoterLimits = 2000
  chainMM9_10 = paste0(genome.dir, "mm9ToMm10.over.chain")
  scaleLimit = 3
  minDiff = 2
  set.seed(0526)
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
  
  # mm10 genes
  {
    mm10Genes = "UCSC.genes.gtf"
    txdb = makeTxDbFromGFF(paste0(genome.dir, mm10Genes), format = "gtf")
    mm10Promoters = data.frame(promoters(txdb, upstream = promoterLimits, downstream = promoterLimits/10, columns = c("tx_name", "gene_id")), stringsAsFactors = FALSE)
    mm10Genes = data.frame(genes(txdb, columns = c("tx_name", "gene_id"), single.strand.genes.only = FALSE), stringsAsFactors = FALSE)
    
    mm10Genes = mm10Genes[grep("_", mm10Genes$seqnames, invert = TRUE), ]
    mm10Promoters = mm10Promoters[grep("_", mm10Promoters$seqnames, invert = TRUE), ]
    mm10Promoters$gene_id = as.character(mm10Promoters$gene_id)
    mm10Promoters = do.call(annotateEnsembles, list(mm10Promoters, "bg"))
    mm10Promoters$fpkm = geneExpression$meanRPKM[match(mm10Promoters$gene_id, geneExpression$Gene)]
    mm10Promoters = do.call(annotateEnsembles, list(mm10Promoters, "bg"))
    
    # Identifying gene aliases
    geneAlias = paste0(working.dir, "../mm10/MGI_MRK_List2.rpt.txt")
    geneAlias = read.table(geneAlias, sep = '\t', header = TRUE, stringsAsFactors = FALSE, quote = "", fill = TRUE)
    geneAlias = geneAlias[, c("Marker.Symbol", "Feature.Type", "Marker.Synonyms..pipe.separated.")]
    myTypes = c("lncRNA gene", "protein coding gene", "antisense lncRNA gene", "snoRNA gene", "lincRNA gene", "gene", 
                "bidirectional promoter lncRNA gene", "miRNA gene", "snRNA gene", "rRNA gene", "scRNA gene")
    geneAlias = geneAlias[geneAlias$Feature.Type %in% myTypes, ]
    colnames(geneAlias) = c("Symbol", "Type", "Alias")
    geneAlias$Alias[geneAlias$Alias == ""] = NA
    
    myPipes = geneAlias$Alias
    myNumber = 1
    for (myPipe in myPipes) {
      tempNum = length(unlist(strsplit(myPipe, "\\|")))
      if (tempNum > myNumber) myNumber = tempNum
      remove(tempNum)
    }
    geneAlias[, paste0("Alias_", seq(1, myNumber))] = NA
    for (myPipe in seq_along(myPipes)) {
      tempPipe = unlist(strsplit(myPipes[myPipe], "\\|"))
      tempNum = sapply(tempPipe, function(X){length(unlist(strsplit(X, "\\ ")))}) == 1
      tempPipe = tempPipe[tempNum]
      tempPipe = tempPipe[grep("kDa", tempPipe, invert = TRUE)]
      for (myNum in seq_along(tempPipe)) {
        geneAlias[myPipe, paste0("Alias_", myNum)] = tempPipe[myNum]
      }
    }
    
    geneAlias = geneAlias %>% 
      dplyr::select(-Alias) %>% 
      gather(key = "aliasNumber", value = "Alias", -Symbol, -Type) %>% 
      dplyr::select(-aliasNumber) %>% 
      dplyr::filter(!is.na(Alias)) %>% 
      as.data.frame()
    remove(myPipes)
    geneAlias$Symbol = gsub("\\ ", "", geneAlias$Symbol)
    geneAlias$Alias = gsub("\\ ", "", geneAlias$Alias)
    
    geneAlias[, c("seqnames", "start", "end", "strand")] = mm10Genes[match(geneAlias$Symbol, mm10Genes$group_name), c("seqnames", "start", "end", "strand")]
    geneAlias[is.na(geneAlias$seqnames), c("seqnames", "start", "end", "strand")] = mm10Genes[match(geneAlias$Alias[is.na(geneAlias$seqnames)], mm10Genes$group_name), c("seqnames", "start", "end", "strand")]
    
    file.name = paste0(working.dir, "../mm10/mart_export.txt")
    bmart = read.table(file.name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    bmart$Chromosome.scaffold.name = paste0("chr", bmart$Chromosome.scaffold.name)
    colnames(bmart) = c("Alias", "Symbol", "seqnames", "start", "end")
    bmart$seqnames[bmart$seqnames == 'chrMT'] = 'chrM'
  }
  
  # Load Biomart genes
  {
    library("biomaRt")
    ensembl = useMart("ensembl")
    ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)
    myAttributes = c("chromosome_name", "start_position", "end_position", "strand", "external_gene_name", 
                     "external_synonym", "hgnc_symbol", "mgi_symbol")
    biomartGene = getBM(attributes = myAttributes, mart = ensembl)
    colnames(biomartGene) = c("seqnames", "start", "end", "strand", "Symbol", "Alias", "HGNC", "MGI")
    biomartGene$seqnames = paste0("chr", biomartGene$seqnames)
    biomartGene$seqnames[biomartGene$seqnames == "chrMT"] = "chrM"
    biomartGene$strand[biomartGene$strand == 1] = '+'
    biomartGene$strand[biomartGene$strand == -1] = '-'
    biomartGene = biomartGene %>% 
      tidyr::gather(key = "geneClass", value = "gene", -c(seqnames, start, end, strand)) %>% 
      dplyr::select(-geneClass) %>% 
      distinct() %>% 
      dplyr::filter(gene != "") %>% 
      as.data.frame()
    biomartGene$fpkm = geneExpression$meanRPKM[match(biomartGene$gene, geneExpression$Gene)]
    biomartGene = biomartGene %>% 
      dplyr::filter(!is.na(fpkm)) %>% 
      as.data.frame()
    biomartGene = do.call(annotateEnsembles, list(biomartGene, "bg"))
  }
  
  # Load plac-seq data
  {
    placReduced = paste0(output.dir, "consolidatedEnsemblesSummary.txt")
    placReduced = read.table(placReduced, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  }
  
  # Load DLX/NKX REs
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


### Main Figure
# 2A - Schematics on Combinatorial TF binding
{
  # created by Rinaldo in Biorender
}

# 2B (proximal) / 2C (distal) - BG clusters
{
  # from external source - deeptools
}

# 2D - No. of TFs distribution by cluster
{
  file.name = paste0(output.dir, "numTFsInClusters.pdf")
  myColumns = colnames(bgPeaks)[grep("normHeight", colnames(bgPeaks))]
  tfDistrib = data.frame()
  for (myCluster in unique(bgPeaks$bgCluster[!is.na(bgPeaks$bgCluster)])) {
    tempObj = bgPeaks %>% 
      select_at(myColumns) %>% 
      dplyr::mutate(numTFs = rowSums(. > 0)) %>% 
      dplyr::rename_with(~ sub("_.*", "", .x)) %>% 
      dplyr::mutate(cluster = bgPeaks$bgCluster) %>% 
      dplyr::filter(cluster == myCluster) %>% 
      dplyr::select(numTFs) %>% 
      unlist() %>% unname()
    tempObj = replicate(nPerm, mean(sample(tempObj, size = samplingRatio * length(tempObj), replace = FALSE)))
    tempObj = data.frame(cluster = myCluster, numTFs = tempObj, stringsAsFactors = FALSE)
    if (nrow(tfDistrib) == 0) {
      tfDistrib = tempObj
    } else {
      tfDistrib = rbind(tfDistrib, tempObj)
    }
    remove(tempObj)
  }
  
  # random sample
  tempObj = randomSample %>% 
    select_at(myColumns) %>% 
    dplyr::mutate(numTFs = rowSums(. > 0)) %>% 
    dplyr::rename_with(~ sub("_.*", "", .x)) %>% 
    dplyr::select(numTFs) %>% 
    unlist() %>% unname()
  tempObj[is.na(tempObj)] = 0
  tempObj = replicate(nPerm, mean(sample(tempObj, size = samplingRatio * length(tempObj), replace = FALSE), na.rm = TRUE))
  tempObj = data.frame(cluster = "random", numTFs = tempObj, stringsAsFactors = FALSE)
  tfDistrib = rbind(tfDistrib, tempObj)
  remove(tempObj)
  
  myLevels = tfDistrib %>% group_by(cluster) %>% dplyr::mutate(mean = mean(numTFs)) %>% arrange(mean) %>% dplyr::select(cluster) %>% unique() %>% unlist()
  tfDistrib$cluster = factor(tfDistrib$cluster, levels = myLevels)
  
  ggplot(data = tfDistrib, aes(x = cluster, y = numTFs, color = cluster)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_classic() +
    scale_y_log10(breaks = seq(0, 10, 1), limits = c(0.9, 11)) +
    xlab("cluster") +
    ylab("mean number of TFs bound") +
    coord_flip() +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          legend.position = "none")
  ggsave(file.name, height = 4, width = 3)
  remove(tfDistrib)
}

# 2E - Distribution of width means
{
  # raw
  {
    # portrait
    {
      file.name = paste0(output.dir, "clusterDistribWidths.portrait.pdf")
      myColumns = colnames(bgPeaks)[grep("normHeight", colnames(bgPeaks))]
      lociWidths = data.frame()
      for (myCluster in unique(bgPeaks$bgCluster[!is.na(bgPeaks$bgCluster)])) {
        tempObj0 = bgPeaks[bgPeaks$bgCluster == myCluster, c("bgCluster", "width", myColumns)] %>% 
          dplyr::filter(!is.na(bgCluster))
        tempObj0$numTFs = rowSums(tempObj0[, myColumns] > 0)
        tempObj0$maxHeight = apply(tempObj0[, myColumns], 1, max)
        tempObj0 = tempObj0 %>% 
          dplyr::select(width) %>% 
          unlist()    
        tempObj0 = replicate(nPerm, mean(sample(tempObj0, size = samplingRatio * length(tempObj0), replace = FALSE)))
        tempObj0 = data.frame(cluster = myCluster, meanWidth = tempObj0, stringsAsFactors = FALSE)
        if (nrow(lociWidths) == 0) {
          lociWidths = tempObj0
        } else {
          lociWidths = rbind(lociWidths, tempObj0)
        }
        remove(tempObj0)
      }
      
      # random sample
      {
        tempObj = randomSample[rowSums(!is.na(randomSample[, -1])) > 0, c("width", myColumns)]
        tempObj = tempObj$width
        tempObj = replicate(nPerm, mean(sample(tempObj, size = samplingRatio * length(tempObj), replace = FALSE), na.rm = TRUE))
        tempObj = data.frame(cluster = "random", meanWidth = tempObj, stringsAsFactors = FALSE)
      }
      
      lociWidths = rbind(lociWidths, tempObj)
      myLevels = lociWidths %>% dplyr::group_by(cluster) %>% dplyr::mutate(mean = mean(meanWidth)) %>% dplyr::arrange(mean) %>% dplyr::select(cluster) %>% unique() %>% unlist()
      lociWidths$cluster = factor(lociWidths$cluster, levels = myLevels)
      ggplot(data = lociWidths, aes(x = cluster, y = meanWidth / 1000, fill = cluster)) +
        geom_boxplot(outlier.size = 0.1) +
        theme_classic() +
        ylab("normal. mean loci width (kb)") +
        xlab("cluster") +
        coord_flip() +
        theme(axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 14),
              legend.position = "none")
      ggsave(file.name, height = 4, width = 3)
      remove(list = c("lociWidths", "tempObj"))
    }
    
    # landscape
    {
      file.name = paste0(output.dir, "clusterDistribWidths.pdf")
      myColumns = colnames(bgPeaks)[grep("normHeight", colnames(bgPeaks))]
      lociWidths = data.frame()
      for (myCluster in unique(bgPeaks$bgCluster[!is.na(bgPeaks$bgCluster)])) {
        tempObj = bgPeaks[bgPeaks$bgCluster == myCluster, c("bgCluster", "width", myColumns)] %>% 
          dplyr::filter(!is.na(bgCluster))
        tempObj$numTFs = rowSums(tempObj[, myColumns] > 0)
        tempObj = tempObj %>% 
          dplyr::select(width) %>% 
          unlist()    
        tempObj = replicate(nPerm, mean(sample(tempObj, size = samplingRatio * length(tempObj), replace = FALSE)))
        tempObj = data.frame(cluster = myCluster, meanWidth = tempObj, stringsAsFactors = FALSE)
        if (nrow(lociWidths) == 0) {
          lociWidths = tempObj
        } else {
          lociWidths = rbind(lociWidths, tempObj)
        }
        remove(tempObj)
      }
      
      # random sample
      {
        tempObj = randomSample[rowSums(!is.na(randomSample[, -1])) > 0, c("width", myColumns)]
        tempObj = replicate(nPerm, mean(sample(tempObj, size = samplingRatio * length(tempObj), replace = FALSE), na.rm = TRUE))
        tempObj = data.frame(cluster = "random", meanWidth = tempObj, stringsAsFactors = FALSE)
      }
      
      lociWidths = rbind(lociWidths, tempObj)
      myLevels = lociWidths %>% dplyr::group_by(cluster) %>% dplyr::mutate(mean = mean(meanWidth)) %>% dplyr::arrange(mean) %>% dplyr::select(cluster) %>% unique() %>% unlist()
      lociWidths$cluster = factor(lociWidths$cluster, levels = myLevels)
      
      ggplot(data = lociWidths, aes(x = cluster, y = meanWidth / 1000, fill = cluster)) +
        geom_boxplot(outlier.size = 0.1) +
        theme_classic() +
        ylab("mean loci width (kb)") +
        xlab("cluster") +
        # coord_flip() +
        theme(axis.text.y = element_text(size = 12),
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 12),
              legend.position = "none")
      ggsave(file.name, height = 3, width = 5.5)
      remove(list = c("lociWidths", "tempObj"))
    }
  }
}

# 2F - Core motifs across clusters
{
  relativeMotifs = data.frame()
  for (myCluster in unique(bgPeaks$bgCluster[!is.na(bgPeaks$bgCluster)])) {
    myMotifs = paste0(motif.dir, "noRepsStripped/final/cluster_", myCluster, ".bed/targeted/", "knownResults.txt")
    myMotifs = read.table(myMotifs, sep = '\t', header = FALSE, skip = 1, stringsAsFactors = FALSE, 
                          col.names = c("motifName", "V2", "PValue", "negLogPValue", "V5", "V6", "pctTarget", "V7", "pctBkgrd"))
    myMotifs = myMotifs[, c("motifName", "PValue", "negLogPValue", "pctTarget", "pctBkgrd")]
    if (nrow(myMotifs) == 0) myMotifs = data.frame(
      motifName = TFs, PValue = 1, negLogPValue = 0, pctTarget = 0, pctBkgrd = 0
    )
    myMotifs$pctTarget = as.numeric(sub('%', '', myMotifs$pctTarget))
    myMotifs$pctBkgrd = as.numeric(sub('%', '', myMotifs$pctBkgrd))
    myMotifs$motifName = sub("\\(.*", "", myMotifs$motifName)
    myMotifs$relatEnrich = myMotifs$pctTarget - myMotifs$pctBkgrd
    if (nrow(relativeMotifs) == 0) {
      relativeMotifs = data.frame(cluster = sub('cluster_', '', myCluster), myMotifs, stringsAsFactors = FALSE)
    } else {
      relativeMotifs = rbind(relativeMotifs, data.frame(cluster = sub('cluster_', '', myCluster), myMotifs, stringsAsFactors = FALSE))
    }
    remove(myMotifs)
  }
  
  relativeMotifs = relativeMotifs %>%
    dplyr::filter(pctTarget > 10,
                  negLogPValue < -100) %>%
    as.data.frame()
  
  file.name = paste0(output.dir, "coreRelativeMotifs.portrait.pdf")
  relativeMotifs = relativeMotifs %>% 
    dplyr::select(cluster, motifName, relatEnrich) %>% 
    tidyr::spread(key = motifName, value = relatEnrich, fill = 0) %>% 
    as.data.frame()
  rownames(relativeMotifs) = relativeMotifs$cluster
  relativeMotifs$cluster = NULL
  pheatmap::pheatmap(relativeMotifs
                     , cluster_rows = TRUE
                     , cluster_cols = TRUE
                     , show_rownames = TRUE
                     , show_colnames = TRUE
                     , color = colorRampPalette(c("black", "yellow", "red"), bias = 0.2)(100)
                     , border_color = "grey10"
                     , angle_col = 45
                     , fontsize_row = 14
                     , fontsize_col = 14
                     , height = 6
                     , width = 5
                     , filename = file.name)  
  
  file.name = paste0(output.dir, "primaryRelativeMotifs.landscape.pdf")
  relativeMotifs = t(relativeMotifs) %>% 
    as.data.frame()
  pheatmap::pheatmap(relativeMotifs
                     , cluster_rows = TRUE
                     , cluster_cols = TRUE
                     , show_rownames = TRUE
                     , show_colnames = TRUE
                     #, color = colorRampPalette(c("black", "yellow", "red"), bias = 0.5)(100)
                     , color = colorRampPalette(c("black", "yellow", "red", "purple"))(50)
                     , border_color = "grey10"
                     , angle_col = 45
                     , fontsize_row = 14
                     , fontsize_col = 14
                     , height = 4.5
                     , width = 8.5
                     , filename = file.name)  
}


### Supplementary
# S2a - Example track
{
  # obtained from UCSC genome browser
}

# S2b - Primary motif average coverage across clusters
{
  file.name = paste0(output.dir, "primaryMotifCoverage.pdf")
  homeobox = c("ARX", "DLX1", "DLX2", "DLX5", "GSX2", "LHX6")
  # TF colors
  {
    tfColors = data.frame(TF = c("ARX", "ASCL1", "DLX1", "DLX2", "DLX5", "GSX2", 
                                 "LHX6", "NKX2.1", "NR2F1", "OTX2", "PBX1", "SP9", 
                                 "homeobox"), 
                          tfcolor = c("#fa958b", "#e4a550", "#c5b450", "#96bf52", "#5ac968", "#5bcda3", 
                                      "#59ccd0", "#4ec5f3", "#7eb2ff", "#d299ff", "#f889e9", "#ff88c0", 
                                      "#96bf52"), 
                          stringsAsFactors = FALSE)
    tfColors$colorcode = tfColors$tfcolor
    myColor = tfColors$colorcode
    names(myColor) = tfColors$TF
  }
  
  totalMotifs = data.frame()
  for (myTF in bgTFs) {
    for (myCluster in unique(bgPeaks$bgCluster[!is.na(bgPeaks$bgCluster)])) {
      myMotif = paste0(motif.dir, "noRepsStripped/final/cluster_", myCluster, ".bed/targeted/cluster_", myCluster, ".bed.hist.txt")
      myMotif = read.table(myMotif, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
      colnames(myMotif)[1] = "distance"
      myMotif = myMotif %>% 
        dplyr::select(distance, grep("total.sites", colnames(myMotif)) & grep(myTF, colnames(myMotif))) 
      colnames(myMotif)[2] = "coverage"
      myMotif = data.frame(TF = myTF, cluster = myCluster, myMotif, stringsAsFactors = FALSE)
      if (nrow(totalMotifs) == 0) {
        totalMotifs = myMotif
      } else {
        totalMotifs = rbind(totalMotifs, myMotif)
      }
    }
  }
  totalMotifs$class = totalMotifs$TF
  totalMotifs$class[totalMotifs$class %in% homeobox] = "homeobox"
  totalMotifs = totalMotifs %>% 
    group_by(class, distance, cluster) %>% 
    mutate(mean = mean(coverage)) %>% 
    as.data.frame()
  totalMotifs$smoothed = NA
  for (myClass in unique(totalMotifs$class)) {
    for (myCluster in unique(totalMotifs$cluster)) {
      tempObj = totalMotifs[totalMotifs$class == myClass & totalMotifs$cluster == myCluster, ]
      if (nrow(tempObj) > 0) {
        tempObj$smoothed = as.numeric(sma(tempObj$mean, order = 3, h = 2)$fitted)
        totalMotifs$smoothed[totalMotifs$cluster == myCluster & totalMotifs$class == myClass] = tempObj$smoothed[match(totalMotifs$distance, tempObj$distance)]
      }
      remove(tempObj)
    }
  }
  myLevels = c(
    gtools::mixedsort(unique(totalMotifs$cluster[grep("_D", totalMotifs$cluster)])),
    gtools::mixedsort(unique(totalMotifs$cluster[grep("_P", totalMotifs$cluster)]))
  )
  totalMotifs$cluster = factor(totalMotifs$cluster, levels = myLevels)
  
  maxScale = 0.015
  totalMotifs %>% 
    dplyr::select(class, cluster, distance, smoothed) %>% 
    dplyr::distinct() %>% 
    ggplot(aes(x = distance, y = smoothed, color = class)) +
    geom_line() + 
    coord_cartesian(ylim = c(0, maxScale)) +
    scale_y_continuous(breaks = c(0, round(maxScale/2, 3), maxScale)) +
    xlab("distance from peak center (bp)") +
    ylab("") +
    scale_color_manual(values = myColor) +
    theme_classic() +
    theme(legend.position = "none") +
    facet_grid(cluster ~ class, scales = "fixed") +
    theme(strip.background = element_rect(fill = "lightyellow"), 
          strip.text.y = element_text(angle = 0))
  ggsave(file.name, height = 15, width = 7)
}

