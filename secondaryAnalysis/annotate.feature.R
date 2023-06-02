annotate.feature = function(base.peak.set, feature.types) {
  print(noquote(" "))
  txdb = makeTxDbFromGFF(paste0(genome.dir, gene.file), format = "gtf")
  promoter = promoters(txdb, upstream = promoter.limits, downstream = promoter.limits/10, columns = c("tx_name", "gene_id"))
  TSS = promoters(txdb, upstream = 0, downstream = 1, columns = c("tx_name", "gene_id"))
  transcript = transcripts(txdb, columns = c("tx_name", "gene_id"))
  genes = GenomicFeatures::genes(txdb, columns = c("tx_name", "gene_id"), single.strand.genes.only = FALSE)
  genebody = cds(txdb, columns = c("tx_name", "gene_id"))
  exons = exons(txdb, columns = c("tx_name", "gene_id"))
  
  threeUTR = unlist(threeUTRsByTranscript(txdb, use.names=TRUE))
  mcols(threeUTR)$tx_name = names(threeUTR)
  mcols(threeUTR)$gene_id = mcols(transcript)$gene_id[match(mcols(threeUTR)$tx_name, mcols(transcript)$tx_name)]
  
  fiveUTR = unlist(fiveUTRsByTranscript(txdb, use.names=TRUE))
  mcols(fiveUTR)$tx_name = names(fiveUTR)
  mcols(fiveUTR)$gene_id = mcols(transcript)$gene_id[match(mcols(fiveUTR)$tx_name, mcols(transcript)$tx_name)]
  
  introns = unlist(intronsByTranscript(txdb, use.names=TRUE))
  mcols(introns)$tx_name = names(introns)
  mcols(introns)$gene_id = mcols(transcript)$gene_id[match(mcols(introns)$tx_name, mcols(transcript)$tx_name)]
  
  for (feature.type in feature.types) {
    print(noquote(" "))
    if (feature.type == "genomic_features") {
      print(noquote(paste0(" ... annotating ", "Genomic Features")))
      for (feature in c("promoter", "genebody", "fiveUTR", "threeUTR", "exons", "introns")) {
          print(noquote(paste0("   ... ", feature)))
          
          temp.nearest = nearest(base.peak.set, subject = get(feature), ignore.strand = TRUE)
          mcols(base.peak.set)[, paste0("nearest.", feature)] = mcols(get(feature))$gene_id[temp.nearest]
          temp.distance = suppressWarnings(distanceToNearest(base.peak.set, subject = get(feature), ignore.strand = TRUE))
          mcols(base.peak.set)[, paste0("dist.nearest.", feature)] = as.data.frame(temp.distance)[,3]
      }
      
      mcols(base.peak.set)$Feature = "Intergenic"
      mcols(base.peak.set)$Feature[mcols(base.peak.set)$dist.nearest.genebody == 0] = "GeneBody"
      mcols(base.peak.set)$Feature[mcols(base.peak.set)$dist.nearest.exons == 0] = "Exon"
      mcols(base.peak.set)$Feature[mcols(base.peak.set)$dist.nearest.threeUTR == 0 &
                                     mcols(base.peak.set)$dist.nearest.exons == 0] = "threeUTR"
      mcols(base.peak.set)$Feature[mcols(base.peak.set)$dist.nearest.fiveUTR == 0 &
                                     mcols(base.peak.set)$dist.nearest.exons == 0] = "fiveUTR"
      mcols(base.peak.set)$Feature[unlist(mcols(base.peak.set)$nearest.promoter) == unlist(mcols(base.peak.set)$nearest.introns) & 
                                     mcols(base.peak.set)$dist.nearest.introns == 0] = "Intronic"
      mcols(base.peak.set)$Feature[mcols(base.peak.set)$dist.nearest.promoter == 0] = "Promoter"
      
      # Create a Gene column, based on nearest gene TSS
      temp.nearest = nearest(base.peak.set, subject = TSS, ignore.strand = TRUE)
      mcols(base.peak.set)$Gene = mcols(TSS)$gene_id[temp.nearest]
      mcols(base.peak.set)$PeakID = paste0("Merged.", seq(1, length(base.peak.set)))
    }
    
    if (feature.type == "conservancy") {
      print(noquote(paste0(" ... annotating ", "Evolutionary Conservancy")))
      for (feature in c("VertebratePhastcons", "PlacentalPhastcons", "EuarchPhastcons", "GERP")) {
        if (feature == "EuarchPhastcons") {
          print(noquote(paste0("   ... ", "Euarchontoglires (supraprimates) Conserved Elements")))
        } else {
          print(noquote(paste0("   ... ", feature, " Conserved Elements")))
        }
        
        if (feature == "GERP") {
          file.name = paste0(genome.dir, "mm10.GERP.txt")
          if (file.exists(file.name)) {
            temp.feature = GRanges(read.table(file.name, sep='\t', header=FALSE,
                                              comment.char="#", as.is=TRUE, stringsAsFactors=FALSE,
                                              col.names = c("chr", "start", "end", "width", "strand", "score", "V5", "PValue")))
          } else {
            library(rtracklayer)
            temp.feature = GRanges(read.table(paste0("/Volumes/GoogleDrive/Team Drives/Nord Lab - Personal Folders/Rinaldo/mm9/mm9.GERP.txt"), sep='\t', header=FALSE,
                                              comment.char="#", as.is=TRUE, stringsAsFactors=FALSE,
                                              col.names = c("chr", "start", "end", "score", "V5", "PValue")))
            chain = import.chain(paste0(genome.dir, "mm9ToMm10.over.chain"))
            seqlevelsStyle(temp.feature) = "UCSC"
            temp.feature = unlist(liftOver(temp.feature, chain))
            write.table(temp.feature, file.name, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
          }
        } else {
          if (feature == "VertebratePhastcons") {
            file.name = paste0(genome.dir, "phastCons60way.txt")
          } else if (feature == "PlacentalPhastcons") {
            file.name = paste0(genome.dir, "phastCons60wayPlacental.txt")
          } else {
            file.name = paste0(genome.dir, "phastCons60wayEuarch.txt")
          }
          temp.feature = GRanges(read.table(file.name, sep='\t', header=FALSE,
                                            comment.char="#", as.is=TRUE, stringsAsFactors=FALSE,
                                            col.names = c("bin", "chr", "start", "end", "lod", "score")))
        }
        
        ovlp = as.data.frame(findOverlaps(base.peak.set, temp.feature))
        ovlp = ovlp %>% 
          mutate(origValue = temp.feature$score[subjectHits]) %>% 
          group_by(queryHits) %>% 
          mutate(myValue = max(temp.feature$score[subjectHits])) %>% 
          filter(!duplicated(queryHits)) %>% 
          dplyr::select(queryHits, myValue) %>% 
          as.data.frame()
        ovlp$sequence = seq(1, nrow(ovlp))
        mcols(base.peak.set)[, feature] = NA
        mcols(base.peak.set)[ovlp$queryHits, feature] = ovlp$myValue[ovlp$sequence]
        # for (myInterval in unique(ovlp$queryHits)) {
        #   myValue = max(temp.feature$score[ovlp$subjectHits[ovlp$queryHits == myInterval]])
        #   mcols(base.peak.set)[ovlp$queryHits, feature] = myValue
        # }
        remove(temp.feature)
      }
    }
  }
  
  print(noquote(" "))
  print(noquote(" ... finished annotating features"))
  return(base.peak.set)
}