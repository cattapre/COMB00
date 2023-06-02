# Dataset definition
source("COMB00.datasets.R")

allReps = FALSE

if (region == "BG") {
  dataset = dataset.BG
} else if (region == "CX") {
  dataset = dataset.CX
} else if (region == "combined") {
  dataset = c(dataset.BG, dataset.CX)
}

# 
histones = c("H3K27ac", "H3K27me3", "H3K4me3", "H3K4me1")
TFs = unique(gsub("_.*", "", dataset))

BGReps = data.frame(TF = c("ARX", "ASCL1", "DLX1", "DLX2", "DLX5", "GSX2", "LHX6", "NKX2.1", "NR2F1", 
                           "OTX2", "PBX1", "RBBP4", "SP9"), 
                    rep = c(3, 2, "N4", "N1", "O2", 1, 3, 2, 1, 
                            3, 1, "N1", "N1")) 
CXReps = data.frame(TF = c("EMX2", "LHX2", "NR2F1", "NR2F2", "PAX6", "PBX1"), 
                    rep = c(2, 1, 1, 1, 1, 1))

myClusters = data.frame(region = c("BG", "CX", "BG", "CX"), 
                        allReps = c(TRUE, TRUE, FALSE, FALSE), 
                        NoOfReps = c(11, 8, 20, 8))

if (exists("strippedVersion")) {
  if (strippedVersion) {
    myClusters = data.frame(region = c("BG", "CX", "CX", "BG", "BG", "BG"), 
                            allReps = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), 
                            NoOfReps = c(11, 8, 8, 12, 7, 8), 
                            location = c("all", "all", "all", "all", "Distal", "Proximal")
                            )
  }
}

use.histones = TRUE
diff.histones = FALSE
valid.enhancers = FALSE
DEX = FALSE
atac = FALSE
plac = FALSE
chromatin = FALSE

filter.peaks = TRUE

filter.peaks.by.peakHeight = FALSE   # either peakheight or PValue
pkHeight.thrsh = 0

filter.peaks.by.PValue = TRUE   # either FC or PValue
PV.thrsh = 0

DE.FDR.thrsh = 0.1
DE.PValue.thrsh = 0.05
too.broad.peak.thrsh = 4500
remove.too.broad.from.bed = as.integer(TRUE)

gene.file = "UCSC.genes.gtf"
enhancer.file = "VISTA.enhancers.20180710.mm10.txt"
enhancer.CX.gradient.file = "VISTA.Forebrain.New.txt"
placseq.CX.file.prefx = "plac.annotations_Iros_Jan2019_mm9_"
placseq.BG.file.prefx = ""

genome.used = "mm10"
promoter.limits = 2000    # bp from TSS
scopes.dataset = c("TF.narrow", "histone.broad", "histone.diff")

peak.buffer = 1000
peak.width.thrsh = 300


# merge distance criteria for merged peaks
merge.distance = 150     # <<<<<<<<  was 500
filter.proportion = 0.1 # set to proportion overlap with filtering feature for flag

use.export.to.cluster = TRUE
cluster.to.use = "barbera"    # options are crick or barbera



