base.dir = paste0("/Users/cattapre/Library/CloudStorage/GoogleDrive-rpreta@ucdavis.edu/Shared drives/")
working.dir = paste0(base.dir, "NordLabRinaldo/COMB00/")

# Input data
atac.diff.dir = paste0(working.dir, "atac.diff/")
atac.seq.dir = paste0(working.dir, "atac.seq/")
broad.peak.dir = paste0(working.dir, "broad.peaks/")
diff.hist.dir = paste0(working.dir, "diff.peaks/")
enhancer.dir = paste0(base.dir, "NordLabRinaldo/VISTA_enhancers/")
gene.grad.dir = paste0(working.dir, "external.info/")
genome.dir = paste0(base.dir, "NordLabRinaldo/mm10/")
homer.dir = paste0(working.dir, "homer/")
motif.dir = paste0(working.dir, "motifs/")
narrow.peak.dir = paste0(working.dir, "narrow.peaks/")
plac.seq.dir = paste0(working.dir, "plac.seq/")

# Output data
output.dir = paste0(working.dir, "output/")

if (allReps) {
  cluster.dir = paste0(working.dir, "external.info/", region, ".allReps/")
  plot.dir = paste0(output.dir, "plots/", region, "/allReps/")
} else {
  cluster.dir = paste0(working.dir, "external.info/", region, ".noReps/")
  plot.dir = paste0(output.dir, "plots/", region, "/noReps/")
}

# Redefine base
base.dir = paste0(working.dir, "COMB00/")
