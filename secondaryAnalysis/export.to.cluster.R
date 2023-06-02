export.to.cluster = function(bedfiles.to.transfer) {
  bedfiles.to.transfer = gsub(" ", "\\\\ ", bedfiles.to.transfer)
  if (cluster.to.use == "crick") {
    system(paste0("rsync -auv ", bedfiles.to.transfer,
                  " cattapre@crick.cse.ucdavis.edu:/group/nordlab/users/rinaldo/COMB00/ChIPseq/data/bedfiles/"))
  } else if (cluster.to.use == "barbera") {
    system(paste0("rsync -auv ", bedfiles.to.transfer,
                  " rinaldo@barbera.genomecenter.ucdavis.edu:/share/nordlab/users/rinaldo/COMB00/ChIPseq/data/bedfiles/"))
  } else {
    print(noquote(" "))
    print(noquote(" "))
    print(noquote("    >>>>>>>> This is not a valid cluster  <<<<<<<<<"))
    print(noquote(" "))
    print(noquote(" "))
  }
}
