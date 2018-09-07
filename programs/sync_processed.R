rm(list=ls())
setwd("/Volumes/DATA_LOCAL/Projects/open_ms_data/cross_sectional/processed")
ids = list.dirs(recursive = FALSE, path = getwd(),
                full.names = FALSE)
id = ids[2]
structural = "/dcl01/smart/data/structural"
jdenig = "jmuschel@jhpce02.jhsph.edu"
# jdenig = "jmuschel@transfer01.jhpce.jhu.edu"
regspace = "FLAIR"
xsuffix = paste0("_N4_noneck_reduced_winsor_regto", 
                 regspace, "_brain_N4_regtoMNI")
mods = c("T1", "T2", "FLAIR", "T1POST")


for (id in ids) {
  local_dir = file.path(id)
  dir.create(local_dir, showWarnings = FALSE)
  cluster_dir = file.path(structural, 
                          "open_ms", "open_ms_data", "cross_sectional", 
                          "raw",
                          id, "prenorm")
  hdr = paste0(jdenig, ":", cluster_dir)
  
  print(id)
  isuffix = "_quantile"
  
  suffixes =  c("", "_quantile", "_trimmedz")
  check = function(x) all(file.exists(x))
  cmds = sapply(suffixes, function(isuffix){
    suffix = paste0(xsuffix, isuffix, ".nii.gz")
    paste0('"', hdr, '"/*', suffix)
  })
  cmds = c(cmds)
  cmds = paste(cmds, collapse = " ")
  
  infiles = sapply(suffixes, function(isuffix){
    suffix = paste0(xsuffix, isuffix, ".nii.gz")
    files = paste0(mods, suffix)
    file.path(local_dir, files)
  })
  
  if (!check(infiles)) {
    rsync_head = paste0('/usr/local/bin/rsync --progress -av ', cmds, " ", id)
    system(rsync_head)
  }
  gs_fname = paste0(
    'GOLD_STANDARD_N4_noneck_reduced_winsor_regto', regspace, 
    '_regtoMNI.nii.gz')
  
  if (!file.exists(file.path(id, gs_fname))) {
    cmd = paste0('/usr/local/bin/rsync --progress "', hdr, '"/', 
                 gs_fname, " ",
                 id)
    system(cmd)
  }
}
