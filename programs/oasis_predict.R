#!/usr/bin/env RScript
#######################
# Processing Cross-sectional data
#######################
rm(list = ls())
# library(dcmtk)
# library(dcm2niir)
library(fslr)
# library(divest)
library(dplyr)
library(matrixStats)
library(here)
library(tidyr)
library(readr)
library(oasis)

rootdir = here::here()
proc_types = c("none", "trimmedz")
proc_type = proc_types[2]
models = list(
  # ALL = oasis::oasis_model,
  NO_T2_PD = oasis::not2_nopd_oasis_model,
  # NO_T2 = oasis::not2_oasis_model,
  NO_PD = oasis::nopd_oasis_model)

for (proc_type in proc_types) {

  pred_dir = here("cross_sectional", 
                  "predictions", "OASIS", 
                  "untrained", proc_type)
  dir.create(pred_dir, showWarnings = FALSE,
           recursive = TRUE)

  pred_dirs = paste0(pred_dir, "_",
                     names(models))
  sapply(pred_dirs, dir.create, showWarnings = FALSE)
  names(pred_dirs) = names(models)
  
  
  normalize = switch(proc_type,
                     trimmedz = FALSE,
                     none = TRUE)

  proc_dir = here("cross_sectional", 
                "atlases", proc_type)
  imgs = list.files(
    pattern = ".nii.gz",
    path = proc_dir, full.names = TRUE)
  df = data_frame(img = imgs)
  df = df %>% 
    mutate(
      stub = nii.stub(img, bn = TRUE)) %>% 
    separate(stub, into = c("pid", "mod"))
  df = df %>% 
    spread(mod, value = img) %>% 
    mutate(
      id = sub("atlas", "", pid),
      id = as.numeric(id)) %>% 
    arrange(id)
  
  isub = as.numeric(
    Sys.getenv("SGE_TASK_ID")
  )
  # idir <- as.numeric(
  #     Sys.getenv("LSB_JOBINDEX"))
  if (is.na(isub) || isub < 1) {
    isub = 1
  }
  
  idf = df[isub, c("FL", "T1", "T2", 
                   "T1POST", "mask", "pid")]
  pid = idf$pid
  idf$pid = NULL
  maskfile = idf$mask
  if (is.na(maskfile)) {
    maskfile = NULL
  }
  idf$mask = NULL
  
  fname = file.path(
    pred_dir,
    paste0(pid, "_oasis_list.rds"))
  
  imgs = as.list(idf)
  if (proc_type == "none") {
    imgs = lapply(imgs, readnii)
    imgs = lapply(imgs, function(x) {
      x[ x < 0 ] = 0
      x
    })
  }
  
  t1 = check_nifti(imgs$T1)
  brain_mask = t1 > 0
  
  
  if (!file.exists(fname)) {
    message("processing in data")
    oasis_list = oasis_train_dataframe(
      flair = imgs$FL,
      t1 = imgs$T1, 
      t2 = imgs$T2,
      pd = imgs$PD,
      gold_standard = maskfile,
      brain_mask = brain_mask, 
      preproc = FALSE,
      normalize = normalize,
      return_preproc = FALSE
    )
    saveRDS(oasis_list, 
            file = fname, compress = "xz")
  } else {
    message(paste0("reading in data: ", fname))
    oasis_list = readRDS(fname)
  }
  
  flair = check_nifti(imgs$FL)
  
  oasis_dataframe = oasis_list$oasis_dataframe
  voxel_selection = oasis_list$voxel_selection
  
  
  
  results = lapply(models, function(model) {
    x = oasis_predict(
      oasis_dataframe = oasis_dataframe,
      voxel_selection = voxel_selection,
      brain_mask = brain_mask,
      model = model,
      binary = TRUE
    )
    x
  })
  
  stopifnot(all(names(models) == names(results)))
  x = names(models)[1]
  
  outfiles = lapply(names(models), function(x) {
    pdir = pred_dirs[x]
    fnames = file.path(
      pdir,
      paste0(pid, c("_probability", "_mask"), 
             ".nii.gz"))
    
    img = results[[x]]$oasis_map
    write_nifti(img, fnames[1])
    
    img = results[[x]]$binary_map
    write_nifti(img, fnames[2])	
    return(fnames)
  })
  names(outfiles) = names(models)
  
}
