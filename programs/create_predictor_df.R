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
library(smri.process)
library(here)
library(readr)
library(RNifti)

rootdir = here::here()
proc_dir = file.path(rootdir,
    "cross_sectional", "raw")
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")


data_df = read_rds(df_file)
df = data_df

df = df %>% 
  filter(file.exists(outfile))

iid = 1

for (iid in seq(nrow(df))) {

  idf = df[iid,]
  ofile = idf$out_df

  if (!file.exists(ofile)) {
    pred_df = read_rds(idf$outfile)
    mask = readNifti(idf$mask_file)
    check_mask_fail(mask)
    mask_vec = c(mask) > 0

    if (!is.na(idf$gs_file)) {
      Y = readNifti(idf$gs_file)
      check_mask_fail(Y)
      Y = c(Y)
      # check all inside the mask
      mask_vec = mask_vec | Y > 0
      # stopifnot(!any(mask_vec == 0 & Y > 0))
      Y = Y[mask_vec]
    } else {
      Y = rep(NA, length = sum(mask_vec))
    }


    res = pbapply::pbsapply(pred_df$preds,
      function(x) {
      c(readNifti(x))
    })

    res = res[ mask_vec, ]
    res = as_data_frame(res)
    colnames(res) = pred_df$type
    res$Y = Y

    write_rds(res, path = ofile)
    rm(res)
    for (i in 1:10) gc()
  }
  print(iid)
}
