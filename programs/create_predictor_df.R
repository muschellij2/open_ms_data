#!/usr/bin/env RScript
#######################
# Processing Cross-sectional data
#######################
rm(list = ls())
# library(dcmtk)
# library(dcm2niir)
library(RNifti)
library(fslr)
# library(divest)
library(dplyr)
library(matrixStats)
library(smri.process)
library(here)
library(readr)

rootdir = here::here()
proc_dir = file.path(rootdir,
    "cross_sectional", "raw")

iid = as.numeric(
  Sys.getenv("SGE_TASK_ID")
)
# idir <- as.numeric(
#     Sys.getenv("LSB_JOBINDEX"))
if (is.na(iid) || iid < 1) {
  iid = 18
}


itemplate = "MNI"
# itemplate = "MNI"
  # for (itemplate in c("none", "MNI", "Eve")) {
    print(itemplate)
    pre = switch(itemplate,
                 none = "",
                 MNI = "MNI_",
                 Eve = "Eve_")

    df_file = here(
      "cross_sectional", "raw",
      paste0(pre, "filename_df.rds"))
    df = read_rds(df_file)

    idf = df[iid,]
    ofile = idf$out_df
    print(file.exists(ofile))

# df = df %>% 
#   filter(file.exists(outfile))

# for (iid in seq(nrow(df))) {

  if (!file.exists(ofile)) {

    pred_df = read_rds(idf$outfile)
    mask = readNifti(idf$mask_file)
    check_mask_fail(mask)
    mask_vec = c(mask) > 0

    if (!is.na(idf$gs_file) & 
      file.exists(idf$gs_file)) {
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

    stopifnot(!anyNA(mask_vec))
    for (icol in colnames(res)) {
      res[ is.nan(res[, icol, drop = TRUE]), icol] = 0
      res[ is.infinite(res[, icol, drop = TRUE]), icol] = 0
    }
    stopifnot(!anyNA(res))    
    
    res = res[ mask_vec, ]
    res = as_data_frame(res)
    colnames(res) = pred_df$type
    res$Y = Y

    write_rds(res, path = ofile)
    rm(res)
    for (i in 1:10) gc()
  }
  print(iid)
# }
