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
  print(iid)
  idf = df[iid,]

    if (!is.na(idf$gs_file)) {
      Y = readNifti(idf$gs_file)
      check_mask_fail(Y)
      sum_y = sum(Y)
      vol_y = sum_y * prod(pixdim(Y))/1000
    } 
    if (!is.na(idf$non_res_gs_file)) {
      nonres = readNifti(idf$non_res_gs_file)
      check_mask_fail(nonres)
      sum_check_y = sum(nonres)
      vol_check_y = sum_check_y * prod(pixdim(nonres))/1000
      rat = (vol_check_y/vol_y)
      if (
        (sum_y == 0 & sum_check_y > 0)
        ) {
        message(paste0("ID: ", idf$id, " is bad"))
      } else if (rat > 1.02 || rat < 0.98) {
        message(paste0("Ratio: ", round(rat, 3),
        " id: ", idf$id, " vol_y: ", vol_y))
      }
    }     
}
