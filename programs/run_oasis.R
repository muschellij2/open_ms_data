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
library(oasis)
library(here)
library(readr)

rootdir = here::here()
proc_dir = here("cross_sectional", "raw")

mods = c("FLAIR" = "FLAIR",
  "T1" = "T1W",
  "T2" = "T2W",
  "GOLD_STANDARD" = 
  "consensus_gt")
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")

df = read_rds(df_file)
df$oasis_outfile = file.path(df$proc_dir,
  "oasis_train_list.rds")

isub = as.numeric(
  Sys.getenv("SGE_TASK_ID")
  )
# idir <- as.numeric(
#     Sys.getenv("LSB_JOBINDEX"))
if (is.na(isub) || isub < 1) {
  isub = 4
}
print(isub)

idf = df[isub,]

if (!file.exists(idf$oasis_outfile)) {
  
  suffix = paste0("_N4_noneck_reduced_winsor_", 
    "regtoFLAIR_brain_N4_resampled")

  files = file.path(idf$proc_dir,
    paste0(names(mods), suffix, ".nii.gz"))

  brain_mask_file = file.path(idf$proc_dir,
    "Brain_Mask_resampled.nii.gz")
  names(files) = names(mods)
  files["GOLD_STANDARD"] = sub("brain_N4_", "",
      files["GOLD_STANDARD"])  

  result = oasis_train_dataframe(
    flair = files["FLAIR"],
    t1 = files["T1"],
    t2 = files["T2"],
    pd = NULL,
    gold_standard = files["GOLD_STANDARD"],
    brain_mask = brain_mask_file,
    preproc = FALSE, 
    normalize = TRUE)
  result$preproc = NULL

  write_rds(result, path = idf$oasis_outfile)

}
