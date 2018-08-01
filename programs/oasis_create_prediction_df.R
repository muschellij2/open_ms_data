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
library(readr)
library(RNifti)

rootdir = here::here()
proc_dir = file.path(rootdir,
    "cross_sectional", "raw")
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")

eg = expand.grid(cv = c(FALSE, TRUE),
  trained = c(TRUE, FALSE),
  stringsAsFactors = FALSE)
eg = eg[ !(!eg$trained & eg$cv), ]
iscen = as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iscen)) {
  iscen = 3
}
print(iscen)
cv = eg$cv[iscen]
trained = eg$trained[iscen]

data_df = read_rds(df_file)
df = data_df

df = df %>% 
  filter(file.exists(outfile))

df$prob_file = file.path(df$proc_dir, 
  paste0(df$id, "_",
    "oasis",
    ifelse(cv, "_cv", ""),
    ifelse(!trained, "_untrained", ""),
    "_phat.nii.gz"))
df$sm_prob_file = file.path(df$proc_dir, 
  paste0(df$id, "_",
    "oasis",
    ifelse(cv, "_cv", ""),
    ifelse(!trained, "_untrained", ""),
    "_smoothed",    
    "_phat.nii.gz") )

df = df %>% 
  mutate(pred_df = file.path(proc_dir,
    paste0("oasis",
    ifelse(cv, "_cv", ""),
    ifelse(!trained, "_untrained", ""),
    "_predictions.rds")
    )
  )  

iid = 1

for (iid in seq(nrow(df))) {

  idf = df[iid,]
  ofile = idf$pred_df

  if (!file.exists(ofile)) {

    mask = readNifti(idf$mask_file)
    check_mask_fail(mask)
    if (!is.na(idf$gs_file)) {
      Y = readnii(idf$gs_file)
      # check all inside the mask
      mask = mask | Y > 0
      Y = c(Y)
    } else {
      Y = rep(NA, length = prod(dim(mask)))
    }
    mask_vec = c(mask) > 0
    Y = Y[mask_vec]

    prob_img = readNifti(idf$prob_file)
    prob_img = c(prob_img)
    prob_img = prob_img[ mask_vec]

    sm_prob_img = readNifti(idf$sm_prob_file)
    sm_prob_img = c(sm_prob_img)
    sm_prob_img = sm_prob_img[ mask_vec]

    struct_img = readNifti(idf$struct_file)
    struct_img = c(struct_img)
    struct_img = struct_img[ mask_vec]

    proc_df = data_frame(
      Y = Y,
      sm_p = sm_prob_img,
      p = prob_img,
      struct = struct_img)

    write_rds(proc_df, path = ofile)
    rm(proc_df)
    for (i in 1:10) gc()
  }
  print(iid)
}
