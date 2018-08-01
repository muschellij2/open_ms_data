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


data_df = read_rds(df_file)
df = data_df

df = df %>% 
  filter(file.exists(outfile))

isub = as.numeric(
  Sys.getenv("SGE_TASK_ID")
  )
if (is.na(isub) || isub < 1) {
  isub = 3
}
remove_t1_post = c(FALSE, TRUE)
filtered = c(FALSE, TRUE)
eg = expand.grid(
  remove_t1_post = remove_t1_post,
  filtered = filtered,
  stringsAsFactors = FALSE)
eg$model_file = here("cross_sectional", "model", 
    paste0("train_model", 
    ifelse(eg$remove_t1_post, "_nopost", ""),
    ifelse(eg$filtered, "_filtered", ""),
    ".rds"))

ieg = eg[isub,]

remove_t1_post = ieg$remove_t1_post
filtered = ieg$filtered
model_file = ieg$model_file

df$prob_file = file.path(df$proc_dir, 
  paste0(df$id, "_",
    "ranger",
    ifelse(remove_t1_post, "_nopost", ""),
    ifelse(filtered, "_filtered", ""),        
    "_phat.nii.gz"))
df$sm_prob_file = file.path(df$proc_dir, 
  paste0(df$id, "_",
    "ranger",
    ifelse(remove_t1_post, "_nopost", ""),
    ifelse(filtered, "_filtered", ""),        
    "_smoothed",    
    "_phat.nii.gz") )


df = df %>% 
  mutate(pred_df = file.path(proc_dir,
    paste0("predictions", 
      ifelse(remove_t1_post, "_nopost", ""),
      ifelse(filtered, "_filtered", ""),        
      ".rds"))
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

