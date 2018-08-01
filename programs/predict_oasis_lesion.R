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
library(here)
library(oasis)
library(readr)
library(RNifti)
library(caret)

set.seed(20180410)
rootdir = here::here()
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")
mod_dir = here("cross_sectional", "model")

df = read_rds(df_file)

eg = expand.grid(cv = c(FALSE, TRUE),
  trained = c(TRUE, FALSE),
  stringsAsFactors = FALSE)
eg = eg[ !(!eg$trained & eg$cv), ]
iscen = as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iscen)) {
  iscen = 2
}
print(iscen)
cv = eg$cv[iscen]
trained = eg$trained[iscen]

if (trained) {
  if (!cv) {
    model_file = here("cross_sectional", "model", 
      "oasis_model.rds")
  } else {
    model_file = here("cross_sectional", "model", 
      "cv_oasis_model.rds")
  }
}

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
df$oasis_outfile = file.path(df$proc_dir,
  "oasis_train_list.rds")

if (trained) {
  model = read_rds(model_file)
} else {
  model = nopd_oasis_model
}

n_ids = nrow(df)
all_df = vector(mode = "list",
  length = nrow(df))
names(all_df) = df$id

# read in the data
iid = 1

for (iid in seq(nrow(df))) {
  
  idf = df[iid,]
  ofile = idf$oasis_outfile
  prob_file = idf$prob_file
  sm_prob_file = idf$sm_prob_file

  if (file.exists(ofile)) {

    if (!file.exists(prob_file)) {

      proc_df = read_rds(ofile)
      voxsel_mask = proc_df$voxel_selection
      mask = proc_df$brain_mask
      proc_df = proc_df$oasis_dataframe

      stopifnot(nrow(proc_df) == sum(voxsel_mask))
      # need to remove for the cases where missing
      # values in Y but still want a prediction
      if (any(is.na(proc_df))) {
        print(paste0(iid, " is bad"))
        next; 
      }
      if (cv) {
        p = predict(model, 
          newdata = proc_df, 
          type = "prob")
        p = p[, "lesion"]
      } else {
        p = predict(model, 
          newdata = proc_df, 
          type = "response")
      }
      proc_df$p = p
      # putting Y back
      proc_df$Y = proc_df$GoldStandard    
      proc_df = proc_df %>% 
        dplyr::select(p, Y)

      prob_img = remake_img(proc_df$p, 
        voxsel_mask, voxsel_mask)
      writenii(prob_img, prob_file) 
      rm(proc_df)    
    }

    if (!file.exists(sm_prob_file)) {
      if (!exists("prob_img")) {
        prob_img = readnii(prob_file)
      }
      sm.pimg = ichseg::mean_image(prob_img, 
        nvoxels = 1)
      sm.pimg[
        abs(sm.pimg) < .Machine$double.eps^0.5
        ] = 0
      sm.pimg = niftiarr(prob_img, sm.pimg)
      writenii(sm.pimg, sm_prob_file)
      rm(sm.pimg)
      rm(prob_img)
    }
  }
  print(iid)
}

