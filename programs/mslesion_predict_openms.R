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
library(here)
library(readr)
library(ranger)
library(caret)


set.seed(20180410)
rootdir = here::here()
mod_dir = Sys.getenv("structural")
mod_dir = file.path(mod_dir, 
  "ms_lesion_challenge", "v4", "model")

templates = "MNI"
filtered = c(FALSE, TRUE)
remove_t1_post = c(FALSE, TRUE)
runner = c("ranger", "caret")
structures = c(FALSE, TRUE)

templates = "MNI"
itemplate = "MNI"
filtered = c(FALSE, TRUE)
# can only do nopd
nopd = TRUE
runner = c("ranger", "caret")
structures = c(FALSE, TRUE)
eg = expand.grid(
  filtered = filtered,
  itemplate = templates,
  nopd = nopd,
  structures = structures, 
  runner = runner,
  stringsAsFactors = FALSE)
eg$model_file = file.path(
  mod_dir, 
  paste0(
    ifelse(eg$runner == "ranger", "ranger_", ""),
    "train_model", 
    ifelse(eg$nopd, "_nopd", ""),
    ifelse(eg$filtered, "_filtered", ""),
    ifelse(eg$structures, "_withStruct", ""),
    ".rds"))

eg = eg[ file.exists(eg$model_file), ]

isub = as.numeric(
  Sys.getenv("SGE_TASK_ID")
  )
if (is.na(isub) || isub < 1) {
  isub = 1
}

ieg = eg[isub,]

nopd = ieg$nopd
filtered = ieg$filtered
model_file = ieg$model_file
runner = ieg$runner
itemplate = ieg$itemplate
structures = ieg$structures

first_model_file = file.path(
  mod_dir,
  paste0("first_train_model",
        ifelse(structures, "_withStruct", ""),
        ".rds"))

print(itemplate)
pre = switch(itemplate,
             none = "",
             MNI = "MNI_",
             Eve = "Eve_")

df_file = here("cross_sectional", "raw", 
               paste0(pre, "filename_df.rds"))

df = read_rds(df_file)

df$prob_file = file.path(df$proc_dir, 
  paste0(df$id, "_",
    runner, "_", itemplate,
    ifelse(filtered, "_filtered", ""),
    ifelse(structures, "_withStruct", ""),    
    "_ISBI",
    "_phat.nii.gz"))
df$sm_prob_file = sub("_phat", "_smoothed_phat", 
  df$prob_file)

model = read_rds(model_file)
keep_structures = attributes(model)$keep_structures

if (filtered) {
  first_model = read_rds(first_model_file)  
}

n_ids = nrow(df)
all_df = vector(mode = "list",
  length = nrow(df))
names(all_df) = df$id

# read in the data
iid = 30

for (iid in seq(nrow(df))) {
  
  idf = df[iid,]
  ofile = idf$out_df
  prob_file = idf$prob_file
  sm_prob_file = idf$sm_prob_file

  if (file.exists(ofile)) {

    if (!file.exists(prob_file)) {
      
      message("Reading in Data")
      proc_df = read_rds(ofile)
      Y = proc_df$Y
      # need to remove for the cases where missing
      # values in Y but still want a prediction
      proc_df = proc_df %>% 
        select(-Y)
      if (any(is.na(proc_df))) {
        print(paste0(iid, " is bad"))
        next; 
      }

      if (filtered) {
        message("Filtering")
        first_cut = predict(first_model$model,
          type = "response",
          newdata = proc_df)
        names(first_cut) = NULL        
        first_cut = first_cut > 
          first_model$prob_cutoff
      }

      keep = rep(TRUE, nrow(proc_df))
      if (structures) {
        pp = rep(0, nrow(proc_df))
        keep = proc_df$STRUCTURES %in% keep_structures
      }

      message("Predicting")
      if (runner == "caret") {
        p = predict(model, newdata = proc_df[keep, ], 
          type = "prob")
        p = p[, "lesion"] 
      } else {
        p = predict(model, data = proc_df[keep, ])
        if (model$treetype == "Classification") {
          p = p$predictions == "lesion"
        } else {
          p = p$predictions[, "lesion"]
        }
      }
      if (structures) {
        pp[keep] = p
        p = pp
      }

      proc_df$p = p
      # putting Y back
      proc_df$Y = Y    
      proc_df = proc_df %>% 
        dplyr::select(p, Y)

      # need this because may include y > 0 but 
      # not in mask
      mask = readnii(idf$mask_file)
      brain_mask = mask
      if (!is.na(idf$gs_file) 
        & file.exists(idf$gs_file)) {
        Y = readnii(idf$gs_file)
        # check all inside the mask
        mask = mask | Y > 0
        # need to do this because that's dim of proc_df$p
      }

      prob_img = remake_img(proc_df$p, mask, mask)
      
      if (filtered) {
        filtered_mask = remake_img(first_cut,
          mask, mask)
        prob_img = mask_img(prob_img, filtered_mask)
      }
      prob_img = mask_img(prob_img, brain_mask)
      message("Writing")
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

