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
library(readr)
library(RNifti)
library(caret)

reduce_train_object = function(x) {
  x$control$index= NULL
  x$control$indexOut = NULL
  x$trainingData = NULL
  x$finalModel$predictions = NULL
  x
}

set.seed(20180410)
rootdir = here::here()
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")
mod_dir = here("cross_sectional", "model")
dir.create(mod_dir, showWarnings = FALSE)
model_file = here("cross_sectional", "model", 
  "train_model.rds")


isub = as.numeric(
  Sys.getenv("SGE_TASK_ID")
  )
if (is.na(isub) || isub < 1) {
  isub = 1
}

first_model_file = here(
  "cross_sectional", "model", 
    "first_train_model.rds")

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

df = read_rds(df_file)

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

model = read_rds(model_file)

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
        first_cut = predict(first_model$model,
          type = "response",
          newdata = proc_df)
        first_cut = first_cut > 
          first_model$prob_cutoff
      }

      p = predict(model, newdata = proc_df, 
        type = "prob")
      p = p[, "lesion"]
      proc_df$p = p
      # putting Y back
      proc_df$Y = Y    
      proc_df = proc_df %>% 
        dplyr::select(p, Y)

      # need this because may include y > 0 but 
      # not in mask
      mask = readnii(idf$mask_file)
      if (!is.na(idf$gs_file)) {
        Y = readnii(idf$gs_file)
        # check all inside the mask
        mask = mask | Y > 0
      }

      prob_img = remake_img(proc_df$p, mask, mask)
      
      if (filtered) {
        filtered_mask = remake_img(first_cut,
          mask, mask)
        prob_img = mask_img(prob_img, filtered_mask)
      }

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

