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

rootdir = here::here()
proc_dir = here("cross_sectional", "raw")

mods = c("FLAIR" = "FLAIR",
         "T1" = "T1W",
         "T2" = "T2W",
         "T1POST" = "T1WKS",
         "GOLD_STANDARD" = 
           "consensus_gt")
df_file = here("cross_sectional", "raw", 
               "filename_df.rds")

df = read_rds(df_file)

isub = as.numeric(
  Sys.getenv("SGE_TASK_ID")
)
# idir <- as.numeric(
#     Sys.getenv("LSB_JOBINDEX"))
if (is.na(isub) || isub < 1) {
  isub = 1
}

idf = df[isub,]

#template = "none",
itemplate = "none"
for (itemplate in c("none", "MNI", "Eve")) {
  idf = df[isub,]
  pre = switch(itemplate,
         none = "",
         MNI = "MNI_",
         Eve = "Eve_")
  outfile = idf[, paste0(pre, "outfile")]
  outfile = unlist(outfile)
  stopifnot(!is.na(outfile))
  file.exists(outfile)
  if (!file.exists(outfile)) {
    
    dir.create(idf$proc_dir,
               showWarnings = FALSE)
    dir.create(idf$reg_dir,
               showWarnings = FALSE)
    dir.create(idf$brain_malf_dir,
               showWarnings = FALSE)
    dir.create(idf$malf_dir,
               showWarnings = FALSE)
    fmods = setdiff(names(mods),
                    "GOLD_STANDARD")
    files = as.character(idf[, fmods])
    names(files) = fmods
    
    gold_standard = idf$GOLD_STANDARD
    
    print(isub)
    print(idf$id)
    
    outprefix = file.path(
      idf$brain_malf_dir,
      "FLAIR_")
    
    processed = smri_prenormalize(
      x = files,
      outdir = idf$proc_dir,
      gold_standard = gold_standard,
      gs_space = "FLAIR",
      reg_space = "FLAIR",
      malf_transform = "SyN",
      verbose = 2,
      outprefix = outprefix,  
      probs = c(0, 0.995),
      num_templates = 15,
      force_registration = FALSE)
    
    outprefix = file.path(
      idf$malf_dir,
      "T1")
    
    all_resampled = seg_normalize(
      prenormalize = processed,
      #template = "none",
      template = itemplate,
      verbose = TRUE,
      force_registration = FALSE,
      outprefix = outprefix
    )
    normalized = all_resampled$normalized
    
    pred = norm_predictors(
      normalized = normalized,
      normalization = "trimmed_z"
    )
    download_template_img_data()
    
    tdf = template_img_data()
    tdf = tdf %>% 
      mutate(
        modality = case_when(
          (modality == "T1_Pre") ~ "T1",
          (modality == "T1_Post") ~ "T1POST",
          TRUE ~ modality
        ),
        modality = factor(modality,
                          levels = names(files))
      ) %>% 
      filter(!is.na(modality)) %>% 
      arrange(modality)
    
    tdf = split(tdf, tdf$statistic)
    
    mean_imgs = tdf$Mean$file
    sd_imgs = tdf$SD$file
    
    add_suffix = paste0(normalized$suffix, "_", 
                        "trimmed_z", "_ztemp")
    
    ztemp = template_z_score(
      pred$intensity_normalized,
      mask = pred$normalized$brain_mask,
      template = "Eve",
      mean_imgs = mean_imgs,
      sd_imgs = sd_imgs,
      outdir = idf$proc_dir,
      verbose = TRUE,
      remask = TRUE,
      interpolator = "lanczosWindowedSinc",
      suffix = "_ztemp")
    
    names(ztemp) = paste0(names(ztemp), "_ztemp")
    pred$normalized$z_to_template = ztemp
    
    predictors = gather_predictors(pred)
    
    xdf = data_frame(preds = predictors,
                     id = idf$id,
                     type = names(predictors))
    
    write_rds(xdf, path = outfile)
    ############################
    # Running the same thing
    # but with Eve
    ############################
    # all_resampled = seg_normalize(
    #   prenormalize = processed,
    #   norm_outdir = idf$reg_dir,
    #   template = "Eve",
    #   verbose = TRUE
    # )
    # normalized = all_resampled$normalized
    
    # pred = norm_predictors(
    #   normalized = normalized,
    #   normalization = "trimmed_z"
    # )
    
    # predictors = gather_predictors(pred)
  }
}
