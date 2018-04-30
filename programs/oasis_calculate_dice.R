#!/usr/bin/env RScript
#######################
# Processing Cross-sectional data
#######################
rm(list = ls())
# library(dcmtk)
# library(dcm2niir)
# library(divest)
library(dplyr)
library(here)
library(readr)
library(ROCR)
library(neurobase)

rootdir = here::here()
prog_dir = file.path(rootdir, "programs")
source(file.path(prog_dir, "helper_functions.R"))

proc_dir = file.path(rootdir,
    "cross_sectional", "raw")
mod_dir = file.path(rootdir,
    "cross_sectional", "model")
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")

mods = c("FLAIR" = "FLAIR",
  "T1" = "T1W",
  "T2" = "T2W",
  "GOLD_STANDARD" = 
  "consensus_gt")

eg = expand.grid(cv = c(FALSE, TRUE),
  trained = c(TRUE, FALSE),
  stringsAsFactors = FALSE)
eg = eg[ !(!eg$trained & eg$cv), ]
iscen = 1

for (iscen in 1:3) {
  print(iscen)
  cv = eg$cv[iscen]
  trained = eg$trained[iscen]

  run_groups = c("train", "test")
  run_group = run_groups[2]

  for (run_group in run_groups) {
    print(run_group)
    df = read_rds(df_file)

    df = df %>% 
      filter(file.exists(outfile)) %>% 
      filter(group == run_group)

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
    df$pred_df = file.path(df$proc_dir,
      "oasis_train_list.rds")

    iid = 1

    all_df = vector(mode = "list",
      length = nrow(df))
    names(all_df) = df$id
    iid = 1

    outfile = file.path(mod_dir, 
      paste0(run_group, "_", 
        "oasis_",
        ifelse(cv, "cv_", ""),
        ifelse(!trained, "untrained_", ""),      
        "roc_information.rds"))

    data_outfile = file.path(mod_dir, 
      paste0(run_group, "_", 
        "oasis_",      
        "full_prediction_data.rds"))

    for (iid in seq(nrow(df))) {

      print(iid)
      idf = df[iid,]
      ofile = idf$pred_df
      suffix = paste0("_N4_noneck_reduced_winsor_", 
        "regtoFLAIR_brain_N4_resampled")

      files = file.path(idf$proc_dir,
        paste0(names(mods), suffix, ".nii.gz"))

      brain_mask_file = file.path(idf$proc_dir,
        "Brain_Mask_resampled.nii.gz")
      names(files) = names(mods)
      files["GOLD_STANDARD"] = sub("brain_N4_", "",
          files["GOLD_STANDARD"])  

      if (file.exists(ofile)) {
        vsel = read_rds(ofile)$voxel_selection
        p = readnii(idf$prob_file)
        sm = readnii(idf$sm_prob_file)
        y = readnii(files["GOLD_STANDARD"])      
        mask = readnii(brain_mask_file)
        proc_df = data_frame(
          Y = mask_vals(y, mask),
          voxsel = mask_vals(vsel, mask),
          p = mask_vals(p, mask),
          sm_p = mask_vals(sm, mask),
          )
        all_df[[iid]] = proc_df
        rm(proc_df)
        for (i in 1:10) gc()
      }
    }

    full_df = bind_rows(all_df, .id = "id")

    write_rds(full_df,
      path = data_outfile)
    roc = run_roc(full_df$p, full_df$Y)
    tab = dice(
      full_df$p > roc$dice_cutoff,
      full_df$Y)
    sm_roc = run_roc(full_df$sm_p, full_df$Y)
    sm_tab = dice(
      full_df$sm_p > sm_roc$dice_cutoff,
      full_df$Y)

    L = list(
      non_smoothed = roc,
      smoothed = sm_roc
      )
    write_rds(L, path = outfile)
    rm(full_df)
  }
}
