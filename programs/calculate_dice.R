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

rootdir = here::here()
prog_dir = file.path(rootdir, "programs")
source(file.path(prog_dir, "helper_functions.R"))

proc_dir = file.path(rootdir,
    "cross_sectional", "raw")
mod_dir = file.path(rootdir,
    "cross_sectional", "model")
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")

model = c("ranger", "oasis")
model = model[1]
iscen = 3

if (model == "ranger") {

  eg = expand.grid(
    remove_t1_post = c(FALSE, TRUE),
    filtered = c(FALSE, TRUE),
    stringsAsFactors = FALSE)
  eg$fname = paste0("predictions", 
      ifelse(eg$remove_t1_post, "_nopost", ""),
      ifelse(eg$filtered, "_filtered", ""),
      ".rds")
  eg$suffix = paste0("ranger", 
    ifelse(eg$remove_t1_post, "_nopost", ""),
    ifelse(eg$filtered, "_filtered", "")
    )

} else {
  eg = expand.grid(
    cv = c(FALSE, TRUE),
    trained = c(TRUE, FALSE),
    nopost = c(FALSE, TRUE),
    stringsAsFactors = FALSE)
  eg = eg[ !(!eg$trained & eg$cv), ]
  eg$suffix = paste0("oasis",
    ifelse(eg$cv, "_cv", ""),
    ifelse(!eg$trained, "_untrained", ""))
  eg$fname = paste0(eg$suffix,
    "_predictions.rds")
}

for (iscen in seq(nrow(eg))) {

  fname = eg$fname[iscen]
  suffix = eg$suffix[iscen]

  run_groups = c("train", "test")
  run_group = run_groups[1]

  for (run_group in run_groups) {
    
    print(run_group)
    df = read_rds(df_file)

    df = df %>% 
      filter(file.exists(outfile)) %>% 
      filter(group == run_group)

    df = df %>% 
      mutate(pred_df = 
        file.path(proc_dir, fname))

    iid = 1

    all_df = vector(mode = "list",
      length = nrow(df))
    names(all_df) = df$id
    iid = 1

    outfile = file.path(mod_dir, 
      paste0(run_group, "_", 
        suffix, "_",
        "roc_information.rds"))

    data_outfile = file.path(mod_dir, 
      paste0(run_group, "_", 
        suffix, "_",
        "full_prediction_data.rds"))

    for (iid in seq(nrow(df))) {

      print(iid)
      idf = df[iid,]
      ofile = idf$pred_df

      if (file.exists(ofile)) {
        proc_df = read_rds(ofile)
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

