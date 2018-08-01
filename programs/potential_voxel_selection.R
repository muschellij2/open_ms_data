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
library(tidyr)
library(ggplot2)

rootdir = here::here()
prog_dir = file.path(rootdir, "programs")
source(file.path(prog_dir, "helper_functions.R"))

proc_dir = file.path(rootdir,
    "cross_sectional", "raw")
mod_dir = file.path(rootdir,
    "cross_sectional", "model")
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")
df = read_rds(df_file)

all_eg = NULL
models = c("ranger")
groups = c("train")
for (model in models) {
  if (model == "ranger") {
    eg = expand.grid(
      remove_t1_post = c(FALSE, TRUE),
      group = groups,
      stringsAsFactors = FALSE)
    eg$suffix = paste0("ranger", 
      ifelse(eg$remove_t1_post, "_nopost", ""))

  } else {
    eg = expand.grid(
      cv = c(FALSE, TRUE),
      trained = c(TRUE, FALSE),
      nopost = c(FALSE, TRUE),
      group = groups,
      stringsAsFactors = FALSE)
    eg = eg[ !(!eg$trained & eg$cv), ]
    eg$suffix = paste0("oasis",
      ifelse(eg$cv, "_cv", ""),
      ifelse(!eg$trained, "_untrained", ""))
  }
  all_eg = bind_rows(all_eg, eg)
}
eg = all_eg
eg$full_group = paste0(eg$group, "_", eg$suffix)
eg = eg %>% 
  select(group, suffix, full_group) %>% 
  distinct()

eg$outfile = file.path(mod_dir, 
      paste0(eg$group, "_", 
        eg$suffix, "_",
        "roc_information.rds"))
eg$data_outfile = file.path(mod_dir, 
      paste0(eg$group, "_", 
        eg$suffix, "_",
        "full_prediction_data.rds"))



df = df %>% 
  filter(group %in% groups)

all_patient_df = vector(mode = "list",
  length = nrow(df))
names(all_patient_df) = df$id
iid =1 
for (iid in seq(nrow(df))) {
  print(iid)
  x = read_rds(df$out_df[iid])
  x = x %>% select(FLAIR_quantile, WM, WM.mn)
  all_patient_df[[iid]] = x
}

patient_df = bind_rows(all_patient_df)
rm(all_patient_df)

n_ids = nrow(eg)
all_roc = vector(mode = "list",
  length = nrow(eg))
all_df = all_roc
names(all_df) = eg$full_group

ieg = 1
# for (ieg in seq(nrow(eg))) {
  print(ieg)
  outfile = eg$data_outfile[ieg]
  x = read_rds(outfile)
  x$Y = x$Y > 0
  x = bind_cols(x, patient_df)
  all_df[[ieg]] = x
# }


