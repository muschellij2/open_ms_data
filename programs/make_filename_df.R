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

set.seed(20180410)
rootdir = here::here()
proc_dir = file.path(rootdir,
                     "cross_sectional", "raw")

training_ids = readLines(here("programs",
    "training.txt"))
training_ids = as.numeric(training_ids)
training_ids = sprintf("patient%02.0f", training_ids)

load(here("cs_demog.rda"))

ids = list.dirs(proc_dir, 
                recursive = FALSE)
df = data_frame(
  id = basename(ids),
  id_dir = ids)
df = full_join(cs_demog, df)


df = df %>% mutate(
  raw_dir = file.path(id_dir, "raw"),
  proc_dir = file.path(id_dir, 
                       "prenorm"),
  reg_dir = file.path(id_dir, 
                      "registered"),
  brain_malf_dir = file.path(id_dir, 
                             "brain_malf"),
  malf_dir = file.path(id_dir, 
                       "malf")  
)
mods = c("FLAIR" = "FLAIR",
         "T1" = "T1W",
         "T2" = "T2W",
         "T1POST" = "T1WKS",
         "GOLD_STANDARD" = 
           "consensus_gt")
mod_df = sapply(mods, function(mod) {
  fname = file.path(df$id_dir, 
                    paste0(mod, ".nii.gz"))
})
fe = array(
  file.exists(mod_df),
  dim = dim(mod_df))
mod_df[ !fe ] = NA
mod_df = as_data_frame(mod_df)
mod_df = mod_df %>% 
  mutate(have_data = rowAlls(fe))
df = cbind(df, mod_df)
df = as_data_frame(df)
df = df %>% 
  filter(have_data) %>% 
  select(-have_data)

df = df %>% 
  mutate(
    outfile = file.path(proc_dir,
                        "predictor_names.rds"),
    MNI_outfile = file.path(proc_dir,
                            "MNI_predictor_names.rds"),
    Eve_outfile = file.path(proc_dir,
                            "EVE_predictor_names.rds"),    
    out_df = file.path(proc_dir,
                       "predictors.rds"),
    MNI_out_df = file.path(proc_dir,
                       "MNI_predictors.rds"),
    Eve_out_df = file.path(proc_dir,
                       "Eve_predictors.rds")                           
    
  )  
df$n_voxels = pbapply::pbsapply(df$GOLD_STANDARD, 
                                fslsum)
df$voxres = pbapply::pbsapply(df$GOLD_STANDARD, 
                              voxres)
df$volume = df$voxres * df$n_voxels / 1000

df$volume_cut = cut(df$volume,
                    breaks = quantile(df$volume,
                                      probs = seq(0, 1, by = 0.1)),
                    include.lowest = TRUE
)

# df$volume_cut = df$volume > median(df$volume)
df$age_cut = df$age > median(df$age)

groups = c("train", "test")

# df = df %>% 
#   group_by(sex, volume_cut, age_cut) %>% 
#   mutate(group = sample(groups, size = n(),
#                         replace = TRUE)) %>% 
#   ungroup()

df = df %>% 
  mutate(group = ifelse(id %in% training_ids, 
    "train", "test"))

# ggplot(aes(x = volume, fill = group), data = df) + 
# geom_histogram()

df = df %>% 
  select(-age_cut, -volume_cut)

df = df %>% 
  mutate(mask_file = file.path(df$proc_dir, 
                               "Brain_Mask_resampled.nii.gz"))
df = df %>% 
  mutate(
    gs_file = file.path(df$proc_dir,
                        paste0("GOLD_STANDARD", 
                               "_N4_noneck_reduced_winsor", 
                               "_regtoFLAIR_resampled", 
                               ".nii.gz"))
  )

df = df %>% 
  mutate(
    non_res_gs_file = file.path(df$proc_dir,
                                paste0("GOLD_STANDARD", 
                                       "_N4_noneck_reduced_winsor", 
                                       "_regtoFLAIR", 
                                       ".nii.gz"))
  )

df = df %>% 
  mutate(
    struct_file = file.path(df$proc_dir,
                            "STRUCTURES_resampled.nii.gz")
  )  
df = df %>% 
  mutate(
    image_file = file.path(df$proc_dir,
                           paste0("FLAIR", 
                                  "_N4_noneck_reduced_winsor", 
                                  "_regtoFLAIR_brain_N4_resampled", 
                                  ".nii.gz"))
  )  

xdf = df
outfile = here("cross_sectional", "raw", 
               "filename_df.rds")
write_rds(df, path = outfile)


df = xdf
outfile = here("cross_sectional", "raw",
               "MNI_filename_df.rds")

ren = function(x) {
  sub("_resampled", "_regtoMNI", x)
}
df = df %>% 
  mutate(
    image_file = ren(image_file),
    struct_file = ren(struct_file),
    gs_file = ren(gs_file),
    mask_file = ren(mask_file),
    outfile = MNI_outfile,
    out_df = MNI_out_df
  )

write_rds(df, path = outfile)


########################
# Eve
########################
df = xdf
outfile = here("cross_sectional", "raw",
               "Eve_filename_df.rds")

ren = function(x) {
  sub("_resampled", "_regtoEve", x)
}
df = df %>% 
  mutate(
    image_file = ren(image_file),
    struct_file = ren(struct_file),
    gs_file = ren(gs_file),
    mask_file = ren(mask_file),
    outfile = Eve_outfile,
    out_df = Eve_out_df
  )

write_rds(df, path = outfile)
