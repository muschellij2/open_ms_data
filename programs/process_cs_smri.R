#!/usr/bin/env RScript
#######################
# Processing Cross-sectional data
#######################
rm(list = ls())
# library(dcmtk)
# library(dcm2niir)
# library(fslr)
# library(divest)
library(dplyr)
library(matrixStats)
library(smri.process)
library(here)

rootdir = here::here()
proc_dir = file.path(rootdir,
    "cross_sectional", "raw")

ids = list.dirs(proc_dir, 
	recursive = FALSE)
df = data_frame(id = basename(ids),
  id_dir = ids)
ids = unique(df$id)

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

isub = as.numeric(
  Sys.getenv("SGE_TASK_ID")
  )
# idir <- as.numeric(
#     Sys.getenv("LSB_JOBINDEX"))
if (is.na(isub) || isub < 1) {
  isub = 1
}

idf = df[isub,]
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
  num_templates = 35)


all_resampled = seg_normalize(
  prenormalize = processed,
  template = "none",
  verbose = TRUE
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

mean_imgs = tdf$Mean$file
sd_imgs = tdf$SD$file

z = template_z_score(
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

predictors = gather_predictors(pred)


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

