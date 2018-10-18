rm(list = ls())
library(here)
library(dplyr)
library(neurobase)
library(tidyr)
proc_dir = here("cross_sectional", "processed")
atlas_dir = here("cross_sectional", "atlases")

template = "MNI"
x = list.files(pattern = ".nii.gz", 
           recursive = TRUE,
           path = proc_dir)
df = data_frame(img = x)
df = df %>% 
  mutate(id = dirname(img))
gs = df %>% 
  filter(grepl("GOLD", img))
df = df %>% 
  filter(!grepl("GOLD", img))
df = df %>% 
  mutate(stub = nii.stub(img, bn = TRUE),
         norm = sub(".*MNI", "", stub))
df = df %>% 
  mutate(mod = sub("_.*", "", stub)) %>% 
  select(-stub) 
suffixes =  c( none = "", quantile = "_quantile", trimmedz = "_trimmedz")

outdir = file.path(atlas_dir, names(suffixes))
x = sapply(outdir, dir.create, showWarnings = FALSE)

suff_df = data_frame(norm = suffixes,
                     suff_dir =  file.path(atlas_dir, names(suffixes)))

pat_df = df %>% 
  select(id) %>% 
  distinct() %>% 
  arrange(id)
pat_df = pat_df %>% 
  mutate(atlas_id = paste0("atlas", seq(nrow(pat_df))))

stopifnot(all(df$norm %in% suffixes))
df = left_join(df, suff_df)
df = left_join(df, pat_df)

gs = left_join(gs, pat_df)
gs = gs %>% 
  mutate(out_fname = paste0(atlas_id, "_mask.nii.gz"),
         full_img = file.path(proc_dir, img))
# copy to all 
for (ioutdir in outdir) {
  file.copy(gs$full_img, file.path(ioutdir, gs$out_fname))
}

df = df %>% 
  mutate(mod = ifelse(mod == "FLAIR", "FL", mod))
df = df %>% 
  mutate(out_fname = paste0(atlas_id, "_", mod, ".nii.gz"),
         out_fname = file.path(suff_dir, out_fname),
         full_img = file.path(proc_dir, img))

file.copy(df$full_img, df$out_fname)
# df = df %>% 
#   spread(key = "norm", value = "img")

