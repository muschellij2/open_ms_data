rm(list = ls())
library(flexconn)
library(flexconnr)
library(dplyr)
library(keras)
library(tidyr)
library(neurobase)
library(here)
user = Sys.getenv("USER")
if (user == "johnmuschelli") {
  reticulate::use_python(paste0(
    "/Library/Frameworks/Python.framework/Versions/3.5/bin/python3"))
}

rootdir = here::here()
proc_dir = here("cross_sectional", "processed")
outdir = here(
  "cross_sectional",
  "predictions",
  "FLEXCONN",
  "untrained")
dir.create(outdir, showWarnings = FALSE)

mods = c("FLAIR" = "FLAIR",
         "T1" = "T1W",
         "T2" = "T2W",
         "T1POST" = "T1WKS",
         "GOLD_STANDARD" = 
           "consensus_gt")

files = list.files(path = proc_dir, 
                   pattern = ".nii.gz",
                   recursive = TRUE, 
                   full.names = TRUE)

df = data_frame(img = files,
                stub = nii.stub(img, bn = TRUE))
df = df %>% 
  mutate(id = basename(dirname(img)),
         modality = sub("_N4.*", "", stub)) %>% 
  select(-stub)
df = df %>% 
  spread(key = modality, value = img)
back = df %>% gather(key = "var", "value", -id)
any(is.na(back$value))



num_atlases = "61"
atlas_dir = file.path(
  outdir,
  paste0(num_atlases, "atlases"))
dir.create(atlas_dir, showWarnings = FALSE)
df = df %>% 
  mutate(
    out_prob = file.path(
      atlas_dir,
      paste0(id, "_", "LesionMembership.nii.gz")
    ),
    out_pred = sub("_LesionMembership", "_LesionMask", 
                   out_prob)
  )

iid = 1
for (iid in seq(nrow(df))) {
  idf = df[iid,]
  
  t1 = idf$T1
  flair = idf$FLAIR
  outfiles = c(idf$out_prob, idf$out_pred)
  
  if (!all(file.exists(outfiles))) {
    f_res = flexconnr::predict_flexconn(
      t1 = t1,
      flair = flair,
      num_atlases = num_atlases)
    file.copy(f_res, to = outfiles, overwrite = TRUE)
  }
  
}
