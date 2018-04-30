# Rnosave run_predict.R -N PRED \
# -l mem_free=40G,h_vmem=41G -hold_jid MODEL
library(here)
library(dplyr)
library(readr)
library(tidyr)
library(neurobase)
library(pbapply)

rootdir = here::here()
prog_dir = file.path(rootdir, "programs")
source(file.path(prog_dir, "helper_functions.R"))
res_dir = here("cross_sectional", "results")

proc_dir = file.path(rootdir,
    "cross_sectional", "raw")
mod_dir = file.path(rootdir,
    "cross_sectional", "model")
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")
df = read_rds(df_file)

all_eg = NULL
models = c("ranger", "oasis")
groups = c("train")
for (model in models) {
  if (model == "ranger") {
    eg = expand.grid(
      remove_t1_post = c(FALSE, TRUE),
      filtered = c(FALSE, TRUE),
      group = groups,
      stringsAsFactors = FALSE)
    eg$fname = paste0("predictions", 
        ifelse(eg$remove_t1_post, "_nopost", ""),
        ifelse(eg$filtered, "_filtered", ""),
        ".rds")
    eg$mod = paste0("ranger", 
      ifelse(eg$remove_t1_post, "_nopost", ""),
      ifelse(eg$filtered, "_filtered", "")
      )
  } else {
    eg = expand.grid(
      cv = c(FALSE, TRUE),
      trained = c(TRUE, FALSE),
      nopost = c(FALSE, TRUE),
      group = groups,
      stringsAsFactors = FALSE)
    eg = eg[ !(!eg$trained & eg$cv), ]
    eg$mod = paste0("oasis",
      ifelse(eg$cv, "_cv", ""),
      ifelse(!eg$trained, "_untrained", ""))
  }
  all_eg = bind_rows(all_eg, eg)
}
eg = all_eg
eg$full_group = paste0(eg$group, "_", eg$mod)
eg = eg %>% 
  select(group, mod, full_group) %>% 
  distinct()
run_group = "train"
eg$outfile = file.path(mod_dir, 
      paste0(run_group, "_", 
        eg$mod, "_",
        "roc_information.rds"))

extract_cutoffs = function(outfile) {
  x = read_rds(outfile)
  data_frame(outfile = outfile,
    smoothed = x$smoothed$dice_cutoff,
    non_smoothed = x$non_smoothed$dice_cutoff)
}
cuts = pblapply(eg$outfile, extract_cutoffs)
cuts = bind_rows(cuts)
eg = left_join(eg, cuts, by = "outfile")

cutoff_df = gather(eg, smooth, cutoff, 
  smoothed, non_smoothed)
cutoff_df$smooth = !grepl("non", cutoff_df$smooth)


phat = list.files(path = proc_dir, 
  pattern = "phat[.]nii[.]gz$",
  recursive = TRUE)

img_df = data_frame(
  file = phat)
img_df = img_df %>% 
  mutate(fname = basename(file),
    id = sapply(strsplit(file, "/"), dplyr::first),
    mod = sub(".*\\d_(.*)_phat.*", "\\1", fname),
    smooth = grepl("smoothed", fname),
    smod = mod,
    mod = sub("_smoothed", "", mod)
    )
img_df$file = file.path(proc_dir, img_df$file)
eg_mods = sort(unique(cutoff_df$mod))
img_mods = sort(unique(img_df$mod))
stopifnot(all(img_mods %in% eg_mods))

img_df = left_join(img_df, cutoff_df)

df = df %>% 
  select(gs_file, id, image_file)
df = left_join(img_df, df)

df$image_plot_file = file.path(res_dir, 
  paste0(df$id, "_image.png"))
df$roi_plot_file = file.path(res_dir, 
  paste0(df$id, "_image_roi.png"))
df$plot_file = file.path(res_dir,
  paste0(nii.stub(df$file, bn = TRUE), ".png"))
iid = 1

for (iid in seq(nrow(df))) {
  
  idf = df[iid,]
  print(iid)
  
  prob_file = idf$file
  roi_fname = idf$gs_file
  
  img_fname = idf$image_file
  fnames = c(idf$plot_file,
             idf$image_plot_file,
             idf$roi_plot_file)
  run_cutoff = idf$cutoff
  
  if (!all(file.exists(fnames))) {
    
    img = readnii(img_fname)    
    roi = readnii(roi_fname)    
    check_mask_fail(roi)    
    # brain_mask = readnii(df$reg_brain_mask[iid])
    # check_mask_fail(brain_mask)    
    
    prob_img = readnii(prob_file) 
    
    img = robust_window(img)  
    plot_xyz = xyz(roi)
    ################################
    # Plot
    ################################
    if (!file.exists(idf$image_plot_file)) {
      png(idf$image_plot_file, 
          res = 600, width = 5, height=5,
          units = "in", type = "cairo")
      ortho2(img, xyz = plot_xyz,
        text = "Image", useRaster = FALSE)
      dev.off()
    }

    if (!file.exists(idf$roi_plot_file)) {
      png(idf$roi_plot_file, 
          res = 600, width = 5, height=5,
          units = "in", type = "cairo")
      ortho2(img, roi, xyz = plot_xyz,
        text = "ROI")
      dev.off()
    }    
    
    pred = prob_img > run_cutoff
    fname = idf$plot_file
    
    cols = c("#56B4E9", "#D55E00", "#009E73")
    levels = c("False Negative", 
               "False Positive", 
               "True Positive")  
    
    txt = idf$smod

    if (!file.exists(fname)) {
      
      diff = roi + pred*2
      diff[roi == 0 & pred == 0] = NA
      
      png(fname, 
          res = 600, width = 5, height=5,
          units = "in", type = "cairo")
      ortho2(x = img, 
             y = diff, 
             # don't do alpha blending
             col.y = cols,
             xyz = plot_xyz, 
             addlegend = TRUE,
             legend = levels, 
             leg.col = cols, 
             leg.cex = 1.5,
             ybreaks = c(0, 1.1, 2.1, 3.1),
             leg.title = txt,
             leg.x = 15,
             leg.y = 45)
      dev.off()
    }
  }
}

