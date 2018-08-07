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
devtools::source_gist("524eade46135f6348140",
  filename = "ggplot_smooth_func.R")

rootdir = here::here()

transparent_legend =  theme(
  legend.background = element_rect(
    fill = "transparent"),
  legend.key = element_rect(fill = "transparent", 
                            color = "transparent") )

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
groups = c("train", "test")
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
    eg$suffix = paste0("ranger", 
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


n_ids = nrow(eg)
all_roc = vector(mode = "list",
  length = nrow(eg))
all_df = all_roc
names(all_df) = eg$full_group
names(all_roc) = eg$full_group

ieg = 1
for (ieg in seq(nrow(eg))) {
  print(ieg)
  outfile = eg$outfile[ieg]
  x = read_rds(outfile)
  all_roc[[ieg]] = x
  outfile = eg$data_outfile[ieg]
  x = read_rds(outfile)
  x$Y = x$Y > 0
  x = gather(x, key = "smooth", 
    value = "prob", p, sm_p)  
  x$smooth = x$smooth == "sm_p"
  all_df[[ieg]] = x
}


##################################
# Getting some that good ROC action
##################################
perfs = lapply(all_roc,
  function(x) {
    df = data_frame(
      y = x$smoothed$perf@y.values[[1]],
      x = x$smoothed$perf@x.values[[1]], 
      smooth = TRUE)
    df2 = data_frame(
      y = x$non_smoothed$perf@y.values[[1]],
      x = x$non_smoothed$perf@x.values[[1]], 
      smooth = FALSE)    
    df = bind_rows(df, df2)
  })
perfs = bind_rows(perfs, .id = "full_group")

perfs = perfs %>% 
  filter(x <= 0.10)

eg = eg %>% 
  select(group, suffix, full_group)

###################################
# Combine All the data
###################################
paucs = lapply(all_roc,
  function(x) {
    data_frame(
      pauc = 
      c(x$smoothed$pauc, 
        x$non_smoothed$pauc),
      smooth = c(TRUE, FALSE))        
  })
paucs = bind_rows(paucs, .id = "full_group")

dice_cutoffs = lapply(all_roc,
  function(x) {
    data_frame(
      dice_cutoff = 
      c(x$smoothed$dice_cutoff, 
        x$non_smoothed$dice_cutoff),
      smooth = c(TRUE, FALSE))    
  })
dice_cutoffs = bind_rows(dice_cutoffs, 
  .id = "full_group")

cutoffs = lapply(all_roc,
  function(x) {
    data_frame(
      cutoff = 
      c(x$smoothed$cutoff, x$non_smoothed$cutoff),
      smooth = c(TRUE, FALSE))
  })
cutoffs = bind_rows(cutoffs, 
  .id = "full_group")

res = left_join(paucs, dice_cutoffs)
res = left_join(res, cutoffs)
rm(paucs, dice_cutoffs, cutoffs)
rm(all_roc)


# ss = strsplit(perfs$full_group, "_")
# perfs$group = sapply(ss, 
#   dplyr::first)
# perfs$run_model = sapply(
#   ss, function(x) paste(x[2:length(x)], 
#     collapse = "_"))

pngname = file.path(res_dir, 
  "pROC_Curves.png")
png(pngname, height = 8, width = 8, 
  res = 600, units = "in")
perfs %>% 
  filter(x <= 0.05) %>% 
  ggplot(aes(x = x, y = y, colour = run_model)) +
  geom_line() +
  facet_wrap(group ~ smooth)
dev.off()

pngname = file.path(res_dir, 
  "ranger_pROC_Curves.png")
png(pngname, height = 8, width = 8, 
  res = 600, units = "in")
perfs %>% 
  filter(x <= 0.05) %>% 
  filter(grepl("ranger", run_model)) %>% 
  filter(smooth) %>% 
  ggplot(aes(x = x, y = y, colour = run_model)) +
  geom_line() +
  facet_wrap(~ group)
dev.off()

rr = split(res, res$full_group)
# reordering so always using train cutoffs
full_groups = names(all_df)
rename_groups = sub("test", "train", full_groups)
rr = lapply(rename_groups, function(x) rr[[x]])
names(rr) = rename_groups
mapply(function(r, g) {
  stopifnot(all(r$full_group == g))
}, rr, rename_groups)


dice_res = mapply(function(df, cutoff) {
  cutoff = cutoff %>% 
    select(smooth, dice_cutoff)
  df = left_join(df, cutoff)
  df$pred = df$prob > df$dice_cutoff
  all_dice = df %>% 
    group_by(smooth) %>% 
    summarize(dice = dice(pred, Y),
      est_n_voxels = sum(pred),
      n_voxels = sum(Y)
      )
  all_dice$id = "overall"
  result = df %>% 
    group_by(id, smooth) %>% 
    summarize(dice = dice(pred, Y),
      est_n_voxels = sum(pred),
      n_voxels = sum(Y))
  result = bind_rows(all_dice, result)
  result
}, all_df, rr, SIMPLIFY = FALSE)
names(dice_res) = names(all_df)


dice_res = bind_rows(dice_res, 
  .id = "full_group")
dice_res$group = sapply(
  strsplit(dice_res$full_group, "_"), 
  dplyr::first)
dice_res$run_model = sapply(
  strsplit(dice_res$full_group, "_"), 
  function(x) paste(x[2:length(x)], 
    collapse = "_"))
dice_res = dice_res %>% 
  arrange(group, desc(dice))

overall = dice_res %>% 
  filter(id == "overall") 

patient_data = dice_res %>% 
  filter(id != "overall")  


patient_data = patient_data %>% 
  left_join(df[, c("id", "volume", "voxres")], by = "id")  

patient_data = patient_data %>% 
  mutate(est_volume = est_n_voxels * voxres / 1000,
    true_volume = n_voxels * voxres / 1000)



pngname = file.path(res_dir, 
  "volume_results.png")
png(pngname, height = 4, width = 8, 
  res = 600, units = "in")

vol_plot = patient_data %>% 
  filter(run_model == "ranger_nopost") %>% 
  filter(smooth) %>% 
  ggplot(aes(
    x = est_volume, 
    y = true_volume)) + 
  facet_wrap(~ group) +  
  geom_point() +
  geom_abline(intercept = 0, slope = 1) + 
  stat_smooth_func(geom="text",
      method="lm",
      hjust=0,parse=TRUE) +
  geom_smooth(method = "lm", se = FALSE) +
  ylab("Automatic Lesion Volume (mL)") +
  xlab("Manual Lesion Volume (mL)")  +
  theme(text = element_text(size = 20)) 
print(vol_plot)
dev.off()

pngname = file.path(res_dir, 
  "volume_results_oasis.png")
png(pngname, height = 4, width = 8, 
  res = 600, units = "in")
res = patient_data %>% 
  filter(run_model == "oasis") %>% 
  filter(smooth) 
print({ vol_plot %+% res  })
dev.off()



pngname = file.path(res_dir, 
  "dice_results.png")
png(pngname, height = 8, width = 8, 
  res = 600, units = "in")
  g = patient_data %>% 
    ggplot(aes(x = run_model, 
      y = dice, colour = smooth)) + 
    facet_wrap(~ group) +
    geom_boxplot() + 
    theme(text = element_text(size = 20))
  gg = g + 
    geom_point(data = overall,
      position = position_dodge(width = 0.75),
      shape = 2) + coord_flip()
  print(gg)
dev.off()


pngname = file.path(res_dir, 
  "ranger_dice_results.png")
png(pngname, height = 8, width = 8, 
  res = 600, units = "in")
d = patient_data %>% 
  filter(grepl("ranger", run_model)) 
gg = g %+% d
gg = gg + 
  geom_point(data = overall %>% 
      filter(grepl("ranger", run_model)),
    position = position_dodge(width = 0.75),
    shape = 2) + coord_flip()
print(gg)
dev.off()



pngname = file.path(res_dir, 
  "compare_dice_results.png")
png(pngname, height = 8, width = 8, 
  res = 600, units = "in")

  d = patient_data %>% 
    mutate(type = ifelse(
      grepl("ranger", run_model),
      "Random Forest", 
      "OASIS")
    ) %>% 
    filter(smooth) %>% 
    filter(
      run_model %in% c("oasis", "oasis_untrained", 
        "ranger", "ranger_nopost")) %>% 
    mutate(
      run_model = recode(run_model,
        oasis = "OASIS Trained",
        oasis_untrained = "OASIS",
        ranger_nopost = "RF, no T1Post",
        ranger = "RF")
        )
  o = overall %>% 
    mutate(type = ifelse(
      grepl("ranger", run_model),
      "Random Forest", 
      "OASIS")
    )   %>% filter(smooth) %>% 
    filter(
      run_model %in% c("oasis", "oasis_untrained", 
        "ranger", "ranger_nopost")) %>% 
    mutate(
      run_model = recode(run_model,
        oasis = "OASIS Trained",
        oasis_untrained = "OASIS",
        ranger_nopost = "RF, no T1Post",
        ranger = "RF")
        )

  g = d %>% 
    ggplot(aes(x = run_model, 
      y = dice, colour = type)) + 
    facet_wrap(~ group) +
    geom_boxplot() + 
    theme(text = element_text(size = 20)) + 
    theme(legend.position = c(0.7, 0.6)) +
    ylab("Dice Coefficient") + xlab("") +
    transparent_legend + 
    guides(colour = guide_legend(title = "Model"))
  gg = g + 
    geom_point(data = o,
      position = position_dodge(width = 0.75),
      shape = 2) +  coord_flip()

  print(gg)
dev.off()



pngname = file.path(res_dir, 
  "ranger_tapas_compare.png")
png(pngname, height = 8, width = 8, 
  res = 600, units = "in")
patient_data %>% 
  filter(run_model == "ranger_nopost") %>% 
  ggplot(aes(x = volume, 
    y = dice, colour = smooth)) + 
  facet_wrap(smooth ~ group) +  
  geom_point()
dev.off()

pngname = file.path(res_dir, 
  "ranger_tapas.png")
png(pngname, height = 8, width = 8, 
  res = 600, units = "in")
patient_data %>% 
  filter(run_model == "ranger_nopost") %>% 
  filter(smooth) %>% 
  ggplot(aes(x = volume, 
    y = dice)) + 
  facet_wrap(~ group) +  
  geom_point(size = 3) + 
  theme(text = element_text(size = 20)) + 
  ylab("Dice Coefficient") + 
  xlab("Manual Lesion Volume (mL)")
dev.off()


rm(all_df)
