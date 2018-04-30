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
library(mgcv)
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

run_groups = c("train", "test")
run_group = run_groups[2]

all_df = vector(mode = "list",
  length = length(run_groups))
names(all_df) = run_groups
all_roc = all_df

for (run_group in run_groups) {
  print(run_group)
  outfile = file.path(mod_dir, 
    paste0(run_group, "_", "roc_information.rds"))
  df = read_rds(outfile)
  all_roc[[run_group]] = df
  data_outfile = file.path(mod_dir, 
    paste0(run_group, "_", 
      "full_prediction_data.rds"))
  df = read_rds(data_outfile)
  all_df[[run_group]] = df  
}

cutoff = all_roc$train$non_smoothed$dice_cutoff
sm_cutoff = all_roc$train$smoothed$dice_cutoff

ddf = data_frame(
  cutoff = cutoff,
  sm_cutoff = sm_cutoff)

all_df$train$group = "train"
all_df$test$group = "test"
df = bind_rows(all_df$train, all_df$test)
df$pred = df$p > cutoff
df$spred = df$sm_p > sm_cutoff
res = df %>% 
  group_by(group, id) %>% 
  summarize(
    dice = dice(pred, Y),
    sm_dice = dice(spred, Y),
    n_y = sum(Y),
    dice_cutoff = run_roc(p, Y)$dice_cutoff,
    sm_dice_cutoff = run_roc(sm_p, Y)$dice_cutoff,
    n_pred = sum(pred),
    n_spred = sum(spred)
    )

train_df = res %>% 
  filter(group == "train")
mod = gam(dice_cutoff ~ s(n_pred), 
  data = train_df)

gm_mod = gam(dice_cutoff ~ s(n_pred), 
  data = train_df, family = quasipoisson(link = log))


test_df = res %>% 
  filter(group == "test") 
test_df$new_cutoff = 
  predict(mod, newdata = test_df)
test_df$new_gm_cutoff = 
  exp(predict(gm_mod, newdata = test_df, 
    link = "response")) 
test_df = left_join(
  test_df[, c("id", "new_cutoff", "new_gm_cutoff")], 
  all_df$test)
test_df$pred = test_df$p > test_df$new_cutoff

new_res = test_df %>% 
  group_by(group, id) %>% 
  summarize(
    dice = dice(pred, Y),
    n_pred = sum(pred)
    )

res %>% 
  ggplot(aes(y = dice_cutoff, x = n_pred, 
    colour = group)) +
  geom_point() + geom_smooth(se = FALSE) + 
  geom_hline(aes(yintercept = cutoff), 
    data = ddf)  

  res %>% 
    ggplot(aes(y = sm_dice_cutoff, x = n_spred,
      colour = group)) +
    geom_point() + geom_smooth(se = FALSE) + 
    geom_hline(aes(yintercept = sm_cutoff), 
      data = ddf)  


df = all_df$test
reg_dice = dice(
  df$p > cutoff,
  df$Y)
sm_dice = dice(
  df$sm_p > sm_cutoff,
  df$Y)
