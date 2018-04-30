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
library(smri.process)
library(here)
library(readr)
library(RNifti)
library(caret)

reduce_train_object = function(x) {
  x$control$index= NULL
  x$control$indexOut = NULL
  x$trainingData = NULL
  x$finalModel$predictions = NULL
  x
}

set.seed(20180410)
rootdir = here::here()
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")
prog_dir = file.path(rootdir, "programs")
source(file.path(prog_dir, "helper_functions.R"))
mod_dir = here("cross_sectional", "model")
dir.create(mod_dir, showWarnings = FALSE)



first_model_file = here(
  "cross_sectional", "model", 
    "first_train_model.rds")

remove_t1_post = c(FALSE, TRUE)
filtered = c(FALSE, TRUE)
eg = expand.grid(
  remove_t1_post = remove_t1_post,
  filtered = filtered,
  stringsAsFactors = FALSE)
eg$model_file = here("cross_sectional", "model", 
    paste0("train_model", 
    ifelse(eg$remove_t1_post, "_nopost", ""),
    ifelse(eg$filtered, "_filtered", ""),
    ".rds"))
isub = as.numeric(
  Sys.getenv("SGE_TASK_ID")
  )
if (is.na(isub) || isub < 1) {
  isub = 3
}
ieg = eg[isub,]

remove_t1_post = ieg$remove_t1_post
filtered = ieg$filtered
model_file = ieg$model_file


data_df = read_rds(df_file)

df = data_df %>% 
  filter(group == "train")


n_ids = nrow(df)
all_df = vector(mode = "list",
  length = nrow(df))
names(all_df) = df$id

# read in the data
iid = 1
for (iid in seq(nrow(df))) {
  idf = df[iid,]
  ofile = idf$out_df

  if (file.exists(ofile)) {
    res = read_rds(ofile)
    all_df[[iid]] = res
    rm(res); gc()
  }
  print(iid)
}

full_df = bind_rows(all_df, .id = "id")
full_df$y = ifelse(full_df$Y > 0,
  "lesion", "non_lesion")
full_df$y = factor(full_df$y,
  levels = c("non_lesion", "lesion"))
n_y = sum(full_df$Y > 0)
n_not_y = sum(full_df$Y == 0)
nr = nrow(full_df)
rm(all_df); gc()

if (filtered) {
  if (!file.exists(first_model_file)) {

    first_df = full_df %>% 
      select(y, FLAIR_quantile, 
        WM, WM.mn, GM.mn)

    first_mod = glm(
      y ~ FLAIR_quantile + WM + WM.mn + GM.mn,
      data = first_df,
      control= list(trace = TRUE),
      family = binomial())

    full_df$p = predict(first_mod, 
      type = "link")
    first_mod = reduce_glm_mod(first_mod)
    quantile_used = 0.005
    q01 = quantile(full_df$p[full_df$Y > 0], 
      probs = quantile_used)
    prob_cutoff = exp(q01)

    L = list(
      model = first_mod,
      prob_cutoff = prob_cutoff,
      quantile_used = quantile_used)

    write_rds(L, path = first_model_file)
  } else {
    L = read_rds(first_model_file)
    first_df = full_df %>% 
      select(y, FLAIR_quantile, 
        WM, WM.mn, GM.mn)    
    full_df$p = predict(L$model, 
      newdata = first_df, type = "response")
    prob_cutoff = L$prob_cutoff
  }
}

if (filtered) {
  full_df = full_df %>% 
    filter(p >= prob_cutoff) 
  full_df = full_df %>% 
    select(-p)
}
n_lesion = sum(full_df$Y > 0)
n_non_lesion = sum(full_df$Y == 0)

samp = full_df %>% 
    group_by(id, y) %>% 
    sample_frac(size = 0.1) %>% 
    ungroup()
print(nrow(samp))
rm(full_df); gc()

cn = colnames(samp)
keep = !cn %in% c("mask", "id", "Y")
samp = samp[, keep]

gc(); gc(); gc();

y = samp$y
samp = samp %>% 
  select(-y)
  
if (remove_t1_post) {
  cn = colnames(samp)
  keep = cn[!grepl("T1POST", cn)]
  samp = samp[, keep]
}
myControl <- trainControl(
    method = "cv", number = 5,
    summaryFunction = twoClassSummary,
    classProbs = TRUE, # IMPORTANT!
    verboseIter = TRUE
  )

if (!file.exists(model_file)) {
  model <- train(
    x = samp,
    y = y,
    tuneLength = 1,
    method = "ranger",
    trControl = myControl,
    importance = "permutation"
  )
  model = reduce_train_object(model)    
  write_rds(model, path = model_file)
}