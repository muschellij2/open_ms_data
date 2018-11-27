#!/usr/bin/env RScript
#######################
# Processing Cross-sectional data
#######################
rm(list = ls())
# library(dcmtk)
# library(dcm2niir)
library(RNifti)
library(fslr)
# library(divest)
library(dplyr)
library(matrixStats)
library(smri.process)
library(here)
library(readr)
library(ranger)
library(caret)

rerun = FALSE
reduce_train_object = function(x) {
  x$control$index= NULL
  x$control$indexOut = NULL
  x$trainingData = NULL
  x$finalModel$predictions = NULL
  x
}
reduce_ranger = function(x) {
  x$predictions = NULL
  x
}

set.seed(20180410)
rootdir = here::here()
# proc_dir = here("images")

# templates = c("none", "MNI", "Eve")
templates = "MNI"
itemplate = "MNI"
filtered = c(FALSE, TRUE)
remove_t1_post = c(FALSE, TRUE)
runner = c("ranger", "caret")
training = c("first15", "random")
# training = "first15"
eg = expand.grid(
  filtered = filtered,
  itemplate = templates,
  remove_t1_post = remove_t1_post,
  training = training,
  runner = runner,
  stringsAsFactors = FALSE)
eg$model_file = here(
  "model", 
  paste0(
    ifelse(eg$runner == "ranger", "ranger_", ""),
    "train_model", 
    ifelse(eg$remove_t1_post, "_not1post", ""),
    ifelse(eg$filtered, "_filtered", ""),
    ifelse(eg$training == "random", "", "_first15"),
    ".rds"))
isub = as.numeric(
  Sys.getenv("SGE_TASK_ID")
)
if (is.na(isub) || isub < 1) {
  isub = 4
}

ieg = eg[isub,]
model_file = ieg$model_file
itemplate = ieg$itemplate
filtered = ieg$filtered
runner = ieg$runner
remove_t1_post = ieg$remove_t1_post
training = ieg$training

# for (iid in seq(nrow(xdf))) {

if (!file.exists(model_file) | rerun) {

  # itemplate = "MNI"
  # for (itemplate in c("none", "MNI", "Eve")) {
  print(itemplate)
  pre = switch(itemplate,
               none = "",
               MNI = "MNI_",
               Eve = "Eve_")

  df_file = here("cross_sectional", "raw", 
                 paste0(pre, "filename_df.rds"))

  prog_dir = file.path(rootdir, "programs")
  source(file.path(prog_dir, "helper_functions.R"))
  mod_dir = here("model")
  dir.create(mod_dir, showWarnings = FALSE)

  first_model_file = here(
    "model", 
    paste0("first_train_model",
          ifelse(training == "random", "", "_first15"),
          ".rds"))

  data_df = read_rds(df_file)

  if (training == "random") {
    df = data_df %>% 
      filter(group == "train")
  } else {
    df = data_df %>% 
      filter(id %in% sprintf("patient%02.0f", 1:15))
  }


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
      if (remove_t1_post) {
        cn = colnames(res)
        cn = cn[ !grepl("T1POST", cn)]
        res = res[, cn, drop = FALSE]
      }   
      for (icol in colnames(res)) {
        res[ is.nan(res[, icol, drop = TRUE]), 
          icol] = 0
        res[ is.infinite(res[, icol, drop = TRUE]), 
          icol] = 0
      }
      stopifnot(!anyNA(res))  
      all_df[[iid]] = res
      rm(res); gc()
    } else {
      stop(paste(ofile, "does not exist"))
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
    # don't need to worry about T1POST here
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
                          newdata = first_df, 
                          type = "response")
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

  if (runner == "caret") {
    y = samp$y
    samp = samp %>% 
      select(-y)
    
    myControl <- trainControl(
      method = "cv", number = 5,
      summaryFunction = twoClassSummary,
      classProbs = TRUE, # IMPORTANT!
      verboseIter = TRUE
    )
    
      model <- train(
        x = samp,
        y = y,
        tuneLength = 1,
        method = "ranger",
        trControl = myControl,
        importance = "permutation",
        num.threads = 10
      )
      model = reduce_train_object(model)    
      write_rds(model, path = model_file)
  }
  if (runner == "ranger") {
      model = ranger(
        formula = y ~ .,
        data = samp,
        importance = "permutation",
        num.threads = 10,
        probability = TRUE
      )
      model = reduce_ranger(model)
      write_rds(model, path = model_file)
  }
}