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
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")
prog_dir = file.path(rootdir, "programs")
source(file.path(prog_dir, "helper_functions.R"))
mod_dir = here("cross_sectional", "model")
dir.create(mod_dir, showWarnings = FALSE)



model_file = here("cross_sectional", "model", 
    "first_train_model.rds")

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
    res = res %>% 
      select(Y, FLAIR_quantile, 
        WM, WM.mn, GM.mn)
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
rm(all_df); gc()

first_mod = glm(
  y ~ FLAIR_quantile + WM + WM.mn + GM.mn,
  data = full_df,
  control= list(trace = TRUE),
  family = binomial())

full_df$p = predict(first_mod)

first_mod = reduce_glm_mod(first_mod)

q01 = quantile(full_df$p[full_df$Y>0], 
  probs = 0.001)
q_prob01 = exp(q01)

# ggplot(aes(x = Y > 0, y = p), data = full_df) + 
#   geom_boxplot()


L = list(model = first_mod,
  prob_cutoff = q_prob01)
write_rds(L, path = model_file)

# samp = full_df %>% 
#     group_by(id, y) %>% 
#     sample_frac(size = 0.1) %>% 
#     ungroup()
# print(nrow(samp))
# rm(full_df); gc()

# cn = colnames(samp)
# keep = !cn %in% c("mask", "id", "Y")
# samp = samp[, keep]

# gc(); gc(); gc();

# y = samp$y
# samp = samp %>% 
#   select(-y)
  
# if (remove_t1_post) {
#   cn = colnames(samp)
#   keep = cn[!grepl("T1POST", cn)]
#   samp = samp[, keep]
# }
# myControl <- trainControl(
#     method = "cv", number = 5,
#     summaryFunction = twoClassSummary,
#     classProbs = TRUE, # IMPORTANT!
#     verboseIter = TRUE
#   )

# if (!file.exists(model_file)) {
#   model <- train(
#     x = samp,
#     y = y,
#     tuneLength = 1,
#     method = "ranger",
#     trControl = myControl,
#     importance = "permutation"
#   )
#   model = reduce_train_object(model)    
#   write_rds(model, path = model_file)
# }