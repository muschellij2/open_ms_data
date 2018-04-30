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
library(oasis)
library(here)
library(readr)
library(RNifti)
library(caret)


set.seed(20180410)
rootdir = here::here()
df_file = here("cross_sectional", "raw", 
  "filename_df.rds")
prog_dir = file.path(rootdir, "programs")
source(file.path(prog_dir, "helper_functions.R"))
mod_dir = here("cross_sectional", "model")
dir.create(mod_dir, showWarnings = FALSE)

df = read_rds(df_file)
df = df %>% 
  filter(group == "train")
df$oasis_outfile = file.path(df$proc_dir,
  "oasis_train_list.rds")

model_file = here("cross_sectional", "model", 
    "oasis_model.rds")

cv_model_file = here("cross_sectional", "model", 
    "cv_oasis_model.rds")

n_ids = nrow(df)
all_df = vector(mode = "list",
  length = nrow(df))
names(all_df) = df$id

# read in the data
iid = 4

for (iid in seq(nrow(df))) {
  idf = df[iid,]
  ofile = idf$oasis_outfile

  if (file.exists(ofile)) {
    res = read_rds(ofile)
    all_df[[iid]] = res$oasis_dataframe
    # rm(res); gc()
  }
  print(iid)
}

full_df = bind_rows(all_df, .id = "id")


full_df$y = ifelse(full_df$GoldStandard > 0,
  "lesion", "non_lesion")
full_df$GoldStandard = factor(full_df$y,
  levels = c("non_lesion", "lesion"))
full_df$y = NULL


form = formals(oasis_training)$formula
form_terms = as.character(form)[[3]]
form_terms = strsplit(form_terms, "+", 
  fixed = TRUE)[[1]]
form_terms = form_terms[!grepl("PD", form_terms)]
form_terms = paste(form_terms, collapse= "+")
form_terms = paste0("GoldStandard ~ ", form_terms)
form = as.formula(form_terms)

model = glm(
  formula = form,
  data = full_df,
  family = binomial,
  control= list(trace = TRUE))
model = reduce_glm_mod(model)

write_rds(model, path = model_file)

myControl <- trainControl(
    method = "cv", number = 5,
    summaryFunction = twoClassSummary,
    classProbs = TRUE, # IMPORTANT!
    verboseIter = TRUE
  )

model <- train(
  form = form,
  data = full_df,
  family = binomial(),    
  method = "glm",
  trControl = myControl
)
model = reduce_train_object(model)    
write_rds(model, path = cv_model_file)

