reduce_glm_mod = function(model){
  model$y = c()
  model$model = c()
  model$residuals = c()
  model$fitted.values = c()
  model$effects = c()
  model$qr$qr = c()  
  model$linear.predictors = c()
  model$weights = c()
  model$prior.weights = c()
  model$data = c()
  attr(model$terms,".Environment") = c()
  attr(model$formula,".Environment") = c()
  model
}

opt.cut = function(perf){
    cut.ind = mapply(
      FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        ind = min(ind)
        c(sensitivity = y[[ind]], 
          specificity = 1-x[[ind]], 
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, 
    perf@alpha.values)
}

tab_dice = function(tab) {
  stopifnot(all(dim(tab) == c(2,2)))
  2*tab[2,2] / (2*tab[2,2] + tab[1,2] + tab[2,1])
}

my.tab = function(
  x, 
  y, 
  dnames=c("x", "y")) {
  x = as.numeric(x)
  y = as.numeric(y)
  stopifnot(all(unique(c(x,y)) %in% c(0, 1, NA)))
  tt = sum(x * y)
  t1 = sum(x)
  t2 = sum(y)
  tab = matrix(c(length(x) - t1 - t2 + tt,  t1 - tt, t2 - tt, tt), 2, 2)
  n = list(c("FALSE", "TRUE"), c("FALSE", "TRUE"))
  names(n) = dnames
  dimnames(tab) = n
  tab = as.table(tab)
  return(tab) 
}

dice = function(x, y) {
  tab = my.tab(c(x), c(y))
  tab_dice(tab)
}

opt.dice = function(pred){
    cut.ind = mapply(
      FUN=function(tp, fp, fn, p) {
        dice = 2*tp / (2*tp + fn + fp)
        ind = which(dice == max(dice))
        ind = min(ind)
        c(dice = dice[[ind]], 
          cutoff = p[[ind]])        
    }, pred@tp, pred@fp, pred@fn, 
    pred@cutoffs)
}


run_roc = function(x, y, fpr.stop = 0.01) {
  pred = prediction(x, y)
  auc = performance(pred, "auc", 
    fpr.stop = fpr.stop)
  pauc = unlist(auc@y.values)/fpr.stop
  print(pauc)
  perf = performance(pred, "tpr", "fpr")
  cutoff = opt.cut(perf)
  dice_cutoff = opt.dice(pred)

  cutoff = cutoff["cutoff",]
  dice_cutoff = dice_cutoff["cutoff",]
  return(list(pauc = pauc,
    perf = perf,
    cutoff = cutoff,
    dice_cutoff = dice_cutoff))
}


reduce_train_object = function(x) {
  x$control$index= NULL
  x$control$indexOut = NULL
  x$trainingData = NULL
  x$finalModel$predictions = NULL
  x
}



pred_to_df = function(df) {
  df = df %>% 
    dplyr::mutate(
      dice = (2 * tp) / (2 * tp + fp + fn),
      n.pos = tp + fn,
      n.neg = tn + fp,
      n.pos.pred = tp + fp,
      n.neg.pred = tn + fn      
      ) %>% 
    dplyr::mutate(      
      vol = n.pos,
      pred_vol = n.pos.pred)
  return(df)
}

condense_prediction = function(x,
  y,
  digits = 3,
  by = 10^{-digits}
  ) {

  library(magrittr)

  mult = 10^digits
  x = round(x * mult)
  # tab = table(x, y)
  # df = dplyr::as_data_frame(tab)
  # df$x = as.integer(df$x)ff
  # df$x = df$x / mult
  vals = seq(0, 1, by = by) * mult
  # vals = seq(by, 1-by, by = by) * mult
  vals = round(vals)
  tabs = pbapply::pbsapply(vals, function(val) {
    tab = my.tab(x > val, y)
    tab = c(tab)
  })
  tabs = t(tabs)
  colnames(tabs) = c("tn", "fp", "fn", "tp")
  tabs = dplyr::as_data_frame(tabs)
  tabs$cutoffs = vals / mult

  pred_df = pred_to_df(tabs)
  pred_df = pred_df %>% 
    dplyr::arrange(desc(cutoffs))
  pred_df
}



run_df_roc = function(df, fpr.stop = 0.01) {

  pred = new("prediction", 
    predictions = list(NA), labels = list(NA), 
      cutoffs = list(df$cutoffs), 
      fp = list(df$fp), 
      tp = list(df$tp), 
      fn = list(df$fn), 
      tn = list(df$tn),
      n.pos.pred = list(df$n.pos.pred),
      n.neg.pred = list(df$n.neg.pred),
      n.pos = list(df$n.pos),
      n.neg = list(df$n.neg)        
    )

  auc = performance(pred, "auc", 
    fpr.stop = fpr.stop)
  pauc = unlist(auc@y.values)/fpr.stop
  print(pauc)
  perf = performance(pred, "tpr", "fpr")
  cutoff = opt.cut(perf)
  dice_cutoff = opt.dice(pred)

  cutoff = cutoff["cutoff",]
  dice_cutoff = dice_cutoff["cutoff",]
  return(list(pauc = pauc,
    perf = perf,
    cutoff = cutoff,
    dice_cutoff = dice_cutoff))
}
