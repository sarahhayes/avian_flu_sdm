# In this script we train a BART model to classify sites for avian flu
# presence/absence based on environmental and species abundance factors.
# This script includes cross-validation of BART parameters.

SAVE_FITS <- TRUE

# Global storing optimised k value for BART - which may be updated
K_OPT <- 2

# Set path to folder containing data, and where output will be stored
PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(dplyr)
library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

#### EMBARCADERO FUNCTIONS ADJUSTED TO ALLOW CUSTOM BART PARAMETERS: ####

bart.flex <- function(x.data, y.data,
                      y.name = NULL,
                      n.trees = 200,
                      k = 2.0, power = 2.0, base = 0.95) {
  
  train <- cbind(y.data, x.data) 
  if(!is.null(y.name)) {colnames(train)[1] <- y.name}
  train <- na.omit(train)
  model <- bart(y.train = train[,1], 
                x.train = train[,2:ncol(train)], 
                k = k,
                power = power,
                base = base,
                ntree = n.trees, keeptrees=TRUE)
  return(model)
}

varimp.diag <- function(x.data, y.data, iter=50, quiet=FALSE,
                        k = 2.0, power = 2.0, base = 0.95) {
  
  nvars <- ncol(x.data)
  varnums <- c(1:nvars)
  varlist <- colnames(x.data)
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY 4 THIS CODE :) 
  
  ###############
  
  # auto-drops 
  
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               n.trees = 200,
                               k = k,
                               power = power,
                               base = base))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  
  if(length(dropnames)==0) {} else{
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames,collapse = ' '), ' \n')
  }
  
  x.data %>% dplyr::select(-dropnames) -> x.data  
  
  ###############
  
  for (n.trees in c(10, 20, 50, 100, 150, 200)) {
    
    cat(paste('\n', n.trees, 'tree models:', iter, 'iterations\n'))
    if(!quiet){pb <- txtProgressBar(min = 0, max = iter, style = 3)}
    
    for(index in 1:iter) {
      quietly(model.j <- bart.flex(x.data = x.data[,varnums], y.data = y.data, 
                                   n.trees = n.trees,
                                   k = k,
                                   power = power,
                                   base = base))
      
      vi.j <- varimp(model.j)
      if(index==1) {
        vi.j.df <- vi.j
      } else {
        vi.j.df[,index+1] <- vi.j[,2]
      }
      if(!quiet){setTxtProgressBar(pb, index)}
    }
    vi.j <- data.frame(vi.j.df[,1],
                       rowMeans(vi.j.df[,-1]))
    
    if(n.trees==10) { vi <- vi.j } else {  vi <- cbind(vi,vi.j[,2])  }
  }
  
  colnames(vi) <- c('variable','10','20','50','100','150','200')
  vi <- reshape::melt(vi, "variable")
  colnames(vi) <- c('variable','trees','imp')
  
  vi %>% group_by(variable) %>% summarise(max = max(imp)) -> vi.fac
  vi.fac <- vi.fac[order(-vi.fac$max),]
  
  vi$names <- factor(vi$variable, levels=vi.fac$variable)
  
  g1 <- ggplot2::ggplot(vi, aes(y=imp, x=names, group=trees)) +
    geom_line(aes(color=trees)) + geom_point(size=3) + theme_classic() +
    ylab("Relative contribution\n") + xlab("\nVariables dropped") +
    ggpubr::rotate_x_text(angle = 35) + 
    theme(axis.text = element_text(size=10),
          axis.title = element_text(size=14,face="bold")); print(g1)
  
}

variable.step <- function(x.data, y.data, n.trees=10, iter=50, quiet=FALSE,
                          k = 2.0, power = 2.0, base = 0.95) {
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY
  
  comp <- complete.cases(x.data)
  
  if(length(comp) < (nrow(x.data))) {
    message("Some rows with NA's have been automatically dropped. \n")
  }
  x.data <- x.data[comp,]
  y.data <- y.data[comp]
  
  ###############
  
  # auto-drops 
  
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               n.trees = 200,
                               k = k,
                               power = power,
                               base = base))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  
  if(length(dropnames) > 0) {
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames,collapse = ' '), ' \n')
  }
  
  x.data %>% dplyr::select(-any_of(dropnames)) -> x.data  
  
  ###############
  
  nvars <- ncol(x.data)
  varnums <- c(1:nvars)
  varlist.orig <- varlist <- colnames(x.data)
  
  rmses <- data.frame(Variable.number=c(),RMSE=c())
  dropped.varlist <- c()
  
  for(var.j in c(nvars:3)) {
    
    print(noquote(paste("Number of variables included:",var.j)))
    print(noquote("Dropped:"))
    print(if(length(dropped.varlist)==0) {noquote("")} else {noquote(dropped.varlist)})
    
    rmse.list <- c()
    
    if(!quiet){pb <- txtProgressBar(min = 0, max = iter, style = 3)}
    for(index in 1:iter) {
      quietly(model.j <- bart.flex(x.data = x.data[,varnums], y.data = y.data, 
                                   n.trees = n.trees,
                                   k = k,
                                   power = power,
                                   base = base))
      
      quietly(vi.j <- varimp(model.j))
      if(index==1) {
        vi.j.df <- vi.j
      } else {
        vi.j.df[,index+1] <- vi.j[,2]
      }
      
      pred.p <- colMeans(pnorm(model.j$yhat.train))[y.data==1]
      pred.a <- colMeans(pnorm(model.j$yhat.train))[y.data==0]
      #e <- evaluate(p=pred.p,
      #              a=pred.a)
      #aucs <- rbind(aucs,c(var.j,e@auc)); colnames(aucs) <- c('Vars','AUC')
      
      pred.c <- c(pred.p, pred.a)
      true.c <- c(rep(1,length(pred.p)), rep(0,length(pred.a)))
      rmsej.i <- Metrics::rmse(true.c,pred.c)
      rmse.list <- c(rmse.list,rmsej.i)
      if(!quiet){setTxtProgressBar(pb, index)}
    }
    
    vi.j <- data.frame(vi.j.df[,1],
                       rowMeans(vi.j.df[,-1]))
    vi.j <- vi.j[order(vi.j[,2]),]
    
    drop.var <- vi.j[1,1]
    dropped.varlist <- c(dropped.varlist,as.character(drop.var))
    
    rmsej <- mean(rmse.list)
    
    rmses <- rbind(rmses,c(nvars-var.j,rmsej)); colnames(rmses) <- c('VarsDropped','RMSE')
    
    varnums <- varnums[!(varnums==which(varlist.orig==drop.var))]
    varlist <- varlist.orig[varnums]
    print(noquote("---------------------------------------"))
  }
  
  g1 <- ggplot2::ggplot(rmses, aes(y=RMSE, x=VarsDropped)) +
    geom_line(color="black") + geom_point(size=3) + theme_bw() +
    ylab("RMSE of model\n") + xlab("\nVariables dropped") +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14,face="bold")) +
    scale_x_discrete(limits=c(0:(nrow(rmses)))); print(g1)
  
  print(noquote("---------------------------------------"))
  print(noquote("Final recommended variable list"))
  varlist.final <- varlist.orig[!(varlist.orig %in% dropped.varlist[0:(which(rmses$RMSE==min(rmses$RMSE))-1)])]
  print(noquote(varlist.final))
  invisible(varlist.final)
}

bart.step <- function(x.data, y.data,
                      iter.step=100, tree.step=10,
                      iter.plot=100,
                      k = 2.0, power = 2.0, base = 0.95,
                      full=FALSE,
                      quiet=FALSE) {
  
  ###############
  
  # auto-drops 
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY
  
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               n.trees = 200,
                               k = k,
                               power = power,
                               base = base))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  
  if(length(dropnames)==0) {} else{
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames,collapse = ' '), ' \n')
  }
  
  x.data %>% dplyr::select(-dropnames) -> x.data  
  
  ###############
  
  quiet2 <- quiet
  if(full==TRUE){varimp.diag(x.data, y.data, iter=iter.plot, quiet=quiet2)}
  vs <- variable.step(x.data,
                      y.data,
                      n.trees=tree.step,
                      iter=iter.step,
                      quiet=quiet2)
  
  invisible(best.model <- bart.flex(x.data = x.data[,vs], y.data = y.data, 
                                    n.trees=200,
                                    k = k,
                                    power = power,
                                    base = base))
  if(full==TRUE){varimp(best.model, plots=TRUE)}
  if(full==TRUE) {p <- summary(best.model, plots=TRUE)
  print(p)} else 
  {p <- summary(best.model, plots=FALSE)
  print(p)}
  invisible(best.model)
}

get_threshold <- function(object){
  fitobj <- object$fit
  
  true.vector <- fitobj$data@y 
  
  pred <- prediction(colMeans(pnorm(object$yhat.train)), true.vector)
  
  perf.tss <- performance(pred,"sens","spec")
  tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
  tss.df <- data.frame(alpha=perf.tss@alpha.values[[1]],tss=tss.list)
  thresh <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))])
  return(thresh)
}

get_sens_and_spec <- function(sdm, xtest, ytest, ri, cutoff){
  pred <- xtest[which(complete.cases(xtest)), ] %>%
    stats::predict(object=sdm, type = "bart", group.by=ri) %>%
    pnorm() %>%
    colMeans() %>%
    prediction(labels = ytest[which(complete.cases(xtest))])
  perf <-  performance(pred, measure = "sens", x.measure = "spec")
  tss_list <- (perf@x.values[[1]] + perf@y.values[[1]] - 1)
  tss_df <- data.frame(alpha=perf@alpha.values[[1]],tss=tss_list)
  # cutoff <- min(tss_df$alpha[which(tss_df$tss==max(tss_df$tss))])
  sens <- perf@x.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  spec <- perf@y.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  tss <- tss_df[which.min(abs(tss_df$alpha-cutoff)),'tss']
  auc <- performance(pred,"auc")@y.values[[1]]
  return(pairlist("pred"=pred,
                  "perf"=perf,
                  "sens"=sens,
                  "spec"=spec,
                  "tss"=tss,
                  "auc"=auc))
}



#### Period A Q1 ####

training_data <- read.csv("training_sets/training_data_A_Q1.csv")
xtrain <- training_data %>% dplyr::select(!("y"|"ri"))
ytrain <- training_data$y


# Fold construction

fold_ids <- caret::createFolds(training_data$y, k = 5)

folds <- lapply(1:length(fold_ids),
                FUN = function(i){
                  training_data[-fold_ids[[i]],]})
antifolds <- lapply(1:length(fold_ids),
                    FUN = function(i){
                      training_data[fold_ids[[i]],]})

## Check training representation
# lapply(folds, function(x) with(x, table(ri,y)))


k_vals = c(1, 2, 3)
power_vals = c(1.6, 1.8, 2)
base_vals = c(0.75, 0.85, 0.95)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         auc1=numeric(),
                         auc2=numeric(),
                         auc3=numeric(),
                         auc4=numeric(),
                         auc5=numeric())
cv_results[1:kl*pl*bl, ] <- 0
for (i in 1:length(k_vals)){
  k_val <- k_vals[i]
  for (j in 1:length(power_vals)){
    power_val <- power_vals[j]
    for (m in 1:length(base_vals)){
      base_val <- base_vals[m]
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      cat("idx=", idx, "\n")
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){
        model <- bart2(y ~ . - ri,
                          folds[[fold_no]],
                          test = antifolds[[fold_no]],
                          k = k_val,
                          power = power_val,
                          base = base_val,
                          n.trees = 200,
                          n.chains = 1L,
                          n.threads = 1L,
                          keepTrees = TRUE,
                          verbose = FALSE)
        antifold_x <- subset(antifolds[[fold_no]], select=-c(y, ri))
        antifold_y <- antifolds[[fold_no]]$y
        antifold_ri <- antifolds[[fold_no]]$ri
        cutoff <- get_threshold(model)
        auc <- get_sens_and_spec(model, antifold_x, antifold_y, antifold_ri, cutoff)$auc
        cv_results[idx, 3+fold_no] <- auc
        rm(model)
      }
    }
  }
}

cv_results <- cv_results %>%
  rowwise() %>%
  mutate(mean_err = mean(auc1,
                         auc2,
                         auc3,
                         auc4,
                         auc5))
argmin <- which.max(cv_results$mean_err)
K_OPT <- cv_results$k[argmin]
power_opt <- cv_results$power[argmin]
base_opt <- cv_results$base[argmin]


test_data <- read.csv("training_sets/test_data_A_Q1.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y

# Initialise model
basic_model <- bart( xtrain,
                         ytrain,
                     x.test = xtest,
                         k = K_OPT,
                         power = power_opt,
                         base = base_opt,
                     keeptrees = TRUE)
invisible(basic_model$fit$state)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_A_Q1.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = K_OPT,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)

if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_A_Q1.rds", sep = ""))
}




#### Period A Q2 ####

training_data <- read.csv("training_sets/training_data_A_Q2.csv")
xtrain <- training_data %>% dplyr::select(!("y"|"ri"))
ytrain <- training_data$y


# Fold construction

fold_ids <- caret::createFolds(training_data$y, k = 5)

folds <- lapply(1:length(fold_ids),
                FUN = function(i){
                  training_data[-fold_ids[[i]],]})
antifolds <- lapply(1:length(fold_ids),
                    FUN = function(i){
                      training_data[fold_ids[[i]],]})

## Check training representation
# lapply(folds, function(x) with(x, table(ri,y)))


k_vals = c(1, 2, 3)
power_vals = c(1.6, 1.8, 2)
base_vals = c(0.75, 0.85, 0.95)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         auc1=numeric(),
                         auc2=numeric(),
                         auc3=numeric(),
                         auc4=numeric(),
                         auc5=numeric())
cv_results[1:kl*pl*bl, ] <- 0
for (i in 1:length(k_vals)){
  k_val <- k_vals[i]
  for (j in 1:length(power_vals)){
    power_val <- power_vals[j]
    for (m in 1:length(base_vals)){
      base_val <- base_vals[m]
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      cat("idx=", idx, "\n")
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){
        model <- bart2(y ~ . - ri,
                       folds[[fold_no]],
                       test = antifolds[[fold_no]],
                       k = k_val,
                       power = power_val,
                       base = base_val,
                       n.trees = 200,
                       n.chains = 1L,
                       n.threads = 1L,
                       keepTrees = TRUE,
                       verbose = FALSE)
        antifold_x <- subset(antifolds[[fold_no]], select=-c(y, ri))
        antifold_y <- antifolds[[fold_no]]$y
        antifold_ri <- antifolds[[fold_no]]$ri
        cutoff <- get_threshold(model)
        auc <- get_sens_and_spec(model, antifold_x, antifold_y, antifold_ri, cutoff)$auc
        cv_results[idx, 3+fold_no] <- auc
        rm(model)
      }
    }
  }
}

cv_results <- cv_results %>%
  rowwise() %>%
  mutate(mean_err = mean(auc1,
                         auc2,
                         auc3,
                         auc4,
                         auc5))
argmin <- which.max(cv_results$mean_err)
K_OPT <- cv_results$k[argmin]
power_opt <- cv_results$power[argmin]
base_opt <- cv_results$base[argmin]


test_data <- read.csv("training_sets/test_data_A_Q2.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    k = K_OPT,
                         power = power_opt,
                         base = base_opt)
invisible(basic_model$fit$state)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_A_Q2.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = K_OPT,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)

if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_A_Q2.rds", sep = ""))
}



#### Period A Q3 ####

training_data <- read.csv("training_sets/training_data_A_Q3.csv")
xtrain <- training_data %>% dplyr::select(!("y"|"ri"))
ytrain <- training_data$y


# Fold construction

fold_ids <- caret::createFolds(training_data$y, k = 5)

folds <- lapply(1:length(fold_ids),
                FUN = function(i){
                  training_data[-fold_ids[[i]],]})
antifolds <- lapply(1:length(fold_ids),
                    FUN = function(i){
                      training_data[fold_ids[[i]],]})

## Check training representation
# lapply(folds, function(x) with(x, table(ri,y)))


k_vals = c(1, 2, 3)
power_vals = c(1.6, 1.8, 2)
base_vals = c(0.75, 0.85, 0.95)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         auc1=numeric(),
                         auc2=numeric(),
                         auc3=numeric(),
                         auc4=numeric(),
                         auc5=numeric())
cv_results[1:kl*pl*bl, ] <- 0
for (i in 1:length(k_vals)){
  k_val <- k_vals[i]
  for (j in 1:length(power_vals)){
    power_val <- power_vals[j]
    for (m in 1:length(base_vals)){
      base_val <- base_vals[m]
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      cat("idx=", idx, "\n")
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){
        model <- bart2(y ~ . - ri,
                       folds[[fold_no]],
                       test = antifolds[[fold_no]],
                       k = k_val,
                       power = power_val,
                       base = base_val,
                       n.trees = 200,
                       n.chains = 1L,
                       n.threads = 1L,
                       keepTrees = TRUE,
                       verbose = FALSE)
        antifold_x <- subset(antifolds[[fold_no]], select=-c(y, ri))
        antifold_y <- antifolds[[fold_no]]$y
        antifold_ri <- antifolds[[fold_no]]$ri
        cutoff <- get_threshold(model)
        auc <- get_sens_and_spec(model, antifold_x, antifold_y, antifold_ri, cutoff)$auc
        cv_results[idx, 3+fold_no] <- auc
        rm(model)
      }
    }
  }
}

cv_results <- cv_results %>%
  rowwise() %>%
  mutate(mean_err = mean(auc1,
                         auc2,
                         auc3,
                         auc4,
                         auc5))
argmin <- which.max(cv_results$mean_err)
K_OPT <- cv_results$k[argmin]
power_opt <- cv_results$power[argmin]
base_opt <- cv_results$base[argmin]


test_data <- read.csv("training_sets/test_data_A_Q3.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    k = K_OPT,
                         power = power_opt,
                         base = base_opt)
invisible(basic_model$fit$state)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_A_Q3.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = K_OPT,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)

if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_A_Q3.rds", sep = ""))
}



#### Period A Q4 ####

training_data <- read.csv("training_sets/training_data_A_Q4.csv")
xtrain <- training_data %>% dplyr::select(!("y"|"ri"))
ytrain <- training_data$y


# Fold construction

fold_ids <- caret::createFolds(training_data$y, k = 5)

folds <- lapply(1:length(fold_ids),
                FUN = function(i){
                  training_data[-fold_ids[[i]],]})
antifolds <- lapply(1:length(fold_ids),
                    FUN = function(i){
                      training_data[fold_ids[[i]],]})

## Check training representation
# lapply(folds, function(x) with(x, table(ri,y)))


k_vals = c(1, 2, 3)
power_vals = c(1.6, 1.8, 2)
base_vals = c(0.75, 0.85, 0.95)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         auc1=numeric(),
                         auc2=numeric(),
                         auc3=numeric(),
                         auc4=numeric(),
                         auc5=numeric())
cv_results[1:kl*pl*bl, ] <- 0
for (i in 1:length(k_vals)){
  k_val <- k_vals[i]
  for (j in 1:length(power_vals)){
    power_val <- power_vals[j]
    for (m in 1:length(base_vals)){
      base_val <- base_vals[m]
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      cat("idx=", idx, "\n")
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){
        model <- bart2(y ~ . - ri,
                       folds[[fold_no]],
                       test = antifolds[[fold_no]],
                       k = k_val,
                       power = power_val,
                       base = base_val,
                       n.trees = 200,
                       n.chains = 1L,
                       n.threads = 1L,
                       keepTrees = TRUE,
                       verbose = FALSE)
        antifold_x <- subset(antifolds[[fold_no]], select=-c(y, ri))
        antifold_y <- antifolds[[fold_no]]$y
        antifold_ri <- antifolds[[fold_no]]$ri
        cutoff <- get_threshold(model)
        auc <- get_sens_and_spec(model, antifold_x, antifold_y, antifold_ri, cutoff)$auc
        cv_results[idx, 3+fold_no] <- auc
        rm(model)
      }
    }
  }
}

cv_results <- cv_results %>%
  rowwise() %>%
  mutate(mean_err = mean(auc1,
                         auc2,
                         auc3,
                         auc4,
                         auc5))
argmin <- which.max(cv_results$mean_err)
K_OPT <- cv_results$k[argmin]
power_opt <- cv_results$power[argmin]
base_opt <- cv_results$base[argmin]


test_data <- read.csv("training_sets/test_data_A_Q4.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    k = K_OPT,
                         power = power_opt,
                         base = base_opt)
invisible(basic_model$fit$state)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_A_Q4.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = K_OPT,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)

if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_A_Q4.rds", sep = ""))
}



#### Period B Q1 ####

training_data <- read.csv("training_sets/training_data_B_Q1.csv")
xtrain <- training_data %>% dplyr::select(!("y"|"ri"))
ytrain <- training_data$y


# Fold construction

fold_ids <- caret::createFolds(training_data$y, k = 5)

folds <- lapply(1:length(fold_ids),
                FUN = function(i){
                  training_data[-fold_ids[[i]],]})
antifolds <- lapply(1:length(fold_ids),
                    FUN = function(i){
                      training_data[fold_ids[[i]],]})

## Check training representation
# lapply(folds, function(x) with(x, table(ri,y)))


k_vals = c(1, 2, 3)
power_vals = c(1.6, 1.8, 2)
base_vals = c(0.75, 0.85, 0.95)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         auc1=numeric(),
                         auc2=numeric(),
                         auc3=numeric(),
                         auc4=numeric(),
                         auc5=numeric())
cv_results[1:kl*pl*bl, ] <- 0
for (i in 1:length(k_vals)){
  k_val <- k_vals[i]
  for (j in 1:length(power_vals)){
    power_val <- power_vals[j]
    for (m in 1:length(base_vals)){
      base_val <- base_vals[m]
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      cat("idx=", idx, "\n")
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){
        model <- bart2(y ~ . - ri,
                       folds[[fold_no]],
                       test = antifolds[[fold_no]],
                       k = k_val,
                       power = power_val,
                       base = base_val,
                       n.trees = 200,
                       n.chains = 1L,
                       n.threads = 1L,
                       keepTrees = TRUE,
                       verbose = FALSE)
        antifold_x <- subset(antifolds[[fold_no]], select=-c(y, ri))
        antifold_y <- antifolds[[fold_no]]$y
        antifold_ri <- antifolds[[fold_no]]$ri
        cutoff <- get_threshold(model)
        auc <- get_sens_and_spec(model, antifold_x, antifold_y, antifold_ri, cutoff)$auc
        cv_results[idx, 3+fold_no] <- auc
        rm(model)
      }
    }
  }
}

cv_results <- cv_results %>%
  rowwise() %>%
  mutate(mean_err = mean(auc1,
                         auc2,
                         auc3,
                         auc4,
                         auc5))
argmin <- which.max(cv_results$mean_err)
K_OPT <- cv_results$k[argmin]
power_opt <- cv_results$power[argmin]
base_opt <- cv_results$base[argmin]

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    k = K_OPT,
                         power = power_opt,
                         base = base_opt)
invisible(basic_model$fit$state)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_B_Q1.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = K_OPT,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)

if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_B_Q1.rds", sep = ""))
}


#### Period B Q2 ####

training_data <- read.csv("training_sets/training_data_B_Q2.csv")
xtrain <- training_data %>% dplyr::select(!("y"|"ri"))
ytrain <- training_data$y


# Fold construction

fold_ids <- caret::createFolds(training_data$y, k = 5)

folds <- lapply(1:length(fold_ids),
                FUN = function(i){
                  training_data[-fold_ids[[i]],]})
antifolds <- lapply(1:length(fold_ids),
                    FUN = function(i){
                      training_data[fold_ids[[i]],]})

## Check training representation
# lapply(folds, function(x) with(x, table(ri,y)))


k_vals = c(1, 2, 3)
power_vals = c(1.6, 1.8, 2)
base_vals = c(0.75, 0.85, 0.95)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         auc1=numeric(),
                         auc2=numeric(),
                         auc3=numeric(),
                         auc4=numeric(),
                         auc5=numeric())
cv_results[1:kl*pl*bl, ] <- 0
for (i in 1:length(k_vals)){
  k_val <- k_vals[i]
  for (j in 1:length(power_vals)){
    power_val <- power_vals[j]
    for (m in 1:length(base_vals)){
      base_val <- base_vals[m]
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      cat("idx=", idx, "\n")
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){
        model <- bart2(y ~ . - ri,
                       folds[[fold_no]],
                       test = antifolds[[fold_no]],
                       k = k_val,
                       power = power_val,
                       base = base_val,
                       n.trees = 200,
                       n.chains = 1L,
                       n.threads = 1L,
                       keepTrees = TRUE,
                       verbose = FALSE)
        antifold_x <- subset(antifolds[[fold_no]], select=-c(y, ri))
        antifold_y <- antifolds[[fold_no]]$y
        antifold_ri <- antifolds[[fold_no]]$ri
        cutoff <- get_threshold(model)
        auc <- get_sens_and_spec(model, antifold_x, antifold_y, antifold_ri, cutoff)$auc
        cv_results[idx, 3+fold_no] <- auc
        rm(model)
      }
    }
  }
}

cv_results <- cv_results %>%
  rowwise() %>%
  mutate(mean_err = mean(auc1,
                         auc2,
                         auc3,
                         auc4,
                         auc5))
argmin <- which.max(cv_results$mean_err)
K_OPT <- cv_results$k[argmin]
power_opt <- cv_results$power[argmin]
base_opt <- cv_results$base[argmin]

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    k = K_OPT,
                         power = power_opt,
                         base = base_opt)
invisible(basic_model$fit$state)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_B_Q2.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = K_OPT,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)

if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_B_Q2.rds", sep = ""))
}


#### Period B Q3 ####

training_data <- read.csv("training_sets/training_data_B_Q3.csv")
xtrain <- training_data %>% dplyr::select(!("y"|"ri"))
ytrain <- training_data$y


# Fold construction

fold_ids <- caret::createFolds(training_data$y, k = 5)

folds <- lapply(1:length(fold_ids),
                FUN = function(i){
                  training_data[-fold_ids[[i]],]})
antifolds <- lapply(1:length(fold_ids),
                    FUN = function(i){
                      training_data[fold_ids[[i]],]})

## Check training representation
# lapply(folds, function(x) with(x, table(ri,y)))


k_vals = c(1, 2, 3)
power_vals = c(1.6, 1.8, 2)
base_vals = c(0.75, 0.85, 0.95)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         auc1=numeric(),
                         auc2=numeric(),
                         auc3=numeric(),
                         auc4=numeric(),
                         auc5=numeric())
cv_results[1:kl*pl*bl, ] <- 0
for (i in 1:length(k_vals)){
  k_val <- k_vals[i]
  for (j in 1:length(power_vals)){
    power_val <- power_vals[j]
    for (m in 1:length(base_vals)){
      base_val <- base_vals[m]
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      cat("idx=", idx, "\n")
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){
        model <- bart2(y ~ . - ri,
                       folds[[fold_no]],
                       test = antifolds[[fold_no]],
                       k = k_val,
                       power = power_val,
                       base = base_val,
                       n.trees = 200,
                       n.chains = 1L,
                       n.threads = 1L,
                       keepTrees = TRUE,
                       verbose = FALSE)
        antifold_x <- subset(antifolds[[fold_no]], select=-c(y, ri))
        antifold_y <- antifolds[[fold_no]]$y
        antifold_ri <- antifolds[[fold_no]]$ri
        cutoff <- get_threshold(model)
        auc <- get_sens_and_spec(model, antifold_x, antifold_y, antifold_ri, cutoff)$auc
        cv_results[idx, 3+fold_no] <- auc
        rm(model)
      }
    }
  }
}

cv_results <- cv_results %>%
  rowwise() %>%
  mutate(mean_err = mean(auc1,
                         auc2,
                         auc3,
                         auc4,
                         auc5))
argmin <- which.max(cv_results$mean_err)
K_OPT <- cv_results$k[argmin]
power_opt <- cv_results$power[argmin]
base_opt <- cv_results$base[argmin]

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    k = K_OPT,
                         power = power_opt,
                         base = base_opt)
invisible(basic_model$fit$state)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_B_Q3.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = K_OPT,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)

if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_B_Q3.rds", sep = ""))
}


#### Period B Q4 ####

training_data <- read.csv("training_sets/training_data_B_Q4.csv")
xtrain <- training_data %>% dplyr::select(!("y"|"ri"))
ytrain <- training_data$y


# Fold construction

fold_ids <- caret::createFolds(training_data$y, k = 5)

folds <- lapply(1:length(fold_ids),
                FUN = function(i){
                  training_data[-fold_ids[[i]],]})
antifolds <- lapply(1:length(fold_ids),
                    FUN = function(i){
                      training_data[fold_ids[[i]],]})

## Check training representation
# lapply(folds, function(x) with(x, table(ri,y)))


k_vals = c(1, 2, 3)
power_vals = c(1.6, 1.8, 2)
base_vals = c(0.75, 0.85, 0.95)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         auc1=numeric(),
                         auc2=numeric(),
                         auc3=numeric(),
                         auc4=numeric(),
                         auc5=numeric())
cv_results[1:kl*pl*bl, ] <- 0
for (i in 1:length(k_vals)){
  k_val <- k_vals[i]
  for (j in 1:length(power_vals)){
    power_val <- power_vals[j]
    for (m in 1:length(base_vals)){
      base_val <- base_vals[m]
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      cat("idx=", idx, "\n")
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){
        model <- bart2(y ~ . - ri,
                       folds[[fold_no]],
                       test = antifolds[[fold_no]],
                       k = k_val,
                       power = power_val,
                       base = base_val,
                       n.trees = 200,
                       n.chains = 1L,
                       n.threads = 1L,
                       keepTrees = TRUE,
                       verbose = FALSE)
        antifold_x <- subset(antifolds[[fold_no]], select=-c(y, ri))
        antifold_y <- antifolds[[fold_no]]$y
        antifold_ri <- antifolds[[fold_no]]$ri
        cutoff <- get_threshold(model)
        auc <- get_sens_and_spec(model, antifold_x, antifold_y, antifold_ri, cutoff)$auc
        cv_results[idx, 3+fold_no] <- auc
        rm(model)
      }
    }
  }
}

cv_results <- cv_results %>%
  rowwise() %>%
  mutate(mean_err = mean(auc1,
                         auc2,
                         auc3,
                         auc4,
                         auc5))
argmin <- which.max(cv_results$mean_err)
K_OPT <- cv_results$k[argmin]
power_opt <- cv_results$power[argmin]
base_opt <- cv_results$base[argmin]

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    k = K_OPT,
                         power = power_opt,
                         base = base_opt)
invisible(basic_model$fit$state)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_B_Q4.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = K_OPT,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)

if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_B_Q4.rds", sep = ""))
}