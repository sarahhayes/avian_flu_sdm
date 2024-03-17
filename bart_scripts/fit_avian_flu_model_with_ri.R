 # In this script we train a BART model to classify sites for avian flu
# presence/absence based on environmental and species abundance factors.
# This script includes cross-validation of BART parameters.

SAVE_FITS <- FALSE
SAVE_PLOTS <- FALSE

# Set path to folder containing data, and where output will be stored
PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(dplyr)
library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

#### EMBARCADERO FUNCTIONS ADJUSTED TO ALLOW CUSTOM BART PARAMETERS: ####

bart.flex <- function(x.data, y.data, ri.data = NULL,
                      y.name = NULL, ri.name = NULL,
                      n.trees = 200,
                      k = 2.0, power = 2.0, base = 0.95) {
  
  if(is.null(ri.data)) {
    train <- cbind(y.data, x.data) 
    if(!is.null(y.name)) {colnames(train)[1] <- y.name}
    train <- na.omit(train)
    model <- bart(y.train = train[,1], 
                  x.train = train[,2:ncol(train)], 
                  ntree = n.trees, keeptrees=TRUE)
  } else { 
    train <- cbind(y.data, x.data, ri.data) 
    if(!is.null(y.name)) {colnames(train)[1] <- y.name}
    if(!is.null(ri.name)) {colnames(train)[ncol(train)] <- ri.name}
    f <- as.formula(paste(paste(colnames(train)[1],paste(colnames(train)[2:(ncol(train)-1)], 
                                                         collapse=' + '), sep = ' ~ '), 
                          colnames(train)[ncol(train)], sep=' - '))
    
    train <- na.omit(train) 
    model <- rbart_vi(f, 
                      group.by=train[,ncol(train)],
                      data=train,
                      n.samples=1000,
                      n.burn=100,
                      n.chains=1,
                      n.threads=1,
                      n.trees = n.trees,
                      keepTrees = TRUE,
                      k = k,
                      power = power,
                      base = base) 
  }
  return(model)
}

varimp.diag <- function(x.data, y.data, ri.data=NULL, iter=50, quiet=FALSE,
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
                               ri.data = ri.data,
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
                                   ri.data = ri.data,
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

variable.step <- function(x.data, y.data, ri.data=NULL, n.trees=10, iter=50, quiet=FALSE,
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
                               ri.data = ri.data,
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
                                   ri.data = ri.data,
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

bart.step <- function(x.data, y.data, ri.data=NULL,
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
                               ri.data = ri.data,
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
  if(full==TRUE){varimp.diag(x.data, y.data, ri.data, iter=iter.plot, quiet=quiet2)}
  vs <- variable.step(x.data,
                      y.data,
                      ri.data,
                      n.trees=tree.step,
                      iter=iter.step,
                      quiet=quiet2)
  
  invisible(best.model <- bart.flex(x.data = x.data[,vs], y.data = y.data, 
                                    ri.data = ri.data, n.trees=200,
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



#### Period A Q1 ####

training_data <- read.csv("training_sets/training_data_A_Q1.csv")

#### Sketch of fold construction ####
fold_order <- sample(1:nrow(training_data), nrow(training_data))
reordered_training_data <- training_data[fold_order,]
split_pts <- c(0,
               as.integer(floor(nrow(training_data)/5)),
               as.integer(floor(2*nrow(training_data)/5)),
               as.integer(floor(3*nrow(training_data)/5)),
               as.integer(floor(4*nrow(training_data)/5)),
               nrow(training_data))
folds <- lapply(1:5,
                FUN = function(i){
                  reordered_training_data[(split_pts[i]+1):split_pts[i+1],]})
antifolds <- lapply(1:length(folds),
                    FUN = function(i){
                      bind_rows(folds[c(-i)])
                    })


k_vals = c(1, 2, 3)
power_vals = c(1.5,1.7,1.9)
base_vals = c(0.7,0.8,0.9)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         err1=numeric(),
                         err2=numeric(),
                         err3=numeric(),
                         err4=numeric(),
                         err5=numeric())
cv_results[1:kl*pl*bl, ] <- 0
for (i in 1:length(k_vals)){
  k_val <- eval(k_vals[i])
  for (j in 1:length(power_vals)){
    power_val <- eval(power_vals[j])
    for (m in 1:length(base_vals)){
      base_val <- eval(base_vals[m])
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      cat("idx=", idx, "\n")
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){
        model <- rbart_vi(y ~ . - ri,
                          folds[[fold_no]],
                          group.by = ri,
                          group.by.test = ri,
                          test = antifolds[[fold_no]],
                          k = k_val,
                          power = power_val,
                          base = base_val,
                          n.chains = 1L,
                          n.threads = 1L,
                          keepTrees = TRUE,
                          verbose = FALSE)
        mean_err <- model %>% residuals %>% abs %>% mean
        cv_results[idx, 3+fold_no] <- mean_err
      }
    }
  }
}

cv_results <- cv_results %>%
              rowwise() %>%
              mutate(mean_err = mean(err1,
                                     err2,
                                     err3,
                                     err4,
                                     err5))
argmin <- which.min(cv_results$mean_err)
k_opt <- cv_results$k[argmin]
power_opt <- cv_results$power[argmin]
base_opt <- cv_results$base[argmin]

########

xtrain <- training_data %>% dplyr::select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_A_Q1.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                    training_data,
                    group.by = ri,
                    group.by.test = ri,
                    test = test_data,
                    k = k_opt,
                    power = power_opt,
                    base = base_opt,
                    n.chains = 1L,
                    n.threads = 1L,
                    keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_A_Q1.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 k = k_opt,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q1.rds", sep = ""))
}


#### Period A Q2 ####

training_data <- read.csv("training_sets/training_data_A_Q2.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_A_Q2.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_A_Q2.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 k = k_opt,
                 power = power_opt,
                 base = base_opt,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q2.rds", sep = ""))
}

#### Period A Q3 ####

training_data <- read.csv("training_sets/training_data_A_Q3.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_A_Q3.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_A_Q3.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q3.rds", sep = ""))
}

#### Period A Q4 ####

training_data <- read.csv("training_sets/training_data_A_Q4.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_A_Q4.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_A_Q4.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q4.rds", sep = ""))
}

#### Period B Q1 ####

training_data <- read.csv("training_sets/training_data_B_Q1.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_B_Q1.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_B_Q1.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_B_Q1.rds", sep = ""))
}



#### Period B Q2 ####

training_data <- read.csv("training_sets/training_data_B_Q2.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_B_Q2.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_B_Q2.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_B_Q2.rds", sep = ""))
}


#### Period B Q3 ####

training_data <- read.csv("training_sets/training_data_B_Q3.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_B_Q3.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_B_Q3.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_B_Q3.rds", sep = ""))
}


#### Period B Q4 ####

training_data <- read.csv("training_sets/training_data_B_Q4.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_B_Q4.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_B_Q4.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_B_Q4.rds", sep = ""))
}