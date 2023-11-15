# In this script we define a function for doing both cross validating models
# with variable selection


# ADAPTED FROM EMBARCADERO:

variable.step.cv <- function(x.data,
                          y.data,
                          n.trees=10,
                          iter=50,
                          k = 2.0,
                          power = 2.0,
                          base = 0.95,
                          quiet=FALSE) {
  
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
  
  quietly(model.0 <- bart(x.train = x.data, y.train = y.data,
                          k = k,
                          power = power,
                          base = base,
                          ntree = 200,
                          keeptrees = TRUE))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  
  # dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  # 
  # if(length(dropnames) > 0) {
  #   message("Some of your variables have been automatically dropped by dbarts.")
  #   message("(This could be because they're characters, homogenous, etc.)")
  #   message("It is strongly recommended that you remove these from the raw data:")
  #   message(paste(dropnames,collapse = ' '), ' \n')
  # }
  # 
  # x.data %>% dplyr::select(-any_of(dropnames)) -> x.data  
  
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
      quietly(model.j <- bart(x.train = x.data[,varnums], y.train = y.data, 
                              k = k,
                              power = power,
                              base = base,
                              ntree = n.trees,
                              keeptrees = TRUE))
      
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

bart.step.cv <- function(x.data, y.data, ri.data=NULL,
                      k = 2.0,
                      power = 2.0,
                      base = 0.95,
                      iter.step=100, tree.step=10,
                      iter.plot=100,
                      full=FALSE,
                      quiet=FALSE) {
  
  ###############
  
  # auto-drops 
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY
  
  quietly(model.0 <- bart(x.train = x.data, y.train = y.data,
                          k = k,
                          power = power,
                          base = base,
                          ntree = 200,
                          keeptrees = TRUE))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  
  # dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  # 
  # if(length(dropnames)==0) {} else{
  #   message("Some of your variables have been automatically dropped by dbarts.")
  #   message("(This could be because they're characters, homogenous, etc.)")
  #   message("It is strongly recommended that you remove these from the raw data:")
  #   message(paste(dropnames,collapse = ' '), ' \n')
  # }
  # 
  # x.data %>% dplyr::select(-dropnames) -> x.data  
  
  ###############
  
  quiet2 <- quiet
  if(full==TRUE){varimp.diag(x.data, y.data, ri.data, iter=iter.plot, quiet=quiet2)}
  vs <- variable.step.cv(x.data, y.data,
                         k = k,
                         power = power,
                         base = base, n.trees=tree.step, iter=iter.step, quiet=quiet2)
  
  invisible(best.model <- bart(x.train = x.data[,vs], y.train = y.data, 
                               k = k,
                               power = power,
                               base = base,
                               ntree=200,
                               keeptrees = TRUE))
  if(full==TRUE){varimp(best.model, plots=TRUE)}
  if(full==TRUE) {p <- summary(best.model, plots=TRUE)
  print(p)} else 
  {p <- summary(best.model, plots=FALSE)
  print(p)}
  invisible(best.model)
}  


cv_with_vs <- function(xtrain, ytrain, 
                       k = 2.0,
                       power = 2.0,
                       base = 0.95, n.reps, n_folds){
  loss_vect <- zeros(n.reps)
  for (n in 1:n.reps){
    folds <- split(sample(1:nrow(xtrain)), 1:n_folds)
    for (i in 1:n_folds){
      xfold <- xtrain[-unlist(folds[i]), ]
      yfold <- ytrain[-unlist(folds[i])]
      model <- bart.step.cv(xfold, yfold, 
                            k,
                            power,
                            base)
      
      pred.p <- colMeans(pnorm(model$yhat.train))[ytrain==1]
      pred.a <- colMeans(pnorm(model$yhat.train))[ytrain==0]
      
      pred.c <- c(pred.p, pred.a)
      true.c <- c(rep(1,length(pred.p)), rep(0,length(pred.a)))
      rmse <- Metrics::rmse(true.c,pred.c)
      loss_vect[i] <- loss_vect[i] + (1/n_folds) * rmse
    }
  }
  return(loss_vect)
}

retune.vs <- function(xtrain, ytrain, reps = 10) {
  
  # auto-drops 
  
  #quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
  #                             ri.data = ri.data,
  #                             n.trees = 200))
  
  #dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(model.0$fit$data@x,"drop"))==FALSE)))]
  
  #  if(length(dropnames)==0) {} else{
  #    message("Some of your variables have been automatically dropped by dbarts.")
  #    message("(This could be because they're characters, homogenous, etc.)")
  #    message("It is strongly recommended that you remove these from the raw data:")
  #    message(paste(dropnames,collapse = ' '), ' \n')
  #   }
  
  # x.data %>% select(-dropnames) -> x.data  
  
  ####
  
  k_vals <- c(1,2,3)
  power_vals <- c(1.5, 1.6, 1.7, 1.8, 1.9, 2)
  base_vals <- c(0.75, 0.8, 0.85, 0.9, 0.95)
  
  loss_array <- array(0, dim=c(length(k_vals), length(power_vals), length(base_vals)))
  
  for (i in 1:length(k_vals)){
    for (j in 1:length(power_vals)){
      for (k in 1:length(base_vals)){
        x <- cv_with_vs(xtrain,
                        ytrain,
                        n_folds = 5, n.reps = reps,
                        k = k_vals[i],
                        power = power_vals[j],
                        base = base_vals[k])
        loss_array[i] <- x
      }
    }
  }
  
  priors <- which(loss_array==min(loss_array), arr.ind = TRUE)
  
  model <- bart.step.cv(object$fit$data@x,
                        object$fit$data@y,
                        k = c(1,2,3)[priors[1]],
                        power = c(1.5, 1.6, 1.7, 1.8, 1.9, 2)[priors[2]],
                        base = c(0.75, 0.8, 0.85, 0.9, 0.95)[priors[3]],)
  return(model)
  
}

