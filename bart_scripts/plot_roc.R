# In this script we plot risk projections using the outputs from
# fit_avian_flu_model.R.

SKIP_AQ3 <- TRUE # Set to TRUE to skip small period A quarter 3 dataset

B_TEST_DATA_AVAILABLE <- TRUE

# Optional command line arguments, must be passed as strings:
args <- commandArgs(trailingOnly = T)
if (length(args)<4){
  INCLUDE_CROSSTERMS <- "with-crossterms-multichain" # Set to "no-crossterms" to do model without crossterms or "with-crossterms" to do model with crossterms
}else{
  INCLUDE_CROSSTERMS <- args[4]
}
if (length(args)<3){
  CV_OR_RI <- "cv" # Set to "cv" to do crossvalidated model or "ri" to do crossvalidation + random intercept model
}else{
  CV_OR_RI <- args[3]
}
if (length(args)<1){
  # Set path to folder containing data, and where output will be stored
  PATH_TO_OUTPUTS <- ""
}else{
  PATH_TO_OUTPUTS <- args[1]
}

# Decide whether to use models with ecological season boundaries
ECO_SEASONS <- TRUE
if (ECO_SEASONS){
  PATH_TO_MODELS <- paste(PATH_TO_OUTPUTS,
                          "",
                          sep="")
}else{
  PATH_TO_MODELS <- PATH_TO_OUTPUTS
}

library(embarcadero)
library(raster)
library(terra)
library(viridis)
set.seed(12345)

preds <- list()

if (!SKIP_AQ3){
  plt_idx <- 1:4
}else{
  plt_idx <- c(1, 2, 4) %>% as.integer()
}

for (idx in plt_idx){
  # Slightly hacky way to load in predictions and give it an arbitrary name.
  # This is needed because readRDS appears to be deprecated and the load
  # function brings in the original variable name. If you do x<-load(x.rds) then
  # x just gives you the variable name; when we pipe with get() we can get the
  # value of the thing with that variable name, provided it's been loaded in.
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_metrics_A_Q",
                    idx,
                    ".rds",
                    sep = ""))
  
}

{
  cat("Dataset A Q1, test metrics:\n sens =",
      signif(metrics_A_Q1$sens, digits=2),
      " (",
      signif(metric_cis_A_Q1$sens, digits=2),
      "),\n spec =",
      signif(metrics_A_Q1$spec, digits=2),
      " (",
      signif(metric_cis_A_Q1$spec, digits=2),
      "),\n AUC =",
      signif(metrics_A_Q1$auc, digits=2),
      " (",
      signif(metric_cis_A_Q1$auc, digits=2),
      ")\n")
  cat("Dataset A Q2, test metrics:\n sens =",
      signif(metrics_A_Q2$sens, digits=2),
      " (",
      signif(metric_cis_A_Q2$sens, digits=2),
      "),\n spec =",
      signif(metrics_A_Q2$spec, digits=2),
      " (",
      signif(metric_cis_A_Q2$spec, digits=2),
      "),\n AUC =",
      signif(metrics_A_Q2$auc, digits=2),
      " (",
      signif(metric_cis_A_Q2$auc, digits=2),
      ")\n")
  if (!SKIP_AQ3){
    cat("Dataset A Q3, test metrics:\n sens =",
        signif(metrics_A_Q3$sens, digits=2),
        " (",
        signif(metric_cis_A_Q3$sens, digits=2),
        "),\n spec =",
        signif(metrics_A_Q3$spec, digits=2),
        " (",
        signif(metric_cis_A_Q3$spec, digits=2),
        "),\n AUC =",
        signif(metrics_A_Q3$auc, digits=2),
        " (",
        signif(metric_cis_A_Q3$auc, digits=2),
        ")\n")
  }
  cat("Dataset A Q4, test metrics:\n sens =",
      signif(metrics_A_Q4$sens, digits=2),
      " (",
      signif(metric_cis_A_Q4$sens, digits=2),
      "),\n spec =",
      signif(metrics_A_Q4$spec, digits=2),
      " (",
      signif(metric_cis_A_Q4$spec, digits=2),
      "),\n AUC =",
      signif(metrics_A_Q4$auc, digits=2),
      " (",
      signif(metric_cis_A_Q4$auc, digits=2),
      ")\n")
}

Q1_x <- performance(metrics_A_Q1$pred, "tpr", "fpr")
Q2_x <- performance(metrics_A_Q2$pred, "tpr", "fpr")
Q4_x <- performance(metrics_A_Q4$pred, "tpr", "fpr")
Q1_rocdf <- data.frame(fpr=Q1_x@x.values[[1]],
                       tpr=Q1_x@y.values[[1]])
Q2_rocdf <- data.frame(fpr=Q2_x@x.values[[1]],
                       tpr=Q2_x@y.values[[1]])
Q4_rocdf <- data.frame(fpr=Q4_x@x.values[[1]],
                       tpr=Q4_x@y.values[[1]])
pal <- c("#2271B2",
         "#F748A5",
         "#359B73",
         "#e69f00")
roc_plot <- ggplot() + 
  geom_line(aes(x=Q1_rocdf$fpr,y=Q1_rocdf$tpr, colour = "A Q1")) + 
  geom_line(aes(x=Q2_rocdf$fpr,y=Q2_rocdf$tpr, colour = "A Q2")) + 
  geom_line(aes(x=Q4_rocdf$fpr,y=Q4_rocdf$tpr, colour = "A Q4")) + 
  scale_colour_manual("",
                      breaks = c("A Q1", "A Q2", "A Q4"),
                      values = pal) +
  ggtitle('Receiver-operating characteristic') + 
  xlab('False positive rate') + 
  ylab('True positive rate') + 
  geom_abline(intercept=0,slope=1,col='black')

ggsave(paste("plots/",
       INCLUDE_CROSSTERMS,
       "_",
       CV_OR_RI,
       "_RO_curve_A.png",
       sep=""),
       plot = roc_plot,
       width = 6.5,
       height = 4.5)

################################################################################
# Make AUC plot and get numbers for table for dataset B
if (B_TEST_DATA_AVAILABLE){
  for (idx in 1:4){
    # Slightly hacky way to load in predictions and give it an arbitrary name.
    # This is needed because readRDS appears to be deprecated and the load
    # function brings in the original variable name. If you do x<-load(x.rds) then
    # x just gives you the variable name; when we pipe with get() we can get the
    # value of the thing with that variable name, provided it's been loaded in.
    load(file = paste(PATH_TO_OUTPUTS,
                      "fitted-BART-models-",
                      INCLUDE_CROSSTERMS,
                      "/",
                      CV_OR_RI,
                      "_metrics_B_Q",
                      idx,
                      ".rds",
                      sep = ""))
    
  }
  
  cat("Dataset B Q1, test metrics:\n sens =",
      signif(metrics_B_Q1$sens, digits=2),
      " (",
      signif(metric_cis_B_Q1$sens, digits=2),
      "),\n spec =",
      signif(metrics_B_Q1$spec, digits=2),
      " (",
      signif(metric_cis_B_Q1$spec, digits=2),
      "),\n AUC =",
      signif(metrics_B_Q1$auc, digits=2),
      " (",
      signif(metric_cis_B_Q1$auc, digits=2),
      ")\n")
  cat("Dataset B Q2, test metrics:\n sens =",
      signif(metrics_B_Q2$sens, digits=2),
      " (",
      signif(metric_cis_B_Q2$sens, digits=2),
      "),\n spec =",
      signif(metrics_B_Q2$spec, digits=2),
      " (",
      signif(metric_cis_B_Q2$spec, digits=2),
      "),\n AUC =",
      signif(metrics_B_Q2$auc, digits=2),
      " (",
      signif(metric_cis_B_Q2$auc, digits=2),
      ")\n")
  cat("Dataset B Q3, test metrics:\n sens =",
      signif(metrics_B_Q3$sens, digits=2),
      " (",
      signif(metric_cis_B_Q3$sens, digits=2),
      "),\n spec =",
      signif(metrics_B_Q3$spec, digits=2),
      " (",
      signif(metric_cis_B_Q3$spec, digits=2),
      "),\n AUC =",
      signif(metrics_B_Q3$auc, digits=2),
      " (",
      signif(metric_cis_B_Q3$auc, digits=2),
      ")\n")
  cat("Dataset B Q4, test metrics:\n sens =",
      signif(metrics_B_Q4$sens, digits=2),
      " (",
      signif(metric_cis_B_Q4$sens, digits=2),
      "),\n spec =",
      signif(metrics_B_Q4$spec, digits=2),
      " (",
      signif(metric_cis_B_Q4$spec, digits=2),
      "),\n AUC =",
      signif(metrics_B_Q4$auc, digits=2),
      " (",
      signif(metric_cis_B_Q4$auc, digits=2),
      ")\n")
  
  Q1_x <- performance(metrics_B_Q1$pred, "tpr", "fpr")
  Q2_x <- performance(metrics_B_Q2$pred, "tpr", "fpr")
  Q3_x <- performance(metrics_B_Q3$pred, "tpr", "fpr")
  Q4_x <- performance(metrics_B_Q4$pred, "tpr", "fpr")
  Q1_rocdf <- data.frame(fpr=Q1_x@x.values[[1]],
                         tpr=Q1_x@y.values[[1]])
  Q2_rocdf <- data.frame(fpr=Q2_x@x.values[[1]],
                         tpr=Q2_x@y.values[[1]])
  Q3_rocdf <- data.frame(fpr=Q3_x@x.values[[1]],
                         tpr=Q3_x@y.values[[1]])
  Q4_rocdf <- data.frame(fpr=Q4_x@x.values[[1]],
                         tpr=Q4_x@y.values[[1]])
  pal <- c("#2271B2",
           "#F748A5",
           "#359B73",
           "#e69f00")
  roc_plot <- ggplot() + 
    geom_line(aes(x=Q1_rocdf$fpr,y=Q1_rocdf$tpr, colour = "B Q1")) + 
    geom_line(aes(x=Q2_rocdf$fpr,y=Q2_rocdf$tpr, colour = "B Q2")) + 
    geom_line(aes(x=Q3_rocdf$fpr,y=Q3_rocdf$tpr, colour = "B Q3")) + 
    geom_line(aes(x=Q4_rocdf$fpr,y=Q4_rocdf$tpr, colour = "B Q4")) + 
    scale_colour_manual("",
                        breaks = c("B Q1", "B Q2", "B Q3", "B Q4"),
                        values = pal) +
    ggtitle('Receiver-operating characteristic') + 
    xlab('False positive rate') + 
    ylab('True positive rate') + 
    geom_abline(intercept=0,slope=1,col='black')
  
  ggsave(paste("plots/",
               INCLUDE_CROSSTERMS,
               "_",
               CV_OR_RI,
               "_RO_curve_B.png",
               sep=""),
         plot = roc_plot,
         width = 6.5,
         height = 4.5)
}
