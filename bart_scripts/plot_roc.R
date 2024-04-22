# In this script we load in the outputs from fit_avian_flu_model.R.

# Optional command line arguments, must be passed as strings:
args <- commandArgs(trailingOnly = T)
if (length(args)<4){
  INCLUDE_CROSSTERMS <- "no-crossterms" # Set to "no-crossterms" to do model without crossterms or "with-crossterms" to do model with crossterms
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
  PATH_TO_OUTPUTS <- "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/"
}else{
  PATH_TO_OUTPUTS <- args[1]
}

B_TEST_DATA_AVAILABLE <- FALSE # Set to TRUE if we have test data for dataset B, otherwise this should be FALSE

library(embarcadero)
library(ggplot2)
library(raster)
library(RColorBrewer)
library(terra)
set.seed(12345)

################################################################################
# Make AUC plot and get numbers for table for dataset A

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
                    "_metrics_A_Q",
                    idx,
                    ".rds",
                    sep = ""))
  metrics <- load(file = paste(PATH_TO_OUTPUTS,
                               "fitted-BART-models-",
                               INCLUDE_CROSSTERMS,
                               "/",
                               CV_OR_RI,
                               "_metrics_A_Q",
                               idx,
                               ".rds",
                               sep = "")) %>% get
  cat("Dataset A Q",
      idx,
      "\nsensitivity =",
      metrics$sens,
      "\nspecificity =",
      metrics$spec,
      "\nTSS =",
      metrics$tss,
      "\nAUC=",
      metrics$auc,
      "\n\n"
  )
  
}

Q1_x <- performance(metrics_A_Q1$pred, "tpr", "fpr")
Q2_x <- performance(metrics_A_Q2$pred, "tpr", "fpr")
Q3_x <- performance(metrics_A_Q3$pred, "tpr", "fpr")
Q4_x <- performance(metrics_A_Q4$pred, "tpr", "fpr")
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
  geom_line(aes(x=Q1_rocdf$fpr,y=Q1_rocdf$tpr, colour = "A Q1")) + 
  geom_line(aes(x=Q2_rocdf$fpr,y=Q2_rocdf$tpr, colour = "A Q2")) + 
  geom_line(aes(x=Q3_rocdf$fpr,y=Q3_rocdf$tpr, colour = "A Q3")) + 
  geom_line(aes(x=Q4_rocdf$fpr,y=Q4_rocdf$tpr, colour = "A Q4")) + 
  scale_colour_manual("",
                      breaks = c("A Q1", "A Q2", "A Q3", "A Q4"),
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
if (B_TEST_DATA_AVAILABLE)
  {for (idx in 1:4){
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
    metrics <- load(file = paste(PATH_TO_OUTPUTS,
                                 "fitted-BART-models-",
                                 INCLUDE_CROSSTERMS,
                                 "/",
                                 CV_OR_RI,
                                 "_metrics_B_Q",
                                 idx,
                                 ".rds",
                                 sep = "")) %>% get
    cat("Dataset B Q",
        idx,
        "\nsensitivity =",
        metrics$sens,
        "\nspecificity =",
        metrics$spec,
        "\nTSS =",
        metrics$tss,
        "\nAUC=",
        metrics$auc,
        "\n\n"
    )
    
  }
  
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
               "_RO_curve_A.png",
               sep=""),
         plot = roc_plot,
         width = 6.5,
         height = 4.5)
}