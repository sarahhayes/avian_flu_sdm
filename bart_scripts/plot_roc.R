# In this script we load in the outputs from fit_avian_flu_model.R.

library(embarcadero)
library(ggplot2)
library(raster)
library(RColorBrewer)
library(terra)
set.seed(12345)

# Function for directly getting optimal cutoff
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

################################################################################
# Lookup table for renaming landcover layers
lc_lookup <- pairlist("lc_1" = "Water_bodies",
                      "lc_2" = "Evergreen_Needleleaf_Forests",
                      "lc_3" = "Evergreen_Broadleaf_Forests",
                      "lc_4" = "Deciduous_Needleleaf_Forests",
                      "lc_5" = "Deciduous_Broadleaf_Forests",
                      "lc_6" = "Mixed_Forests",
                      "lc_7" = "Closed_Shrublands",
                      "lc_8" = "Open_Shrublands",
                      "lc_9" = "Woody_Savannas",
                      "lc_10" = "Savannas",
                      "lc_11" = "Grasslands",
                      "lc_12" = "Permanent_Wetlands",
                      "lc_13" = "Croplands",
                      "lc_14" = "Urband_and_Built-up_Lands",
                      "lc_15" = "Cropland/Natural_Vegetation_Mosaics",
                      "lc_16" = "Non-Vegetated_Lands",
                      "lc_17" = "Unclassified")

################################################################################
# Make AUC plot and getting numbers for table

# Load models:
load("output/fitted-BART-models/sdm_Q1.rds")
q1_sdm <- sdm
load("output/fitted-BART-models/sdm_Q2.rds")
q2_sdm <- sdm
load("output/fitted-BART-models/sdm_Q3.rds")
q3_sdm <- sdm
load("output/fitted-BART-models/sdm_Q4.rds")
q4_sdm <- sdm
rm(sdm)

load(file = "q1_train_test_data.rds")
q1_xtrain <- xtrain
q1_ytrain <- ytrain
q1_xtest <- xtest
q1_ytest <- ytest
load(file = "q2_train_test_data.rds")
q2_xtrain <- xtrain
q2_ytrain <- ytrain
q2_xtest <- xtest
q2_ytest <- ytest
load(file = "q3_train_test_data.rds")
q3_xtrain <- xtrain
q3_ytrain <- ytrain
q3_xtest <- xtest
q3_ytest <- ytest
load(file = "q4_train_test_data.rds")
q4_xtrain <- xtrain
q4_ytrain <- ytrain
q4_xtest <- xtest
q4_ytest <- ytest
rm(xtrain, ytrain, xtest, ytest)

get_sens_and_spec <- function(sdm, xtest, ytest){
  pred <- xtest[which(complete.cases(xtest)), ] %>%
    stats::predict(object=sdm) %>%
    pnorm() %>%
    colMeans() %>%
    prediction(labels = ytest[which(complete.cases(xtest))])
  perf <-  performance(pred, measure = "sens", x.measure = "spec")
  tss_list <- (perf@x.values[[1]] + perf@y.values[[1]] - 1)
  tss_df <- data.frame(alpha=perf@alpha.values[[1]],tss=tss_list)
  cutoff <- min(tss_df$alpha[which(tss_df$tss==max(tss_df$tss))])
  sens <- perf@x.values[[1]][which(perf@alpha.values[[1]]==cutoff)]
  spec <- perf@y.values[[1]][which(perf@alpha.values[[1]]==cutoff)]
  tss <- tss_df[which(tss_df$alpha==cutoff),'tss']
  return(pairlist("pred"=pred,
                  "perf"=perf,
                  "sens"=sens,
                  "spec"=spec,
                  "tss"=tss,
                  "cutoff"=cutoff))
}

# Try getting sensitivity and specificity:
q1_cutoff <- get_threshold(q1_sdm)
q2_cutoff <- get_threshold(q2_sdm)
q3_cutoff <- get_threshold(q3_sdm)
q4_cutoff <- get_threshold(q4_sdm)

q1_sens_spec <- get_sens_and_spec(q1_sdm, q1_xtest, q1_ytest)
q2_sens_spec <- get_sens_and_spec(q2_sdm, q2_xtest, q2_ytest)
q3_sens_spec <- get_sens_and_spec(q3_sdm, q3_xtest, q3_ytest)
q4_sens_spec <- get_sens_and_spec(q4_sdm, q4_xtest, q4_ytest)

cat("Q1 sensitivity =",
    q1_sens_spec$sens,
    ", specificity =",
    q1_sens_spec$spec,
    ", TSS =",
    q1_sens_spec$tss
    )
cat("Q2 sensitivity =",
    q2_sens_spec$sens,
    ", specificity =",
    q2_sens_spec$spec,
    ", TSS =",
    q2_sens_spec$tss
)
cat("Q3 sensitivity =",
    q3_sens_spec$sens,
    ", specificity =",
    q3_sens_spec$spec,
    ", TSS =",
    q3_sens_spec$tss
)
cat("Q4 sensitivity =",
    q4_sens_spec$sens,
    ", specificity =",
    q4_sens_spec$spec,
    ", TSS =",
    q4_sens_spec$tss
)

q1_x <- performance(q1_sens_spec$pred, "tpr", "fpr")
q2_x <- performance(q2_sens_spec$pred, "tpr", "fpr")
q3_x <- performance(q3_sens_spec$pred, "tpr", "fpr")
q4_x <- performance(q4_sens_spec$pred, "tpr", "fpr")
q1_rocdf <- data.frame(fpr=q1_x@x.values[[1]],
                    tpr=q1_x@y.values[[1]])
q2_rocdf <- data.frame(fpr=q2_x@x.values[[1]],
                       tpr=q2_x@y.values[[1]])
q3_rocdf <- data.frame(fpr=q3_x@x.values[[1]],
                       tpr=q3_x@y.values[[1]])
q4_rocdf <- data.frame(fpr=q4_x@x.values[[1]],
                       tpr=q4_x@y.values[[1]])
pal <- c("#2271B2",
         "#F748A5",
         "#359B73",
         "#e69f00")
roc_plot <- ggplot() + 
  geom_line(aes(x=q1_rocdf$fpr,y=q1_rocdf$tpr, colour = "Q1")) + 
  geom_line(aes(x=q2_rocdf$fpr,y=q2_rocdf$tpr, colour = "Q2")) + 
  geom_line(aes(x=q3_rocdf$fpr,y=q3_rocdf$tpr, colour = "Q3")) + 
  geom_line(aes(x=q4_rocdf$fpr,y=q4_rocdf$tpr, colour = "Q4")) + 
  scale_colour_manual("",
                      breaks = c("Q1", "Q2", "Q3", "Q4"),
                      values = pal) +
  ggtitle('Receiver-operator curve') + 
  xlab('False positive rate') + 
  ylab('True positive rate') + 
  geom_abline(intercept=0,slope=1,col='black')

ggsave("plots/RO_curve.png", plot = roc_plot, width = 6.5, height = 4.5)
