# In this script we load in the outputs from fit_avian_flu_model.R.

library(embarcadero)
library(raster)
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
# Example: load in Q1 stuff

# If you want to get covariates as a raster, run this with the appropriate address:
# covstack <- raster::stack("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/q1_covs.tif")

# Load fitted and optimised model:
load("output/fitted-BART-models/sdm_Q1.rds")
xtrain <- sdm$fit$data@x
colnames(xtrain) <- lapply(colnames(xtrain),
                                  FUN = function(n){
                                    if (n %in% names(lc_lookup)){
                                      n <- lc_lookup[[n]]}else{
                                        n <- n}
                                      })
ytrain <- sdm$fit$data@y
xtest <- sdm$fit$data@x.test

# Try a few bits of analysis to check everything's loaded properly:
summary(sdm)

# Try getting sensitivity and specificity:
cutoff <- get_threshold(sdm)
perf <- sdm$yhat.train %>% pnorm() %>%
  colMeans() %>%
  prediction(labels = sdm$fit$data@y) %>%
  performance(measure = "sens", x.measure = "spec")
sens <- perf@x.values[[1]][which(perf@alpha.values[[1]]==cutoff)]
spec <- perf@y.values[[1]][which(perf@alpha.values[[1]]==cutoff)]

# Load in the prediction map:
load("output/fitted-BART-models/prediction_Q1.rds")
plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Mean prediction, Q1')
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5th percentile, Q1')
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5th percentile, Q1')
