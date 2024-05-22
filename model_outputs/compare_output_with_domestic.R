# 22nd May 2024
# Comparing the model results to the domestic cases. 

rm(list = ls())
library(embarcadero)
library(raster)
library(terra)
library(viridis)
set.seed(12345)


PLOT_RAW_CASE_DATA <- FALSE 

# Below is adapted from Joe's code in plot_projections

## For outputs from model B 
PATH_TO_OUTPUTSB <- "model_outputs/preds_B/"

predsB <- list()

for (idx in 1:4){
  # Slightly hacky way to load in predictions and give it an arbitrary name.
  # This is needed because readRDS appears to be deprecated and the load
  # function brings in the original variable name. If you do x<-load(x.rds) then
  # x just gives you the variable name; when we pipe with get() we can get the
  # value of the thing with that variable name, provided it's been loaded in.
  load(file = paste(PATH_TO_OUTPUTSB,
                    "cv_predictions_B_Q",
                    idx,
                    ".rds",
                    sep = ""))
}


pred_layers_B_Q1_mean <- pred_layers_B_Q1$Mean
pred_layers_B_Q2_mean <- pred_layers_B_Q2$Mean
pred_layers_B_Q3_mean <- pred_layers_B_Q3$Mean
pred_layers_B_Q4_mean <- pred_layers_B_Q4$Mean

## we want these as spatRasters for use later
crs <- "epsg:3035"

pred_layers_B_Q1_mean <- terra::rast(pred_layers_B_Q1_mean)
set.crs(pred_layers_B_Q1_mean, value = crs)
pred_layers_B_Q2_mean <- terra::rast(pred_layers_B_Q2_mean)
set.crs(pred_layers_B_Q2_mean, value = crs)
pred_layers_B_Q3_mean <- terra::rast(pred_layers_B_Q3_mean)
set.crs(pred_layers_B_Q3_mean, value = crs)
pred_layers_B_Q4_mean <- terra::rast(pred_layers_B_Q4_mean)
set.crs(pred_layers_B_Q4_mean, value = crs)

# read in the domestic cases

dom_bq1 <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q1_pres_abs.tif")
dom_bq2 <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q2_pres_abs.tif")
dom_bq3 <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q3_pres_abs.tif")
dom_bq4 <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q4_pres_abs.tif")

# dom_bq1_counts <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q1_counts.tif")
# dom_bq2_counts <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q2_counts.tif")
# dom_bq3_counts <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q3_counts.tif")
# dom_bq4_counts <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q4_counts.tif")

png("plots/dom_over_preds_b1.png")
plot(pred_layers_B_Q1_mean)
plot(dom_bq1, add = T, col = "black", legend = F)
dev.off()

png("plots/dom_over_preds_b2.png")
plot(pred_layers_B_Q2_mean)
plot(dom_bq2, add = T, col = "black", legend = F)
dev.off()

png("plots/dom_over_preds_b3.png")
plot(pred_layers_B_Q3_mean)
plot(dom_bq3, add = T, col = "black", legend = F)
dev.off()

png("plots/dom_over_preds_b4.png")
plot(pred_layers_B_Q4_mean)
plot(dom_bq4, add = T, col = "black", legend = F)
dev.off()


# be good to produce multi-panel blot. Want to trim the white space around the map though
dev.off()
png("plots/dom_over_predsB_allq.png")
par(mfrow = c(2,2))
plot(pred_layers_B_Q1_mean, axes = F, box = F, mar = c(1,1,1,1), legend = F)
plot(dom_bq1, add = T, col = "black", legend = F)
plot(pred_layers_B_Q2_mean, axes = F, box = F, mar = c(1,1,1,1), legend = T)
plot(dom_bq2, add = T, col = "black", legend = F)
plot(pred_layers_B_Q3_mean, axes = F, box = F, mar = c(1,1,1,1), legend = F)
plot(dom_bq3, add = T, col = "black", legend = F)
plot(pred_layers_B_Q4_mean, axes = F, box = F, mar = c(1,1,1,1), legend = T)
plot(dom_bq4, add = T, col = "black",legend = F)
dev.off()

# look at trying to outline the area of the rasters that have values and overlaying this on the map.
# Try chat GPT to get the code for this??

### Look at producing some histograms

dom_probs_bq1 <- terra::mask(pred_layers_B_Q1_mean, dom_bq1)
dom_probs_bq2 <- terra::mask(pred_layers_B_Q2_mean, dom_bq2)
dom_probs_bq3 <- terra::mask(pred_layers_B_Q3_mean, dom_bq3)
dom_probs_bq4 <- terra::mask(pred_layers_B_Q4_mean, dom_bq4)

dev.off()
png("plots/hists_of_probs_domestic.png")
par(mfrow = c(2,2))
hist(values(dom_probs_bq1), col = "#6699CC", main = "B - Q1", xlab = "Predicted probability")
hist(values(dom_probs_bq2), col = "#339966", main = "B - Q2", xlab = "Predicted probability")
hist(values(dom_probs_bq3), col = "#FF9933", main = "B - Q3", xlab = "Predicted probability")
hist(values(dom_probs_bq4), col = "#99CCFF", main = "B - Q4", xlab = "Predicted probability")
dev.off()
       

##################################################################################################

# Repeat for A


## For outputs from model B 
PATH_TO_OUTPUTSA <- "model_outputs/preds_A/"

for (idx in 1:4){
 load(file = paste(PATH_TO_OUTPUTSA,
                    "cv_predictions_A_Q",
                    idx,
                    ".rds",
                    sep = ""))
}


pred_layers_A_Q1_mean <- pred_layers_A_Q1$Mean
pred_layers_A_Q2_mean <- pred_layers_A_Q2$Mean
pred_layers_A_Q3_mean <- pred_layers_A_Q3$Mean
pred_layers_A_Q4_mean <- pred_layers_A_Q4$Mean

pred_layers_A_Q1_mean <- terra::rast(pred_layers_A_Q1_mean)
set.crs(pred_layers_A_Q1_mean, value = crs)
pred_layers_A_Q2_mean <- terra::rast(pred_layers_A_Q2_mean)
set.crs(pred_layers_A_Q2_mean, value = crs)
pred_layers_A_Q3_mean <- terra::rast(pred_layers_A_Q3_mean)
set.crs(pred_layers_A_Q3_mean, value = crs)
pred_layers_A_Q4_mean <- terra::rast(pred_layers_A_Q4_mean)
set.crs(pred_layers_A_Q4_mean, value = crs)

# read in the domestic cases

dom_aq1 <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q1_pres_abs.tif")
dom_aq2 <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q2_pres_abs.tif")
dom_aq3 <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q3_pres_abs.tif")
dom_aq4 <- terra::rast("data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q4_pres_abs.tif")


png("plots/dom_over_preds_a1.png")
plot(pred_layers_A_Q1_mean)
plot(dom_aq1, add = T, col = "black", legend = F)
dev.off()

png("plots/dom_over_preds_a2.png")
plot(pred_layers_A_Q2_mean)
plot(dom_aq2, add = T, col = "black", legend = F)
dev.off()

png("plots/dom_over_preds_a3.png")
plot(pred_layers_A_Q3_mean)
plot(dom_aq3, add = T, col = "black", legend = F)
dev.off()

png("plots/dom_over_preds_a4.png")
plot(pred_layers_A_Q4_mean)
plot(dom_aq4, add = T, col = "black", legend = F)
dev.off()

dev.off()
png("plots/dom_over_predsA_allq.png")
par(mfrow = c(2,2))
plot(pred_layers_A_Q1_mean, axes = F, box = F, mar = c(1,1,1,1), legend = F)
plot(dom_aq1, add = T, col = "black", legend = F)
plot(pred_layers_A_Q2_mean, axes = F, box = F, mar = c(1,1,1,1), legend = T)
plot(dom_aq2, add = T, col = "black", legend = F)
plot(pred_layers_A_Q3_mean, axes = F, box = F, mar = c(1,1,1,1), legend = F)
plot(dom_aq3, add = T, col = "black", legend = F)
plot(pred_layers_A_Q4_mean, axes = F, box = F, mar = c(1,1,1,1), legend = T)
plot(dom_aq4, add = T, col = "black",legend = F)
dev.off()

## Look at producing some histograms

dom_probs_aq1 <- terra::mask(pred_layers_A_Q1_mean, dom_aq1)
dom_probs_aq2 <- terra::mask(pred_layers_A_Q2_mean, dom_aq2)
dom_probs_aq3 <- terra::mask(pred_layers_A_Q3_mean, dom_aq3)
dom_probs_aq4 <- terra::mask(pred_layers_A_Q4_mean, dom_aq4)

dev.off()
png("plots/hists_of_probs_domestic_A.png")
par(mfrow = c(2,2))
hist(values(dom_probs_aq1), col = "#6699CC", main = "A - Q1", xlab = "Predicted probability")
hist(values(dom_probs_aq2), col = "#339966", main = "A - Q2", xlab = "Predicted probability")
hist(values(dom_probs_aq3), col = "#FF9933", main = "A - Q3", xlab = "Predicted probability")
hist(values(dom_probs_aq4), col = "#99CCFF", main = "A - Q4", xlab = "Predicted probability")
dev.off()
