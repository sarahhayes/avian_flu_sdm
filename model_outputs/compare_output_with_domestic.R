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

# png("plots/dom_over_preds_b1.png")
# plot(pred_layers_B_Q1_mean)
# plot(dom_bq1, add = T, col = "black", legend = F)
# dev.off()
# 
# png("plots/dom_over_preds_b2.png")
# plot(pred_layers_B_Q2_mean)
# plot(dom_bq2, add = T, col = "black", legend = F)
# dev.off()
# 
# png("plots/dom_over_preds_b3.png")
# plot(pred_layers_B_Q3_mean)
# plot(dom_bq3, add = T, col = "black", legend = F)
# dev.off()
# 
# png("plots/dom_over_preds_b4.png")
# plot(pred_layers_B_Q4_mean)
# plot(dom_bq4, add = T, col = "black", legend = F)
# dev.off()


# be good to produce multi-panel blot. Want to trim the white space around the map though
# dev.off()
# png("plots/dom_over_predsB_allq.png")
# par(mfrow = c(2,2))
# plot(pred_layers_B_Q1_mean, axes = F, box = F, mar = c(1,1,1,1), legend = F)
# plot(dom_bq1, add = T, col = "black", legend = F)
# plot(pred_layers_B_Q2_mean, axes = F, box = F, mar = c(1,1,1,1), legend = T)
# plot(dom_bq2, add = T, col = "black", legend = F)
# plot(pred_layers_B_Q3_mean, axes = F, box = F, mar = c(1,1,1,1), legend = F)
# plot(dom_bq3, add = T, col = "black", legend = F)
# plot(pred_layers_B_Q4_mean, axes = F, box = F, mar = c(1,1,1,1), legend = T)
# plot(dom_bq4, add = T, col = "black",legend = F)
# dev.off()

# look at trying to outline the area of the rasters that have values and overlaying this on the map.
# Try chat GPT to get the code for this??

# 
euro_rast <- terra::rast("output/euro_rast_10k.tif")
euro_map <- terra::vect("output/euro_map.shp")

masked_bq1 <- terra::mask(pred_layers_B_Q1_mean, dom_bq1) 
masked_bq2 <- terra::mask(pred_layers_B_Q2_mean, dom_bq2) 
masked_bq3 <- terra::mask(pred_layers_B_Q3_mean, dom_bq3) 
masked_bq4 <- terra::mask(pred_layers_B_Q4_mean, dom_bq4) 


# four panel plot of these maps
#set the colour palette
library(viridis)

#dev.off()
png("plots/dom_over_predsB_masked.png", width = 700, height = 600)
#pdf("plots/dom_over_predsB_masked.pdf", width = 7, height = 6)
par(mfrow = c(2,2))
par(mar = c(0,0,0,0))
plot(euro_rast, col = "black", alpha = 0.4, legend = F)
#plot(masked_bq1, add = T, col = turbo(100))
plot(masked_bq1, add = T, col = viridis(100))
title("B-Q1", adj = 0)

plot(euro_rast, col = "black", alpha = 0.4, legend = F)
#plot(masked_bq2, add = T, col = turbo(100))
plot(masked_bq2, add = T, col = viridis(100))
title("B-Q2", adj = 0)

plot(euro_rast, col = "black", alpha = 0.4, legend = F)
#plot(masked_bq3, add = T, col = turbo(100))
plot(masked_bq3, add = T, col = viridis(100))
title("B-Q3", adj = 0)

plot(euro_rast, col = "black", alpha = 0.4, legend = F)
#plot(masked_bq4, add = T, col = turbo(100))
plot(masked_bq4, add = T, col = viridis(100))
title("B-Q4", adj = 0)
dev.off()


### Look at producing some histograms

# The ones below don't account for the rest of the model predictions or the background level of domestic birds. 

# dom_probs_bq1 <- terra::mask(pred_layers_B_Q1_mean, dom_bq1)
# dom_probs_bq2 <- terra::mask(pred_layers_B_Q2_mean, dom_bq2)
# dom_probs_bq3 <- terra::mask(pred_layers_B_Q3_mean, dom_bq3)
# dom_probs_bq4 <- terra::mask(pred_layers_B_Q4_mean, dom_bq4)

# dev.off()
# png("plots/hists_of_probs_domestic.png")
# par(mfrow = c(2,2))
# hist(values(dom_probs_bq1), col = "#6699CC", main = "B - Q1", xlab = "Predicted probability")
# hist(values(dom_probs_bq2), col = "#339966", main = "B - Q2", xlab = "Predicted probability")
# hist(values(dom_probs_bq3), col = "#FF9933", main = "B - Q3", xlab = "Predicted probability")
# hist(values(dom_probs_bq4), col = "#99CCFF", main = "B - Q4", xlab = "Predicted probability")
# dev.off()

# read in the chicken and duck density layers 
ducks_2010 <- terra::rast("variable_manipulation/variable_outputs/duck_density_2010_10kres.tif")
plot(ducks_2010)

chucks_2010 <- terra::rast("variable_manipulation/variable_outputs/chicken_density_2010_10kres.tif")
plot(chucks_2010)

par(mfrow = c(1,2))
plot(chucks_2010); plot(ducks_2010)

ducks_n_chucks <- chucks_2010 + ducks_2010 
  #terra::mosaic(ducks_2010, chucks_2010, fun = "sum")
plot(ducks_n_chucks)
ducks_n_chucks
names(ducks_n_chucks) <- "chicken_duck_density"

range(ducks_2010)
range(chucks_2010)
range(ducks_n_chucks)

par(mfrow = c(2,2))
plot(chucks_2010); plot(ducks_2010); plot(ducks_n_chucks) # the combined one looks like the chicken layer but
# I think this is just because the numbers are so much larger with the chickens.
chucks_2010 == ducks_n_chucks

dev.off()
hist(ducks_n_chucks) # values are too spread out
hist(log(ducks_n_chucks), breaks = seq(-5,15,0.5), main = "Histogram of combined chicken and duck density", 
     xlab = "log of density")

## trial some different values for the cutoff

dnc_1000 <- ducks_n_chucks
dnc_1000[dnc_1000 <= 1000] <- NA
plot(dnc_1000)
dnc_1000
#table(values(dnc_1000))

dnc_5000 <- ducks_n_chucks
dnc_5000[dnc_5000 <= 5000] <- NA
plot(dnc_5000)


# extract the cells from the probability raster that have values for the dnc layer
bq1_dnc <- terra::mask(pred_layers_B_Q1_mean, dnc_1000)
bq1_dnc
plot(bq1_dnc)
## from this layer, make another layer with the 0-10% probabilities
bq1_dnc_0010 <- bq1_dnc
bq1_dnc_0010[bq1_dnc_0010 >= 0.1] <- NA
#plot(bq1_dnc_0010)
#plot(dom_bq1, add = T , col = "red")
pos_doms_try <- terra::mask(bq1_dnc_0010, dom_bq1)
plot(pos_doms_try, col = "red")
inf_cells <- sum(table(values(pos_doms_try))) # this is the number of cells with domestic cases
non_inf_cells <- sum(table(values(bq1_dnc_0010))) - inf_cells

odds_0010 <- inf_cells/non_inf_cells

# this needs to go in to a function otherwise will take forever! 

# bq1_dnc_0120 <- bq1_dnc
# bq1_dnc_0120[bq1_dnc_0120 < 0.1 | bq1_dnc_0120 >= 0.2] <- NA
# range(bq1_dnc_0120)

odds_fun <- function(prediction_rast, cd_rast, min_prob, max_prob, dom_rast){
  pred_masked <- terra::mask(prediction_rast, cd_rast)
  prob_rast <- pred_masked
  prob_rast[prob_rast < min_prob | prob_rast >= max_prob] <- NA
  pos_doms <- terra::mask(prob_rast, dom_rast)
  inf_cells <- sum(table(values(pos_doms))) 
  non_inf_cells <- sum(table(values(prob_rast))) - inf_cells
  odds_out <- inf_cells/non_inf_cells
  res_out <- c(inf_cells, non_inf_cells, odds_out)
  return(res_out)
}

## check to see if get same results for 0 - 0.1 probs using function as doing manually

odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0, max_prob = 0.1, dom_rast = dom_bq1)
odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0.1, max_prob = 0.2, dom_rast = dom_bq1)
odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0.2, max_prob = 0.3, dom_rast = dom_bq1)
odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0.3, max_prob = 0.4, dom_rast = dom_bq1)
odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0.4, max_prob = 0.5, dom_rast = dom_bq1)
odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0.5, max_prob = 0.6, dom_rast = dom_bq1)
odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0.6, max_prob = 0.7, dom_rast = dom_bq1)
odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0.7, max_prob = 0.8, dom_rast = dom_bq1)
odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0.8, max_prob = 0.9, dom_rast = dom_bq1)
odds_fun(prediction_rast =  pred_layers_B_Q1_mean, cd_rast = dnc_1000, min_prob = 0.9, max_prob = 1.00001, dom_rast = dom_bq1)

# looks good

min_prob_vect <- seq(0, 0.9, 0.1)
max_prob_vect <- seq(0.1, 1, 0.1)

BQ1_inf_cells <- c()
BQ1_non_inf_cells <- c()
BQ1_odds <- c()

for (i in 1:length(min_prob_vect)) {
  ressy <- odds_fun(prediction_rast =  pred_layers_B_Q1_mean,
                    cd_rast = dnc_1000,
                    min_prob = min_prob_vect[i],
                    max_prob = max_prob_vect[i],
                    dom_rast = dom_bq1)
  BQ1_inf_cells[i] <- ressy[1]
  BQ1_non_inf_cells[i] <- ressy[2]
  BQ1_odds[i] <- ressy[3]
}

sum(BQ1_inf_cells) + sum(BQ1_non_inf_cells)
sum(!is.na(values(dnc_1000)))
sum(is.na(values(dnc_1000)))
sum(!is.na(values(dnc_1000))) + sum(is.na(values(dnc_1000)))
490*440

# something not quite adding up. 
# next try looking at the number of values in each band in the masked predictions rast. USe hist and values
pred_masked <- terra::mask(pred_layers_B_Q1_mean, dnc_1000)
sum(hist(pred_masked, breaks = seq(0,1,0.1))$counts)
sum(!is.na(values(pred_layers_B_Q1_mean)))

## are there cells that are NA in prediction that are not NA in 

if (!compareGeom(raster1, raster2, stopOnError = FALSE)) {
  stop("The rasters do not have the same extent and resolution")
}

# Find the cells in raster1 that are NA
na_raster1 <- is.na(pred_layers_B_Q1_mean)

# Find the cells in raster2 that are not NA
not_na_raster2 <- !is.na(dnc_1000)

# Combine the conditions to find cells that are NA in raster1 and not NA in raster2
na_in_raster1_not_in_raster2 <- na_raster1 & not_na_raster2

# Count the number of such cells
na_count <- sum(values(na_in_raster1_not_in_raster2))
plot(na_in_raster1_not_in_raster2)

# so it seems that there are cells that have a value for chicken and ducks but don't have a prediction, 
# which is why they don't match
# this is probably where the missing values are. 
# so produce the histograms from the numbers that I have,

bq1_hist_dat <- data.frame(
  Odds = BQ1_odds,
  Predicted_prob = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", 
                      "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
)

# Create a bar plot
q1_hist <- ggplot(bq1_hist_dat, aes(x = Predicted_prob, y = BQ1_odds)) +
  geom_bar(stat = "identity", fill = "#2271B2",) +
  labs(title = "B-Q1", x = "Predicted probability", y = "Odds") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
q1_hist

## Rpt for the other quarters
BQ2_inf_cells <- c()
BQ2_non_inf_cells <- c()
BQ2_odds <- c()

for (i in 1:length(min_prob_vect)) {
  ressy <- odds_fun(prediction_rast =  pred_layers_B_Q2_mean,
                    cd_rast = dnc_1000,
                    min_prob = min_prob_vect[i],
                    max_prob = max_prob_vect[i],
                    dom_rast = dom_bq2)
  BQ2_inf_cells[i] <- ressy[1]
  BQ2_non_inf_cells[i] <- ressy[2]
  BQ2_odds[i] <- ressy[3]
}

bq2_hist_dat <- data.frame(
  Odds = BQ2_odds,
  Predicted_prob = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", 
                     "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
)

# Create a bar plot
q2_hist <- ggplot(bq2_hist_dat, aes(x = Predicted_prob, y = BQ2_odds)) +
  geom_bar(stat = "identity", fill = "#F748A5") +
  labs(title = "B-Q2", x = "Predicted probability", y = "Odds") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
q2_hist

# Q3

BQ3_inf_cells <- c()
BQ3_non_inf_cells <- c()
BQ3_odds <- c()

for (i in 1:length(min_prob_vect)) {
  ressy <- odds_fun(prediction_rast =  pred_layers_B_Q3_mean,
                    cd_rast = dnc_1000,
                    min_prob = min_prob_vect[i],
                    max_prob = max_prob_vect[i],
                    dom_rast = dom_bq3)
  BQ3_inf_cells[i] <- ressy[1]
  BQ3_non_inf_cells[i] <- ressy[2]
  BQ3_odds[i] <- ressy[3]
}

bq3_hist_dat <- data.frame(
  Odds = BQ3_odds,
  Predicted_prob = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", 
                     "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
)

# Create a bar plot
q3_hist <- ggplot(bq3_hist_dat, aes(x = Predicted_prob, y = BQ3_odds)) +
  geom_bar(stat = "identity", fill = "#359B73") +
  labs(title = "B-Q3", x = "Predicted probability", y = "Odds") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
q3_hist

# Q4
BQ4_inf_cells <- c()
BQ4_non_inf_cells <- c()
BQ4_odds <- c()

for (i in 1:length(min_prob_vect)) {
  ressy <- odds_fun(prediction_rast =  pred_layers_B_Q4_mean,
                    cd_rast = dnc_1000,
                    min_prob = min_prob_vect[i],
                    max_prob = max_prob_vect[i],
                    dom_rast = dom_bq4)
  BQ4_inf_cells[i] <- ressy[1]
  BQ4_non_inf_cells[i] <- ressy[2]
  BQ4_odds[i] <- ressy[3]
}

bq4_hist_dat <- data.frame(
  Odds = BQ4_odds,
  Predicted_prob = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", 
                     "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
)

# Create a bar plot
q4_hist <- ggplot(bq4_hist_dat, aes(x = Predicted_prob, y = BQ4_odds)) +
  geom_bar(stat = "identity", fill = "#e69f00") +
  labs(title = "B-Q4", x = "Predicted probability", y = "Odds") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
q4_hist


combo_hists_B <- gridExtra::grid.arrange(q1_hist, q2_hist, q3_hist, q4_hist, ncol = 2, nrow = 2)
#ggsave("plots/histograms_for_B_with_cnd.png", combo_hists_B)


##################################################################################################

# Repeat for A


## For outputs from model A
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


masked_aq1 <- terra::mask(pred_layers_A_Q1_mean, dom_aq1) 
masked_aq2 <- terra::mask(pred_layers_A_Q2_mean, dom_aq2) 
masked_aq3 <- terra::mask(pred_layers_A_Q3_mean, dom_aq3) 
masked_aq4 <- terra::mask(pred_layers_A_Q4_mean, dom_aq4) 


# four panel plot of these maps

#dev.off()
png("plots/dom_over_predsA_masked.png", width = 700, height = 600)
#pdf("plots/dom_over_predsA_masked.pdf", width = 7, height = 6)
par(mfrow = c(2,2))
par(mar = c(0,0,0,0))
plot(euro_rast, col = "black", alpha = 0.4, legend = F)
plot(masked_aq1, add = T, col = viridis(100))
title("A-Q1", adj = 0)

plot(euro_rast, col = "black", alpha = 0.4, legend = F)
plot(masked_aq2, add = T, col = viridis(100))
title("A-Q2", adj = 0)

plot(euro_rast, col = "black", alpha = 0.4, legend = F)
plot(masked_aq3, add = T, col = viridis(100))
title("A-Q3", adj = 0)

plot(euro_rast, col = "black", alpha = 0.4, legend = F)
plot(masked_aq4, add = T, col = viridis(100))
title("A-Q4", adj = 0)
dev.off()


### Look at producing some histograms


AQ1_inf_cells <- c()
AQ1_non_inf_cells <- c()
AQ1_odds <- c()

for (i in 1:length(min_prob_vect)) {
  ressy <- odds_fun(prediction_rast =  pred_layers_A_Q1_mean,
                    cd_rast = dnc_1000,
                    min_prob = min_prob_vect[i],
                    max_prob = max_prob_vect[i],
                    dom_rast = dom_aq1)
  AQ1_inf_cells[i] <- ressy[1]
  AQ1_non_inf_cells[i] <- ressy[2]
  AQ1_odds[i] <- ressy[3]
}

sum(AQ1_inf_cells) + sum(AQ1_non_inf_cells)

aq1_hist_dat <- data.frame(
  Odds = AQ1_odds,
  Predicted_prob = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", 
                     "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
)

# Create a bar plot
aq1_hist <- ggplot(aq1_hist_dat, aes(x = Predicted_prob, y = AQ1_odds)) +
  geom_bar(stat = "identity", fill = "#2271B2") +
  labs(title = "A-Q1", x = "Predicted probability", y = "Odds") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
aq1_hist

## Rpt for the other quarters
AQ2_inf_cells <- c()
AQ2_non_inf_cells <- c()
AQ2_odds <- c()

for (i in 1:length(min_prob_vect)) {
  ressy <- odds_fun(prediction_rast =  pred_layers_A_Q2_mean,
                    cd_rast = dnc_1000,
                    min_prob = min_prob_vect[i],
                    max_prob = max_prob_vect[i],
                    dom_rast = dom_aq2)
  AQ2_inf_cells[i] <- ressy[1]
  AQ2_non_inf_cells[i] <- ressy[2]
  AQ2_odds[i] <- ressy[3]
}

aq2_hist_dat <- data.frame(
  Odds = AQ2_odds,
  Predicted_prob = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", 
                     "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
)

# Create a bar plot
aq2_hist <- ggplot(aq2_hist_dat, aes(x = Predicted_prob, y = AQ2_odds)) +
  geom_bar(stat = "identity", fill = "#F748A5") +
  labs(title = "A-Q2", x = "Predicted probability", y = "Odds") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
aq2_hist

# Q3

AQ3_inf_cells <- c()
AQ3_non_inf_cells <- c()
AQ3_odds <- c()

for (i in 1:length(min_prob_vect)) {
  ressy <- odds_fun(prediction_rast =  pred_layers_A_Q3_mean,
                    cd_rast = dnc_1000,
                    min_prob = min_prob_vect[i],
                    max_prob = max_prob_vect[i],
                    dom_rast = dom_aq3)
  AQ3_inf_cells[i] <- ressy[1]
  AQ3_non_inf_cells[i] <- ressy[2]
  AQ3_odds[i] <- ressy[3]
}

aq3_hist_dat <- data.frame(
  Odds = AQ3_odds,
  Predicted_prob = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", 
                     "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
)

# Create a bar plot
aq3_hist <- ggplot(aq3_hist_dat, aes(x = Predicted_prob, y = AQ3_odds)) +
  geom_bar(stat = "identity", fill = "#359B73") +
  labs(title = "A-Q3", x = "Predicted probability", y = "Odds") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
aq3_hist

# Q4
AQ4_inf_cells <- c()
AQ4_non_inf_cells <- c()
AQ4_odds <- c()

for (i in 1:length(min_prob_vect)) {
  ressy <- odds_fun(prediction_rast =  pred_layers_A_Q4_mean,
                    cd_rast = dnc_1000,
                    min_prob = min_prob_vect[i],
                    max_prob = max_prob_vect[i],
                    dom_rast = dom_aq4)
  AQ4_inf_cells[i] <- ressy[1]
  AQ4_non_inf_cells[i] <- ressy[2]
  AQ4_odds[i] <- ressy[3]
}

aq4_hist_dat <- data.frame(
  Odds = AQ4_odds,
  Predicted_prob = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", 
                     "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
)

# Create a bar plot
aq4_hist <- ggplot(aq4_hist_dat, aes(x = Predicted_prob, y = AQ4_odds)) +
  geom_bar(stat = "identity", fill = "#e69f00") +
  labs(title = "A-Q4", x = "Predicted probability", y = "Odds") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
aq4_hist


combo_hists_A <- gridExtra::grid.arrange(aq1_hist, aq2_hist, aq3_hist, aq4_hist, ncol = 2, nrow = 2)
#ggsave("plots/histograms_for_A_with_cnd.png", combo_hists_A)





# 
# png("plots/dom_over_preds_a1.png")
# plot(pred_layers_A_Q1_mean)
# plot(dom_aq1, add = T, col = "black", legend = F)
# dev.off()
# 
# png("plots/dom_over_preds_a2.png")
# plot(pred_layers_A_Q2_mean)
# plot(dom_aq2, add = T, col = "black", legend = F)
# dev.off()
# 
# png("plots/dom_over_preds_a3.png")
# plot(pred_layers_A_Q3_mean)
# plot(dom_aq3, add = T, col = "black", legend = F)
# dev.off()
# 
# png("plots/dom_over_preds_a4.png")
# plot(pred_layers_A_Q4_mean)
# plot(dom_aq4, add = T, col = "black", legend = F)
# dev.off()
# 
# dev.off()
# png("plots/dom_over_predsA_allq.png")
# par(mfrow = c(2,2))
# plot(pred_layers_A_Q1_mean, axes = F, box = F, mar = c(1,1,1,1), legend = F)
# plot(dom_aq1, add = T, col = "black", legend = F)
# plot(pred_layers_A_Q2_mean, axes = F, box = F, mar = c(1,1,1,1), legend = T)
# plot(dom_aq2, add = T, col = "black", legend = F)
# plot(pred_layers_A_Q3_mean, axes = F, box = F, mar = c(1,1,1,1), legend = F)
# plot(dom_aq3, add = T, col = "black", legend = F)
# plot(pred_layers_A_Q4_mean, axes = F, box = F, mar = c(1,1,1,1), legend = T)
# plot(dom_aq4, add = T, col = "black",legend = F)
# dev.off()

masked_aq1 <- terra::mask(pred_layers_A_Q1_mean, dom_bq1) 
masked_aq2 <- terra::mask(pred_layers_A_Q2_mean, dom_bq2) 
masked_aq3 <- terra::mask(pred_layers_A_Q3_mean, dom_bq3) 
masked_aq4 <- terra::mask(pred_layers_A_Q4_mean, dom_bq4) 





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


