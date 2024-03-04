# In this script we use pre-assembled covariates and testing and training data
# coordinates to construct training and testing datasets which we save as csv's.

PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

################################################################################
# Dataset A Q1

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q1.RDS")
ytrain <- training_coords$pos
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtrain)))
training_data <- xtrain[-bad_rows, ]
training_data$y <- ytrain[-bad_rows]
write.csv(training_data, "training_sets/training_data_A_Q1.csv", row.names = FALSE)

test_coords <- readRDS("training_sets/test_coords_A_Q1.RDS")
ytest <- test_coords$const
xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtest)))
test_data <- xtest[-bad_rows, ]
test_data$y <- ytest[-bad_rows]
write.csv(test_data, "training_sets/test_data_A_Q1.csv", row.names = FALSE)

################################################################################
# Dataset A Q2

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q2.RDS")
ytrain <- training_coords$pos
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtrain)))
training_data <- xtrain[-bad_rows, ]
training_data$y <- ytrain[-bad_rows]
write.csv(training_data, "training_sets/training_data_A_Q2.csv", row.names = FALSE)

test_coords <- readRDS("training_sets/test_coords_A_Q2.RDS")
ytest <- test_coords$const
xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtest)))
test_data <- xtest[-bad_rows, ]
test_data$y <- ytest[-bad_rows]
write.csv(test_data, "training_sets/test_data_A_Q2.csv", row.names = FALSE)

################################################################################
# Dataset A Q3

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q3.RDS")
ytrain <- training_coords$pos
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtrain)))
training_data <- xtrain[-bad_rows, ]
training_data$y <- ytrain[-bad_rows]
write.csv(training_data, "training_sets/training_data_A_Q3.csv", row.names = FALSE)

test_coords <- readRDS("training_sets/test_coords_A_Q3.RDS")
ytest <- test_coords$const
xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtest)))
test_data <- xtest[-bad_rows, ]
test_data$y <- ytest[-bad_rows]
write.csv(test_data, "training_sets/test_data_A_Q3.csv", row.names = FALSE)

################################################################################
# Dataset A Q4

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q4.RDS")
ytrain <- training_coords$pos
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtrain)))
training_data <- xtrain[-bad_rows, ]
training_data$y <- ytrain[-bad_rows]
write.csv(training_data, "training_sets/training_data_A_Q4.csv", row.names = FALSE)

test_coords <- readRDS("training_sets/test_coords_A_Q4.RDS")
ytest <- test_coords$const
xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtest)))
test_data <- xtest[-bad_rows, ]
test_data$y <- ytest[-bad_rows]
write.csv(test_data, "training_sets/test_data_A_Q4.csv", row.names = FALSE)

################################################################################
# Dataset B Q1

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_B_Q1.RDS")
ytrain <- training_coords$pos
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtrain)))
training_data <- xtrain[-bad_rows, ]
training_data$y <- ytrain[-bad_rows]
write.csv(training_data, "training_sets/training_data_B_Q1.csv", row.names = FALSE)

################################################################################
# Dataset B Q2

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_B_Q2.RDS")
ytrain <- training_coords$pos
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtrain)))
training_data <- xtrain[-bad_rows, ]
training_data$y <- ytrain[-bad_rows]
write.csv(training_data, "training_sets/training_data_B_Q2.csv", row.names = FALSE)

################################################################################
# Dataset B Q3

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_B_Q3.RDS")
ytrain <- training_coords$pos
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtrain)))
training_data <- xtrain[-bad_rows, ]
training_data$y <- ytrain[-bad_rows]
write.csv(training_data, "training_sets/training_data_B_Q3.csv", row.names = FALSE)

################################################################################
# Dataset B Q4

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_B_Q4.RDS")
ytrain <- training_coords$pos
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtrain)))
training_data <- xtrain[-bad_rows, ]
training_data$y <- ytrain[-bad_rows]
write.csv(training_data, "training_sets/training_data_B_Q4.csv", row.names = FALSE)