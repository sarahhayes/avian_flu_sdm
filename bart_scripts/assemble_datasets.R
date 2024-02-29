# In this script we use pre-assembled covariates and testing and training data
# coordinates to construct training and testing datasets which we save as csv's.

################################################################################


covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q1.RDS")
ytrain <- training_coords$pos
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtrain)))
xtrain <- xtrain[-bad_rows, ]
ytrain <- ytrain[-bad_rows]

test_coords <- readRDS("training_sets/test_coords_A_Q2.RDS")
ytest <- test_coords$const
xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))
bad_rows <- which(is.na(rowSums(xtest)))
xtest <- xtest[-bad_rows, ]
ytest <- ytest[-bad_rows]