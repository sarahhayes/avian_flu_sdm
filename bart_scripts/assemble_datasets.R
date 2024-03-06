# In this script we use pre-assembled covariates and testing and training data
# coordinates to construct training and testing datasets which we save as csv's.

PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

################################################################################
# Try making raster of country idx's, adapted from
# https://gis.stackexchange.com/questions/44139/convert-a-spatialpolygonsdataframe-to-raster-using-rasterize-function


# Create a blank raster with appropriate projection and extent
blank_3035 <- terra::rast("output/euro_rast_10k.tif")
euro_ext <- ext(blank_3035)
crs <- "epsg:3035"

countriesSP <- getMap(resolution='high') %>% spTransform(CRS(crs))
countriesSP$Grd_ranks <- rank(countriesSP$ADMIN)
country_lookup <- data.frame(country=countriesSP$ADMIN,
                             val=rank(countriesSP$Grd_ranks))

r <- rast(nlyrs=1,
          crs=crs,
          extent=extent(countriesSP),
          res=res(blank_3035)) %>% raster()

country_rast <- rasterize(countriesSP, r, field="Grd_ranks", fun="first") %>%
  terra::rast() %>%
  crop(euro_ext)


################################################################################
# Dataset A Q1

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q1.RDS")

country_df <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
mismatch_samples <- which(is.na(rowSums(country_df)))
mismatch_pts <- terra::vect(training_coords[mismatch_samples, ], geom=c("X", "Y"),
                            crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
crs <- "epsg:3035"
euro_ext <- extent(covstack[[1]])

# # change projection and extent. 
# # using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs)
euro_map_crop <- terra::crop(euro_map, euro_ext)
plot(euro_map_crop,
     col = "white",
     background = "azure2",
     axes = FALSE,
     buffer = FALSE,
     mar = c(0, 0, 0, 0))
plot(mismatch_pts, add = T, col = "red", pch = 16, cex = .3)
# 
# # Remove NA's from data - although we should think about how to fix this
training_coords <- training_coords[-mismatch_samples,]

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       mar = c(0, 0, 0, 0))
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
    plot(pts_pos, add = T, col = "red", pch = 16, cex = .3)
    this_country <- country_lookup$country[which(country_lookup$val==country_df$layer[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    cat("This point is in",
        this_country,
        ".\n")
    Sys.sleep(1)
  }
}

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(cov_df)))
bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                       crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
plot(euro_map_crop,
     col = "white",
     background = "azure2",
     axes = FALSE,
     buffer = FALSE,
     mar = c(0, 0, 0, 0))
plot(bad_pts, add = T, col = "red", pch = 16, cex = .3)

# Remove bad points
training_coords <- training_coords[-bad_rows, ]
# Just redraw countries:
country_df <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))

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