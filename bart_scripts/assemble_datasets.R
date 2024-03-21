# In this script we use pre-assembled covariates and testing and training data
# coordinates to construct training and testing datasets which we save as csv's.

PLOT_COUNTRY_VALIDATION <- FALSE
PLOT_PROBLEM_SAMPLES <- FALSE

PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(embarcadero)
library(raster)
library(sp)
library(terra)
library(rworldmap)
set.seed(12345)

################################################################################
# Try making raster of country idx's, adapted from
# https://gis.stackexchange.com/questions/44139/convert-a-spatialpolygonsdataframe-to-raster-using-rasterize-function


# Create a blank raster with appropriate projection and extent
blank_3035 <- terra::rast("output/euro_rast_10k.tif")
euro_ext <- ext(blank_3035)
crs <- "epsg:3035"

countriesSP <- getMap(resolution='low') %>% spTransform(CRS(crs))
countriesSP$Grd_ranks <- rank(countriesSP$ADMIN)
country_lookup <- data.frame(country=countriesSP$ADMIN,
                             val=rank(countriesSP$Grd_ranks))

r <- rast(nlyrs=1,
          crs=crs,
          extent=extent(countriesSP),
          res=res(blank_3035)) %>% raster()

country_rast <- raster::rasterize(countriesSP, r, field="Grd_ranks", fun="first")


################################################################################
# Dataset A Q1

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q1.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_train <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
names(ri_train) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_train)))
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
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, A Q1 training")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_train
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, training_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, A Q1 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  training_coords <- training_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_train <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, A Q1 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, A Q1 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(training_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtrain)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, A Q1 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtrain_shifted <- xtrain
  for (i in bad_rows){
    bad_fields <- which(is.na(xtrain[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], training_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtrain_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtrain_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, A Q1 training")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtrain_shifted <- xtrain_shifted[-bad_rows, ]
    training_coords <- training_coords[-bad_rows, ]
    ri_train <- ri_train[-bad_rows, , drop=FALSE]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtrain <- xtrain_shifted
}
# Assemble training data
training_data <- xtrain
training_data$y <- training_coords$pos
training_data$ri <- sapply(ri_train$country, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(training_data, "training_sets/training_data_A_Q1.csv", row.names = FALSE)

#### Now do testing ####
# Load testing data
test_coords <- readRDS("training_sets/test_coords_A_Q1.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_test <- data.frame(raster::extract(country_rast, test_coords[, 1:2]))
names(ri_test) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_test)))
mismatch_pts <- terra::vect(test_coords[mismatch_samples, ], geom=c("X", "Y"),
                            crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
crs <- "epsg:3035"
euro_ext <- extent(covstack[[1]])

# # change projection and extent. 
# # using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs)
euro_map_crop <- terra::crop(euro_map, euro_ext)
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, A Q1 test")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_test
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, test_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, A Q1 test")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  test_coords <- test_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_test <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, A Q1 test")
  pts_to_plot <- sample(1:nrow(test_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(test_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_test$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, A Q1 test")
  pts_to_plot <- sample(1:nrow(test_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(test_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_test$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtest)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, A Q1 test")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtest_shifted <- xtest
  for (i in bad_rows){
    bad_fields <- which(is.na(xtest[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], test_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtest_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtest_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, A Q1 test")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtest_shifted <- xtest_shifted[-bad_rows, ]
    test_coords <- test_coords[-bad_rows, ]
    ri_test <- ri_test[-bad_rows, ]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtest <- xtest_shifted
}

# Assemble test data
test_data <- xtest
test_data$y <- test_coords$pos
test_data$ri <- sapply(ri_test, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(test_data, "training_sets/test_data_A_Q1.csv", row.names = FALSE)

################################################################################
# Dataset A Q2

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q2.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_train <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
names(ri_train) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_train)))
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
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, A Q2 training")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_train
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, training_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, A Q2 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  training_coords <- training_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_train <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, A Q2 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, A Q2 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(training_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtrain)))
if (length(bad_rows>0)){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, A Q2 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtrain_shifted <- xtrain
  for (i in bad_rows){
    bad_fields <- which(is.na(xtrain[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], training_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtrain_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtrain_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, A Q2 training")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtrain_shifted <- xtrain_shifted[-bad_rows, ]
    training_coords <- training_coords[-bad_rows, ]
    ri_train <- ri_train[-bad_rows, , drop=FALSE]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtrain <- xtrain_shifted
}
# Assemble training data
training_data <- xtrain
training_data$y <- training_coords$pos
training_data$ri <- sapply(ri_train$country, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(training_data, "training_sets/training_data_A_Q2.csv", row.names = FALSE)

#### Now do testing ####
# Load testing data
test_coords <- readRDS("training_sets/test_coords_A_Q2.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_test <- data.frame(raster::extract(country_rast, test_coords[, 1:2]))
names(ri_test) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_test)))
mismatch_pts <- terra::vect(test_coords[mismatch_samples, ], geom=c("X", "Y"),
                            crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
crs <- "epsg:3035"
euro_ext <- extent(covstack[[1]])

# # change projection and extent. 
# # using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs)
euro_map_crop <- terra::crop(euro_map, euro_ext)
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, A Q2 test")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_test
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, test_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, A Q2 test")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  test_coords <- test_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_test <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, A Q2 test")
  pts_to_plot <- sample(1:nrow(test_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(test_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_test$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, A Q2 test")
  pts_to_plot <- sample(1:nrow(test_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(test_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_test$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtest)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, A Q2 test")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtest_shifted <- xtest
  for (i in bad_rows){
    bad_fields <- which(is.na(xtest[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], test_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtest_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtest_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, A Q2 test")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtest_shifted <- xtest_shifted[-bad_rows, ]
    test_coords <- test_coords[-bad_rows, ]
    ri_test <- ri_test[-bad_rows, ]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtest <- xtest_shifted
}
# Assemble test data
test_data <- xtest
test_data$y <- test_coords$pos
test_data$ri <- sapply(ri_test, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(test_data, "training_sets/test_data_A_Q2.csv", row.names = FALSE)

################################################################################
# Dataset A Q3

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q3.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_train <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
names(ri_train) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_train)))
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
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, A Q3 training")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_train
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, training_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, A Q3 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  training_coords <- training_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_train <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, A Q3 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, A Q3 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(training_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtrain)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, A Q3 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtrain_shifted <- xtrain
  for (i in bad_rows){
    bad_fields <- which(is.na(xtrain[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], training_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtrain_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtrain_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, A Q3 training")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtrain_shifted <- xtrain_shifted[-bad_rows, ]
    training_coords <- training_coords[-bad_rows, ]
    ri_train <- ri_train[-bad_rows, , drop=FALSE]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtrain <- xtrain_shifted
}
# Assemble training data
training_data <- xtrain
training_data$y <- training_coords$pos
training_data$ri <- sapply(ri_train$country, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(training_data, "training_sets/training_data_A_Q3.csv", row.names = FALSE)

#### Now do testing ####
# Load testing data
test_coords <- readRDS("training_sets/test_coords_A_Q3.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_test <- data.frame(raster::extract(country_rast, test_coords[, 1:2]))
names(ri_test) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_test)))
mismatch_pts <- terra::vect(test_coords[mismatch_samples, ], geom=c("X", "Y"),
                            crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
crs <- "epsg:3035"
euro_ext <- extent(covstack[[1]])

# # change projection and extent. 
# # using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs)
euro_map_crop <- terra::crop(euro_map, euro_ext)
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, A Q3 test")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_test
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, test_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, A Q3 test")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  test_coords <- test_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_test <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, A Q3 test")
  pts_to_plot <- sample(1:nrow(test_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(test_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_test$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, A Q3 test")
  pts_to_plot <- sample(1:nrow(test_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(test_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_test$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtest)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, A Q3 test")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtest_shifted <- xtest
  for (i in bad_rows){
    bad_fields <- which(is.na(xtest[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], test_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtest_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtest_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, A Q3 test")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtest_shifted <- xtest_shifted[-bad_rows, ]
    test_coords <- test_coords[-bad_rows, ]
    ri_test <- ri_test[-bad_rows, ]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtest <- xtest_shifted
}
# Assemble test data
test_data <- xtest
test_data$y <- test_coords$pos
test_data$ri <- sapply(ri_test, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(test_data, "training_sets/test_data_A_Q3.csv", row.names = FALSE)

################################################################################
# Dataset A Q4

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_A_Q4.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_train <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
names(ri_train) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_train)))
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
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, A Q4 training")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_train
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, training_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, A Q4 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  training_coords <- training_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_train <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, A Q4 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, A Q4 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(training_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtrain)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, A Q4 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtrain_shifted <- xtrain
  for (i in bad_rows){
    bad_fields <- which(is.na(xtrain[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], training_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtrain_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtrain_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, A Q4 training")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtrain_shifted <- xtrain_shifted[-bad_rows, ]
    training_coords <- training_coords[-bad_rows, ]
    ri_train <- ri_train[-bad_rows, , drop=FALSE]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtrain <- xtrain_shifted
}
# Assemble training data
training_data <- xtrain
training_data$y <- training_coords$pos
training_data$ri <- sapply(ri_train$country, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(training_data, "training_sets/training_data_A_Q4.csv", row.names = FALSE)

#### Now do testing ####
# Load testing data
test_coords <- readRDS("training_sets/test_coords_A_Q4.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_test <- data.frame(raster::extract(country_rast, test_coords[, 1:2]))
names(ri_test) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_test)))
mismatch_pts <- terra::vect(test_coords[mismatch_samples, ], geom=c("X", "Y"),
                            crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
crs <- "epsg:3035"
euro_ext <- extent(covstack[[1]])

# # change projection and extent. 
# # using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs)
euro_map_crop <- terra::crop(euro_map, euro_ext)
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, A Q4 test")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_test
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, test_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, A Q4 test")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  test_coords <- test_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_test <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, A Q4 test")
  pts_to_plot <- sample(1:nrow(test_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(test_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_test$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, A Q4 test")
  pts_to_plot <- sample(1:nrow(test_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(test_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_test$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtest)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, A Q4 test")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtest_shifted <- xtest
  for (i in bad_rows){
    bad_fields <- which(is.na(xtest[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], test_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtest_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtest_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(test_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, A Q4 test")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtest_shifted <- xtest_shifted[-bad_rows, ]
    test_coords <- test_coords[-bad_rows, ]
    ri_test <- ri_test[-bad_rows, ]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtest <- xtest_shifted
}
# Assemble test data
test_data <- xtest
test_data$y <- test_coords$pos
test_data$ri <- sapply(ri_test, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(test_data, "training_sets/test_data_A_Q4.csv", row.names = FALSE)

################################################################################
# Dataset B Q1

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_B_Q1.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_train <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
names(ri_train) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_train)))
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
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, B Q1 training")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_train
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, training_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, B Q1 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  training_coords <- training_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_train <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, B Q1 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, B Q1 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(training_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtrain)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, B Q1 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtrain_shifted <- xtrain
  for (i in bad_rows){
    bad_fields <- which(is.na(xtrain[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], training_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtrain_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtrain_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, B Q1 training")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtrain_shifted <- xtrain_shifted[-bad_rows, ]
    training_coords <- training_coords[-bad_rows, ]
    ri_train <- ri_train[-bad_rows, , drop=FALSE]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtrain <- xtrain_shifted
}
# Assemble training data
training_data <- xtrain
training_data$y <- training_coords$pos
training_data$ri <- sapply(ri_train$country, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(training_data, "training_sets/training_data_B_Q1.csv", row.names = FALSE)

################################################################################
# Dataset B Q2

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_B_Q2.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_train <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
names(ri_train) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_train)))
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
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, B Q2 training")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_train
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, training_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, B Q2 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  training_coords <- training_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_train <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, B Q2 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, B Q2 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(training_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtrain)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, B Q2 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtrain_shifted <- xtrain
  for (i in bad_rows){
    bad_fields <- which(is.na(xtrain[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], training_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtrain_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtrain_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, B Q2 training")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtrain_shifted <- xtrain_shifted[-bad_rows, ]
    training_coords <- training_coords[-bad_rows, ]
    ri_train <- ri_train[-bad_rows, , drop=FALSE]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtrain <- xtrain_shifted
}
# Assemble training data
training_data <- xtrain
training_data$y <- training_coords$pos
training_data$ri <- sapply(ri_train$country, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(training_data, "training_sets/training_data_B_Q2.csv", row.names = FALSE)

################################################################################
# Dataset B Q3

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_B_Q3.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_train <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
names(ri_train) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_train)))
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
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, B Q3 training")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_train
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, training_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, B Q3 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  training_coords <- training_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_train <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, B Q3 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, B Q3 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(training_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtrain)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, B Q3 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtrain_shifted <- xtrain
  for (i in bad_rows){
    bad_fields <- which(is.na(xtrain[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], training_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtrain_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtrain_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, B Q3 training")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtrain_shifted <- xtrain_shifted[-bad_rows, ]
    training_coords <- training_coords[-bad_rows, ]
    ri_train <- ri_train[-bad_rows, , drop=FALSE]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtrain <- xtrain_shifted
}
# Assemble training data
training_data <- xtrain
training_data$y <- training_coords$pos
training_data$ri <- sapply(ri_train$country, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(training_data, "training_sets/training_data_B_Q3.csv", row.names = FALSE)

################################################################################
# Dataset B Q4

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_B_Q4.RDS")
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
ri_train <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
names(ri_train) <- c("country")
mismatch_samples <- which(is.na(rowSums(ri_train)))
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
if (PLOT_PROBLEM_SAMPLES){
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       main = "Points with missing country data, B Q4 training")
  plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
}
# Try replacing NA pixels with nearest non-NA
country_shift <- ri_train
for (i in mismatch_samples){
  # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
  candidate_pts <- country_rast[which.min(replace(distanceFromPoints(country_rast, training_coords[i, 1:2]), is.na(country_rast), NA))]
  bad_candidates <- which(is.na(candidate_pts))
  if (length(bad_candidates>0)){
    candidate_pts <- candidate_pts[-bad_candidates, ]
  }
  if (length(candidate_pts)>0){
    country_shift[i, 1] <- candidate_pts[1]
  }
}

# Check if this has worked:
bad_rows <- which(is.na(rowSums(country_shift)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points missing country data after imputation, B Q4 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
  country_shift <- country_shift[-bad_rows, ]
  training_coords <- training_coords[-bad_rows, ]
}else{
  print("Shifting to nearest pixel removed all NA values")
}

# Replace mismatches with shifted versions:
ri_train <- country_shift

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Visual inspection of country assignments, B Q4 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    
  }
  
  # Also specifically check mismatches:
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       main = "Country assignments obtained through imputation, B Q4 training")
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in mismatch_samples){
    pts_pos <- terra::vect(training_coords[i, ], geom=c("X", "Y"),
                           crs =  crs)
    plot(pts_pos, add = T, col = "red", pch = 16, cex = 1.)
    this_country <- country_lookup$country[which(country_lookup$val==ri_train$country[i])]
    text(pts_pos, labels=this_country)
    
  }
}

# Now do covariates
xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(xtrain)))
if (length(bad_rows)>0){
  bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                         crs =  crs)
  if (PLOT_PROBLEM_SAMPLES){
    plot(euro_map_crop,
         col = "white",
         background = "azure2",
         axes = FALSE,
         buffer = FALSE,
         main = "Points with missing covariates, B Q4 training")
    plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
  }
  
  # Try replacing NA pixels with nearest non-NA
  xtrain_shifted <- xtrain
  for (i in bad_rows){
    bad_fields <- which(is.na(xtrain[i, ]))
    # Adapted from https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
    candidate_pts <- covstack[which.min(replace(distanceFromPoints(covstack[[bad_fields[1]]], training_coords[i, 1:2]), is.na(covstack[[bad_fields[1]]]), NA))]
    non_na_candidates <- which(!is.na(rowSums(candidate_pts)))
    if (length(non_na_candidates>0)){
      xtrain_shifted[i, ] <- candidate_pts[non_na_candidates[1], ]
    }
  }
  
  # Check if this has worked:
  bad_rows <- which(is.na(rowSums(xtrain_shifted)))
  if (length(bad_rows)>0){
    bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                           crs =  crs)
    if (PLOT_PROBLEM_SAMPLES){
      plot(euro_map_crop,
           col = "white",
           background = "azure2",
           axes = FALSE,
           buffer = FALSE,
           main = "Points with missing covariates after imputation, B Q4 training")
      plot(bad_pts, add = T, col = "red", pch = 16, cex = 1.)
    }
    
    cat("Removing ", length(bad_rows), " datapoints with NA covariates.\n")
    xtrain_shifted <- xtrain_shifted[-bad_rows, ]
    training_coords <- training_coords[-bad_rows, ]
    ri_train <- ri_train[-bad_rows, , drop=FALSE]
  }else{
    print("Shifting to nearest pixel removed all NA values")
  }
  
  # Replace with shifted covariates
  xtrain <- xtrain_shifted
}
# Assemble training data
training_data <- xtrain
training_data$y <- training_coords$pos
training_data$ri <- sapply(ri_train$country, FUN=function(i){country_lookup$country[country_lookup$val==i]})
write.csv(training_data, "training_sets/training_data_B_Q4.csv", row.names = FALSE)