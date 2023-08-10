############
# Startups #
############

rm(list=ls())

library(tidyverse)
library(magrittr)
library(flexsdm)
library(terra)

setwd("E:\\Working\\avian_flu_sdm")

# Read in positive flu site data 

# pos_sites <- read.csv("data_offline\\Avian flu data\\pos_points_proj_area_all_sources.csv") %>%
#   rename(date = X.1, dataset = observation.date) %>% # fix names
#   select(date, dataset, X, Y) %>%
#   mutate(date = as.Date(date),
#          Q = case_when(months(date) %in% month.name[1:3] ~ "Q1",
#                        months(date) %in% month.name[4:6] ~ "Q2",
#                        months(date) %in% month.name[7:9] ~ "Q3",
#                        months(date) %in% month.name[10:12] ~ "Q4"))

pos_sites <- read.csv("data_offline\\Avian flu data\\pos_points_proj_area_all_sources_duplicates_removed.csv") %>%
  rename(date = X.1, dataset = observation.date) %>% # fix names
  select(date, dataset, X, Y) %>%
  mutate(date = as.Date(date),
         Q = case_when(months(date) %in% month.name[1:3] ~ "Q1",
                       months(date) %in% month.name[4:6] ~ "Q2",
                       months(date) %in% month.name[7:9] ~ "Q3",
                       months(date) %in% month.name[10:12] ~ "Q4"))

pos_sites %>% pull(Q) %>% table

# Read in extent of Europe mapped
base_map <- terra::rast("output\\euro_rast.tif")

#######################
# One-time processing #
#######################

# Read in and assemble environmental predictor layers from .csv, save as .tif for faster read-in
# cov_coast <- read.csv("data_offline\\Environmental variable csvs\\dist_to_coast_output.csv") %>% rename(x = X, y = Y) %>% select(-X.1)
# cov_water <- read.csv("data_offline\\Environmental variable csvs\\dist_to_water_output.csv") %>% select(-ID)
# cov_alt <- read.csv("data_offline\\Environmental variable csvs\\elevation_outputs.csv") %>% select(-ID, -X)
# 
# purrr::reduce(
#   list(
#     data.table::fread("data_offline\\Environmental variable csvs\\climate_output.csv",
#                       select = c("x", "y", "mean_prec_first_quart", "mean_diff_first_quart", "mean_tmin_first_quart", "mean_tmax_first_quart")),
#     cov_coast, cov_water, cov_alt),
#   dplyr::left_join,
#   by = c("x","y")) %>%
#   relocate(x, y) %>%
#   terra::rast(type = "xyz", crs = "epsg:3035") %>%
#   writeRaster("data_offline\\rasters\\env_vars_Q1.tif", overwrite=TRUE)
# 
# purrr::reduce(
#   list(
#     data.table::fread("data_offline\\Environmental variable csvs\\climate_output.csv",
#                       select = c("x", "y", "mean_prec_second_quart", "mean_diff_second_quart", "mean_tmin_second_quart", "mean_tmax_second_quart")),
#     cov_coast, cov_water, cov_alt),
#   dplyr::left_join,
#   by = c("x","y")) %>%
#   relocate(x, y) %>%
#   terra::rast(type = "xyz", crs = "epsg:3035") %>%
#   writeRaster("data_offline\\rasters\\env_vars_Q2.tif", overwrite=TRUE)
# 
# purrr::reduce(
#   list(
#     data.table::fread("data_offline\\Environmental variable csvs\\climate_output.csv",
#                       select = c("x", "y", "mean_prec_third_quart", "mean_diff_third_quart", "mean_tmin_third_quart", "mean_tmax_third_quart")),
#     cov_coast, cov_water, cov_alt),
#   dplyr::left_join,
#   by = c("x","y")) %>%
#   relocate(x, y) %>%
#   terra::rast(type = "xyz", crs = "epsg:3035") %>%
#   writeRaster("data_offline\\rasters\\env_vars_Q3.tif", overwrite=TRUE)
# 
# purrr::reduce(
#   list(
#     data.table::fread("data_offline\\Environmental variable csvs\\climate_output.csv",
#                       select = c("x", "y", "mean_prec_fourth_quart", "mean_diff_fourth_quart", "mean_tmin_fourth_quart", "mean_tmax_fourth_quart")),
#     cov_coast, cov_water, cov_alt),
#   dplyr::left_join,
#   by = c("x","y")) %>%
#   relocate(x, y) %>%
#   terra::rast(type = "xyz", crs = "epsg:3035") %>%
#   writeRaster("data_offline\\rasters\\env_vars_Q4.tif", overwrite=TRUE)
#
# gc()

# # Read in and assemble environmental predictor layers from .tif
# env_vars_Q1 <- terra::rast("data_offline\\rasters\\env_vars_Q1.tif")
# env_vars_Q2 <- terra::rast("data_offline\\rasters\\env_vars_Q2.tif")
# env_vars_Q3 <- terra::rast("data_offline\\rasters\\env_vars_Q3.tif")
# env_vars_Q4 <- terra::rast("data_offline\\rasters\\env_vars_Q4.tif")

#########################
# Sample pseudoabsences #
#########################

set.seed(1047)
for (i in 1:4){
  
  # Define buffer area around positive points as the spatial area to sample pseudoabsences from
  ca <- calib_area(
    data = pos_sites %>% filter(Q == paste0("Q",i)),
    x = "X",
    y = "Y",
    method = c("buffer", width = 75000), 
    crs = crs(base_map)
  )
  
  # Sample pseudoabsences randomly
  pseudoabs_rand <- sample_pseudoabs(
    data = pos_sites %>% filter(Q == paste0("Q",i)),
    x = "X",
    y = "Y",
    n = nrow(pos_sites %>% filter(Q == paste0("Q",i))), # Sample pseudoabsences at 1:1 ratio with presences
    method = "random",
    rlayer = base_map,
    calibarea = ca # Use calibration area as the valid total sampling area for pseudoabsences
  )
  
  # Sample pseudoabsences based on a geographical buffer
  pseudoabs_geo <- sample_pseudoabs(
    data = pos_sites %>% filter(Q == paste0("Q",i)),
    x = "X",
    y = "Y",
    n = nrow(pos_sites %>% filter(Q == paste0("Q",i))), # Sample pseudoabsences at 1:1 ratio with presences
    method=c('geo_const', width='25000'),  
    rlayer = base_map,
    calibarea = ca # Use calibration area as the valid total sampling area for pseudoabsences
  )
  
  # # Sample pseudoabsences based on three-level procedure: spatial buffer from presences, constrained environmentally (like Negative Dissimilarity Based Selection), and distribute pseudoabsences across environmental space with k-means cluster analysis
  # # Takes 9h with i5-8500 CPU
  # pseudoabs_thr <- sample_pseudoabs(
  #   data = pos_sites %>% filter(Q == paste0("Q",i)),
  #   x = "X",
  #   y = "Y",
  #   n = nrow(pos_sites %>% filter(Q == paste0("Q",i))), # Sample pseudoabsences at 1:1 ratio with presences
  #   method=c('geo_env_km_const', 
  #            width='25000', 
  #            env = terra::rast(paste0("data_offline\\rasters\\env_vars_Q",i,".tif"))[[c(1:6,8)]]), # 
  #   rlayer = base_map,
  #   calibarea = ca # Use calibration area as the valid total sampling area for pseudoabsences
  # )
  
  png(paste0("plots//resampling//pseudoabs_rand_Q",i,".png"), width = 10, height = 10, units = "in", res = 600)
  plot(base_map, col = "gray95",  main = "Randomly sampled pseudoabsences")
  plot(ca, add = TRUE)
  points(pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(X),
         pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(Y), 
         pch=16, cex=0.1)
  points(pseudoabs_rand, pch=16, cex=0.1, col = "red")
  dev.off()
  
  saveRDS(pseudoabs_rand, paste0("plots//resampling//pseudoabs_rand_Q",i,".RDS"))
  
  png(paste0("plots//resampling//pseudoabs_geo_Q",i,".png"), width = 10, height = 10, units = "in", res = 600)
  plot(base_map, col = "gray95",  main = "Geographically sampled pseudoabsences (25km buffer)")
  plot(ca, add = TRUE)
  points(pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(X),
         pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(Y), 
         pch=16, cex=0.1)
  points(pseudoabs_geo, pch=16, cex=0.1, col = "purple")
  dev.off()
  
  saveRDS(pseudoabs_geo, paste0("plots//resampling//pseudoabs_geo_Q",i,".RDS"))
  
  # png(paste0("plots//resampling//pseudoabs_thr_Q",i,".png"), width = 10, height = 10, units = "in", res = 600)
  # plot(base_map, col = "gray95", main = "Three-level sampled pseudoabsences")
  # plot(ca, add = TRUE)
  # points(pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(X),
  #        pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(Y), 
  #        pch=16, cex=0.1)
  # points(pseudoabs_thr, pch=16, cex=0.1, col = "blue")
  # dev.off()
  # 
  # saveRDS(pseudoabs_thr, paste0("plots//resampling//pseudoabs_thr_Q",i,".RDS"))
  
  gc()
}


###########################
# Data filtering/thinning #
###########################

## Positive records

set.seed(1533)

# Thin data randomly

pos_thinned_rand <- pos_sites %>% filter(Q == paste0("Q",i)) %>% slice_sample(n = 1000)

# Thin data inversely proportional to human population density

#pos_thinned_popdens <- pos_sites %>% filter(Q == paste0("Q",i)) %>% slice_sample(n = 1000, weight_by = 1/density)

# Thin data based on geographic distance

# pos_thinned_geo_dist <- occfilt_geo(
#   data = pos_sites %>% filter(Q == paste0("Q",i)) %>% mutate(id = row_number()),
#   x = "X",
#   y = "Y",
#   method = c("defined", d = "25"),
#   env_layer = terra::rast(paste0("data_offline\\rasters\\env_vars_Q",i,".tif"))[[c(1:9)]],
#   prj = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # EPSG 3035 as defined by https://epsg.io/3035
# )

# pos_thinned_geo_moran <- occfilt_geo(
#   data = pos_sites %>% filter(Q == paste0("Q",i)) %>% mutate(id = row_number()),
#   x = "X",
#   y = "Y",
#   method = "moran",
#   env_layer = terra::rast(paste0("data_offline\\rasters\\env_vars_Q",i,".tif"))[[c(1:9)]],
#   prj = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # EPSG 3035 as defined by https://epsg.io/3035
# )

png(paste0("plots//resampling//thinning//pos_thinned_rand_Q",i,".png"), width = 10, height = 10, units = "in", res = 600)
plot(base_map, col = "gray95", main = "Data thinning: random")
points(pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(X),
       pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(Y), 
       pch=16, cex=0.1, col = "gray70")
points(pos_thinned_rand %>% pull(X),
       pos_thinned_rand %>% pull(Y), 
       pch=16, cex=0.1, col = "red")
dev.off()

# png(paste0("plots//resampling//thinning//pos_thinned_geo_dist_Q",i,".png"), width = 10, height = 10, units = "in", res = 600)
# plot(base_map, col = "gray95", main = "Data thinning: geographic distance")
# points(pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(X),
#        pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(Y), 
#        pch=16, cex=0.1, col = "gray70")
# points(pos_thinned_geo_dist %>% pull(X),
#        pos_thinned_geo_dist %>% pull(Y), 
#        pch=16, cex=0.1, col = "red")
# dev.off()

# png(paste0("plots//resampling//thinning//pos_thinned_moran_Q",i,".png"), width = 10, height = 10, units = "in", res = 600)
# plot(base_map, main = "Data thinning: Moran's i distance")
# points(pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(X),
#        pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(Y),  
#        pch=16, cex=0.1, col = "gray60")
# points(pos_thinned_geo_moran %>% pull(X),
#        pos_thinned_geo_moran %>% pull(Y), 
#        pch=16, cex=0.1)
# dev.off()

# Thin data based on stratifying multidimensional environmental space
# Computation time (and filtered data) scales with number of bins, as you're randomly sampling each section of a stratified multidimensional space
# "Therefore, it is recommended to use a small number of bins: between 2-5 if more than ten variables are used."

pos_thinned_env_9_7 <- occfilt_env(
  data = pos_sites %>% filter(Q == paste0("Q",i)) %>% mutate(id = row_number()),
  x = "X",
  y = "Y",
  id = "id",
  env_layer = terra::rast(paste0("data_offline\\rasters\\env_vars_Q",i,".tif"))[[c(1:6,8)]],
  nbins = 9
)

png(paste0("plots//resampling//thinning//pos_thinned_env_9_7_Q",i,".png"), width = 10, height = 10, units = "in", res = 600)
plot(base_map, col = "gray95", main = "Data thinning: environmental space (7 variables, 9 bins)")
points(pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(X),
       pos_sites %>% filter(Q == paste0("Q",i)) %>% pull(Y), 
       pch=16, cex=0.1, col = "gray70")
points(pos_thinned_env_9_7 %>% pull(X),
       pos_thinned_env_9_7 %>% pull(Y), 
       pch=16, cex=0.1, col = "red")
dev.off()

## Pseudoabsences

# Select geographically-sampled pseudoabsences
pseudoabs <- readRDS(paste0("plots//resampling//pseudoabs_geo_Q",i,".RDS"))

set.seed(1447)

# Thin data randomly

pseudoabs_thinned_rand <- pseudoabs %>% slice_sample(n = 1000)

# Thin data based on stratifying multidimensional environmental space

pseudoabs_thinned_env_9_7 <- occfilt_env(
  data = pseudoabs %>% mutate(id = row_number()),
  x = "X",
  y = "Y",
  id = "id",
  env_layer = terra::rast(paste0("data_offline\\rasters\\env_vars_Q",i,".tif"))[[c(1:6,8)]],
  nbins = 9
)

png(paste0("plots//resampling//thinning//pseudoabs_thinned_rand_Q",i,".png"), width = 10, height = 10, units = "in", res = 600)
plot(base_map, col = "gray95", main = "Data thinning: random")
points(pseudoabs %>% pull(X),
       pseudoabs %>% pull(Y), 
       pch=16, cex=0.1, col = "gray70")
points(pseudoabs_thinned_rand %>% pull(X),
       pseudoabs_thinned_rand %>% pull(Y), 
       pch=16, cex=0.1, col = "blue")
dev.off()

png(paste0("plots//resampling//thinning//pseudoabs_thinned_env_9_7_Q",i,".png"), width = 10, height = 10, units = "in", res = 600)
plot(base_map, col = "gray95", main = "Data thinning: environmental space (7 variables, 9 bins)")
points(pseudoabs %>% pull(X),
       pseudoabs %>% pull(Y), 
       pch=16, cex=0.1, col = "gray70")
points(pseudoabs_thinned_env_9_7 %>% pull(X),
       pseudoabs_thinned_env_9_7 %>% pull(Y), 
       pch=16, cex=0.1, col = "blue")
dev.off()








# Data partitions for k-fold cross-validation
# SEEMS TO BE TAKEN CARE OF BY XBART FUNCTION?

# Extract environmental data from rasters
all_points <- bind_rows(spp[,2:4], pa_3)

ex_spp <- sdm_extract(
  data = all_points,
  x = "x",
  y = "y",
  env_layer = somevar, # Raster with environmental variables
  variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
  filter_na = TRUE
)

ex_spp
