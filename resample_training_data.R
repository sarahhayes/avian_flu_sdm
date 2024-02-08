############
# Startups #
############

rm(list=ls())

library(tidyverse)
library(magrittr)
library(flexsdm)
library(ibis.iSDM)
library(terra)

# Read in extent of Europe mapped
base_map <- terra::rast("output/euro_rast_10k.tif")

# Read in weighting layers for pseudoabsence sampling
# weight_layer <- terra::rast("variable_manipulation/variable_outputs/human_density_10k.tif")
# layername <- "Human density"

weight_layer <- terra::rast("variable_manipulation/variable_outputs/ebird_records_10k.tif")
layername <- "eBird record"

#weight_layer_log <- weight_layer %>% app(function(x) log(x+1))

# Read in positive flu site data 
pos_sites <- read.csv("data_offline\\Avian flu data\\hpai_pos_birds_nobvbrc.csv") %>% 
  rename(date = observation.date) %>%
  mutate(date = as.Date(date),
         Q = case_when(months(date) %in% month.name[1:3] ~ "Q1",
                       months(date) %in% month.name[4:6] ~ "Q2",
                       months(date) %in% month.name[7:9] ~ "Q3",
                       months(date) %in% month.name[10:12] ~ "Q4"))

pos_sites %>% pull(Q) %>% table
pos_sites %>% pull(serotype_HN) %>% table

# Training set A: all AI before 2020/21 H5N8 outbreak (which includes a H5N8 outbreak in 2017)
train_A <- pos_sites %>% filter(date < as.Date("2020-01-01"))

# Training set A: 2020/21 H5N8 outbreak
test_A <- pos_sites %>% filter(date > as.Date("2020-01-01") & serotype_HN == "H5N8")

# Training set B: all AI after 2020/2021 H5N8 outbreak (almost all H5N1)
train_B <- pos_sites %>% filter(date > as.Date("2021-09-01")) 

# How many per quarter
bind_rows(train_A %>% mutate(df = "train_A"),
          test_A %>% mutate(df = "test_A"),
          train_B %>% mutate(df = "train_B")
) %>% with(., table(df, Q))

# # Plot over time
# train_A %>% 
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#00BA38")
# 
# test_A %>% 
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#F8766D")
# 
# train_B %>% 
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#619CFF")

bind_rows(train_A %>% mutate(df = "train_A"),
          test_A %>% mutate(df = "test_A"),
          train_B %>% mutate(df = "train_B")
) %>% 
  ggplot(aes(x = date, fill = df)) +
  geom_histogram()

bind_rows(train_A %>% mutate(df = "train_A"),
          test_A %>% mutate(df = "test_A"),
          train_B %>% mutate(df = "train_B")
) %>% 
  ggplot(aes(x = date, fill = df)) +
  geom_histogram() +
  facet_wrap(~ df, scales = "free")

#######################
# One-time processing #
#######################

# # Read in and assemble environmental predictor layers from individual csvs and rasters
# cov_coast <- read.csv("data_offline\\Environmental variable csvs\\dist_to_coast_output_10kres.csv") %>% rename(x = X, y = Y) %>% select(-X.1)
# cov_water <- read.csv("data_offline\\Environmental variable csvs\\dist_to_water_output_10kres.csv") %>% select(-ID)
# 
# t <- purrr::reduce(
#   list(cov_coast, cov_water),
#   dplyr::left_join,
#   by = c("x", "y")) %>%
#   relocate(x, y) %>%
#   terra::rast(type = "xyz", crs = "epsg:3035") %>%
#   c(.,
#     terra::rast("data_offline\\Environmental rasters\\elevation_max_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_diff_first_quart_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_prec_first_quart_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_relative_humidity_q1_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_temp_q1_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\ndvi_first_quart_2022_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\variation_in_quarterly_mean_temp_q1_10kres.tif"))
# t %>% set.names(c("dist_to_coast_km", "dist_to_water", "elev_max", "diurn_temp", "prec", "humid", "mean_temp", "ndvi", "seas_temp"))
# t %>% writeRaster("data_offline\\combi_rasters\\env_vars_Q1.tif", overwrite=TRUE)
# 
# t <- purrr::reduce(
#   list(cov_coast, cov_water),
#   dplyr::left_join,
#   by = c("x", "y")) %>%
#   relocate(x, y) %>%
#   terra::rast(type = "xyz", crs = "epsg:3035") %>%
#   c(.,
#     terra::rast("data_offline\\Environmental rasters\\elevation_max_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_diff_second_quart_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_prec_second_quart_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_relative_humidity_q2_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_temp_q2_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\ndvi_second_quart_2022_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\variation_in_quarterly_mean_temp_q2_10kres.tif"))
# t %>% set.names(c("dist_to_coast_km", "dist_to_water", "elev_max", "diurn_temp", "prec", "humid", "mean_temp", "ndvi", "seas_temp"))
# t %>% writeRaster("data_offline\\combi_rasters\\env_vars_Q2.tif", overwrite=TRUE)
# 
# t <- purrr::reduce(
#   list(cov_coast, cov_water),
#   dplyr::left_join,
#   by = c("x", "y")) %>%
#   relocate(x, y) %>%
#   terra::rast(type = "xyz", crs = "epsg:3035") %>%
#   c(.,
#     terra::rast("data_offline\\Environmental rasters\\elevation_max_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_diff_third_quart_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_prec_third_quart_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_relative_humidity_q3_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_temp_q3_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\ndvi_third_quart_2022_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\variation_in_quarterly_mean_temp_q3_10kres.tif"))
# t %>% set.names(c("dist_to_coast_km", "dist_to_water", "elev_max", "diurn_temp", "prec", "humid", "mean_temp", "ndvi", "seas_temp"))
# t %>% writeRaster("data_offline\\combi_rasters\\env_vars_Q3.tif", overwrite=TRUE)
# 
# t <- purrr::reduce(
#   list(cov_coast, cov_water),
#   dplyr::left_join,
#   by = c("x", "y")) %>%
#   relocate(x, y) %>%
#   terra::rast(type = "xyz", crs = "epsg:3035") %>%
#   c(.,
#     terra::rast("data_offline\\Environmental rasters\\elevation_max_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_diff_fourth_quart_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_prec_fourth_quart_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_relative_humidity_q4_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\mean_temp_q4_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\ndvi_fourth_quart_2022_10kres.tif"),
#     terra::rast("data_offline\\Environmental rasters\\variation_in_quarterly_mean_temp_q4_10kres.tif"))
# t %>% set.names(c("dist_to_coast_km", "dist_to_water", "elev_max", "diurn_temp", "prec", "humid", "mean_temp", "ndvi", "seas_temp"))
# t %>% writeRaster("data_offline\\combi_rasters\\env_vars_Q4.tif", overwrite=TRUE)
# 
# gc()

#########################
# Sample pseudoabsences #
#########################

pseudoabs_gen <- function (df){
  
  dfname <- deparse(substitute(df))
  
  for (i in 1:4){
    
    # Define buffer area around positive points as the spatial area to sample pseudoabsences from
    ca <- calib_area(
      data = df %>% filter(Q == paste0("Q", i)),
      x = "X",
      y = "Y",
      method = c("buffer", width = 75000), 
      crs = crs(base_map)
    )
    
    pseudoabs_dens <- sample_background(
      data = df %>% filter(Q == paste0("Q", i)),
      x = "X",
      y = "Y",
      n = nrow(df %>% filter(Q == paste0("Q", i))), # Sample pseudoabsences at 1:1 ratio with presences
      method=c('biased'), 
      rlayer = base_map,
      rbias = weight_layer %>%
        app(function(x) (x - minmax(weight_layer)[1])/(minmax(weight_layer)[2] - minmax(weight_layer)[1])),  # RBIAS ARGUMENT - THIS HAS TO BE 0-1 SCALED
      calibarea = ca # Use calibration area as the valid total sampling area for pseudoabsences
    )
    
    pseudoabs_dens_nolim <- sample_background(
      data = df %>% filter(Q == paste0("Q", i)),
      x = "X",
      y = "Y",
      n = nrow(df %>% filter(Q == paste0("Q", i))), # Sample pseudoabsences at 1:1 ratio with presences
      method=c('biased'), 
      rlayer = base_map,
      rbias = weight_layer %>%
        app(function(x) (x - minmax(weight_layer)[1])/(minmax(weight_layer)[2] - minmax(weight_layer)[1])),  # RBIAS ARGUMENT - THIS HAS TO BE 0-1 SCALED
    )
    
    pseudoabs_dens_buff <- add_pseudoabsence(df %>% filter(Q == paste0("Q", i)) %>% select(X,Y) %>% mutate(pr_ab = 1),
                                             field_occurrence = "pr_ab",
                                             template = base_map,
                                             settings = pseudoabs_settings(nrpoints = nrow(df %>% filter(Q == paste0("Q", i))),
                                                                           method = "buffer",
                                                                           buffer_distance = 25000,
                                                                           bias = weight_layer))
    
    
    png(paste0("plots//resampling//weighted//", dfname, "_pseudoabs_", layername %>% gsub(" ", "", .), "_lim_Q", i, ".png"), width = 10, height = 10, units = "in", res = 600)
    plot(base_map, col = "gray95",  main = paste0(layername, "-weighted pseudoabsences (75km limit, no buffer) [flexsdm]"))
    plot(ca, add = TRUE)
    points(df %>% filter(Q == paste0("Q", i)) %>% pull(X),
           df %>% filter(Q == paste0("Q", i)) %>% pull(Y), 
           pch=16, cex=0.2)
    points(pseudoabs_dens %>% pull(X),
           pseudoabs_dens %>% pull(Y),
           pch=16, cex=0.2, col = "magenta3")
    dev.off()
    
    saveRDS(pseudoabs_dens, paste0("plots//resampling//weighted//", dfname, "_pseudoabs_", layername %>% gsub(" ", "", .), "_lim_Q", i, ".RDS"))
    
    
    png(paste0("plots//resampling//weighted//", dfname, "_pseudoabs_", layername %>% gsub(" ", "", .), "_nolim_Q", i, ".png"), width = 10, height = 10, units = "in", res = 600)
    plot(base_map, col = "gray95",  main = paste0(layername, "-weighted pseudoabsences (no limit, no buffer) [flexsdm]"))
    #plot(ca, add = TRUE)
    points(df %>% filter(Q == paste0("Q", i)) %>% pull(X),
           df %>% filter(Q == paste0("Q", i)) %>% pull(Y), 
           pch=16, cex=0.2)
    points(pseudoabs_dens_nolim %>% pull(X),
           pseudoabs_dens_nolim %>% pull(Y),
           pch=16, cex=0.2, col = "magenta3")
    dev.off()
    
    saveRDS(pseudoabs_dens_nolim, paste0("plots//resampling//weighted//", dfname, "_pseudoabs_", layername %>% gsub(" ", "", .), "_nolim_Q", i, ".RDS"))
    
    
    png(paste0("plots//resampling//weighted//", dfname, "_pseudoabs_", layername %>% gsub(" ", "", .), "_buff_Q", i, ".png"), width = 10, height = 10, units = "in", res = 600)
    plot(base_map, col = "gray95",  main = paste0(layername, "-weighted pseudoabsences (no limit, 25km buffer) [ibis]"))
    #plot(ca, add = TRUE)
    points(pseudoabs_dens_buff %>% filter(pr_ab == 1) %>% pull(x),
           pseudoabs_dens_buff %>% filter(pr_ab == 1) %>% pull(y),
           pch=16, cex=0.2)
    points(pseudoabs_dens_buff %>% filter(pr_ab == 0) %>% pull(x),
           pseudoabs_dens_buff %>% filter(pr_ab == 0) %>% pull(y),
           pch=16, cex=0.2, col = "magenta3")
    dev.off()
    
    saveRDS(pseudoabs_dens_buff %>% 
              filter(pr_ab == 0) %>% as.data.frame() %>% select(x, y, pr_ab) %>% rename(X = x, Y = y),
            paste0("plots//resampling//weighted//", dfname, "_pseudoabs_", layername %>% gsub(" ", "", .), "_buff_Q", i, ".RDS"))
    
    
    gc()
    
  }
}

set.seed(1047)
pseudoabs_gen(train_A)
set.seed(1047)
pseudoabs_gen(train_B)

###########################
# Data filtering/thinning #
###########################

# Read in and assemble environmental predictor layers used in data resampling from .tif
env_vars_Q1 <- terra::rast("data_offline\\combi_rasters\\env_vars_Q1.tif")
env_vars_Q2 <- terra::rast("data_offline\\combi_rasters\\env_vars_Q2.tif")
env_vars_Q3 <- terra::rast("data_offline\\combi_rasters\\env_vars_Q3.tif")
env_vars_Q4 <- terra::rast("data_offline\\combi_rasters\\env_vars_Q4.tif")

# Read pseudoabsences back in - selected pseudoabs weighted by eBird records
pseudoabs_A <- bind_rows(readRDS(paste0("plots//resampling//weighted//train_A_pseudoabs_eBirdrecord_buff_Q1.RDS")) %>% mutate(Q = "Q1"),
                         readRDS(paste0("plots//resampling//weighted//train_A_pseudoabs_eBirdrecord_buff_Q2.RDS")) %>% mutate(Q = "Q2"),
                         readRDS(paste0("plots//resampling//weighted//train_A_pseudoabs_eBirdrecord_buff_Q3.RDS")) %>% mutate(Q = "Q3"),
                         readRDS(paste0("plots//resampling//weighted//train_A_pseudoabs_eBirdrecord_buff_Q4.RDS")) %>% mutate(Q = "Q4")
)

pseudoabs_B <- bind_rows(readRDS(paste0("plots//resampling//weighted//train_B_pseudoabs_eBirdrecord_buff_Q1.RDS")) %>% mutate(Q = "Q1"),
                         readRDS(paste0("plots//resampling//weighted//train_B_pseudoabs_eBirdrecord_buff_Q2.RDS")) %>% mutate(Q = "Q2"),
                         readRDS(paste0("plots//resampling//weighted//train_B_pseudoabs_eBirdrecord_buff_Q3.RDS")) %>% mutate(Q = "Q3"),
                         readRDS(paste0("plots//resampling//weighted//train_B_pseudoabs_eBirdrecord_buff_Q4.RDS")) %>% mutate(Q = "Q4")
)

# Conduct thinning
thinner <- function (df, nbins){
  
  dfname <- deparse(substitute(df))
  n_thinned <- matrix(nrow = 4, ncol = 3)
  
  for (i in 1:4){
    
    # Thin data randomly
    
    thinned_rand <- df %>% filter(Q == paste0("Q", i)) %>% 
      slice_sample(n = min(1000, nrow(df %>% filter(Q == paste0("Q", i)))))
    
    ## Thin data based on geographic distance
    #
    # pos_thinned_geo_dist <- occfilt_geo(
    #   data = pos_sites %>% filter(Q == paste0("Q", i)) %>% mutate(id = row_number()),
    #   x = "X",
    #   y = "Y",
    #   method = c("defined", d = "25"),
    #   env_layer = terra::rast(paste0("data_offline\\combi_rasters\\env_vars_Q", i, ".tif"))[[c(1:9)]],
    #   prj = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # EPSG 3035 as defined by https://epsg.io/3035
    # )
    
    # Thin data based on stratifying multidimensional environmental space
    # Computation time (and filtered data) scales with number of bins, as you're randomly sampling each section of a stratified multidimensional space
    # "Therefore, it is recommended to use a small number of bins: between 2-5 if more than ten variables are used."
    
    thinned_env <-  df %>% filter(Q == paste0("Q", i)) %>% mutate(id = row_number()) %>%
      occfilt_env(
        x = "X",
        y = "Y",
        id = "id",
        env_layer = terra::rast(paste0("data_offline\\combi_rasters\\env_vars_Q", i, ".tif")),
        nbins = nbins
      )
    
    png(paste0("plots//resampling//thinning//", dfname, "_thinned_rand_Q", i, ".png"), width = 10, height = 10, units = "in", res = 600)
    plot(base_map, col = "gray95", main = paste0(dfname, " data thinning: random"))
    points(df %>% filter(Q == paste0("Q", i)) %>% pull(X),
           df %>% filter(Q == paste0("Q", i)) %>% pull(Y), 
           pch=16, cex=0.2, col = "gray40")
    points(thinned_rand %>% pull(X),
           thinned_rand %>% pull(Y), 
           pch=16, cex=0.2, col = "red")
    dev.off()
    
    thinned_rand %>% saveRDS(paste0("training_sets//raw//", dfname, "_thinned_rand_Q", i, ".RDS"))
    
    png(paste0("plots//resampling//thinning//", dfname, "_thinned_env_", nbins, "_Q", i, ".png"), width = 10, height = 10, units = "in", res = 600)
    plot(base_map, col = "gray95", main = paste0(dfname, " data thinning: environmental space (9 variables, ", nbins, " bins)"))
    points(df %>% filter(Q == paste0("Q", i)) %>% pull(X),
           df %>% filter(Q == paste0("Q", i)) %>% pull(Y), 
           pch=16, cex=0.1, col = "gray40")
    points(thinned_env %>% pull(X),
           thinned_env %>% pull(Y), 
           pch=16, cex=0.1, col = "red")
    dev.off()
    
    thinned_env %>% saveRDS(paste0("training_sets//raw//", dfname, "_thinned_env_", nbins, "_Q", i, ".RDS"))
    
    n_thinned[i,1] <- nrow(df %>% filter(Q == paste0("Q", i)))
    n_thinned[i,2] <- nrow(thinned_rand)
    n_thinned[i,3] <- nrow(thinned_env)
    
    gc()
    
  }
  
  colnames(n_thinned) <- c("none", "rand", "env")
  rownames(n_thinned) <- c("Q1", "Q2", "Q3", "Q4")
  return(n_thinned)
  
}

for (j in c(4,6,8,10)){
  
  set.seed(1650)
  n_thinned_train_A <- thinner(train_A, nbins = j)
  
  set.seed(1650)
  n_thinned_pseudoabs_A <- thinner(pseudoabs_A, nbins = j)
  
  set.seed(1650)
  n_thinned_train_B <- thinner(train_B, nbins = j)
  
  set.seed(1650)
  n_thinned_pseudoabs_B <- thinner(pseudoabs_B, nbins = j)
  
  sink(paste0("training_sets//n_thinned_", j, "bins.txt"))
  print(list(n_thinned_train_A = n_thinned_train_A, 
       n_thinned_pseudoabs_A = n_thinned_pseudoabs_A, 
       n_thinned_train_B = n_thinned_train_B, 
       n_thinned_pseudoabs_B = n_thinned_pseudoabs_B))
  sink()
  
}

# ## TO DO: Better data wrangling for plotting rather than constructing df manually
# t <- structure(list(set = c("train_A", "train_A", "train_A", "train_A", 
#                             "train_A", "train_A", "train_A", "train_A", "train_A", "train_A", 
#                             "train_A", "train_A", "train_A", "train_A", "train_A", "train_A", 
#                             "train_A", "train_A", "train_A", "train_A", "train_B", "train_B", 
#                             "train_B", "train_B", "train_B", "train_B", "train_B", "train_B", 
#                             "train_B", "train_B", "train_B", "train_B", "train_B", "train_B", 
#                             "train_B", "train_B", "train_B", "train_B", "train_B", "train_B"
# ), Q = c("Q1", "Q2", "Q3", "Q4", "Q1", "Q2", "Q3", "Q4", "Q1", 
#          "Q2", "Q3", "Q4", "Q1", "Q2", "Q3", "Q4", "Q1", "Q2", "Q3", "Q4", 
#          "Q1", "Q2", "Q3", "Q4", "Q1", "Q2", "Q3", "Q4", "Q1", "Q2", "Q3", 
#          "Q4", "Q1", "Q2", "Q3", "Q4", "Q1", "Q2", "Q3", "Q4"), strat = c("N", 
#                                                                           "N", "N", "N", "10", "10", "10", "10", "8", "8", "8", "8", "6", 
#                                                                           "6", "6", "6", "4", "4", "4", "4", "N", "N", "N", "N", "10", 
#                                                                           "10", "10", "10", "8", "8", "8", "8", "6", "6", "6", "6", "4", 
#                                                                           "4", "4", "4"), n = c(1145L, 113L, 53L, 438L, 632L, 82L, 30L, 
#                                                                                                 248L, 586L, 81L, 31L, 237L, 483L, 81L, 30L, 209L, 322L, 70L, 
#                                                                                                 25L, 142L, 2874L, 1314L, 659L, 1699L, 995L, 564L, 279L, 645L, 
#                                                                                                 843L, 466L, 260L, 542L, 633L, 362L, 207L, 417L, 359L, 187L, 122L, 
#                                                                                                 211L)), class = "data.frame", row.names = c(NA, -40L))
# 
# g <- t %>%
#   mutate(strat = factor(strat, levels = c("N", "10", "8", "6", "4"))) %>%
#   ggplot(aes(x = strat, y = n, color = Q, fill = Q, group = Q)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~ set, scales = "free", nrow = 2) +
#   theme_bw()

#####################################
# Combine final training data files #
#####################################

selectedbins <- 6

for (i in 1:4){
  bind_rows(readRDS(paste0("training_sets//raw//pseudoabs_A_thinned_env_", selectedbins, "_Q", i, ".RDS")) %>% mutate(pos = 0) %>% select(-id),
            readRDS(paste0("training_sets//raw//train_A_thinned_env_", selectedbins, "_Q", i, ".RDS")) %>% mutate(pos = 1) %>% select(-id)) %>%
    saveRDS(paste0("training_sets//training_coords_A_Q", i, ".RDS"))
  
  bind_rows(readRDS(paste0("training_sets//raw//pseudoabs_B_thinned_env_", selectedbins, "_Q", i, ".RDS")) %>% mutate(pos = 0) %>% select(-id),
            readRDS(paste0("training_sets//raw//train_B_thinned_env_", selectedbins, "_Q", i, ".RDS")) %>% mutate(pos = 1) %>% select(-id)) %>%
    saveRDS(paste0("training_sets//training_coords_B_Q", i, ".RDS"))
}