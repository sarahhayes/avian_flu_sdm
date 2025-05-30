############
# Startups #
############

rm(list=ls())

library(tidyverse)
library(magrittr)
# library(flexsdm)
# library(ibis.iSDM)
library(terra)
library(sf)
library(patchwork)

# Set ratio of pseudoabsences to positives (x:1)
pseud_ratio <- 1

# Do you need to initialise the environmental rasters for use in environmental thinning?
# If this is the first time this script is run, should be set to TRUE

INIT_ENV_VARS <- FALSE

# Set season plot palette

pal <- c("#2271B2",
         "#F748A5",
         "#359B73",
         "#e69f00")

#############
# Read data #
#############

# Read in extent of Europe mapped
base_map <- terra::rast("output/euro_rast_10k.tif")
euro_map_st <- st_read("output/euro_map.shp") # country borders

# Read in weighting layers for pseudoabsence sampling
# weight_layer <- terra::rast("variable_manipulation/variable_outputs/human_density.tif")
# layername <- "Human density"

weight_layer <- terra::rast("variable_manipulation/variable_outputs/ebird_records.tif")
layername <- "eBird record"

#weight_layer_log <- weight_layer %>% app(function(x) log(x+1))

# Define starts of ecological seasons
year_start <- as.Date("0000-01-01")
q2_start <- as.Date("0000-03-01") # pre-breeding migration
q3_start <- as.Date("0000-06-07") # breeding season
q4_start <- as.Date("0000-08-10") # post-breeding migration
q1_start <- as.Date("0000-11-30") # non-breeding season
year_end <- as.Date("0000-12-31")

# Read in positive flu site data 
pos_sites <- read.csv("hpai_pos_birds_nobvbrc.csv") %>% 
  rename(date = observation.date) %>%
  mutate(date = as.Date(date),
         Q = case_when((format(date, "%m%-%d") >= format(year_start, "%m%-%d")) & (format(date, "%m%-%d") < format(q2_start-1, "%m%-%d)")) ~ "Q1",
                       (format(date, "%m%-%d") >= format(q2_start, "%m%-%d")) & (format(date, "%m%-%d") < format(q3_start-1, "%m%-%d)")) ~ "Q2",
                       (format(date, "%m%-%d") >= format(q3_start, "%m%-%d")) & (format(date, "%m%-%d") < format(q4_start-1, "%m%-%d)")) ~ "Q3",
                       (format(date, "%m%-%d") >= format(q4_start, "%m%-%d")) & (format(date, "%m%-%d") < format(q1_start-1, "%m%-%d)")) ~ "Q4",
                       (format(date, "%m%-%d") >= format(q1_start, "%m%-%d")) & (format(date, "%m%-%d") <= format(year_end, "%m%-%d)")) ~ "Q1"))

pos_sites %>% pull(Q) %>% table
pos_sites %>% pull(serotype_HN) %>% table

# Data set A: 2.3.4.4b H5NX before H5N1  (includes H5N8, H5N6, retain ambiguous H5 or unlabelled HPAI [n = 80])
# Training set A: 2.3.4.4b H5NX in distinct 2016/2017 outbreak
# Test set A: 2.3.4.4b H5NX in distinct 2020/2021 outbreak
df_A <- pos_sites %>% 
  filter(date >= as.Date("2016-08-10") & date < as.Date("2021-08-10") & serotype_HN %in% c("H5N8", "H5N6", "H5", "")) %>%
  mutate(df = case_when(date < as.Date("2020-08-10") ~ "train_A",
                        date >= as.Date("2020-08-10") ~ "test_A"))

# Data set B: 2.3.4.4b H5N1 (retain ambiguous H5 or unlabelled HPAI [n = 64])
# Training set B: 2.3.4.4b H5N1 from Sep 21 - Mar 23 %>% 
# Test set B: H5N1 from Apr 23 - Mar 24 
df_B <- pos_sites %>% 
  filter(date >= as.Date("2021-08-10") & date < as.Date("2024-03-01") & serotype_HN %in% c("H5N1", "H5", "")) %>%
  mutate(df = case_when(date < as.Date("2023-03-01") ~ "train_B",
                        date >= as.Date("2023-03-01") ~ "test_B"))

# How many per quarter
pre_table <- bind_rows(df_A, df_B) %>% with(., table(df, Q))
pre_table

# Plot over time
# df_A %>% filter(grepl("train", df)) %>%
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#00BA38")
# 
# df_A %>% filter(grepl("test", df)) %>%
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#63db87")
# 
# df_B %>% filter(grepl("train", df)) %>%
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#619CFF")
# 
# df_B %>% filter(grepl("test", df)) %>%
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#99bfff")

bind_rows(df_A, df_B) %>% 
  ggplot(aes(x = date, fill = df)) +
  geom_histogram(bins=100)

timeplot <- pos_sites %>% 
  filter(date > as.Date("2016-08-10") & serotype_HN %in% c("H5N8", "H5N6", "H5N1", "H5", "")) %>%
  mutate(serotype_HN = case_when(
    serotype_HN == "H5N1" ~ "H5N1",
    serotype_HN == "H5N8" ~ "H5N8",
    serotype_HN == "H5N6" ~ "H5N6",
    TRUE ~ "H5NX"
  )) %>%
  ggplot(aes(x = date, fill = serotype_HN)) +
  geom_histogram(bins = round(as.numeric((max(pos_sites$date)-min(pos_sites$date))/30)), position = "stack") +
  geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_vline(xintercept = as.Date("2021-08-10")) +
  geom_vline(xintercept = as.Date("2023-03-01")) +
  geom_text(aes(x = as.Date("2016-08-10")+(as.Date("2020-08-10")-as.Date("2016-08-10"))/2, y = 330 ,label = "A. Train", hjust = 0.5)) +
  geom_text(aes(x = as.Date("2020-08-10")+(as.Date("2021-08-10")-as.Date("2020-08-10"))/2, y = 330 ,label = "A. Test", hjust = 0.5)) +
  geom_text(aes(x = as.Date("2021-08-10")+(as.Date("2023-03-01")-as.Date("2021-08-10"))/2, y = 330 ,label = "B. Train", hjust = 0.5)) +
  geom_text(aes(x = as.Date("2024-01-01"), y = 330 ,label = "B. Test", hjust = 0.5)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  ylab("Monthly reports") +
  xlab("Date") +
  labs(fill = "subtype") +
  theme_bw() +
  theme(legend.position = c(0.05,0.77),
        legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5)))

timeplot_seas <- pos_sites %>% 
  filter(date > as.Date("2016-08-10") & serotype_HN %in% c("H5N8", "H5N6", "H5N1", "H5", "")) %>%
  mutate(serotype_HN = case_when(
    serotype_HN == "H5N1" ~ "H5N1",
    serotype_HN == "H5N8" ~ "H5N8",
    serotype_HN == "H5N6" ~ "H5N6",
    TRUE ~ "H5NX"
  )) %>%
  ggplot(aes(x = date)) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2016-08-10"), xmax = as.Date("2016-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2017-08-10"), xmax = as.Date("2017-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2018-08-10"), xmax = as.Date("2018-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2019-08-10"), xmax = as.Date("2019-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2020-08-10"), xmax = as.Date("2020-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2021-08-10"), xmax = as.Date("2021-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2022-08-10"), xmax = as.Date("2022-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2023-08-10"), xmax = as.Date("2023-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2016-11-30"), xmax = as.Date("2017-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2017-11-30"), xmax = as.Date("2018-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2018-11-30"), xmax = as.Date("2019-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2019-11-30"), xmax = as.Date("2020-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2020-11-30"), xmax = as.Date("2021-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2021-11-30"), xmax = as.Date("2022-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2022-11-30"), xmax = as.Date("2023-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2023-11-30"), xmax = as.Date("2024-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2017-03-01"), xmax = as.Date("2017-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2018-03-01"), xmax = as.Date("2018-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2019-03-01"), xmax = as.Date("2019-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2020-03-01"), xmax = as.Date("2020-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2021-03-01"), xmax = as.Date("2021-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2022-03-01"), xmax = as.Date("2022-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2023-03-01"), xmax = as.Date("2023-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2017-06-07"), xmax = as.Date("2017-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2018-06-07"), xmax = as.Date("2018-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2019-06-07"), xmax = as.Date("2019-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2020-06-07"), xmax = as.Date("2020-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2021-06-07"), xmax = as.Date("2021-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2022-06-07"), xmax = as.Date("2022-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2023-06-07"), xmax = as.Date("2023-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_histogram(bins = round(as.numeric((max(pos_sites$date)-min(pos_sites$date))/30)), position = "stack", fill = "grey20") +
  geom_vline(xintercept = as.Date("2021-08-10")) +
  geom_vline(xintercept = as.Date("2023-03-01")) +
  geom_text(aes(x = as.Date("2016-08-10")+(as.Date("2020-08-10")-as.Date("2016-08-10"))/2, y = 330 ,label = "A. Train", hjust = 0.5)) +
  geom_text(aes(x = as.Date("2020-08-10")+(as.Date("2021-08-10")-as.Date("2020-08-10"))/2, y = 330 ,label = "A. Test", hjust = 0.5)) +
  geom_text(aes(x = as.Date("2021-08-10")+(as.Date("2023-03-01")-as.Date("2021-08-10"))/2, y = 330 ,label = "B. Train", hjust = 0.5)) +
  geom_text(aes(x = as.Date("2024-01-01"), y = 330 ,label = "B. Test", hjust = 0.5)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  ylab("Monthly reports") +
  xlab("Date") +
  theme_bw()

ggsave(paste("plots/timeplot.png"),
       plot = timeplot,
       width = 12,
       height = 3)

ggsave(paste("plots/timeplot_seas.png"),
       plot = timeplot_seas,
       width = 12,
       height = 3)

points_a <- st_as_sf(df_A, coords = c("X", "Y"), crs = crs(base_map)) %>% 
  mutate(serotype_HN = case_when(
    serotype_HN == "H5N1" ~ "H5N1",
    serotype_HN == "H5N8" ~ "H5N8",
    serotype_HN == "H5N6" ~ "H5N6",
    TRUE ~ "H5NX"
  ))
points_b <- st_as_sf(df_B, coords = c("X", "Y"), crs = crs(base_map)) %>% 
  mutate(serotype_HN = case_when(
    serotype_HN == "H5N1" ~ "H5N1",
    serotype_HN == "H5N8" ~ "H5N8",
    serotype_HN == "H5N6" ~ "H5N6",
    TRUE ~ "H5NX"
  ))

map_a_data <- 
  ggplot() +
  geom_sf(data = euro_map_st, fill = "grey", color = "black", alpha = 0.4)  +
  # geom_sf(data = points_a, color = "#00BFC4", size = 0.2) +
  geom_sf(data = points_a, aes(color = serotype_HN), size = 0.2, alpha = 0.6, show.legend = F) +
  scale_color_manual(values = c("H5N1" = "#F8766D", "H5N6" = "#7CAE00", "H5N8" = "#00BFC4", "H5NX" = "#C77CFF")) +
  theme_minimal() +
  theme(legend.position="none") +
  labs(x = "Longitude",
       y = "Latitude",
       color = "Serotype")

map_b_data <- 
  ggplot() +
  geom_sf(data = euro_map_st, fill = "grey", color = "black", alpha = 0.4)  +
  geom_sf(data = points_b, aes(color = serotype_HN), size = 0.2, alpha = 0.6, show.legend = F) +
  scale_color_manual(values = c("H5N1" = "#F8766D", "H5N6" = "#7CAE00", "H5N8" = "#00BFC4", "H5NX" = "#C77CFF")) +
  theme_minimal() +
  theme(legend.position="none") +
  labs(x = "Longitude",
       y = "Latitude",
       color = "Serotype")

fig_data_combi <- (map_a_data|map_b_data)/(timeplot) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(4,1.2)) &
  theme(legend.position = 'bottom')

ggsave(paste("plots/data_fig.png"),
       plot = fig_data_combi,
       width = 11,
       height = 8.5)


bind_rows(df_A, df_B)  %>% 
  ggplot(aes(x = date, fill = df)) +
  geom_histogram() +
  facet_wrap(~ df, scales = "free")

# Assign to 10km grid cells by quarter (so that overlapping points can still count for individual quarters)
rast_train_A_list <- list()
rast_test_A_list <- list()
rast_train_B_list <- list()
rast_test_B_list <- list()

for (i in 1:4){
  
  rast_train_A_list[[i]] <- df_A %>% filter(df == "train_A") %>%
    filter(Q == paste0("Q", i)) %>%
    st_as_sf(coords = c("X", "Y"), crs = "EPSG:3035") %>%  # Map to 10km grid
    rasterize(y = base_map, fun = "max") %>%  # Consider positives as 1 no matter how many records in the cell
    as.data.frame(xy = TRUE) %>%  # Convert back to data frame
    filter(complete.cases(.)) %>%  # Select only cells with presences
    set_names(c("X", "Y", "pos")) %>% 
    mutate(Q = paste0("Q",i), df = "train_A")
  
  rast_test_A_list[[i]] <- df_A %>% filter(df == "test_A") %>%
    filter(Q == paste0("Q", i)) %>%
    st_as_sf(coords = c("X", "Y"), crs = "EPSG:3035") %>%  # Map to 10km grid
    rasterize(y = base_map, fun = "max") %>%  # Consider positives as 1 no matter how many records in the cell
    as.data.frame(xy = TRUE) %>%  # Convert back to data frame
    filter(complete.cases(.)) %>%  # Select only cells with presences
    set_names(c("X", "Y", "pos")) %>% 
    mutate(Q = paste0("Q",i), df = "test_A")
  
  rast_train_B_list[[i]] <- df_B %>% filter(df == "train_B") %>%
    filter(Q == paste0("Q", i)) %>%
    st_as_sf(coords = c("X", "Y"), crs = "EPSG:3035") %>%  # Map to 10km grid
    rasterize(y = base_map, fun = "max") %>%  # Consider positives as 1 no matter how many records in the cell
    as.data.frame(xy = TRUE) %>%  # Convert back to data frame
    filter(complete.cases(.)) %>%  # Select only cells with presences
    set_names(c("X", "Y", "pos")) %>% 
    mutate(Q = paste0("Q",i), df = "train_B")
  
  rast_test_B_list[[i]] <- df_B %>% filter(df == "test_B") %>%
    filter(Q == paste0("Q", i)) %>%
    st_as_sf(coords = c("X", "Y"), crs = "EPSG:3035") %>%  # Map to 10km grid
    rasterize(y = base_map, fun = "max") %>%  # Consider positives as 1 no matter how many records in the cell
    as.data.frame(xy = TRUE) %>%  # Convert back to data frame
    filter(complete.cases(.)) %>%  # Select only cells with presences
    set_names(c("X", "Y", "pos")) %>% 
    mutate(Q = paste0("Q",i), df = "test_B")
  
}

df_A <- bind_rows(rast_train_A_list, rast_test_A_list)
df_B <- bind_rows(rast_train_B_list, rast_test_B_list)

# How many per quarter after counting at grid cell level?
post_table <- bind_rows(df_A, df_B) %>% with(., table(df, Q))

pre_table[c(3,1,4,2),]
post_table[c(3,1,4,2),]

#######################
# One-time processing #
#######################

if (INIT_ENV_VARS == TRUE){
  
  # Read in and assemble environmental predictor layers from individual csvs and rasters
  cov_coast <- read.csv("data_offline\\Environmental variable csvs\\dist_to_coast_output_10kres.csv") %>% rename(x = X, y = Y) %>% select(-X.1)
  cov_water <- read.csv("data_offline\\Environmental variable csvs\\dist_to_water_output_10kres.csv") %>% select(-ID)
  
  t <- purrr::reduce(
    list(cov_coast, cov_water),
    dplyr::left_join,
    by = c("x", "y")) %>%
    relocate(x, y) %>%
    terra::rast(type = "xyz", crs = "epsg:3035") %>%
    c(.,
      terra::rast("data_offline\\Environmental rasters\\elevation_max_10kres.tif"),
      terra::rast("data_offline\\Environmental rasters\\isotherm_mean_q1_eco_quarts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_diff_first_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_prec_first_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_relative_humidity_q1_10kres_eco_quarts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_mean_first_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\ndvi_first_quart_2022_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\variation_in_quarterly_mean_temp_q1_eco_rasts.tif"))
  t %>% set.names(c("dist_to_coast_km", "dist_to_water", "elev_max", "isotherm_mean", "diurn_temp", "prec", "humid", "mean_temp", "ndvi", "seas_temp"))
  t %>% writeRaster("data_offline\\combi_rasters\\env_vars_Q1.tif", overwrite=TRUE)
  
  t <- purrr::reduce(
    list(cov_coast, cov_water),
    dplyr::left_join,
    by = c("x", "y")) %>%
    relocate(x, y) %>%
    terra::rast(type = "xyz", crs = "epsg:3035") %>%
    c(.,
      terra::rast("data_offline\\Environmental rasters\\elevation_max_10kres.tif"),
      terra::rast("data_offline\\Environmental rasters\\isotherm_mean_q2_eco_quarts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_diff_second_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_prec_second_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_relative_humidity_q2_10kres_eco_quarts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_mean_second_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\ndvi_second_quart_2022_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\variation_in_quarterly_mean_temp_q2_eco_rasts.tif"))
  t %>% set.names(c("dist_to_coast_km", "dist_to_water", "elev_max", "isotherm_mean", "diurn_temp", "prec", "humid", "mean_temp", "ndvi", "seas_temp"))
  t %>% writeRaster("data_offline\\combi_rasters\\env_vars_Q2.tif", overwrite=TRUE)
  
  t <- purrr::reduce(
    list(cov_coast, cov_water),
    dplyr::left_join,
    by = c("x", "y")) %>%
    relocate(x, y) %>%
    terra::rast(type = "xyz", crs = "epsg:3035") %>%
    c(.,
      terra::rast("data_offline\\Environmental rasters\\elevation_max_10kres.tif"),
      terra::rast("data_offline\\Environmental rasters\\isotherm_mean_q3_eco_quarts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_diff_third_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_prec_third_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_relative_humidity_q3_10kres_eco_quarts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_mean_third_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\ndvi_third_quart_2022_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\variation_in_quarterly_mean_temp_q3_eco_rasts.tif"))
  t %>% set.names(c("dist_to_coast_km", "dist_to_water", "elev_max", "isotherm_mean", "diurn_temp", "prec", "humid", "mean_temp", "ndvi", "seas_temp"))
  t %>% writeRaster("data_offline\\combi_rasters\\env_vars_Q3.tif", overwrite=TRUE)
  
  t <- purrr::reduce(
    list(cov_coast, cov_water),
    dplyr::left_join,
    by = c("x", "y")) %>%
    relocate(x, y) %>%
    terra::rast(type = "xyz", crs = "epsg:3035") %>%
    c(.,
      terra::rast("data_offline\\Environmental rasters\\elevation_max_10kres.tif"),
      terra::rast("data_offline\\Environmental rasters\\isotherm_mean_q4_eco_quarts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_diff_fourth_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_prec_fourth_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_relative_humidity_q4_10kres_eco_quarts.tif"),
      terra::rast("data_offline\\Environmental rasters\\mean_mean_fourth_quart_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\ndvi_fourth_quart_2022_eco_rasts.tif"),
      terra::rast("data_offline\\Environmental rasters\\variation_in_quarterly_mean_temp_q4_eco_rasts.tif"))
  t %>% set.names(c("dist_to_coast_km", "dist_to_water", "elev_max", "isotherm_mean", "diurn_temp", "prec", "humid", "mean_temp", "ndvi", "seas_temp"))
  t %>% writeRaster("data_offline\\combi_rasters\\env_vars_Q4.tif", overwrite=TRUE)
  
  gc()
  
}

#########################
# Sample pseudoabsences #
#########################

pseudoabs_gen <- function (maindf){
  
  dfname <- deparse(substitute(maindf))
  
  for (i in 1:4){
    
    # Select quarter
    df <- maindf %>% filter(Q == paste0("Q",i))
    
    # Define buffer area around positive points as the spatial area to sample pseudoabsences from
    ca <- calib_area(
      data = df,
      x = "X",
      y = "Y",
      method = c("buffer", width = 75000), 
      crs = crs(base_map)
    )
    
    pseudoabs_dens <- sample_background(
      data = df,
      x = "X",
      y = "Y",
      n = round(nrow(df)*pseud_ratio), # Sample pseudoabsences at given ratio
      method=c('biased'), 
      rlayer = base_map,
      rbias = weight_layer %>%
        app(function(x) (x - minmax(weight_layer)[1])/(minmax(weight_layer)[2] - minmax(weight_layer)[1])),  # RBIAS ARGUMENT - THIS HAS TO BE 0-1 SCALED
      calibarea = ca # Use calibration area as the valid total sampling area for pseudoabsences
    )
    
    pseudoabs_dens_nolim <- sample_background(
      data = df,
      x = "X",
      y = "Y",
      n = round(nrow(df)*pseud_ratio), # Sample pseudoabsences at given ratio
      method=c('biased'), 
      rlayer = base_map,
      rbias = weight_layer %>%
        app(function(x) (x - minmax(weight_layer)[1])/(minmax(weight_layer)[2] - minmax(weight_layer)[1])),  # RBIAS ARGUMENT - THIS HAS TO BE 0-1 SCALED
    )
    
    pseudoabs_dens_buff <- add_pseudoabsence(df %>% select(X,Y) %>% mutate(pr_ab = 1),
                                             field_occurrence = "pr_ab",
                                             template = base_map,
                                             settings = pseudoabs_settings(nrpoints = round(nrow(df)*pseud_ratio),
                                                                           method = "buffer",
                                                                           buffer_distance = 25000,
                                                                           bias = weight_layer))
    
    
    png(paste0("plots\\resampling\\weighted_pseudoabs\\", dfname, "_Q", i, "_pseudoabs_", layername %>% gsub(" ", "", .), "_lim.png"), width = 10, height = 10, units = "in", res = 600)
    plot(base_map, col = "gray95",  main = paste0(layername, "-weighted pseudoabsences (75km limit, no buffer) [flexsdm]"))
    plot(ca, add = TRUE)
    points(df %>% pull(X),
           df %>% pull(Y), 
           pch=16, cex=0.2)
    points(pseudoabs_dens %>% pull(X),
           pseudoabs_dens %>% pull(Y),
           pch=16, cex=0.2, col = "magenta3")
    dev.off()
    
    saveRDS(pseudoabs_dens, paste0("plots\\resampling\\weighted_pseudoabs\\", dfname, "_Q", i, "_pseudoabs_", layername %>% gsub(" ", "", .), "_lim.RDS"))
    
    
    png(paste0("plots\\resampling\\weighted_pseudoabs\\", dfname, "_Q", i, "_pseudoabs_", layername %>% gsub(" ", "", .), "_nolim.png"), width = 10, height = 10, units = "in", res = 600)
    plot(base_map, col = "gray95",  main = paste0(layername, "-weighted pseudoabsences (no limit, no buffer) [flexsdm]"))
    #plot(ca, add = TRUE)
    points(df %>% pull(X),
           df %>% pull(Y), 
           pch=16, cex=0.2)
    points(pseudoabs_dens_nolim %>% pull(X),
           pseudoabs_dens_nolim %>% pull(Y),
           pch=16, cex=0.2, col = "magenta3")
    dev.off()
    
    saveRDS(pseudoabs_dens_nolim, paste0("plots\\resampling\\weighted_pseudoabs\\", dfname, "_Q", i, "_pseudoabs_", layername %>% gsub(" ", "", .), "_nolim.RDS"))
    
    png(paste0("plots\\resampling\\weighted_pseudoabs\\", dfname, "_Q", i, "_pseudoabs_", layername %>% gsub(" ", "", .), "_buff.png"), width = 10, height = 10, units = "in", res = 600)
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
              as.data.frame() %>% 
              filter(pr_ab == 0) %>%  
              select(x, y, pr_ab) %>% 
              rename(X = x, Y = y) %>% 
              set_rownames(NULL),
            paste0("plots\\resampling\\weighted_pseudoabs\\", dfname, "_Q", i, "_pseudoabs_", layername %>% gsub(" ", "", .), "_buff.RDS"))
    
    
    gc()
    
  }
}

set.seed(1047)
pseudoabs_gen(df_A)
set.seed(1047)
pseudoabs_gen(df_B)

###########################
# Data filtering/thinning #
###########################

###########################
# If not using a pseudoabs ratio of 1:1, you might need to differentially thin with parameters based on pseudoabs or presence in order to preserve that ratio and not just thin down to the same value??
###########################

# Read in and assemble environmental predictor layers used in data resampling from .tif
env_vars_Q1 <- terra::rast("data_offline\\combi_rasters\\env_vars_Q1.tif")
env_vars_Q2 <- terra::rast("data_offline\\combi_rasters\\env_vars_Q2.tif")
env_vars_Q3 <- terra::rast("data_offline\\combi_rasters\\env_vars_Q3.tif")
env_vars_Q4 <- terra::rast("data_offline\\combi_rasters\\env_vars_Q4.tif")

# Separate training data from test data, each accompanied by same ratio of pseudoabsences
set.seed(1102)

train_A <- df_A %>% filter(df == "train_A")
test_A <- df_A %>% filter(df == "test_A")
train_B <- df_B %>% filter(df == "train_B")
test_B <- df_B %>% filter(df == "test_B")

pseudoabs_train_A <- list()
pseudoabs_test_A <- list()
pseudoabs_train_B <- list()
pseudoabs_test_B <- list()

for (i in 1:4){
  
  # Selected pseudoabs weighted by eBird records
  pseudoabs <- readRDS(paste0("plots\\resampling\\weighted_pseudoabs\\df_A_Q",i,"_pseudoabs_eBirdrecord_buff.RDS")) %>% mutate(Q = paste0("Q",i))
  
  pseudoabs_rows <- pseudoabs %>% nrow
  train_rows <- train_A %>% filter(Q == paste0("Q",i)) %>% nrow
  train_plus_test_rows <- df_A %>% filter(Q == paste0("Q",i)) %>% nrow
  
  train_index <- sample(1:pseudoabs_rows, 
                        size = round(pseudoabs_rows*train_rows/train_plus_test_rows),
                        replace=FALSE)
  
  pseudoabs_train_A[[i]] <- pseudoabs %>% filter(Q == paste0("Q",i)) %>% .[train_index,]
  pseudoabs_test_A[[i]] <- pseudoabs %>% filter(Q == paste0("Q",i)) %>% .[-train_index,]
  
  pseudoabs <- readRDS(paste0("plots\\resampling\\weighted_pseudoabs\\df_B_Q",i,"_pseudoabs_eBirdrecord_buff.RDS")) %>% mutate(Q = paste0("Q",i))
  
  pseudoabs_rows <- pseudoabs %>% filter(Q == paste0("Q",i)) %>% nrow
  train_rows <- train_B %>% filter(Q == paste0("Q",i)) %>% nrow
  train_plus_test_rows <- df_B %>% filter(Q == paste0("Q",i)) %>% nrow
  
  train_index <- sample(1:pseudoabs_rows, 
                        size = round(pseudoabs_rows*train_rows/train_plus_test_rows),
                        replace=FALSE)
  
  pseudoabs_train_B[[i]] <- pseudoabs %>% filter(Q == paste0("Q",i)) %>% .[train_index,]
  pseudoabs_test_B[[i]] <- pseudoabs %>% filter(Q == paste0("Q",i)) %>% .[-train_index,]
  
}

pseudoabs_train_A %<>% bind_rows
pseudoabs_test_A %<>% bind_rows
pseudoabs_train_B %<>% bind_rows
pseudoabs_test_B %<>% bind_rows

# Conduct thinning
thinner <- function (maindf, nbins){
  
  dfname <- deparse(substitute(maindf))
  n_thinned <- matrix(nrow = 4, ncol = 3)
  
  for (i in 1:4){
    
    # Select quarter
    df <- maindf %>% filter(Q == paste0("Q",i))
    
    # Thin data randomly (n = 1000 initially)
    
    thinned_rand <- df %>% 
      slice_sample(n = min(1000, nrow(df)))
    
    ## Thin data based on geographic distance
    #
    # pos_thinned_geo_dist <- occfilt_geo(
    #   data = df %>% mutate(id = row_number()),
    #   x = "X",
    #   y = "Y",
    #   method = c("defined", d = "25"),
    #   env_layer = terra::rast(paste0("data_offline\\combi_rasters\\env_vars_Q", i, ".tif"))[[c(1:9)]],
    #   prj = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # EPSG 3035 as defined by https:\\epsg.io/3035
    # )
    
    # Thin data based on stratifying multidimensional environmental space
    # Computation time (and filtered data) scales with number of bins, as you're randomly sampling each section of a stratified multidimensional space
    # "Therefore, it is recommended to use a small number of bins: between 2-5 if more than ten variables are used."
    
    thinned_env <-  df %>% mutate(id = row_number()) %>%
      occfilt_env(
        x = "X",
        y = "Y",
        id = "id",
        env_layer = terra::rast(paste0("data_offline\\combi_rasters\\env_vars_Q", i, ".tif")),
        nbins = nbins
      )
    
    png(paste0("plots\\resampling\\thinning\\", dfname, "_Q", i, "_thinned_rand.png"), width = 10, height = 10, units = "in", res = 600)
    plot(base_map, col = "gray95", main = paste0(dfname, " data thinning: random"))
    points(df %>% pull(X),
           df %>% pull(Y), 
           pch=16, cex=0.2, col = "gray40")
    points(thinned_rand %>% pull(X),
           thinned_rand %>% pull(Y), 
           pch=16, cex=0.2, col = pal[2])
    dev.off()
    
    thinned_rand %>% saveRDS(paste0("training_sets\\raw\\", dfname, "_Q", i, "_thinned_rand.RDS"))
    
    png(paste0("plots\\resampling\\thinning\\", dfname, "_Q", i, "_thinned_env_", nbins, ".png"), width = 10, height = 10, units = "in", res = 600)
    plot(base_map, col = "gray95", main = paste0(dfname, " data thinning: environmental space (10 variables, ", nbins, " bins)"))
    points(df %>% pull(X),
           df %>% pull(Y), 
           pch=16, cex=0.2, col = "gray40")
    points(thinned_env %>% pull(X),
           thinned_env %>% pull(Y), 
           pch=16, cex=0.2, col = pal[2])
    dev.off()
    
    thinned_env %>% saveRDS(paste0("training_sets\\raw\\", dfname, "_Q", i, "_thinned_env_", nbins, ".RDS"))
    
    n_thinned[i,1] <- nrow(df)
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
  n_thinned_pseudoabs_A <- thinner(pseudoabs_train_A, nbins = j)
  
  set.seed(1650)
  n_thinned_train_B <- thinner(train_B, nbins = j)
  
  set.seed(1650)
  n_thinned_pseudoabs_B <- thinner(pseudoabs_train_B, nbins = j)
  
  bind_rows(n_thinned_train_A %>% reshape2::melt() %>% mutate(set = "train_A"),
            n_thinned_pseudoabs_A %>% reshape2::melt() %>% mutate(set = "pseudoabs_train_A"),
            n_thinned_train_B %>% reshape2::melt() %>% mutate(set = "train_B"),
            n_thinned_pseudoabs_B %>% reshape2::melt() %>% mutate(set = "pseudoabs_train_B")) %>% mutate(nbins = j) %>%
    write.csv(paste0("training_sets\\n_thinned_", j, "bins.csv"))
  
}

thin_df <- list.files("training_sets\\", pattern = "bins.csv", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  select(-X) %>%
  rename(Q = Var1, strat = Var2) %>%
  mutate(nbins = as.character(nbins)) %>%
  mutate(nbins = case_when(strat == "none" ~ "N",
                           TRUE ~ nbins)) %>%
  mutate(nbins = factor(nbins, levels = c("N", "10", "8", "6", "4"))) %>%
  distinct()

thin_df %>% 
  filter(strat != "rand") %>% 
  reshape2::dcast(set+nbins ~ Q) %>% 
  mutate(set = factor(set, levels = c("train_A", "pseudoabs_train_A", "train_B", "pseudoabs_train_B"))) %>%
  filter(grepl("_A", set)) %>%
  arrange(nbins, set) %>%
  write.csv("training_sets\\thin_summary_A.csv")

thin_df %>% 
  filter(strat != "rand") %>% 
  reshape2::dcast(set+nbins ~ Q) %>% 
  mutate(set = factor(set, levels = c("train_A", "pseudoabs_train_A", "train_B", "pseudoabs_train_B"))) %>%
  filter(grepl("_B", set)) %>%
  arrange(nbins, set) %>%
  write.csv("training_sets\\thin_summary_B.csv")

thin_df %>%
  filter(strat != "rand") %>%
  ggplot(aes(x = nbins, y = value, color = Q, fill = Q, group = Q)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ set, scales = "free", nrow = 4) +
  theme_bw()

#####################################
# Combine final training data files #
#####################################

selectedbins <- 6

# Combine thinned presences and thinned pseudoabsences into training sets

for (i in 1:4){
  bind_rows(readRDS(paste0("training_sets\\raw\\pseudoabs_train_A_Q",i,"_thinned_env_",selectedbins,".RDS")) %>% mutate(pos = 0) %>% select(-id),
            readRDS(paste0("training_sets\\raw\\train_A_Q",i,"_thinned_env_",selectedbins,".RDS")) %>% mutate(pos = 1) %>% select(-id)) %>% as.data.frame %>%
    saveRDS(paste0("training_sets\\training_coords_A_Q", i, ".RDS"))
  
  bind_rows(readRDS(paste0("training_sets\\raw\\pseudoabs_train_B_Q",i,"_thinned_env_",selectedbins,".RDS")) %>% mutate(pos = 0) %>% select(-id),
            readRDS(paste0("training_sets\\raw\\train_B_Q",i,"_thinned_env_",selectedbins,".RDS")) %>% mutate(pos = 1) %>% select(-id)) %>% as.data.frame %>%
    saveRDS(paste0("training_sets\\training_coords_B_Q", i, ".RDS"))
}

# Combine presences and pseudoabsences into test sets

for (i in 1:4){
  df <- bind_rows(
    test_A %>% filter(Q == paste0("Q",i)) %>% select(X, Y, pos),
    pseudoabs_test_A %>% filter(Q == paste0("Q",i)) %>% rename(pos = pr_ab) %>% select(X, Y, pos) # Add in test pseudoabsences
  ) 
  
  png(paste0("plots\\resampling\\test\\test_A_Q", i, ".png"), width = 10, height = 10, units = "in", res = 600)
  plot(base_map, col = "gray95", main = paste0("test set A, Q",i," (pos = green, pseudoabs = red)"))
  points(df %>% filter(pos == 1) %>% pull(X),
         df %>% filter(pos == 1) %>% pull(Y), 
         pch=16, cex=0.2, col = "green3")
  points(df %>% filter(pos == 0) %>% pull(X),
         df %>% filter(pos == 0) %>% pull(Y), 
         pch=16, cex=0.2, col = "firebrick1")
  dev.off()
  
  df %>% saveRDS(paste0("training_sets\\test_coords_A_Q", i, ".RDS"))
  
  df <- bind_rows(
    test_B %>% filter(Q == paste0("Q",i)) %>% select(X, Y, pos),
    pseudoabs_test_B %>% filter(Q == paste0("Q",i)) %>% rename(pos = pr_ab) %>% select(X, Y, pos) # Add in test pseudoabsences
  ) 
  
  png(paste0("plots\\resampling\\test\\test_B_Q", i, ".png"), width = 10, height = 10, units = "in", res = 600)
  plot(base_map, col = "gray95", main = paste0("test set B, Q",i," (pos = green, pseudoabs = red)"))
  points(df %>% filter(pos == 1) %>% pull(X),
         df %>% filter(pos == 1) %>% pull(Y), 
         pch=16, cex=0.2, col = "green3")
  points(df %>% filter(pos == 0) %>% pull(X),
         df %>% filter(pos == 0) %>% pull(Y), 
         pch=16, cex=0.2, col = "firebrick1")
  dev.off()
  
  df %>% saveRDS(paste0("training_sets\\test_coords_B_Q", i, ".RDS"))
}

# Plot all in single plot

for (i in 1:4){
  
  tr <- readRDS(paste0("training_sets/training_coords_A_Q", i, ".RDS"))
  ts <- readRDS(paste0("training_sets/test_coords_A_Q", i, ".RDS"))
  
  png(paste0("training_sets/plot_A_Q", i, ".png"), width = 10, height = 10, units = "in", res = 600)
  plot(base_map, col = "gray95", main = paste0("dataset A, Q",i," (pos = green, pseudoabs = red, faded = test)"))
  points(tr %>% filter(pos == 1) %>% pull(X),
         tr %>% filter(pos == 1) %>% pull(Y), 
         pch=16, cex=0.3, col = rgb(red = 16/255, green = 165/255, blue = 53/255, alpha = 1))
  points(tr %>% filter(pos == 0) %>% pull(X),
         tr %>% filter(pos == 0) %>% pull(Y), 
         pch=16, cex=0.3, col = rgb(red = 184/255, green = 22/255, blue = 22/255, alpha = 1))
  points(ts %>% filter(pos == 1) %>% pull(X),
         ts %>% filter(pos == 1) %>% pull(Y), 
         pch=16, cex=0.3, col = rgb(red = 16/255, green = 165/255, blue = 53/255, alpha = 0.4))
  points(ts %>% filter(pos == 0) %>% pull(X),
         ts %>% filter(pos == 0) %>% pull(Y), 
         pch=16, cex=0.3, col = rgb(red = 133/255, green = 81/255, blue = 81/255, alpha = 0.4))
  dev.off()
  
  tr <- readRDS(paste0("training_sets/training_coords_B_Q", i, ".RDS"))
  ts <- readRDS(paste0("training_sets/test_coords_B_Q", i, ".RDS"))
  
  png(paste0("training_sets/plot_B_Q", i, ".png"), width = 10, height = 10, units = "in", res = 600)
  plot(base_map, col = "gray95", main = paste0("dataset B, Q",i," (pos = green, pseudoabs = red, faded = test)"))
  points(tr %>% filter(pos == 1) %>% pull(X),
         tr %>% filter(pos == 1) %>% pull(Y), 
         pch=16, cex=0.3, col = rgb(red = 16/255, green = 165/255, blue = 53/255, alpha = 1))
  points(tr %>% filter(pos == 0) %>% pull(X),
         tr %>% filter(pos == 0) %>% pull(Y), 
         pch=16, cex=0.3, col = rgb(red = 184/255, green = 22/255, blue = 22/255, alpha = 1))
  points(ts %>% filter(pos == 1) %>% pull(X),
         ts %>% filter(pos == 1) %>% pull(Y), 
         pch=16, cex=0.3, col = rgb(red = 16/255, green = 165/255, blue = 53/255, alpha = 0.4))
  points(ts %>% filter(pos == 0) %>% pull(X),
         ts %>% filter(pos == 0) %>% pull(Y), 
         pch=16, cex=0.3, col = rgb(red = 133/255, green = 81/255, blue = 81/255, alpha = 0.4))
  dev.off()
  
}
