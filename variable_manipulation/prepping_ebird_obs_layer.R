##############################################
# Filter the eBird sampling dataset with auk #
##############################################

library(tidyverse)
library(magrittr)
library(auk)
library(sf)
library(terra)

setwd("~/avian_flu_sdm/")
auk_set_ebd_path("data_offline/EBD/sampling")
base_map <- raster::raster("output/euro_rast_latlong.tif")

# Filter to Europe study area; takes up to 10 mins to conduct filtering
auk_sampling("ebd_sampling_relDec-2023.txt") %>%
  auk_bbox(bbox = base_map) %>% # limit to area of flu sampling records
  auk_date(date = c("2005-10-07", "2023-06-30")) %>% # limit to dates of flu sampling records
  auk_filter(file = "output/ebd_sampling_eur.txt")

# Filter to country codes of the Americas
auk_sampling("ebd_sampling_relDec-2023.txt") %>%
  auk_country(country = c("AI","AG","AR","AW","BS","BB","BZ","BM","BO","BR","VG","CA","KY","CL","CO","CR","CU","DM","DO","EC","SV","FK","GF","GL","GD","GP","GT","GY","HT","HN","JM","MQ","MX","MS","NI","PA","PY","PE","PR","KN","LC","PM","VC","SR","TT","TC","US","UY","VE","VI","BV","SH","GS")) %>% 
  auk_date(date = c("2005-10-07", "2023-06-30")) %>% # limit to dates of flu sampling records
  auk_filter(file = "output/ebd_sampling_amer.txt")

# Filter to country codes of Asia
auk_sampling("ebd_sampling_relDec-2023.txt") %>%
  auk_country(country = c("AF","AM","AZ","BH","BD","BT","BN","KH","CN","CX","CC","CY","GE","IN","ID","IR","IQ","IL","JP","JO","KZ","KP","KR","KW","KG","LA","LB","MY","MV","MN","MM","NP","OM","PK","PH","QA","RU","SA","SG","LK","SY","TW","TJ","TH","TR","TM","AE","UZ","VN","YE","HK","MO","IO","TF","HM")) %>% 
  auk_date(date = c("2005-10-07", "2023-06-30")) %>% # limit to dates of flu sampling records
  auk_filter(file = "output/ebd_sampling_asia.txt")


ebird_eur <- read_sampling("data_offline/EBD/sampling/ebd_sampling_eur.txt")
ebird_eur %>% head %>% as.data.frame
ebird_eur %>% nrow
ebird_eur %<>% distinct(latitude, longitude, observation_date, observer_id)
ebird_eur %>% nrow()
ebird_eur %<>% group_by(latitude, longitude) %>% tally

# x %>% distinct(latitude, longitude, observation_date, time_observations_started) %>% nrow()
# x %>% group_by(latitude, longitude, observation_date, time_observations_started) %>% tally %>% arrange(-n)

# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

point_data <- st_as_sf(x = ebird_eur, 
                       coords = c("longitude", "latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform("EPSG:3035")

base_map_10k <- raster::raster("output/euro_rast_10k.tif")

points_sp <- sf::as_Spatial(point_data)

points_rast <- terra::rasterize(x = vect(points_sp), 
                                y = base_map_10k, 
                                field = ebird_eur %>% pull(n), 
                                fun = "sum")

# Assign 0 to any terrestrial cells with no observations
points_rast <- mask(subst(points_rast %>% as("SpatRaster"), NA, 0), 
                    base_map_10k %>% as("SpatRaster"), 
                    maskvalue=NA)

plot(points_rast)

png(paste0("plots//ebird_records.png"), width = 10, height = 10, units = "in", res = 600)
plot(app(points_rast, function(x) log(x+1)), main = "log(eBird records)")
dev.off()

terra::writeRaster(points_rast, "variable_manipulation/variable_outputs/ebird_records.tif", overwrite = TRUE)

########################################
# Plot eBird records for Asia/Americas #
########################################

library(rnaturalearth)
spdf_world <- ne_download(scale = 110, type = "countries")
plot(spdf_world)

americas_map <- spdf_world %>% subset(., REGION_UN == "Americas")
plot(americas_map)
asia_map <- spdf_world %>% subset(., REGION_UN == "Asia")
plot(asia_map)

ebird_asia <- read_sampling("data_offline/EBD/sampling/ebd_sampling_asia.txt")
ebird_asia %>% head %>% as.data.frame
ebird_asia %>% nrow
ebird_asia %<>% distinct(latitude, longitude, observation_date, observer_id)
ebird_asia %>% nrow()
ebird_asia %<>% group_by(latitude, longitude) %>% tally

# x %>% distinct(latitude, longitude, observation_date, time_observations_started) %>% nrow()
# x %>% group_by(latitude, longitude, observation_date, time_observations_started) %>% tally %>% arrange(-n)

# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

point_data <- st_as_sf(x = ebird_asia, 
                       coords = c("longitude", "latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

blank_latlong <- terra::rast(crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", res = 1)

points_sp <- sf::as_Spatial(point_data)

points_rast <- terra::rasterize(x = vect(points_sp), 
                                y = blank_latlong, 
                                field = ebird_asia %>% pull(n), 
                                fun = "sum")

plot(points_rast)

png(paste0("plots//ebird_records_asia.png"), width = 15, height = 10, units = "in", res = 600)
plot(app(points_rast, function(x) log(x+1)), main = "log(eBird records, lat/long)", xlim = c(10,180), ylim = c(-35,80))
plot(spdf_world[1], add = TRUE, color = NA)
dev.off()


ebird_amer <- read_sampling("data_offline/EBD/sampling/ebd_sampling_amer.txt")
ebird_amer %>% head %>% as.data.frame
ebird_amer %>% nrow
ebird_amer %<>% distinct(latitude, longitude, observation_date, observer_id)
ebird_amer %>% nrow()
ebird_amer %<>% group_by(latitude, longitude) %>% tally

# x %>% distinct(latitude, longitude, observation_date, time_observations_started) %>% nrow()
# x %>% group_by(latitude, longitude, observation_date, time_observations_started) %>% tally %>% arrange(-n)

# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

point_data <- st_as_sf(x = ebird_amer, 
                       coords = c("longitude", "latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

blank_latlong <- terra::rast(crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", res = 1)

points_sp <- sf::as_Spatial(point_data)

points_rast <- terra::rasterize(x = vect(points_sp), 
                                y = blank_latlong, 
                                field = ebird_amer %>% pull(n), 
                                fun = "sum")

plot(points_rast)

png(paste0("plots//ebird_records_amer.png"), width = 15, height = 10, units = "in", res = 600)
plot(app(points_rast, function(x) log(x+1)), main = "log(eBird records, lat/long)")
plot(spdf_world[1], add = TRUE, color = NA)
dev.off()

