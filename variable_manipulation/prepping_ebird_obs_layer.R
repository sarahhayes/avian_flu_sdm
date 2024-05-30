###############################################
# Filter the eBird sampling datasets with auk #
###############################################

library(tidyverse)
library(magrittr)
library(auk)
library(sf)
library(terra)

# Prepare background maps for plots
library(rnaturalearth)
spdf_world <- ne_download(scale = 110, type = "countries")
plot(spdf_world)

setwd("~/avian_flu_sdm_fork/")
auk_set_ebd_path("data_offline/EBD/sampling")

# Filter to Europe study area; takes up to 10 mins to conduct filtering
base_map <- terra::rast("output/euro_rast_latlong.tif")
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

# Filter to country codes of Africa
auk_sampling("ebd_sampling_relDec-2023.txt") %>%
  auk_country(country = c("DZ","AO","BJ","BW","BF","BI","CM","CV","CF","TD","KM","CG","CI","DJ","EG","GQ","ER","ET","GA","GM","GH","GN","GW","KE","LS","LR","LY","MG","MW","ML","MR","MU","YT","MA","MZ","NA","NE","NG","RE","RW","ST","SN","SC","SL","SO","ZA","SD","SZ","TZ","TG","TN","UG","EH","ZM","ZW","BS")) %>%
  auk_date(date = c("2005-10-07", "2023-06-30")) %>% # limit to dates of flu sampling records
  auk_filter(file = "output/ebd_sampling_afri.txt")

#################################
# Process data layer for Europe #
#################################

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

base_map_10k <- terra::rast("output/euro_rast_10k.tif")

point_data <- st_as_sf(x = ebird_eur, 
                       coords = c("longitude", "latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform("EPSG:3035")

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

png(paste0("plots//ebird_records.png"), width = 8, height = 8, units = "in", res = 600)
plot(app(points_rast, function(x) log(x+1)), main = NULL, axes = TRUE, box = TRUE)
text(x=7650000, y=3000000, "log(eBird sightings)", srt=90, cex=1, xpd=NA, pos=4)
plot(spdf_world[1], add = TRUE, color = NA)
dev.off()

terra::writeRaster(points_rast, "variable_manipulation/variable_outputs/ebird_records.tif", overwrite = TRUE)

###############################
# Process data layer for Asia #
###############################

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

base_map_10k <- terra::rast("output/asia_russia_rast.tif")

point_data <- st_as_sf(x = ebird_asia, 
                       coords = c("longitude", "latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs(base_map_10k))

points_sp <- sf::as_Spatial(point_data)

points_rast <- terra::rasterize(x = vect(points_sp), 
                                y = base_map_10k, 
                                field = ebird_asia %>% pull(n), 
                                fun = "sum")

# Assign 0 to any terrestrial cells with no observations
points_rast <- mask(subst(points_rast %>% as("SpatRaster"), NA, 0), 
                    base_map_10k %>% as("SpatRaster"), 
                    maskvalue=NA)

plot(points_rast)

png(paste0("plots//ebird_records_asia.png"), width = 15, height = 10, units = "in", res = 600)
plot(app(points_rast, function(x) log(x+1)), main = "log(eBird records, lat/long)")
plot(spdf_world[1], add = TRUE, color = NA)
dev.off()

terra::writeRaster(points_rast, "variable_manipulation/variable_outputs/ebird_records_asia.tif", overwrite = TRUE)

###################################
# Process data layer for Americas #  (takes a long time as ~19GB of data)
###################################

# ebird_amer <- read_sampling("data_offline/EBD/sampling/ebd_sampling_amer.txt")
# ebird_amer %>% head %>% as.data.frame
# ebird_amer %>% nrow
# ebird_amer %<>% distinct(latitude, longitude, observation_date, observer_id)
# ebird_amer %>% nrow()
# ebird_amer %<>% group_by(latitude, longitude) %>% tally

ebird_amer <- readRDS(file = "data_offline/EBD/sampling/ebird_amer.rds")

# x %>% distinct(latitude, longitude, observation_date, time_observations_started) %>% nrow()
# x %>% group_by(latitude, longitude, observation_date, time_observations_started) %>% tally %>% arrange(-n)

# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

base_map_10k <- terra::rast("output/americas_rast.tif")

point_data <- st_as_sf(x = ebird_amer, 
                       coords = c("longitude", "latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs(base_map_10k))

points_sp <- sf::as_Spatial(point_data)

points_rast <- terra::rasterize(x = vect(points_sp), 
                                y = base_map_10k, 
                                field = ebird_amer %>% pull(n), 
                                fun = "sum")

# Assign 0 to any terrestrial cells with no observations
points_rast <- mask(subst(points_rast %>% as("SpatRaster"), NA, 0), 
                    base_map_10k %>% as("SpatRaster"), 
                    maskvalue=NA)

plot(points_rast)

png(paste0("plots//ebird_records_amer.png"), width = 15, height = 10, units = "in", res = 600)
plot(app(points_rast, function(x) log(x+1)), main = "log(eBird records, lat/long)")
plot(spdf_world[1], add = TRUE, color = NA)
dev.off()

terra::writeRaster(points_rast, "variable_manipulation/variable_outputs/ebird_records_amer.tif", overwrite = TRUE)


#################################
# Plot eBird records for Africa #
#################################

ebird_afri <- read_sampling("data_offline/EBD/sampling/ebd_sampling_afri.txt")
ebird_afri %>% head %>% as.data.frame
ebird_afri %>% nrow
ebird_afri %<>% distinct(latitude, longitude, observation_date, observer_id)
ebird_afri %>% nrow()
ebird_afri %<>% group_by(latitude, longitude) %>% tally

# x %>% distinct(latitude, longitude, observation_date, time_observations_started) %>% nrow()
# x %>% group_by(latitude, longitude, observation_date, time_observations_started) %>% tally %>% arrange(-n)

# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

blank_latlong <- terra::rast(crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", res = 1)

point_data <- st_as_sf(x = ebird_afri, 
                       coords = c("longitude", "latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

points_sp <- sf::as_Spatial(point_data)

points_rast <- terra::rasterize(x = vect(points_sp), 
                                y = blank_latlong, 
                                field = ebird_afri %>% pull(n), 
                                fun = "sum")

plot(points_rast)

png(paste0("plots//ebird_records_afri.png"), width = 15, height = 10, units = "in", res = 600)
plot(app(points_rast, function(x) log(x+1)), main = "log(eBird records, lat/long)", xlim = c(-40,80), ylim = c(-50,50))
plot(spdf_world[1], add = TRUE, color = NA)
dev.off()