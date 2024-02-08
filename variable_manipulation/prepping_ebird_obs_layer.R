library(tidyverse)
library(magrittr)
library(auk)
library(sf)
library(terra)

setwd("~/avian_flu_sdm/")
auk_set_ebd_path("data_offline/EBD/sampling")
base_map <- raster::raster("output/euro_rast_latlong.tif")

# Takes up to 10 mins to conduct filtering:
auk_sampling("ebd_sampling_relDec-2023.txt.gz") %>%
  auk_bbox(bbox = base_map) %>% # limit to area of flu sampling records
  auk_date(date = c("2005-10-07", "2023-06-30")) %>% # limit to dates of flu sampling records
  auk_filter(file = "output/ebd_sampling_eur.txt")

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

points_rast <- terra::rasterize(x = points_sp, 
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

