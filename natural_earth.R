## setting up maps as per this information
## https://cornelllabofornithology.github.io/ebird-best-practices/intro.html

library(sf)
library(rnaturalearth)
library(dplyr)

# file to save spatial data
gpkg_dir <- "data"
if (!dir.exists(gpkg_dir)) {
  dir.create(gpkg_dir)
}
f_ne <- file.path(gpkg_dir, "gis-data.gpkg")

# download bcrs (bird conservation range boundaries)

# this is in example but we probably don't need them 

# tmp_dir <- normalizePath(tempdir())
# tmp_bcr <- file.path(tmp_dir, "bcr.zip")
# paste0("https://www.birdscanada.org/research/gislab/download/", 
#        "bcr_terrestrial_shape.zip") %>% 
#   download.file(destfile = tmp_bcr)
# unzip(tmp_bcr, exdir = tmp_dir)
# bcr <- file.path(tmp_dir, "BCR_Terrestrial_master_International.shp") %>% 
#   read_sf() %>% 
#   select(bcr_code = BCR, bcr_name = LABEL) %>% 
#   filter(bcr_code == 27)
# # clean up
# list.files(tmp_dir, "bcr", ignore.case = TRUE, full.names = TRUE) %>% 
#   unlink()

# political boundaries
# land border with lakes removed
ne_land <- ne_download(scale = 50, type = "countries", category = "cultural",
                       returnclass = "sf") %>%
  filter(CONTINENT == "North America") %>%
  st_set_precision(1e6) %>%
  st_union()
# country lines
# downloaded globally then filtered to north america with st_intersect()
ne_country_lines <- ne_download(scale = 50, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") %>% 
  st_geometry()
ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%
  as.logical() %>%
  {ne_country_lines[.]}
# states, north america
ne_state_lines <- ne_download(scale = 50, category = "cultural",
                              type = "admin_1_states_provinces_lines",
                              returnclass = "sf") %>%
  filter(adm0_a3 %in% c("USA", "CAN")) %>%
  mutate(iso_a2 = recode(adm0_a3, USA = "US", CAN = "CAN")) %>% 
  select(country = adm0_name, country_code = iso_a2)

# output
unlink(f_ne)
write_sf(ne_land, f_ne, "ne_land")
write_sf(ne_country_lines, f_ne, "ne_country_lines")
write_sf(ne_state_lines, f_ne, "ne_state_lines")
write_sf(bcr, f_ne, "bcr")