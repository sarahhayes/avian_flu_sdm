# trailling natural earth via this webpage
# https://docs.ropensci.org/rnaturalearth/articles/rnaturalearth.html

library(rnaturalearth)
library(sp)

spdf_world <- ne_download( scale = 110, type = 'countries' )
plot(spdf_world)
sp::plot(ne_countries(type = "countries", scale = "small"))

spdf_europe <- ne_download(scale = 110, type = "countries",  
                           returnclass = "sf") 
plot(spdf_europe)

#%>%
  dplyr::filter(CONTINENT == "Europe")
plot(spdf_europe)


# spdf_europe <- ne_download(scale = 110, type)


eur_land <- ne_download(scale = 50, type = "countries", category = "cultural",
                       returnclass = "sf") %>%
  filter(CONTINENT == "Europe") %>%
  st_set_precision(1e6) %>%
  st_union()

plot(eur_land)
# country lines
# downloaded globally then filtered to north america with st_intersect()
eur_country_lines <- ne_download(scale = 50, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") %>% 
  st_geometry()
plot(eur_country_lines)

ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%
  as.logical() %>%
  {ne_country_lines[.]}