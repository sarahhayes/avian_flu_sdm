# trailling natural earth via this webpage
# https://docs.ropensci.org/rnaturalearth/articles/rnaturalearth.html

library(rnaturalearth)
library(sp)

spdf_world <- ne_download( scale = 110, type = 'countries' )

sp::plot(ne_countries(type = "countries", scale = "small"))
