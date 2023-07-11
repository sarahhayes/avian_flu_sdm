## Script to calculate distance to coast. 
## Based on example that is seen here:
## https://dominicroye.github.io/en/2019/calculating-the-distance-to-the-sea-in-r/

rm(list =ls())
library(tidyverse)
library(terra)
library(sf)

#euromap <- terra::vect("output/euro_map.shp")
euromap <- st_read("output/euro_map.shp")
euromap

##
grid <- st_make_grid(euromap, cellsize = 1000, what = "centers") ## Need to look 
## at this more closely and ensure that our grid will be the same as the blank raster 
## we are using
plot(euromap, max.plot = 1)
plot(eurogrid, add = T)

tictoc::tic()
eurogrid <- st_intersection(grid, euromap) 
tictoc::toc()
# the net should now just be the area of land 
plot(eurogrid)

#transform from polygon shape to line
euroline <- st_cast(euromap, "MULTILINESTRING")
plot(euroline)

#calculation of the distance between the coast and our points
dist <- st_distance(euroline, eurogrid)

#distance with unit in meters
head(dist[1,])


## plotting the distances
df <- data.frame(dist = as.vector(dist)/1000,
                 st_coordinates(eurogrid))

#structure
str(df)

#colors 
col_dist <- brewer.pal(11, "RdGy")


ggplot(df, aes(X, Y, fill = dist))+ #variables
  geom_tile()+ #geometry
  scale_fill_gradientn(colours = rev(col_dist))+ #colors for plotting the distance
  labs(fill = "Distance (km)")+ #legend name
  theme_void()+ #map theme
  theme(legend.position = "bottom") #legend position



## Example from webpage
library(rnaturalearth)
library(sf)
library(raster)
library(tidyverse)
library(RColorBrewer)


iceland <- ne_countries(scale = 10, country = "Iceland", returnclass = "sf")

#info of our spatial vector object
iceland

iceland <- st_transform(iceland, 3055)

grid <- st_make_grid(iceland, cellsize = 5000, what = "centers")

#our fishnet with the extension of Iceland
plot(grid)

grid <- st_intersection(grid, iceland)   

#our fishnet now
plot(grid)

#transform Iceland from polygon shape to line
iceland <- st_cast(iceland, "MULTILINESTRING")

#calculation of the distance between the coast and our points
dist <- st_distance(iceland, grid)

#distance with unit in meters
head(dist[1,])


## plotting the distances
df <- data.frame(dist = as.vector(dist)/1000,
                 st_coordinates(grid))

#structure
str(df)

#colors 
col_dist <- brewer.pal(11, "RdGy")


ggplot(df, aes(X, Y, fill = dist))+ #variables
  geom_tile()+ #geometry
  scale_fill_gradientn(colours = rev(col_dist))+ #colors for plotting the distance
  labs(fill = "Distance (km)")+ #legend name
  theme_void()+ #map theme
  theme(legend.position = "bottom") #legend position
