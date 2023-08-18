## Script to calculate distance to coast. 
## Based on example that is seen here:
## https://dominicroye.github.io/en/2019/calculating-the-distance-to-the-sea-in-r/

rm(list =ls())
library(tidyverse)
library(terra)
library(sf)

#euromap <- terra::vect("output/euro_map.shp")
# euromap <- st_read("output/euro_map.shp")
## need a bigger map to avoid the boundary being classified as the sea

euromap <- st_read("output/large_map.shp")

# look at one/two countries only whilst working it out
# euromap <- euromap[which(euromap$NAME_ENGL %in% c("France", "Spain")),]
# euromap <- euromap[which(euromap$NAME_ENGL %in% c("Denmark")),]

plot(euromap)
 
# euromap
# The information that is printed when this is read in shows that the extent is the
# same as the one we have set for our base raster - see bounding box. 
# (2000000, 6000000, 1000000, 5500000)

cell_size <- c(1000,1000)
crs <- "epsg:3035"

##
grid_centers <- st_make_grid(euromap, cellsize = cell_size, crs = crs, what = "centers") ## Need to look 
grid_centers # the extent is different to our raster but the XY look like they are the 
# midpoints of the boxes, which I think is what we need.

#grid_corners <- st_make_grid(euromap, cellsize = cell_size, what = "corners") ## Need to look 
#grid_corners # this seems to keep the extent the same as we want it to be. 

## at this more closely and ensure that our grid will be the same as the blank raster 
## we are using
plot(euromap, max.plot = 1)
#plot(grid_centers, add = T) # don't run this! If want to check, look at subsection of area.

# below is incredibly slow
#tictoc::tic()
#eurogrid_centers <- st_intersection(grid_centers, euromap) 
#tictoc::toc()
# the net should now just be the area of land 
#plot(eurogrid_centers)

#transform from polygon shape to line
# However, we just want the outline so have to remove the internal lines
class(euromap)

eurounion <- st_union(euromap$geom)
plot(eurounion)
eurounion

euroline <- st_cast(eurounion, "MULTILINESTRING")
euroline
plot(euroline)

## work out the intersection using the line
## use intersects instead of intersection

tictoc::tic()
eurogrid <- st_intersects(grid_centers, eurounion, sparse = F)
tictoc::toc()

#g3 <- grid_centers[which(eurogrid_3[,1] == T)]
#plot(g3, col = "green")

g4 <- grid_centers[eurogrid]
#plot(g4, col = "pink")

#calculation of the distance between the coast and our points
tictoc::tic()
dist <- st_distance(euroline, g4)
tictoc::toc()

#distance with unit in meters
head(dist[1,])


## plotting the distances
df <- data.frame(dist = as.vector(dist)/1000,
                 st_coordinates(g4))

#structure
str(df)

#colors 
library(RColorBrewer)
col_dist <- brewer.pal(11, "RdGy")

#df_small <- df[1000000: 1200000,]

# The larger map
ggplot(df, aes(X, Y, fill = dist)) + #variables
  geom_tile()+ #geometry
  scale_fill_gradientn(colours = rev(col_dist))+ #colors for plotting the distance
  labs(fill = "Distance (km)")+ #legend name
  theme_void()+ #map theme
  theme(legend.position = "bottom") #legend position

# save one that is only shown to the size we want
png("plots/distance_to_coast_map.png")
ggplot(df, aes(X, Y, fill = dist)) + #variables
  geom_tile()+ #geometry
  scale_fill_gradientn(colours = rev(col_dist))+ #colors for plotting the distance
  labs(fill = "Distance (km)")+ #legend name
  theme_void()+ #map theme
  xlim(2600000, 7000000) +
  ylim(1500000, 6400000) +
  theme(legend.position = "bottom") #legend position
dev.off()  
  

#save some of these files as it takes a long time to run

# write.csv(df, "variable_manipulation/variable_outputs/interim_dist_to_coast.csv")


# now we want to extract the points we need

head(df) # we can use these XY coords to do the extraction

# read in the reference raster
shp_for_points <- terra::rast("output/euro_rast.tif")
shp_for_points # shows us the extent of our reference area.

res_dist_to_coast <- df %>%
  dplyr::filter(X > 2600000 & X < 7000000) %>%
  dplyr::filter(Y >1500000 & Y < 6400000)

# The results are the same size as the results dfs from the other variables

head(res_dist_to_coast)

res_dist_to_coast <- rename(res_dist_to_coast, "dist_to_coast_km" = "dist")
range(res_dist_to_coast$dist_to_coast)

# write.csv(res_dist_to_coast, "variable_manipulation/variable_outputs/dist_to_coast_output.csv")



#########################################################################################
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

grid_centers <- st_make_grid(iceland, cellsize = 5000, what = "centers")
#grid_corners <- st_make_grid(iceland, cellsize = 5000, what = "corners")

grid_centers
#grid_corners
#our fishnet with the extension of Iceland
plot(grid_centers)
#plot(grid_corners, col = "blue", add = T)

tictoc::tic()
grid <- st_intersection(grid_centers, iceland)   
tictoc::toc()

# looking for a way to speed things up. 
tictoc::tic()
grid2 <- st_intersects(grid_centers, iceland, sparse = F)
tictoc::toc()

grid3 <- st_intersects(grid_centers, iceland, sparse = T)

grid4 <- grid_centers[grid2]

plot(grid4)

#our fishnet now
plot(grid)
plot(grid2)

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
