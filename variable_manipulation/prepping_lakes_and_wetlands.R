## Exploring and using the GLWD data
## 09/06/2023

rm(list = ls())

library(terra)

glwd_rast <- rast("data/variables/lakes_and_wetlands/glwd_3/w001001x.adf")

glwd_rast
# resolution is supposed to be 30 arc seconds = 1km
# This is in degrees and 0.008 degrees is approx 1km. So all as should be currently.

# change projection and crop 
crs <- "epsg:3035"
# euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) 
# 
# # Create a blank raster with appropriate projection and extent
# blank_raster <- rast(crs=crs, extent=euro_ext, res = 1000) # the unit for epsg is metres
# # thus I'm hoping that this makes it 1km res
# blank_raster

blank_raster <- terra::rast("output/euro_rast.tif") # 1k res
blank_raster <- terra::rast("output/euro_rast_10k.tif") # 1k res
blank_raster
#terra::xyFromCell(blank_3035, 1) # coordinates of the centre of the first cell
plot(blank_raster)

glwd_crop_prj <- terra::project(x = glwd_rast, y = blank_raster, method = "near") #use nearest neighbour as the 
# numbers are categories rather than values
glwd_crop_prj # we can see that this keeps the crs of the blank raster
plot(glwd_crop_prj) 
head(glwd_crop_prj)

#plot over the top of the shapefile

euromap <- terra::vect(x = "output/euro_map.shp")

plot(euromap)
plot(glwd_crop_prj, add = T, axes = F)

# can we project and keep res of underlying glwd raster

# glwd_rast
# tryit <- terra::project(x = glwd_rast, y = crs, method = "near")
# just projecting is crazy slow. 

# blank_raster_glwdres <- rast(crs=crs, extent=euro_ext, res = res(glwd_rast)) # the unit for epsg is metres

# glwd_crop_prj_2 <- terra::project(x = glwd_rast, y = blank_raster_glwdres, method = "near") 
# oddly this is also very slow. 
# glwd_crop_prj_2  

# for now will stick with 1km res raster and change if needed. 


# info on the key is available in the data documentation pdf in the data file
# Not sure we want to differentiate between the different types? Re-code so that they are all 
# either water or not-water?? 

glwd_crop_prj

table(values(glwd_crop_prj))
# This has some values of 9 in it. We don't want to include these as they are intermittent wetland
# We just want to include permanent wetlands.

glwd_crop_prj_combo <- terra::subst(glwd_crop_prj, from = c(1:8), to = c(rep(99,8)),
                                    others = NA)
table(values(glwd_crop_prj_combo))
sum(table(values(glwd_crop_prj_combo)))
sum(table(values(glwd_crop_prj)))
# when using 1k res diff of 9185 between these as there should be
# diff of 99 with 10k res

plot(glwd_crop_prj)
plot(glwd_crop_prj_combo, col = "white", add = T)
# this shows that the level 9 are NA in the new raster as you can still see the 
# green dots from the original before we changed it

plot(euromap)
plot(glwd_crop_prj_combo, add = T, axes = F, col = "blue", legend = F)

# save a copy of this plot
#pdf("plots/inland_water.pdf", height = 5, width = 7)
#plot(euromap)
#plot(glwd_crop_prj_combo, add = T, axes = F, col = "blue", legend = F)
#dev.off()

## Next step is to try and work out distance to nearest water for each cell of the blank raster. 

## We want a small subsection to try this with first. 

small_extent <- terra::ext(5950000, 6000000, 5350000, 5400000) 

# Create a blank raster with appropriate projection and extent
#blank_small <- rast(crs=crs, extent=small_extent, res = 1000) # the unit for epsg is metres
blank_small <- rast(crs=crs, extent=small_extent, res = 10000) # the unit for epsg is metres


glwd_small <- terra::project(x = glwd_rast, y = blank_small, method = "near") 
glwd_small # we can see that this keeps the crs of the blank raster
plot(glwd_small) 
glwd_small <- terra::subst(glwd_small, from = c(1:8), to = c(rep(99,8)),
                           others = NA)
plot(glwd_small, col = "blue")


# using example from here https://gis.stackexchange.com/questions/360516/measure-shortest-distance-between-raster-cell-to-raster-cell-of-another-raster-i

# Get the cell centres of the NA cells in raster one, and the non-NA cells
# in raster two. The only difference in these lines apart from r1 and r2 everywhere is a little ! in the second section:
# The two rasters are blank_small and glwd_small

values(blank_small) <- 1
head(blank_small)

p1 = as.data.frame(blank_small,xy=TRUE)
p1 = p1[,1:2]

p2 = as.data.frame(glwd_small, xy=TRUE)
p2 = p2[!is.na(p2[,3]),1:2]

# plot(blank_small, col="grey")
# plot(glwd_small, add=TRUE, col=c("white","black"),legend=FALSE)
# points(p1$x, p1$y)
# points(p2$x, p2$y, pch=3)

# Now we are set up for knnx.dist.
# install.packages("FNN")

dnear = FNN::knnx.dist(p2, p1, k=1)
# Since k=1 that object is a one-column matrix. We can fill the missing values in r1 with that column:
  
blank_small[] = dnear[,1] # assigning the distances to the raster

# And plot:
  
plot(blank_small) # shows the distances. The 0 seem to be in the right place
points(p2$x, p2$y, pch=3) # add the points from the glwd raster that were non NA values
## this seems to work! 
plot(glwd_small, add = T, col = "blue")


## Now we can do it for the full raster

values(blank_raster) <- 1

p1 = as.data.frame(blank_raster,xy=TRUE)
p1 = p1[,1:2]

p2 = as.data.frame(glwd_crop_prj_combo, xy=TRUE)
p2 = p2[!is.na(p2[,3]),1:2]

# Now we are set up for knnx.dist.
dnear = FNN::knnx.dist(p2, p1, k=1)

blank_raster[] = dnear[,1] # assigning the distances to the raster
#blank_raster[] = dnear[,1]/1000 # assigning the distances to the raster as km

# And plot:
plot(blank_raster) # shows the distances. This is harder to see if correct. 
# Think we want to plot 0 as a particular colour 
# points(p2$x, p2$y, pch=3) # add the points from the glwd raster that were non NA values
max(dnear)

# Because this method uses all NA cells in the GLWD data, it includes all the oceans etc
# read in the ref raster again, but with a different name

#euro_rast <- terra::rast("output/euro_rast.tif")
euro_rast <- terra::rast("output/euro_rast_10k.tif")

blank_raster
euro_rast

masked <- terra::mask(blank_raster, euro_rast )
plot(masked)

#png("plots/distance_to_water_m.png", width = 600, height = 400)
bp_500 <- c(1000,5000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 
            100000, 200000, 500000, 1000000)
bp_500 <- bp_500*10 # if using 10k
plot(masked, breaks = c(0, bp_500), 
     plg = list(title = "Distance in m"),
     col= c("white", viridis::viridis(14)), background = "light blue")
#plot(glwd_crop_prj_combo, add = T, axes = F, col = "red", legend = F)
#dev.off()

# ideally re-do this in km. 
#pdf("plots/distance_to_water_masked.pdf", width = 7, height = 5)
plot(masked, breaks = c(0, bp_500), 
     plg = list(title = "Distance in m"),
     pax=list(side=1:2, cex.axis = 0.8),
     col= c("white", viridis::viridis(14)), background = "grey")
#dev.off()

### Extracting the values 
# make a points object using the centre of each pixel from the reference raster
points_3035 <- terra::as.points(euro_rast)
points_3035

tictoc::tic()
dist_to_inland_water_res <- terra::extract(masked, points_3035, method = "simple", xy = T)
tictoc::toc()


head(dist_to_inland_water_res)
range(dist_to_inland_water_res$layer)
dist_to_inland_water_res$layer[400:440]

dist_to_inland_water_res <- rename(dist_to_inland_water_res, "dist_to_water" = "layer")
head(dist_to_inland_water_res)

# write.csv(dist_to_inland_water_res,
#          "variable_manipulation/variable_outputs/dist_to_water_output.csv",
#          row.names = F)


# write.csv(dist_to_inland_water_res,
#          "variable_manipulation/variable_outputs/dist_to_water_output_10kres.csv",
#          row.names = F)


## looking at how best to plot and exploring options

plot(blank_raster,  col= c("white", grDevices::rainbow(50)))
# try and plot so only zeros are show in white. 
bp <- seq(1, max(dnear), length.out = 10)
plot(blank_raster, breaks = c(0, bp), col= c("white", grDevices::terrain.colors(10)))
points(p2$x, p2$y, pch=3, cex = 0.01) # add the points from the glwd raster that were non NA values

plot(blank_raster, breaks = c(0, bp), col= c("white", viridis::viridis(10)))

bp_km <- seq(1, max(dnear)/1000, length.out = 10)

blank_raster_km <- blank_raster/1000
plot(blank_raster_km, breaks = c(0, bp_km), col= c("white", viridis::viridis(10)))


bp_km_1200 <- c(1,10,20,50,100, 200, 500, 1000, 1200)
plot(blank_raster_km, breaks = c(0, bp_km_1200), 
     plg = list(title = "Distance in km"),
     col= c("white", viridis::magma(8)))
#save this plot

#pdf("plots/distance_to_water.pdf", width = 7, height = 5)
plot(blank_raster_km, breaks = c(0, bp_km_1200), 
     plg = list(title = "Distance in km"),
     col= c("white", viridis::viridis(8)))
#dev.off()
#table(values(blank_raster_km))

#writeRaster(blank_raster, "output/distance_to_inland_water.tif", overwrite=TRUE)

# want to try and make the background white/blue

# make a copy of the raster
trial <- glwd_crop_prj_combo2
#bring in the europe shapefile
euromap <- terra::vect(x = "output/euro_map.shp")

plot(trial)
plot(euromap)
plot(blank_raster_km)
eurorast <- rasterize(euromap, blank_raster_km)
plot(eurorast)


masked <- terra::mask(blank_raster_km, eurorast)
plot(masked)

bp_km_500 <- c(1,5, 10, 15,20, 25, 30, 40, 50,100, 200, 500)
plot(masked, breaks = c(0, bp_km_500), 
     plg = list(title = "Distance in km"),
     col= c("white", viridis::viridis(13)), background = "light blue")

pdf("plots/distance_to_water_masked.pdf", width = 7, height = 5)
plot(masked, breaks = c(0, bp_km_500), 
     plg = list(title = "Distance in km"),
     pax=list(side=1:2, cex.axis = 0.8),
     col= c("white", viridis::viridis(13)), background = "grey")
dev.off()

