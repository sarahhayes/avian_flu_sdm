# Map of abundance for a selected species for use in poster

rm(list = ls())

library(ebirdst)
library(terra)
library(sf)
library(fields)
library(rnaturalearth)
library(raster)
#remotes::install_github("ropensci/rnaturalearthhires")

## look at the species table provided

sp_df <- ebirdst_runs

#################################################################################
# First do an example with one species only to work out the steps. 

# download data - we will look at the herring gull 
# path <- ebirdst_download(species = "hergul")
# path <- ebirdst_download(species = "mallar3")
path <- ebirdst_download(species = "comter")

# load relative abundance raster stack with 52 layers, one for each week
# abd <- load_raster(path = path, resolution = "lr") # currently set as lr = low resolution for speed
 abd <- load_raster(path = path, resolution = "hr")

# set the crs we want to use
crs <- "epsg:3035"

#change the projection of the bird data
abd_prj <- terra::project(x = abd, y = crs, method = "near") # ideally (according to help file) would use our template 
# Spatraster of the area we want as y. 

# Now we want to crop
e <- terra::ext(2000000, 6000000, 1000000, 5500000)
crop_abd <- crop(x = abd_prj, y = e )
plot(crop_abd[[1]], axes = F)

# change the colours
# quantiles of non-zero values
v <- values(crop_abd)
v <- v[!is.na(v) & v > 0]
range(v)
# bins <- quantile(v, seq(0, 1, by = 0.1)) # splits into the 10% quantiles
#bins <- c(1e-10, 1e-6, 1e-3, 1e-2, 1.5e-2, 1e-1, 1.5e-1, 1,100) # bins if low res
bins <- c(1e-3, 1e-2, 1e-1, 1.25e-1, 1.5e-1, 1.75e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 8e-1, 1, 2, 100) # bins if high res

pal <- abundance_palette(length(bins))

winter <- mean(crop_abd[[1:13]])
summer <- mean(crop_abd[[26:39]])
plot(winter, breaks = bins, col = pal, axes = FALSE)
plot(summer, breaks = bins, col = pal, axes = F)

v_wint <- values(winter)
v_wint <- v_wint[!is.na(v_wint) & v_wint > 0]
range(v_wint)

v_aut <- values(summer)
v_aut <- v_aut[!is.na(v_aut) & v_aut > 0]
range(v_aut)


# plot with grey background
plot(crop_abd[[26]], axes = F, col = "#DCDCDC", legend = F)
plot(winter, breaks = bins, col = pal, axes = FALSE, add = T)

plot(crop_abd[[26]], axes = F, col = "#DCDCDC", legend = F)
plot(summer, breaks = bins, col = pal, axes = F, add = T)

#par(mar = c(1,1,2,1))

png("z_old_or_unused_scripts/tern_wint.png", height=480, width=580)
plot(crop_abd[[26]], axes = F, col = "#DCDCDC", legend = F)
plot(winter, breaks = bins, col = pal, axes = FALSE, add = T)
dev.off()

png("z_old_or_unused_scripts/tern_summ.png", height=480, width=580)
plot(crop_abd[[26]], axes = F, col = "#DCDCDC", legend = F)
plot(summer, breaks = bins, col = pal, axes = F, add = T)
dev.off()




