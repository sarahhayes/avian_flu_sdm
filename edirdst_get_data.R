# We want a script that will bring in the ebird weekly abundance layers
# for all the species that are found in Europe. 

# We want to then save these layers 

rm(list = ls())

library(ebirdst)
library(terra)
library(sf)
library(fields)
library(rnaturalearth)

#remotes::install_github("ropensci/rnaturalearthhires")

## look at the species table provided

sp_df <- ebirdst_runs

# download data - we will look at the herring gull 
path <- ebirdst_download(species = "hergul")

# load relative abundance raster stack with 52 layers, one for each week
abd <- load_raster(path = path, resolution = "lr")
