## Creating the avian flu data csv for use within the BART model

rm(list = ls())

library(readxl)
library(tidyverse)
library(terra)

## bring in the different data sets and combine to one large data set. 

## first do fao as it's the biggest
#FAO 
# use the full data as want to include asia as well 
fao_data  <- read.csv("data/flu_data/raw_data/epidemiology-raw-data_202308011345.csv")

table(fao_data$Region)
## filter to Europe and Asia

fao_data <- dplyr::filter(fao_data, Region %in% c("Europe", "Asia"))
colnames(fao_data)

fao_data <-  dplyr::rename(fao_data, observation.date = "Observation.date..dd.mm.yyyy.")
fao_data <-  dplyr::rename(fao_data, report.date = "Report.date..dd.mm.yyyy.")

# There are rows with no entries for observation date so remove these
fao_data <- fao_data[which(fao_data$observation.date != ""),]

# we only want to include wild birds
table(fao_data$Species)

unwanted_sp <- c("Unspecified Mammal", "Stone Marten", 
                 "South American Coati (Nasua Nasua):Procyonidae-Carnivora",
                 "Red Fox", "Polecat", "Polar Fox", "Nyctereutes Viverrinus (Japanese Racoon Dog)",
                 "Mink", "Harbor Seal (Phoca Vitulina):Phocidae-Carnivora", 
                 "Gray Seal (Halichoerus Grypus):Phocidae-Carnivora", "Fox",
                 "European Pine Marten", "Civet", "Caspian Seal (Pusa Caspica)")

fao_data <- fao_data[-which(fao_data$Species %in% unwanted_sp),]

fao_data_trim <- dplyr::select(fao_data, all_of(c("Latitude", "Longitude", "observation.date")))
fao_data_trim$source <- "fao"
fao_data_trim$observation.date <- as.Date(fao_data_trim$observation.date, "%d/%m/%Y")

# WAHIS
wahis <- read_excel("data/flu_data/raw_data/WOAH_july_2023.xlsx", sheet = 2)

# contains all diseases and species so filter these down to the ones we want

inf_data <- wahis[which(wahis$disease_eng %in% 
                          c("High pathogenicity avian influenza viruses (poultry) (Inf. with)",
                            "Influenza A virus (Inf. with)",
                            "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")),]

inf_a_hpai <- dplyr::filter(inf_data, disease_eng == "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)" )
table(inf_a_hpai$region) # all regions represented

inf_a_hpai_trim <- dplyr::select(inf_a_hpai, all_of(c("Latitude", "Longitude", "Outbreak_start_date")))
inf_a_hpai_trim$source <- "woah"


# BVBRC
bvbrc_data <- read.csv("data/flu_data/raw_data/BVBRC_surveillance.csv")

columns_needed <- c("Collection.Date", "Collection.Year", "Collection.Country", 
                    "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                    "Host.Species", "Host.Group", "Host.Natural.State", "Host.Health", 
                    "Subtype", "Strain")

bvbrc_data_slim <- dplyr::select(bvbrc_data, all_of(columns_needed))
table(bvbrc_data_slim$Collection.Country) # still contains global data as required

bvbrc_data_slim <- filter(bvbrc_data_slim, Host.Natural.State == "Wild")

# remove the ones without a year
bvbrc_data_slim <- drop_na(bvbrc_data_slim, Collection.Year)

# There are samples labelled "env" which I assume are environmental so remove these
bvbrc_data_slim <- bvbrc_data_slim[which(bvbrc_data_slim$Host.Species != "Env"),]

# separate positive and negative
pos_data <- filter(bvbrc_data_slim, Pathogen.Test.Result == "Positive")
neg_data <- filter(bvbrc_data_slim, Pathogen.Test.Result == "Negative")

bv_pos_data_trim <- dplyr::select(pos_data, all_of(c("Collection.Date", 
                                                     "Collection.Latitude",
                                                     "Collection.Longitude")))

bv_pos_data_trim$source <- "bvbrc"
bv_pos_data_trim$Collection.Date <- as.Date(bv_pos_data_trim$Collection.Date,
                                            "%Y-%m-%d")


bv_neg_data_trim <- dplyr::select(neg_data, all_of(c("Collection.Date", 
                                                     "Collection.Latitude",
                                                     "Collection.Longitude")))
bv_neg_data_trim$source <- "bvbrc"
bv_neg_data_trim$Collection.Date <- as.Date(bv_neg_data_trim$Collection.Date,
                                            "%Y-%m-%d")

# We now have the trimmed data sets for all the sources. 
# combine them and add a section for pos or neg

# first need all the names to match so we can rbind

colnames(fao_data_trim)
colnames(inf_a_hpai_trim)
colnames(bv_pos_data_trim)
inf_a_hpai_trim <- rename(inf_a_hpai_trim, observation.date = Outbreak_start_date)
bv_pos_data_trim <- rename(bv_pos_data_trim, Latitude = Collection.Latitude,
                           Longitude = Collection.Longitude,
                           observation.date = Collection.Date)

bv_pos_data_trim <- dplyr::select(bv_pos_data_trim, c(Latitude, Longitude, observation.date, source))
colnames(bv_pos_data_trim)

all_ai_data <- rbind(fao_data_trim, inf_a_hpai_trim, bv_pos_data_trim)
all_ai_data$flu <- 1 # denotes a positive


## deal with the negatives

bv_neg_data_trim <- rename(bv_neg_data_trim, Latitude = Collection.Latitude,
                           Longitude = Collection.Longitude,
                           observation.date = Collection.Date)

bv_neg_data_trim <- dplyr::select(bv_neg_data_trim, c(Latitude, Longitude, observation.date, source))
colnames(bv_neg_data_trim)
bv_neg_data_trim$flu <- 0 # denotes negative test

all_ai_data <- rbind(all_ai_data, bv_neg_data_trim)

## check for any NA values
sum(is.na(all_ai_data$Latitude))
sum(is.na(all_ai_data$Longitude))

# remove the rows that are NA for location 
all_ai_data <- drop_na(all_ai_data, Latitude)
sum(is.na(all_ai_data$Latitude))

###################################################################################
## make the points into a raster. 
## Using info from here: 
## # https://gis.stackexchange.com/questions/458522/use-r-to-create-raster-from-data-frame-with-non-uniformly-gridded-points-and-cat

library(sf)
# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

point_data <- st_as_sf(x = all_ai_data, 
                        coords = c("Longitude", "Latitude"),
                        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
class(point_data)

# change the projections
proj_crs <- 3035
point_data <- sf::st_transform(point_data, proj_crs)
#plot(point_data)
point_data # Now a spatial points object in the right projection
class(point_data)
head(point_data)

point_data_df <- point_data %>%
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                     Y = sf::st_coordinates(.)[,2])
head(point_data_df)
st_geometry(point_data)
range(point_data_df$X)
range(point_data_df$Y) # this confirms that we have labelled the columns correctly

## Now select only those rows which are in the area that we want

ai_data_prj_area <- point_data_df %>%
  dplyr::filter(X >= 2600000 & X <= 7000000) %>%
  dplyr::filter(Y >= 1500000 & Y <= 6400000)

# separate pos and neg

pos_ai_data_prj_area <- ai_data_prj_area %>%
  filter(flu == 1)

neg_ai_data_prj_area <- ai_data_prj_area %>%
  filter(flu == 0)

nrow(pos_ai_data_prj_area) + nrow(neg_ai_data_prj_area)
nrow(ai_data_prj_area)

#save these
#write.csv(pos_ai_data_prj_area, "output/avian_flu/pos_points_proj_area_all_sources.csv" )
#write.csv(neg_ai_data_prj_area, "output/avian_flu/neg_points_proj_area_all_sources.csv" )

# remove rows where duplicates of coordinates 
library(data.table)
dt_pos <- data.table(pos_ai_data_prj_area)
#Remove duplicates on specific column
tictoc::tic()
dt_pos <- unique(dt_pos, by = c("X", "Y", "observation.date"))
tictoc::toc()
head(dt_pos)
table(dt_pos$flu)

# repeat for negative
dt_neg <- unique(neg_ai_data_prj_area, by = c("X", "Y", "observation.date"))
head(dt_neg)
table(dt_neg$flu)

#write.csv(dt_pos, "output/avian_flu/pos_points_proj_area_all_sources_duplicates_removed.csv" )
#write.csv(dt_neg, "output/avian_flu/neg_points_proj_area_all_sources_duplicates_removed.csv" )

############################# MAKE RASTER  ############################

# Other manipulations needed to rasterize
table(point_data$flu) # this is pos and neg data

points_sp <- sf::as_Spatial(point_data)
head(points_sp)

points_vect <- terra::vect(points_sp)
head(points_vect)
plot(points_vect) # this is also still the global data set

#read in the raster we are using for the project
euro_rast <- terra::rast("output/euro_rast.tif")
euro_rast <- terra::rast("output/euro_rast_plus.tif")

euro_rast

## use the vector within rasterize
points_rast <- terra::rasterize(points_vect, euro_rast, field = "flu", fun = "max")
# using max because if it has a positive in the cell (1) and a negative (0) we want
# the cell to be classes as a negative
plot(points_rast)
points_rast
head(points_rast)


#read in the map
euro_map <- terra::vect("output/euro_map.shp")

plot(euro_map)
plot(points_rast, add = TRUE, col = "red")

table(terra::values(points_rast))


# Create the csv file for avian flu cases
# make a points object using the centre of each pixel from the ref raster
points_euro_rast <- terra::as.points(euro_rast)
points_euro_rast

tictoc::tic()
flu_res <- terra::extract(points_rast, points_euro_rast, method = "simple", xy = T)
tictoc::toc()

table(flu_res$max)
table(values(points_rast))
euro_rast
nrow(flu_res)
points_rast

# why have we got fewer values in the results table rather than the raster? 
# Looked at it with the 10km raster and wonder if some of them are in the sea/on land too small to be 
# noted on our raster?? 

points_rast
euro_rast
table(values(points_rast))
# we have the higher number in the points raster

points_euro_rast # the points are all within the boundaries of the euro-raster.
# so that's not the issue

# so I think the issue is somewhere in extract or if there were NAs in the euroraster that
# are not in the case data.

euro_nonNA_coords <- as.data.frame(terra::crds(euro_rast, na.rm = T))
cases_nonNA_coords <- as.data.frame(terra::crds(points_rast, na.rm = T))
sum(table(values(points_rast))) - sum(table(flu_res$max))
diff_nas <- setdiff( cases_nonNA_coords, euro_nonNA_coords) # this is the same number
# as are missing between the datasets

# let's turn these into spatial points and plot them 
missing_points <- st_as_sf(x = diff_nas, 
                       coords = c("x", "y"),
                       crs = "3035")
plot(euro_map)
plot(missing_points, add = T, col = "red", pch = 18)





dev.off()
plot(euro_map)
plot(points_rast, xlim= c(4000000, 5000000), ylim = c(3000000, 4000000), 
     col = "hot pink")#, background = "blue")
plot(points_rast, col = "red")
plot(points_rast, col = "red", add = T, axes = F)
plot(euro_rast, col = "White", add = T, axes = F)
#plot(euro_map, add = T, axes = F)

euro_rast
points_rast
euro_rast[1000:2000] # we can see that there are NAs
points_rast[1:100]

# we want to know which points have NA in euro_rast that aren't NA in points rast. 
# if it's a NA in the euro_rast but has a value in the points rast that might explain it

################################################################
## Some trials to aid understanding/visualisation

## Hard to see. Look at with larger raster

#read in the raster we are using for the project
euro_rast_10k <- terra::rast("output/euro_rast_10k.tif")
euro_rast_10k

points_rast_10k <- terra::rasterize(points_vect, euro_rast_10k,"flu", fun = max)
plot(points_rast_10k)
points_rast_10k
head(points_rast_10k)
table(values(points_rast_10k))
# Easier to visualise that it probably is working OK. 
# I think the 1km are just too small to see. 

dev.off()
plot(euro_map)
plot(points_rast_10k, col = "red")
plot(euro_rast_10k, col = "White", add = T, axes = F)



# Create the csv file for avian flu cases
# make a points object using the centre of each pixel from the ref raster
points_euro_rast_10k <- terra::as.points(euro_rast_10k)
points_euro_rast_10k

flu_res_10k <- terra::extract(points_rast_10k, points_euro_rast_10k,
                              method = "simple", xy = T)

table(flu_res_10k$max)
table(values(points_rast_10k))

dev.off()
plot(points_rast_10k, col = "red")
plot(euro_map, add = T, axes = F)



#Perhaps just look at UK? 
GB_ext <- terra::ext(3300000, 3800000, 3100000, 4100000)
GB_crop <- terra::crop(euro_map, GB_ext)
plot(GB_crop)

GB_rast <- terra::crop(points_rast, GB_ext)
plot(GB_crop)
plot(GB_rast, add = T, axes = F)
