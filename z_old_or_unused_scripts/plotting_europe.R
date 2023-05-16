#rm(list = ls())

library(utils)
library(rgdal)

# #unzip("K:/2012 - IPD - Policy Maps/fe_2007_us_zcta500.zip")
# unzip("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_LB_2020_3035.shp.zip")
# unzip("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_LB_2020_3035.shp.zip")
# 
# #zipmap <- readOGR(dsn = "K:/2012 - IPD - Policy Maps/fe_2007_us_zcta500", layer = "myZIPmap" )
# zipmap <- readOGR(dsn = "C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_LB_2020_3035.shp.zip",
#                   layer = "CNTR_LB_2020_3035" )
# 
# plot(zipmap)
# 
# 
# unzip("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_BN_03M_2020_3035.shp.zip")
# zipmap <- readOGR(dsn = "C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_BN_03M_2020_3035.shp.zip",
#                   layer = "CNTR_BN_03M_2020_3035")
# plot(zipmap)
# 
# unzip("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_BN_03M_2020_3857.shp.zip")
# zipmap <- readOGR(dsn = "C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_BN_03M_2020_3857.shp.zip",
#                   layer = "CNTR_BN_03M_2020_3857")
# plot(zipmap)
# 
# unzip("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_BN_03M_2020_4326.shp.zip")
# zipmap <- readOGR(dsn = "C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_BN_03M_2020_4326.shp.zip",
#                   layer = "CNTR_BN_03M_2020_4326")
# plot(zipmap)
# 
# 
# unzip("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_RG_03M_2020_3035.shp.zip")
# zipmap <- readOGR(dsn = "C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_RG_03M_2020_3035.shp.zip",
#                   layer = "CNTR_RG_03M_2020_3035")
# plot(zipmap)
# 
# 
# unzip("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_RG_03M_2020_3857.shp.zip")
# zipmap <- readOGR(dsn = "C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_RG_03M_2020_3857.shp.zip",
#                   layer = "CNTR_RG_03M_2020_3857")
# plot(zipmap)

#unzip("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip")
# can use the above to unzip but that unpacks it into the current location. 
# Alternative is to unzip into the desired location first 
zipmap <- readOGR(dsn = "C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                  layer = "CNTR_RG_03M_2020_4326")
plot(zipmap)


euro_map <- zipmap[zipmap$EU_STAT == "T",]
plot(euro_map)

zipmap$ISO3_CODE
zipmap$CNTR_ID

euro_countries <- c("AL", "AD", "AM", "AT", "AZ", "BY", "BE", "BA", "BG", 
                    "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "GE", "DE", 
                    "GR", "HU", "IC", "IR", "IT", "KZ", "XK", "LV", "LI", 
                    "LT", "LU", "MT", "MD", "MC", "ME", "NL", "MK", "NO", 
                    "PL", "PT", "RO", "RU", "SM", "RS", "SK", "SI", "ES", 
                    "SE", "CH", "UA", "UK", "VA", "TR", "IS") # adding iceland for now

euro_map <- zipmap[zipmap$CNTR_ID %in% euro_countries, ]
plot(euro_map)

euro_countries_no_russia <- c("AL", "AD", "AM", "AT", "AZ", "BY", "BE", "BA", "BG", 
                    "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "GE", "DE", 
                    "GR", "HU", "IC", "IR", "IT", "KZ", "XK", "LV", "LI", 
                    "LT", "LU", "MT", "MD", "MC", "ME", "NL", "MK", "NO", 
                    "PL", "PT", "RO", "SM", "RS", "SK", "SI", "ES", 
                    "SE", "CH", "UA", "UK", "VA", "TR")

euro_map_no_rus <- zipmap[zipmap$CNTR_ID %in% euro_countries_no_russia,]
plot(euro_map_no_rus)

euro_countries_west <- c("AL", "AD", "AT", "BY", "BE", "BA", "BG",
                         "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "DE", 
                         "GR", "HU", "IC", "IE", "IS", "IT", "XK", "LV", "LI", 
                         "LT", "LU", "MT", "MD", "MC", "ME", "NL", "MK", "NO", 
                         "PL", "PT", "RO", "SM", "RS", "SK", "SI", "ES", 
                         "SE", "CH", "UA", "UK", "VA")

euro_map_west <- zipmap[zipmap$CNTR_ID %in% euro_countries_west,]
plot(euro_map_west)
plot(zipmap)

## Crop to the desired extent, then plot
library(raster)
out <- raster::crop(euro_map_west, extent(-25, 40, 20, 70))
plot(out, col="khaki", bg="azure2")

extent(euro_map_west)
extent(out)

par(mar = c(1,1,1,1))
plot(out, col="khaki", bg="azure2")

#xy <- mydf[,c(1,2)]
xy_neg <- neg_data[,c("Collection.Longitude", "Collection.Latitude")]

# neg_data_europe <- neg_data[which(!is.na(neg_data$Collection.Latitude)), 
#                             c("Collection.Longitude", "Collection.Latitude")]
# xy_neg <- neg_data_europe[complete.cases(neg_data_europe),c("Collection.Longitude", "Collection.Latitude")]


#spdf <- SpatialPointsDataFrame(coords = xy, data = mydf,
#                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
spdf_neg <- SpatialPointsDataFrame(coords = xy_neg, data = neg_data_europe, 
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
plot(spdf_neg, add = T)

xy_pos <- pos_data_europe[,c("Collection.Longitude", "Collection.Latitude")]
spdf_pos <- SpatialPointsDataFrame(coords = xy_pos, data = pos_data_europe, 
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# if have run empres i first look
pos_data <- empres_europe_pos

pos_data_europe <- pos_data[which(!is.na(pos_data$Latitude)), 
                            c("Longitude", "Latitude")]
xy_pos <- pos_data_europe[complete.cases(pos_data_europe),c("Longitude", "Latitude")]
spdf_pos <- SpatialPointsDataFrame(coords = xy_pos, data = pos_data_europe, 
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))




# pos_data_europe <- pos_data[which(!is.na(pos_data$Collection.Latitude)), 
#                             c("Collection.Longitude", "Collection.Latitude")]
# xy_pos <- pos_data_europe[complete.cases(pos_data_europe),c("Collection.Longitude", "Collection.Latitude")]
# spdf_pos <- SpatialPointsDataFrame(coords = xy_pos, data = pos_data_europe, 
#                                    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
# 

### 
crop_euro_map <- crop(euro_map, extent(-25, 70, 0 ,70))
plot(crop_euro_map, col = "white", bg = "azure2")
plot(spdf_neg, add = T, pch = 20)
plot(spdf_pos, add = T, col = "red", pch = 16)
