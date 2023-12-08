### Look at the bv-brc data to produce some summaries

rm(list = ls())
library(tidyverse)
library(zoo)

# This is wild and captive wild
pos_data <- read.csv("data/flu_data/prepped_data/bvbrc_pos_europe.csv")
neg_data <- read.csv("data/flu_data/prepped_data/bvbrc_neg_europe.csv")

# look at a few things to check as expected
table(pos_data$Collection.Country)
table(neg_data$Collection.Country)

table(pos_data$Host.Natural.State)
table(neg_data$Host.Natural.State)

table(pos_data$Host.Group)
table(neg_data$Host.Group)

table(pos_data$Pathogen.Test.Result)
table(neg_data$Pathogen.Test.Result)


## Subtypes from the positive cases

bvbrc_subtype <- as.data.frame(table(pos_data$Subtype))

bvbrc_subtype_country <- as.data.frame(table(pos_data$Collection.Country, 
                                             pos_data$Subtype)) %>%
  filter(Freq != 0) %>%
  filter(Var2 != "") %>%
  filter(!Var2 %in% c("HxNx", "HxNy", "HXNX"))

#write.csv(bvbrc_subtype, "output/bvbrc_summaries/bvbrc_subtypes.csv", row.names = F)
#write.csv(bvbrc_subtype_country, "output/bvbrc_summaries/bvbrc_subtypes_countries.csv", row.names = F)

# counts per country

pos_data_counts <- pos_data %>%
  dplyr::count(as.factor(Collection.Country),as.factor(Collection.Year),.drop=FALSE) %>%
  rename(Country = "as.factor(Collection.Country)", Year = "as.factor(Collection.Year)")

neg_data_counts <- neg_data %>%
  count(as.factor(Collection.Country),as.factor(Collection.Year),.drop=FALSE) %>%
  rename(Country = "as.factor(Collection.Country)", Year = "as.factor(Collection.Year)")


ggplot(pos_data_counts, aes(x = Year, y = n, 
       fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "BV-BRC Positive samples from wild birds in Europe")


ggplot(neg_data_counts, aes(x = Year, y = n, 
                                   fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "Negative samples from wild birds in Europe")


## just by country and not by year
pos_data %>%
  dplyr::count(as.factor(Collection.Country), .drop=FALSE) %>%
  rename(Country = "as.factor(Collection.Country)") %>%
  ggplot(., aes(x = Country, y = n, fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", 
       title = "Positive samples from wild birds in Europe 
       (BV-BRC)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
#ggsave("plots/by_country_pos_bvbrc.png" )#, width = 5, height = 5)

neg_data %>%
  dplyr::count(as.factor(Collection.Country), .drop=FALSE) %>%
  rename(Country = "as.factor(Collection.Country)") %>%
  ggplot(., aes(x = Country, y = n, fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "Negative samples from wild  
       birds in Europe (BV-BRC)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
# ggsave("plots/by_country_neg_bvbrc.png")


## plots of time series of the data
## plot a time series of the cases. 
pos_data$date <- strptime(pos_data$Collection.Date,
                                       format = "%Y-%m-%d")

pos_data$month_year <- format(pos_data$date, format = "%m/%Y")
pos_data$month <- format(pos_data$date, format = "%m")

pos_counts <- pos_data %>% group_by(month_year) %>% count()
pos_counts$date <- as.Date(zoo::as.yearmon(pos_counts$month_year, "%m/%Y"))

ggplot(pos_counts, aes(x = date, y = n)) + 
  geom_line( col = "red", linewidth = 1) + 
  labs(title = "Positive cases in wild birds in Europe 
       (BV-BRC)", 
       x = "Date", y = "Count") + 
  theme_bw() + 
  xlim(c(as.Date("2003-01-01"), as.Date("2023-01-01")))

sum(pos_data_counts$n) # matches the number we should have 

# ggsave("plots/timeline_bvbrc_pos.png")

## look at entries per month
month_counts <- as.data.frame(table(pos_data$month))
month_counts$Var1 <- as.character(month_counts$Var1)
ggplot(data = month_counts, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", col = "orange", fill = "orange") + 
  labs(x = "Month", title = "Monthly positive counts bvbrc data") 

# now negative
neg_data$date <- strptime(neg_data$Collection.Date,
                                 format = "%Y-%m-%d")

neg_data$month_year <- format(neg_data$date, format = "%m/%Y")
neg_data$month <- format(neg_data$date, format = "%m")


neg_counts <- neg_data %>% group_by(month_year) %>% count()
neg_counts$date <- as.Date(zoo::as.yearmon(neg_counts$month_year, "%m/%Y"))


ggplot(neg_counts, aes(x = date, y = n)) + 
  geom_line( col = "red", size = 1) + 
  labs(title = "Negative samples in wild birds in Europe 
       (BV-BRC)", 
       x = "Date", y = "Count") + 
  theme_bw() + 
  xlim(c(as.Date("2003-01-01"), as.Date("2023-01-01")))

# ggsave("plots/timeline_bvbrc_neg.png")

# plot of monthly counts of negative cases
month_counts_neg <- as.data.frame(table(neg_data$month))
month_counts_neg$Var1 <- as.character(month_counts_neg$Var1)
ggplot(data = month_counts_neg, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", col = "orangered3", fill = "orangered3", ) + 
  labs(x = "Month", title = "Monthly negative counts bvbrc data")


### Look at the species that are included in the data 

species_table <- as.data.frame(table(pos_data$Host.Species))
colnames(species_table) <- c("species", "freq")

### READ BELOW PARAGRAPH ---------------------------------------------

# to add the order/family/genus need to run the script in 'add_species_info.R'
# some parts need checking as run - so switch to that and run manually rather than using source

#write.csv(species_table, "output/bvbrc_summaries/species_bvbrc_pos.csv", row.names = F)

# with the negatives, we mainly want the counts for the species in which we have positives 
# so that we could work out % positive.

neg_species <- as.data.frame(table(neg_data$Host.Species))
colnames(neg_species) <- c("species", "freq_neg") 
neg_species$species <- as.character(neg_species$species)

pos_neg_sp <- left_join(species_table, neg_species)
# write.csv(pos_neg_sp, "output/bvbrc_summaries/species_bvbrc_inc_neg_counts_of_pos_sp.csv")

# Next part is to plot the data
# Prepare the maps
# (Numerous options of maps to plot. Have looked at others in plotting_europe.R)


library(terra)
zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
plot(zipmap)
crs <- "epsg:3035"
#euro_ext <- terra::ext(2000000, 9000000, 1000000, 9000000)

# set the extent as the same as the euro shapefil
euro_shp <- terra::vect("output/euro_map.shp")
ext(euro_shp)
euro_ext <- ext(euro_shp)

# change projection and extent. 
# using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs) 
euro_map_crop <- terra::crop(euro_map, euro_ext)
plot(euro_map)
plot(euro_map_crop)

# Transform the data to spatial points
pos_data_for_points <- pos_data[,c("Collection.Longitude", "Collection.Latitude")]
pts_pos <- terra::vect(pos_data_for_points, geom=c("Collection.Longitude","Collection.Latitude"),
                   crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
pts_pos <- terra::project(pts_pos,  "epsg:3035")

neg_data_for_points <- neg_data[,c("Collection.Longitude", "Collection.Latitude")]
pts_neg <- terra::vect(neg_data_for_points, geom=c("Collection.Longitude", "Collection.Latitude"),
                       crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
pts_neg <- terra::project(pts_neg,  "epsg:3035")

dev.off()

#pdf("plots/bvbrc_pos_neg_map.pdf", width = 10, height = 6)
par(mfrow = c(1,2))
plot(euro_map_crop, col = "white", background = "azure2", main = "BV-BRC positive")
plot(pts_pos, add = T, col = "red", pch = 18)
plot(euro_map_crop, col = "white", background = "azure2", main = "BV-BRC negative")
plot(pts_neg, add = T, col = "blue", pch = 18)
#dev.off()
