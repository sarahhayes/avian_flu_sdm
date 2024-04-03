## 22/02/2024
## Creating the avian flu database of domestic animal cases
## Plan is to overlay this on our risk map 

rm(list = ls())

library(readxl)
library(tidyverse)
library(terra)
library(RColorBrewer)

## bring in the different data sets and combine to one large data set. 

## first do fao 
## 
## FAO

# when downloading these data, go to the epidemiology page, select the virus, 
# select the date (01/01/9170 - 30/06/2023), select domestic and then download raw data from tab in upper right corner
# use the full data as want to include asia as well 
# I have selected those just listed as domestic and not included 'captive'. I did look at captive on the fao page but these were more
# rare/zoo species. First will focus on just domestic outbreaks. 

fao_data  <- read.csv("data/flu_data/raw_data/fao_domestic_cases_raw_data.csv")

table(fao_data$Region)

colnames(fao_data)

# rename the date columns
fao_data <-  dplyr::rename(fao_data, any_of(c(observation.date = "Observation.date..dd.mm.yyyy.",
                                              report.date = "Report.date..dd.mm.yyyy.")))

# There are 1657 rows with no entries for observation date 
nrow(fao_data[which(fao_data$observation.date == ""),])

no_obs_data <- fao_data[which(fao_data$observation.date == ""),]
table(no_obs_data$Region)
# Very few of these are in europe so will remove them

fao_data <- fao_data[which(fao_data$observation.date != ""),]

# we only want to include domestic species this time
table(fao_data$Species) # selected just those listed as 'domestic' from website. So these should all be domestic. 

# Remove the entries that are mammals or environmental

unwanted_sp_fao <- c("Canine", "Cats", "Dogs", "Donkey", "Ferret", "Goats", "House Crow", 
                     "Swine", "Unspecified Env. Sample")

fao_data <- fao_data %>% filter(!Species %in% unwanted_sp_fao)
table(fao_data$Species)

# Rename some of them with very long names

fao_data[which(fao_data$Species == "Xpeacock"), "Species"] <- "Peacock"
fao_data[which(fao_data$Species == "Phasianus Colchicus (Common Pheasant) - Phasanidae"), "Species"] <- "Pheasant"

#remove bits in brackets to simplify the names
fao_data$Species <- str_replace(fao_data$Species, " \\s*\\([^\\)]+\\)", "")

fao_data[which(fao_data$Species == "Common Quail:coturnix Coturnix(Phasianidae)"), "Species"] <- "Common Quail"
fao_data[which(fao_data$Species == "Muscovy Duck:cairina Moschata(Anatidae)"), "Species"] <- "Muscovy Duck"
fao_data[which(fao_data$Species == "Phasianidae:Phasianidae (Incognita)(Phasianidae)"), "Species"] <- "Phasianidae"
table(fao_data$Species)


fao_data_trim <- dplyr::select(fao_data, all_of(c("Latitude", "Longitude", "observation.date", "Species", "Country", "Serotype")))
fao_data_trim$Epi_unit <- NA # for later merge with WAHIS
fao_data_trim$source <- "fao"
fao_data_trim$observation.date <- as.Date(fao_data_trim$observation.date, "%d/%m/%Y")


####### Now do WAHIS

# WAHIS
wahis <- read_excel("data/flu_data/raw_data/WOAH_july_2023.xlsx", sheet = 2)

# contains all diseases and species so filter these down to the ones we want

inf_data <- wahis[which(wahis$disease_eng %in% 
                          c("High pathogenicity avian influenza viruses (poultry) (Inf. with)",
                            "Influenza A virus (Inf. with)",
                            "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")),]


## first select those labelled as wild = F

inf_a_hpai <- dplyr::filter(inf_data, is_wild == F)
table(inf_a_hpai$disease_eng)

poss_cases <- dplyr::filter(inf_a_hpai, disease_eng == "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")
table(poss_cases$Epi_unit)
table(poss_cases$Species)
table(inf_a_hpai$Epi_unit)
#
#Removing zoo and cage animals. Will save further filtering until later

inf_a_hpai <- dplyr::filter(inf_a_hpai, !Epi_unit %in% c("Zoo", "Cage"))

table(inf_a_hpai$region) # all regions represented

colnames(inf_a_hpai)
inf_a_hpai_trim <- dplyr::select(inf_a_hpai, all_of(c("Latitude", "Longitude", "Outbreak_start_date", "Species", 
                                                           "country", "sero_sub_genotype_eng", "Epi_unit")))
inf_a_hpai_trim$source <- "woah"

# We now have the trimmed data sets for all the sources. 
# first need all the names to match so we can rbind

colnames(fao_data_trim)
colnames(inf_a_hpai_trim)
inf_a_hpai_trim <- rename(inf_a_hpai_trim, observation.date = Outbreak_start_date, 
                          Country = country, Serotype = sero_sub_genotype_eng)



all_ai_data <- rbind(fao_data_trim, inf_a_hpai_trim)

## check for any NA values
sum(is.na(all_ai_data$Latitude))
sum(is.na(all_ai_data$Longitude))

# None without location data

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
range(point_data_df$Y) # checking that we have labelled the columns correctly

## Now select only those rows which are in the area that we want

euro_shp <- terra::vect("output/euro_map.shp")
ext(euro_shp)

ai_data_prj_area <- point_data_df %>%
  dplyr::filter(X >= 2600000 & X <= 7000000) %>%
  dplyr::filter(Y >= 1500000 & Y <= 6400000)


## We need to remove mammals from these data

table(ai_data_prj_area$Species)

unwanted_sp <- c("Cats", "Mustelidae (dom)", "Swine")

nrow(ai_data_prj_area[which(ai_data_prj_area$Species %in% unwanted_sp),])

# 61 mammal entries that we are removing. 

## Remove from the data so just birds

ai_pos_birds <- 
  ai_data_prj_area[-which(ai_data_prj_area$Species %in% unwanted_sp),]

## Below plots are all serotypes - not just HPAI

# png("plots/fao_pos_domestic.png", width = 480, height = 480)
# plot(euro_shp, main = "FAO")
# plot(ai_pos_birds[which(ai_pos_birds$source == "fao"),], add = T, 
#      pch = 18, legend = T, col = "red", cex = 0.5)
# dev.off()
# 
# png("plots/woah_pos_domestic.png", width = 480, height = 480)
# plot(euro_shp, main = "WOAH")
# plot(ai_pos_birds[which(ai_pos_birds$source == "woah"),], add = T, 
#      pch = 18, legend = T, col = "blue", cex = 0.5)
# dev.off()

# Remove any of the data that are labelled LPAI
table(ai_pos_birds$Serotype)
ai_pos_birds$serotype_HN <- ai_pos_birds$Serotype
# some of these state that they are LPAI so remove those. 
ai_pos_birds <- ai_pos_birds %>% filter(!str_detect(Serotype, "LPAI"))
# Now remove extra info
ai_pos_birds$serotype_HN <-  gsub("HPAI","",as.character(ai_pos_birds$serotype_HN))
ai_pos_birds$serotype_HN <-  gsub("LPAI","",as.character(ai_pos_birds$serotype_HN))
ai_pos_birds$serotype_HN<-toupper(ai_pos_birds$serotype_HN) # so they match
ai_pos_birds$serotype_HN <- trimws(ai_pos_birds$serotype_HN) # trim any white space in the entries to aid matching
table(ai_pos_birds$serotype_HN)

table(ai_pos_birds$source, ai_pos_birds$serotype_HN)

# remove the ones where we don't have a serotype as only 7 of them
ai_pos_birds <- dplyr::filter(ai_pos_birds, serotype_HN != "")



## Next I want to look at the duplicates that are present in the data 

library(data.table)
dt_pos <- data.table(ai_pos_birds)

table(dt_pos$Country) # a couple of different spellings for a number of the countries
dt_pos[which(dt_pos$Country == "Faeroe Islands"), "Country"] <- "Faroe Islands"
dt_pos[which(dt_pos$Country == "Moldova  Republic of"), "Country"] <- "Moldova"
dt_pos[which(dt_pos$Country == "Russian Federation"), "Country"] <- "Russia"
dt_pos[which(dt_pos$Country == "Türkiye"), "Country"] <- "Turkey"
dt_pos[which(dt_pos$Country == "U.K. of Great Britain and Northern Ireland"), "Country"] <- "United Kingdom"

no_duplicates_date_and_loc <- unique(dt_pos, by = c("X", "Y", "observation.date"))
no_duplicates_all <- unique(dt_pos, by = c("X", "Y", "observation.date", "source", "Species", "Country"))
no_dups_all_except_source <- unique(dt_pos, by = c("X", "Y", "observation.date", "Species", "Country"))
no_dups_date_loc_serotype <- unique(dt_pos, by = c("X", "Y", "observation.date", "serotype_HN"))

# The middle two are the same, suggesting that there are none that are the same for all details except source. 
# The final one (no_dups_date_loc_serotype) contains 60 more than the one filtered just on date and location.
# On inspection, these are indeed listed as different serotypes in the same location and same date. There are a 
# number of entries in Romania where fao have listed as "H5N1 HPAI" and woah just have "H5".. Imagine these are the same? 

## No duplicates_all will still have some that are the same place and time but different species. However, when looking
## at these, there are instances where we have fao and woah reporting a same location and date but might say, 
## for example, "Birds", "Turkey" or "Unspecified Birds". 
## I think we basically want to know if a location has AI and so no_duplicates_date_and_loc is probably the one to use,
## given the notes above on no_dups_date_loc_serotype. 

# easier to re-classify the dataframe so can use the existing code
ai_pos_birds <- no_duplicates_date_and_loc

table(ai_pos_birds$serotype_HN)

hpai <- ai_pos_birds %>% mutate(geometry = as.character(geometry))
# write.csv(hpai, "data/flu_data/prepped_data/hpai_pos_domestic_birds.csv")

# Add a year and month and year_month to the data

hpai$year <- lubridate::year(hpai$observation.date)
hpai$month <- lubridate::month(hpai$observation.date)
hpai$month_year <-  format(as.Date(hpai$observation.date), "%Y-%m")
hpai$week <- format(hpai$observation.date, "%Y Week %W")
hpai$week_num <- lubridate::isoweek(hpai$observation.date)

range(hpai$observation.date)

# Now to plot the time series by serotype 

serotype_data <- as.data.frame(table(hpai$serotype_HN))

## Need to ensure that have all the weeks represented as currently misses out those with no cases 

## To make sure we include all weeks in a plot, even those without any data, we need a dataframe with 
## all the weeks and months,

#Start on 3rd Jan 2004 as that's a Monday and no data before middle of the year
start.plot <- c("2005-01-03", "2023-12-31")
week_test <- as.character(seq(as.Date(start.plot[1]), as.Date(start.plot[2]), by="weeks"))

all_dates <- as.data.frame(matrix(nrow = length(week_test), ncol = 1))
all_dates[,1]  <- week_test
colnames(all_dates) <- "week_from_start"
all_dates$date <- as.Date(all_dates$week_from_start, "%Y-%m-%d")

all_dates$year <- lubridate::year(all_dates$date)
all_dates$month <- lubridate::month(all_dates$date)
all_dates$month_year <- format(all_dates$date, "%Y-%m")
all_dates$week <- format(all_dates$date, "%Y Week %W")
all_dates$week_num <- lubridate::isoweek(all_dates$date)

# For quite a number of the subtypes there are only a few entries. 
# For the line plot by subtype, only use those that have >10 entries

serotype_data_ten_or_more <- serotype_data[which(serotype_data$Freq > 10),]
serotype_data_ten_or_more <- serotype_data_ten_or_more[which(serotype_data_ten_or_more$Var1 != ""),]

subtype_plot_data <- hpai[which(hpai$serotype_HN %in% serotype_data_ten_or_more$Var1),]
subtype_plot_data <- dplyr::select(subtype_plot_data, c("year", "month", "month_year", "week","week_num", "serotype_HN"))

# by week 
subtype_counts_weekyear <- subtype_plot_data %>%
  group_by(week, serotype_HN) %>% count()

# merge with the weekly counts
all_dates_subtype_data <- left_join(all_dates, subtype_counts_weekyear)
all_dates_subtype_data$year_of_study <- all_dates_subtype_data$year - 2005
all_dates_subtype_data$week_seq <- all_dates_subtype_data$week_num + (52*all_dates_subtype_data$year_of_study)

#write.csv(all_dates_subtype_data, "data/flu_data/prepped_data/time_series_domestic_hpai.csv", row.names = F)

weekbreaks <- all_dates_subtype_data %>% filter(grepl('Week 01', week))
weekbreaks <- unique(weekbreaks[,"week_seq"])
weekbreaks


yearlabs <- c("2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", 
              "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")

par(mar = c(3,3,3,3))

ggplot(data = all_dates_subtype_data, 
       aes(x = week_seq, y = n, col = serotype_HN, group = serotype_HN)) +
  geom_line(lwd = 0.9) + 
  scale_color_brewer(palette="Paired", na.translate = F) + 
  theme(axis.text.x = element_text(angle=90, margin = margin(t = .2, unit = "cm"), 
                                   face = "bold", size = 12, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 12),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 14),
        plot.background = element_rect(fill = "light blue"))+
  #  panel.background = element_rect(fill = "white", colour = "grey50")) + 
  labs(colour = "Subtype", y = "Number of weekly cases", x = "Year", size = 18) +
  scale_x_continuous(breaks = weekbreaks, labels = yearlabs)
#ggsave("plots/subtypes/serotype_line_by_week_hpai_domestic_birds.png")


## If want to repeat by month, will need to do another table to join with 
## so don't miss any months.

subtype_counts_year <- subtype_plot_data %>%
  group_by(year, serotype_HN) %>% count()

table(subtype_counts_year$year) #NB there are no data from 2012 in here. 

# stacked bar plot coloured by serotype

subtype_counts_year$factor_HN <- as.factor(subtype_counts_year$serotype_HN)

ggplot(data = subtype_plot_data, aes(x = year, fill = serotype_HN))+
  geom_bar()
#ggsave("plots/subtypes/serotype_bar_by_year_hpai_domestic_birds.png")

######

### Now want to create a raster at 10km resolution where we label a pixel as positive if it has had a domestic oputbreak
### Then will want to do with counts of numbers of outbreaks 
### and then perhaps look at it in the context of the poultry density? 
### For this latter option, might need to 'cut out' our wild bird risk map based on presence of poultry and overlay? 
### Will need some thought

euro_rast <- terra::rast("output/euro_rast_10k.tif")
plot(euro_rast)

# extract points
p <- cbind(hpai$X, hpai$Y)

plot(euro_rast)
points(p, pch = 18, cex = 0.5)

# rasterize points as a matrix
# points_rast <- rasterize(p, euro_rast, fun=sum)
points_rast <- rasterize(p, euro_rast, value = 1) # this one just marks as present or not. Doesn't give an idea of number. 

plot(points_rast)
points_rast


counts_by_geom <- ai_pos_birds %>% count(geometry)
counts_by_geom_date <- ai_pos_birds %>% 
  count(geometry, observation.date)
  

# from Liam's code prepping eBird data

point_ai_data <- st_as_sf(x = hpai, 
                       coords = c("X", "Y"), 
                       crs = "EPSG:3035")

points_sp <- sf::as_Spatial(point_ai_data)

# bit below not right yet.
points_rast_sum <- terra::rasterize(x = vect(points_sp), 
                                y = euro_rast, 
                                field = hpai %>% pull(n), 
                                fun = "sum")


# points_rast <- terra::rasterize(x = vect(points_sp), 
#                                 y = base_map_10k, 
#                                 field = ebird_eur %>% pull(n), 
#                                 fun = "sum")


# Assign 0 to any terrestrial cells with no observations
points_rast <- mask(subst(points_rast %>% as("SpatRaster"), NA, 0), 
                    base_map_10k %>% as("SpatRaster"), 
                    maskvalue=NA)

plot(points_rast)


# rasterize points as a SpatVector
pv <- vect(p)
xv <- rasterize(pv, r, fun=sum)

