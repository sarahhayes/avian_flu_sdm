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

#fao_data  <- read.csv("data/flu_data/raw_data/fao_domestic_cases_raw_data.csv")
fao_data  <- read.csv("data/flu_data/raw_data/epidemiology-raw-data_domestic_apr_2024.csv") # later download from april 2024

table(fao_data$Region)

colnames(fao_data)
table(fao_data$Diagnosis.status)
table(fao_data$Animal.type)

# the rogue entry is misaligned and is in cats in Korea so remove
fao_data <- fao_data %>% dplyr::filter(Animal.type != "Confirmed")

# rename the date columns
fao_data <-  dplyr::rename(fao_data, any_of(c(observation.date = "Observation.date..dd.mm.yyyy.",
                                              report.date = "Report.date..dd.mm.yyyy.")))

# There are 1684 rows with no entries for observation date 
nrow(fao_data[which(fao_data$observation.date == ""),])

no_obs_data <- fao_data[which(fao_data$observation.date == ""),]
table(no_obs_data$Region)
# Very few of these are in europe so will remove them

fao_data <- fao_data[which(fao_data$observation.date != ""),]

# we only want to include domestic species this time
table(fao_data$Species) # selected just those listed as 'domestic' from website. So these should all be domestic. 

# Remove the entries that are mammals or environmental

unwanted_sp_fao <- c("Canine", "Cats", "Cattle", "Dogs", "Donkey", "Ferret", "Goats", "House Crow", 
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
#wahis <- read_excel("data/flu_data/raw_data/WOAH_july_2023.xlsx", sheet = 2)
wahis <- read.csv("data/flu_data/raw_data/WOAH_april_2024.csv") # later download from april 2024

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
table(inf_a_hpai$disease_eng)

#Removing zoo and cage animals. Will save further filtering until later
inf_a_hpai <- dplyr::filter(inf_a_hpai, !Epi_unit %in% c("Zoo", "Cage"))
table(inf_a_hpai$region) # all regions represented
table(inf_a_hpai$disease_eng)

## Alternative approach would be just to use those labelled as poultry These make up the majority of the records
## However, the definition of poultry which is here (https://www.woah.org/en/what-we-do/standards/codes-and-manuals/terrestrial-code-online-access/?id=169&L=1&htmfile=glossaire.htm) 
## means this would exclude backyard flocks. 

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

unwanted_sp <- c("Cats", "Dogs", "Mustelidae (dom)", "Swine")

nrow(ai_data_prj_area[which(ai_data_prj_area$Species %in% unwanted_sp),])

# 63 mammal entries that we are removing. 

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
#ai_pos_birds$serotype_HN <-  gsub("LPAI","",as.character(ai_pos_birds$serotype_HN))
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
# number of entries in Romania where fao have listed as "H5N1 HPAI" and woah just have "H5". Imagine these are the same? 

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


## To make sure we include all weeks in a plot, even those without any data, we need a dataframe with 
## all the weeks and months,
## # Now to plot the time series by serotype 

## start on 3rd Jan 2005 as that's a Monday
week_calendar <- data.frame(date=seq(as.Date("2005-01-03"), as.Date("2024-04-28"), by="day")) %>% 
  mutate(week_num=isoweek(date),year=year(date)) %>%
  group_by(year,week_num) %>% 
  summarise(weekdate=min(date)) 
week_calendar$week_of_study <- seq(1, nrow(week_calendar), by = 1)

# combine with the data
dates <- merge(hpai, week_calendar)

# plot directly
ggplot(dates, aes(x=weekdate)) + geom_histogram(bins = nrow(week_calendar), col = "purple") + 
  labs(y = "Number of weekly cases", x = "Year", size = 18) +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "grey90"))
ggsave("plots/hist_by_week_hpai_only_domestic.png")

# could also make a dataframe of the counts of the number in each week
count_dates <- dates %>% 
  count(weekdate)

weekly_counts <- left_join(week_calendar, count_dates)
weekly_counts[which(is.na(weekly_counts$n)), "n"] <- 0

## Also then can have as a line plot
ggplot(data = weekly_counts, aes(x = weekdate, y = n)) +
  geom_line(lwd = 0.9, col = "purple")  + 
  labs(y = "Number of weekly cases", x = "Year", size = 18) +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "grey90"))
ggsave("plots/lineplot_by_week_hpai_only_domestic.png")

# save for the time series
#write.csv(weekly_counts, "data/flu_data/prepped_data/time_series_domestic_hpai.csv", row.names = F)

# For quite a number of the subtypes there are only a few entries. 
# For the line plot by subtype, only use those that have >10 entries


subtype_plot_data <- dates
(subtype_plot_data)

yearlabs <- c("2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", 
              "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024")

weekbreaks <- week_calendar[which(week_calendar$week_num == 1), "week_of_study"]
weekbreaks <- weekbreaks[["week_of_study"]]


ggplot(subtype_plot_data, aes(x = week_of_study, fill = serotype_HN)) + 
  geom_bar(position = "stack") + 
  labs(y = "Number of weekly cases", x = "Year", size = 18, fill = "Serotype") +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        panel.background = element_rect(fill = "grey90", colour = "grey90"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "white")) +
  scale_x_continuous(breaks = weekbreaks, labels = yearlabs)
ggsave("plots/subtypes/serotype_by_week_hpai_only_domestic.png")


subtype_counts_year <- subtype_plot_data %>%
  group_by(year, serotype_HN) %>% count()

table(subtype_counts_year$year) #NB there are no data from 2012 in here. 

# stacked bar plot coloured by serotype

subtype_counts_year$factor_HN <- as.factor(subtype_counts_year$serotype_HN)

ggplot(data = subtype_plot_data, aes(x = year, fill = serotype_HN))+
  geom_bar()
 ggsave("plots/subtypes/serotype_bar_by_year_hpai_domestic_birds.png")


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
points_rast <- rasterize(p, euro_rast, value = 1) # this one just marks as present or not. Doesn't give an idea of number. 

plot(points_rast)
points_rast

# We want these split into quarters but also by the period over which we split the data. 
# For data A and B, the date was selected as 1st September 2021

dom_a <- hpai %>% filter(observation.date < "2021-09-01")
dom_b <- hpai %>% filter(observation.date >= "2021-09-01")


# also need to split into quarters of the year
dom_a_q1 <- dom_a %>% dplyr::filter(month %in% c(1,2,3))
dom_a_q2 <- dom_a %>% dplyr::filter(month %in% c(4,5,6))
dom_a_q3 <- dom_a %>% dplyr::filter(month %in% c(7,8,9))
dom_a_q4 <- dom_a %>% dplyr::filter(month %in% c(10,11,12))
# repeat for b
dom_b_q1 <- dom_b %>% dplyr::filter(month %in% c(1,2,3))
dom_b_q2 <- dom_b %>% dplyr::filter(month %in% c(4,5,6))
dom_b_q3 <- dom_b %>% dplyr::filter(month %in% c(7,8,9))
dom_b_q4 <- dom_b %>% dplyr::filter(month %in% c(10,11,12))

# function to make the raster
presence_rast  <- function(pts, ref_rast){
  p <- cbind(pts$X, pts$Y)
#  points(p, pch = 18, cex = 0.5) # this just plots the points 
  points_rast <- rasterize(p, ref_rast, value = 1)
}
  
dom_a_q1_presence_rast <- presence_rast(dom_a_q1, euro_rast)
plot(dom_a_q1_presence_rast, col = "orange")

dom_a_q2_presence_rast <- presence_rast(dom_a_q2, euro_rast)
dom_a_q3_presence_rast <- presence_rast(dom_a_q3, euro_rast)
dom_a_q4_presence_rast <- presence_rast(dom_a_q4, euro_rast)
# repeat for b

dom_b_q1_presence_rast <- presence_rast(dom_b_q1, euro_rast)
dom_b_q2_presence_rast <- presence_rast(dom_b_q2, euro_rast)
dom_b_q3_presence_rast <- presence_rast(dom_b_q3, euro_rast)
dom_b_q4_presence_rast <- presence_rast(dom_b_q4, euro_rast)

# writeRaster(dom_a_q1_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q1_pres_abs.tif")
# writeRaster(dom_a_q2_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q2_pres_abs.tif")
# writeRaster(dom_a_q3_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q3_pres_abs.tif")
# writeRaster(dom_a_q4_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q4_pres_abs.tif")
# 
# writeRaster(dom_b_q1_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q1_pres_abs.tif")
# writeRaster(dom_b_q2_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q2_pres_abs.tif")
# writeRaster(dom_b_q3_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q3_pres_abs.tif")
# writeRaster(dom_b_q4_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q4_pres_abs.tif")

## One issue with this is that you may well just get rasters that overlay the wild bird cases.
## Probably useful to compare the two rasters and see if there are many cells that are independent

## Next part is to count how many unique occurrences there are in each cell.
## We have removed duplicates based on date and location, but possible that we might want to only count an 
## occurrence in a cell once in a day? 

## first will just count them all.

## issue with the code below is that using raster makes the crs tricky. 
# total_count_function <- function(points_obj){
#   p <- cbind(points_obj$X, points_obj$Y)  # make the points object
#   empty_rast <- raster::raster(euro_rast) # create an empty raster
#   empty_rast[] = 0 # set value as 0
#   counts = table(raster::cellFromXY(empty_rast,p)) # make a table of the counts for the cells
#   empty_rast[as.numeric(names(counts))] = counts #add the counts to the raster
#   return(empty_rast)
# }

total_count_function <- function(points_obj){
  p <- cbind(points_obj$X, points_obj$Y)  # make the points object
  empty_rast <- terra::rast(euro_rast ) # create an empty raster
  empty_rast[] = 0 # set value as 0
  counts = table(terra::cellFromXY(empty_rast,p)) # make a table of the counts for the cells
  cells_with_cases <- as.numeric(names(counts))
  cell_vals <- values(empty_rast)
  cell_vals[cells_with_cases] <-counts[] #  counts_vect
  values(empty_rast) <- cell_vals #add the counts to the raster
  return(empty_rast)
}

points_obj <- dom_a_q1

dom_a_q1_rast_counts <- total_count_function(dom_a_q1)
dom_a_q1_rast_counts
plot(dom_a_q1_rast_counts)
table(values(dom_a_q1_rast_counts))

breakpoints <- c(0, 1, 5, 10, 20, 50, 100, 150)
colors <- c("white", rev(RColorBrewer::brewer.pal(7, "Reds")))

plot(dom_a_q1_rast_counts, breaks = breakpoints, col = colors, box = F, axes = F)
plot(euro_shp, add = T)

## create the other rasters
dom_a_q2_rast_counts <- total_count_function(dom_a_q2)
dom_a_q3_rast_counts <- total_count_function(dom_a_q3)
dom_a_q4_rast_counts <- total_count_function(dom_a_q4)

dom_b_q1_rast_counts <- total_count_function(dom_b_q1)
dom_b_q1_rast_counts
dom_b_q2_rast_counts <- total_count_function(dom_b_q2)
dom_b_q3_rast_counts <- total_count_function(dom_b_q3)
dom_b_q4_rast_counts <- total_count_function(dom_b_q4)

# # save the rasters
# writeRaster(dom_a_q1_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q1_counts.tif")
# writeRaster(dom_a_q2_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q2_counts.tif")
# writeRaster(dom_a_q3_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q3_counts.tif")
# writeRaster(dom_a_q4_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q4_counts.tif")
# 
# writeRaster(dom_b_q1_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q1_counts.tif")
# writeRaster(dom_b_q2_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q2_counts.tif")
# writeRaster(dom_b_q3_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q3_counts.tif")
# writeRaster(dom_b_q4_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q4_counts.tif")




#####################################################################################################
# Compare the presence/absence rasters with those of wild bird cases.

# Read in the data

wild <- read.csv("data/flu_data/prepped_data/hpai_pos_birds_nobvbrc_april2024.csv")
head(wild)

## add a month and year to aid splitting
wild$year <- lubridate::year(wild$observation.date)
wild$month <- lubridate::month(wild$observation.date)

# We want these split into quarters but also by the period over which we split the data. 
# For data A and B, the date was selected as 1st September 2021

wild_a <- wild %>% filter(observation.date < "2021-09-01")
wild_b <- wild %>% filter(observation.date >= "2021-09-01")


# also need to split into quarters of the year
wild_a_q1 <- wild_a %>% dplyr::filter(month %in% c(1,2,3))
wild_a_q2 <- wild_a %>% dplyr::filter(month %in% c(4,5,6))
wild_a_q3 <- wild_a %>% dplyr::filter(month %in% c(7,8,9))
wild_a_q4 <- wild_a %>% dplyr::filter(month %in% c(10,11,12))
# repeat for b
wild_b_q1 <- wild_b %>% dplyr::filter(month %in% c(1,2,3))
wild_b_q2 <- wild_b %>% dplyr::filter(month %in% c(4,5,6))
wild_b_q3 <- wild_b %>% dplyr::filter(month %in% c(7,8,9))
wild_b_q4 <- wild_b %>% dplyr::filter(month %in% c(10,11,12))


wild_a_q1_presence_rast <- presence_rast(wild_a_q1, euro_rast)
plot(wild_a_q1_presence_rast, col = "purple")

wild_a_q2_presence_rast <- presence_rast(wild_a_q2, euro_rast)
wild_a_q3_presence_rast <- presence_rast(wild_a_q3, euro_rast)
wild_a_q4_presence_rast <- presence_rast(wild_a_q4, euro_rast)
# repeat for b

wild_b_q1_presence_rast <- presence_rast(wild_b_q1, euro_rast)
wild_b_q2_presence_rast <- presence_rast(wild_b_q2, euro_rast)
wild_b_q3_presence_rast <- presence_rast(wild_b_q3, euro_rast)
wild_b_q4_presence_rast <- presence_rast(wild_b_q4, euro_rast)

# Now we want to compare the domestic and wild rasters and find how many there are with values in
# one of the rasters but not in the other.

# start with q1 a
# set NAs as 0 in both cases

get_raster_diff <- function(rd, rw){
  rd[is.na(rd[])] <- 0 # set the NAs as zero in domestic 
  rw[is.na(rw[])] <- 10 # set the NAs to 10 in wildlife
  rdiff <- rd - rw # subtract one from the other
  return(rdiff)
}


diff_a_q1 <- get_raster_diff(dom_a_q1_presence_rast, wild_a_q1_presence_rast)
table(values(diff_a_q1))

# as done domestic versus wild. 
# pos dom (1) - pos wild(1) = 0 
# neg dom(0) - neg wild(10) = - 10  
# pos dom(1) - neg wild(10) = -9
# neg dom(0) - pos wild = -1

diff_a_q2 <- get_raster_diff(dom_a_q2_presence_rast, wild_a_q2_presence_rast)
table(values(diff_a_q2))

diff_a_q3 <- get_raster_diff(dom_a_q3_presence_rast, wild_a_q3_presence_rast)
table(values(diff_a_q3))

diff_a_q4 <- get_raster_diff(dom_a_q4_presence_rast, wild_a_q4_presence_rast)
table(values(diff_a_q4))

# and for b 
diff_b_q1 <- get_raster_diff(dom_b_q1_presence_rast, wild_b_q1_presence_rast)
table(values(diff_b_q1))
diff_b_q2 <- get_raster_diff(dom_b_q2_presence_rast, wild_b_q2_presence_rast)
table(values(diff_b_q2))
diff_b_q3 <- get_raster_diff(dom_b_q3_presence_rast, wild_b_q3_presence_rast)
table(values(diff_b_q3))
diff_b_q4 <- get_raster_diff(dom_b_q4_presence_rast, wild_b_q4_presence_rast)
table(values(diff_b_q4))
