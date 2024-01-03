## 02/11/2023
## Creating the avian flu data csv for use within the BART model
## Also create some summaries of the individual datasets for the area we are using
## and remove the mammals into a separate set

rm(list = ls())

library(readxl)
library(tidyverse)
library(terra)
library(RColorBrewer)

## bring in the different data sets and combine to one large data set. 

## first do fao as it's the biggest
#FAO 
# use the full data as want to include asia as well 
fao_data  <- read.csv("data/flu_data/raw_data/epidemiology-raw-data_202308011345.csv")

table(fao_data$Region)

colnames(fao_data)
# rename the date columns
fao_data <-  dplyr::rename(fao_data, observation.date = "Observation.date..dd.mm.yyyy.")
fao_data <-  dplyr::rename(fao_data, report.date = "Report.date..dd.mm.yyyy.")

# There are rows with no entries for observation date so remove these
fao_data <- fao_data[which(fao_data$observation.date != ""),]

# we only want to include wild birds
table(fao_data$Species)


fao_data_trim <- dplyr::select(fao_data, all_of(c("Latitude", "Longitude", "observation.date", "Species", "Country", "Serotype")))
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
colnames(inf_a_hpai)
table(inf_a_hpai$is_wild)
table(inf_a_hpai$Epi_unit)
table(inf_a_hpai$Epi_unit, inf_a_hpai$is_wild)

# I think we want to filter by is.wild = T 
# This will leave in some animals that are listed as wild in a zoo - but these could have been wild animals 
# that were found within zoo grounds? Or they could be captive wild... Looking at the numbers, the latter seems more likely

zoo_cases <- inf_a_hpai[which(inf_a_hpai$Epi_unit == "Zoo"),]

table(zoo_cases$wild_type) # this doesn't show the many NAs

# so perhaps we should remove the zoo captives? So select all those that are wild and then remove those that
# are captive 
inf_a_hpai_wild <- inf_a_hpai[which(inf_a_hpai$is_wild == T),]
inf_a_hpai_wild <- inf_a_hpai_wild[which(inf_a_hpai_wild$wild_type != "captive"),]
table(inf_a_hpai_wild$Epi_unit, inf_a_hpai_wild$is_wild)

table(inf_a_hpai_wild$wild_type)

colnames(inf_a_hpai_wild)
inf_a_hpai_trim <- dplyr::select(inf_a_hpai_wild, all_of(c("Latitude", "Longitude", "Outbreak_start_date", "Species", 
                                                           "country", "sero_sub_genotype_eng")))
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
                                                     "Collection.Longitude",
                                                     "Host.Species",
                                                     "Collection.Country",
                                                     "Subtype")))

bv_pos_data_trim$source <- "bvbrc"
bv_pos_data_trim$Collection.Date <- as.Date(bv_pos_data_trim$Collection.Date,
                                            "%Y-%m-%d")


bv_neg_data_trim <- dplyr::select(neg_data, all_of(c("Collection.Date", 
                                                     "Collection.Latitude",
                                                     "Collection.Longitude",
                                                     "Host.Species",
                                                     "Collection.Country",
                                                     "Subtype")))
bv_neg_data_trim$source <- "bvbrc"
bv_neg_data_trim$Collection.Date <- as.Date(bv_neg_data_trim$Collection.Date,
                                            "%Y-%m-%d")

# We now have the trimmed data sets for all the sources. 
# combine them and add a section for pos or neg

# first need all the names to match so we can rbind

colnames(fao_data_trim)
colnames(inf_a_hpai_trim)
colnames(bv_pos_data_trim)
inf_a_hpai_trim <- rename(inf_a_hpai_trim, observation.date = Outbreak_start_date, 
                          Country = country, Serotype = sero_sub_genotype_eng)
bv_pos_data_trim <- rename(bv_pos_data_trim, Latitude = Collection.Latitude,
                           Longitude = Collection.Longitude,
                           observation.date = Collection.Date,
                           Species  = Host.Species,
                           Country = Collection.Country,
                           Serotype = Subtype)

bv_pos_data_trim <- dplyr::select(bv_pos_data_trim, c(Latitude, Longitude, observation.date, Species,
                                                      Country, Serotype, source))
colnames(bv_pos_data_trim)

all_ai_data <- rbind(fao_data_trim, inf_a_hpai_trim, bv_pos_data_trim)
all_ai_data$flu <- 1 # denotes a positive

###### ** Not using negatives so no longer need the below script
## deal with the negatives

# bv_neg_data_trim <- rename(bv_neg_data_trim, Latitude = Collection.Latitude,
#                            Longitude = Collection.Longitude,
#                            observation.date = Collection.Date,
#                            Species  = Host.Species,
#                            Country = Collection.Country,
#                            Serotype = Subtype)
# 
# bv_neg_data_trim <- dplyr::select(bv_neg_data_trim, c(Latitude, Longitude, observation.date, Species,
#                                                       Country, Serotype, source))
# colnames(bv_neg_data_trim)
# bv_neg_data_trim$flu <- 0 # denotes negative test
# 
# all_ai_data <- rbind(all_ai_data, bv_neg_data_trim)

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

euro_shp <- terra::vect("output/euro_map.shp")
ext(euro_shp)

ai_data_prj_area <- point_data_df %>%
  dplyr::filter(X >= 2600000 & X <= 7000000) %>%
  dplyr::filter(Y >= 1500000 & Y <= 6400000)

# select positives - should be all data now not using negative

pos_ai_data_prj_area <- ai_data_prj_area %>%
  filter(flu == 1)

#neg_ai_data_prj_area <- ai_data_prj_area %>%
#  filter(flu == 0)

nrow(pos_ai_data_prj_area)
nrow(ai_data_prj_area)

#save these
#write.csv(pos_ai_data_prj_area, "output/avian_flu/pos_points_proj_area_all_sources_all_wild.csv" )
#write.csv(neg_ai_data_prj_area, "output/avian_flu/neg_points_proj_area_all_sources_all_wild.csv" )

# ## We need to remove mammals from these data 
# ### Next task is to remove the non-birds from the data... 
# 
# # use Taxize to see how many can be identified automatically 
# # If no updates to data can skip this part and just read in the csv below. 
# 
# #install.packages("taxizedb")
# #library(taxizedb)
# library(taxize)
#  
#  species_list <- as.data.frame(table(ai_data_prj_area$Species))
#  
# colnames(species_list) <- c("species", "freq")
# species_list$species <- as.character(species_list$species)
#  
# for (i in 1:nrow(species_list)) {
#   species_list[i,"class"] <- tax_name(species_list[i,"species"],
#                                        get = "class",
#                                        db = "ncbi")$class
#   print(i)
# }
# 
# not_sp <- species_list[which(is.na(species_list$class)),]
# # manual inspection suggests some of these may not be found due to things like having 'incognita' added to name
# 
# not_sp$species_edit <- not_sp$species
# not_sp$species_edit <- sub("\\(.*", "", not_sp$species_edit)
# not_sp$species_edit <- sub("\\:.*", "", not_sp$species_edit)
# 
# for (i in 1:nrow(not_sp)) {
#   not_sp[i,"class"] <- tax_name(not_sp[i,"species_edit"],
#                                       get = "class",
#                                       db = "ncbi")$class
#   print(i)
# }
# 
# 
# # still quite a few without a class which may have to manually look through
# 
# not_sp <- rename(not_sp, class_2 = class)
# species_list <- left_join(species_list, not_sp)
# species_list[which(is.na(species_list$class)), "class"] <- species_list[which(is.na(species_list$class)), "class_2"]
# species_list <- dplyr::select(species_list, c("species", "freq", "class"))
# 
# # Peacock has been mislabelled as an insect so change to Aves
# 
# species_list[which(species_list$species == "Peacock"),"class"] <- "Aves"

#select the ones labelled "Mammalia" for removal
# save the species list so don't have to generate every time'

#write.csv(species_list, "avian_flu_scripts/detecting_mammals.csv", row.names = F)

species_list <- read.csv("avian_flu_scripts/detecting_mammals.csv")

unwanted_sp <- species_list[which(species_list$class == "Mammalia"), "species"]

## manually extract those not automatically done. 
unwanted_sp_manual_ext <- c("Al indeterminatum fau",
                            "Fox",
                            "Nyctereutes Viverrinus (Japanese Racoon Dog)",
                            "Polar Fox",
                            "Polecat",
                            "Stone Marten",
                            "Unspecified Mammal")


unwanted_sp <- c(unwanted_sp, unwanted_sp_manual_ext)

nrow(ai_data_prj_area[which(ai_data_prj_area$Species %in% unwanted_sp),])

# 123 mammal entries that we are removing. 

# unwanted_sp <- c("Unspecified Mammal", "Stone Marten", 
#                  "South American Coati (Nasua Nasua):Procyonidae-Carnivora",
#                  "Red Fox", "Polecat", "Polar Fox", "Nyctereutes Viverrinus (Japanese Racoon Dog)",
#                  "Mink", "Harbor Seal (Phoca Vitulina):Phocidae-Carnivora", 
#                  "Gray Seal (Halichoerus Grypus):Phocidae-Carnivora", "Fox",
#                  "European Pine Marten", "Civet", "Caspian Seal (Pusa Caspica)",
#                  "Halichoerus grypus", "Phoca vitulina", "Mustela putorius")

## Remove from the data so just birds

ai_pos_birds <- 
  ai_data_prj_area[-which(ai_data_prj_area$Species %in% unwanted_sp),]

png("plots/fao_data_pos.png", width = 480, height = 480)
plot(euro_shp, main = "FAO")
plot(ai_pos_birds[which(ai_pos_birds$source == "fao"),], add = T, 
     pch = 18, legend = T, col = "red", cex = 0.5)
dev.off()

png("plots/woah_data_pos.png", width = 480, height = 480)
plot(euro_shp, main = "WOAH")
plot(ai_pos_birds[which(ai_pos_birds$source == "woah"),], add = T, 
     pch = 18, legend = T, col = "blue", cex = 0.5)
dev.off()

png("plots/bvbrc_data_pos.png", width = 480, height = 480)
plot(euro_shp, main = "BV-BRC")
plot(ai_pos_birds[which(ai_pos_birds$source == "bvbrc"),], add = T, 
     pch = 18, legend = T, col = "orange", cex = 0.5)
dev.off()

nrow(ai_pos_birds[which(ai_pos_birds$source == "bvbrc"),])

## Before remove duplicates, try and count number of species/families etc.
## To do this, run the script in add_species_info.R

species_list <- 
  species_list[-which(species_list$species %in% unwanted_sp),]

table(species_list$family)
unique(species_list$family)

## 44 bird families in total 

## Save this table as it takes a while to run
# write.csv(species_list, "avian_flu_scripts/species_list_with_order_and_family.csv")

## Next I want to remove the duplicates that are present in the data 

library(data.table)
dt_pos <- data.table(ai_pos_birds)

no_duplicates_date_and_loc <- unique(dt_pos, by = c("X", "Y", "observation.date"))
no_duplicates_all <- unique(dt_pos, by = c("X", "Y", "observation.date", "source", "Species", "Country"))
no_dups_all_except_source <- unique(dt_pos, by = c("X", "Y", "observation.date", "Species", "Country"))

# The last two are the same, suggesting that there are none that are the same for all details except source. 

## No duplicates_all will still have some that are the same place and time but different species. However, when looking
## at these, there are instances where we have fao and woah reporting a same location and date but one has the common name
## for the bird and one has the scientific name. I think we basically want to know if a location has AI in wild birds
## so we don't mind which species. 
## As such, no_duplicates_date_and_time is probably the one to use. 

# However, if we want to look at HPAI only, it might be worth filtering on this first? ]
# Just in case we might otherwise remove HPAI and leave in LPAI. 

# Let's look at serotypes

# easier to re-classify the dataframe so can use the existing code
ai_pos_birds <- no_duplicates_date_and_loc

table(ai_pos_birds$Serotype)
ai_pos_birds$serotype_HN <- ai_pos_birds$Serotype
ai_pos_birds$serotype_HN <-  gsub("HPAI","",as.character(ai_pos_birds$serotype_HN))
ai_pos_birds$serotype_HN <-  gsub("LPAI","",as.character(ai_pos_birds$serotype_HN))
ai_pos_birds$serotype_HN<-toupper(ai_pos_birds$serotype_HN) # so they match
ai_pos_birds$serotype_HN <- trimws(ai_pos_birds$serotype_HN) # trim any white space in the entries to aid matching
table(ai_pos_birds$serotype_HN)

# manually change the one entry that has the serotype entered twice
ai_pos_birds[(which(ai_pos_birds$serotype_HN == "H7N2,H7N2")), "serotype_HN"] <- "H7N2"

# Difficult to split by LPAI and HPAI as not specified directly. 
# FAO and WOAH only report HPAI so the ones listed in this data base should all be HPAI

table(ai_pos_birds$source, ai_pos_birds$serotype_HN)

serotype_source_df <- as.data.frame(table(ai_pos_birds$source, ai_pos_birds$serotype_HN))
serotype_source_df <- pivot_wider(serotype_source_df, names_from = Var1, values_from = "Freq")

sum(colSums(serotype_source_df[2:4]))

# if we filter by all the entries that don't have an entry in either FAO or WOAH then this
# may well leave only the LPAI? 
# Alternatively just select all the H5 and H7

hpai <- ai_pos_birds %>% filter(str_detect(serotype_HN,"H5") | str_detect(serotype_HN, "H7"))
hpai <- rename(hpai, species = Species)

## To do - need to get species info for the birds that are HPAI only? 

#write.csv(hpai, "data/flu_data/prepped_data/hpai_pos_birds")

# Add a year and month and year_month to the data

hpai$year <- lubridate::year(hpai$observation.date)
hpai$month <- lubridate::month(hpai$observation.date)
hpai$month_year <-  format(as.Date(hpai$observation.date), "%Y-%m")
hpai$week <- format(hpai$observation.date, "%Y Week %W")
hpai$week_num <- lubridate::isoweek(hpai$observation.date)

range(hpai$observation.date)

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

ggsave("plots/subtypes/serotype_line_by_week_hpai_only.png")


## If want to repeat by month, will need to do another table to join with 
## so don't miss any months.

subtype_counts_year <- subtype_plot_data %>%
  group_by(year, serotype_HN) %>% count()

table(subtype_counts_year$year) #NB there are no data from 2012 in here. 

# stacked bar plot coloured by serotype

subtype_counts_year$factor_HN <- as.factor(subtype_counts_year$serotype_HN)

ggplot(data = subtype_plot_data, aes(x = year, fill = serotype_HN))+
  geom_bar()

# Year seems to work OK. suspect this is because year is a numeric variable. 
# The lines below are likely to miss out some dates as they are character variables. 
# This seems to be confirmed by the fact that there are no gaps in the plots. 
# As such, don't use below without creating new data set where specify the weeks/months over the period as we did above. 

# 
# ggplot(data = subtype_plot_data, aes(x = month_year, fill = serotype_HN)) + 
#   geom_bar()
# ggsave("plots/subtypes/serotype_bar_by_month.png")
# 
# ggplot(data = subtype_plot_data, aes(x = week, fill = serotype_HN)) + 
#   geom_bar()
# ggsave("plots/subtypes/serotype_bar_by_week.png")



#Plots of quarters – for each day in the year per quarter do a stacked bar plot coloured by year –
# to visualise which years we are primarily getting data from for each quarter.

q1_sub <- subtype_plot_data[which(subtype_plot_data$month %in% c(1,2,3)),]
q2_sub <- subtype_plot_data[which(subtype_plot_data$month %in% c(4,5,6)),]
q3_sub <- subtype_plot_data[which(subtype_plot_data$month %in% c(7,8,9)),]
q4_sub <- subtype_plot_data[which(subtype_plot_data$month %in% c(10,11,12)),]

q1_sub$year_fact <- as.factor(q1_sub$year)

# set a palette to plot
colourCount = length(unique(q1_sub$year_fact))
colourCount = 16
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

year_cols <- c("2005" = "grey", 
               "2006" = "pink", 
               "2007" = "green", 
               "2008" = "brown",
               "2009" = "#003399",
               "2010" = "#FF3300",
               "2011" = "#339966",
               "2012" = "#FF66FF",
               "2013" = "#FF3366",
               "2014" = "#009999",
               "2015" = "#CCCFFF",
               "2016" = "#CC6633",
               "2017" = "#FFCC33",
               "2018" = "#9900FF",
               "2019" = "#33CC00",
               "2020" = "#CC6699",
               "2021" = "#3399FF",
               "2022" = "#FF0000",
               "2023" = "#00CC33")

#trial run
ggplot(q1_sub, aes(x=week_num, fill = year_fact)) +
   geom_bar() +
  scale_fill_manual(values = c(year_cols)) + 
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"))
  
# these data contain some labelled as week 52 and 53 due to variations in calendar. 
# Not sure of the best way to deal with these. Possibly label as 0 and -1 to aid display? 

q1_sub[which(q1_sub$week_num == 53), "week_num"] <- 0
q1_sub[which(q1_sub$week_num == 52), "week_num"] <- -1

q1plot <- ggplot(q1_sub, aes(x=week_num, fill = year_fact)) +
  geom_bar() +
  scale_fill_manual(values = c(year_cols))+ 
  ggtitle("Q1") + 
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"))

q1plot

q2_sub$year_fact <- as.factor(q2_sub$year)
q2plot <- ggplot(q2_sub, aes(x=week_num, fill = year_fact)) +
  geom_bar() +
  scale_fill_manual(values = c(year_cols))+ 
  ggtitle("Q2") + 
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"))

q2plot

q3_sub$year_fact <- as.factor(q3_sub$year)
q3plot <- ggplot(q3_sub, aes(x=week_num, fill = year_fact)) +
  geom_bar() +
  scale_fill_manual(values = c(year_cols))+
  ggtitle("Q3") + 
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"))


q4_sub$year_fact <- as.factor(q4_sub$year)
q4plot <- ggplot(q4_sub, aes(x=week_num, fill = year_fact)) +
  geom_bar() +
  scale_fill_manual(values = c(year_cols))+
  ggtitle("Q4")+ 
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"))

ggpubr::ggarrange(q1plot, q2plot, q3plot, q4plot, ncol = 2, nrow = 2)

ggsave("plots/counts_by_week_per_quart_coloured_by_year.png")


## Would also be useful to see the spatial distribution of the hpai cases. 
plot(euro_shp)
plot(hpai$geometry, add = T, pch = 18, cex = 0.6, col = "red")


# Alternatively, plot the H5N8 and H5N1
plot(euro_shp, main = "H5N1")
plot(hpai[which(hpai$serotype_HN == "H5N1"), "geometry"], add = T, col = "blue", pch = 18, cex = 0.5)

plot(euro_shp, main = "H5N8")
plot(hpai[which(hpai$serotype_HN == "H5N8"), "geometry"], add = T, col = "orange", pch = 18, cex = 0.5)

# for the H5N8, also split by the first and second waves. 
plot(euro_shp, main = "H5N8 pre-2020")
plot(hpai[which(hpai$serotype_HN == "H5N8" & hpai$year < "2020"), "geometry"], add = T, col = "orange", pch = 18, cex = 0.5)

plot(euro_shp, main = "H5N8 2020 onwards")
plot(hpai[which(hpai$serotype_HN == "H5N8" & hpai$year >= "2020"), "geometry"], add = T, col = "dark green", pch = 18, cex = 0.5)

