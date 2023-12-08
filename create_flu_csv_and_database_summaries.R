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


## deal with the negatives

bv_neg_data_trim <- rename(bv_neg_data_trim, Latitude = Collection.Latitude,
                           Longitude = Collection.Longitude,
                           observation.date = Collection.Date,
                           Species  = Host.Species,
                           Country = Collection.Country,
                           Serotype = Subtype)

bv_neg_data_trim <- dplyr::select(bv_neg_data_trim, c(Latitude, Longitude, observation.date, Species,
                                                      Country, Serotype, source))
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

euro_shp <- terra::vect("output/euro_map.shp")
ext(euro_shp)

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
#write.csv(pos_ai_data_prj_area, "output/avian_flu/pos_points_proj_area_all_sources_all_wild.csv" )
#write.csv(neg_ai_data_prj_area, "output/avian_flu/neg_points_proj_area_all_sources_all_wild.csv" )

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

#write.csv(dt_pos, "output/avian_flu/pos_points_proj_area_all_sources_all_wild_duplicates_removed.csv" )
#write.csv(dt_neg, "output/avian_flu/neg_points_proj_area_all_sources_all_wild_duplicates_removed.csv" )

plot(euro_shp)
plot(pos_ai_data_prj_area, add = T, pch = 18, legend = T)

### Next task is to remove the non-birds from the data... 

# use Taxize to see how many can be identified automatically 
# If no updates to data can skip this part and just read in the csv. 

#install.packages("taxizedb")
#library(taxizedb)
# library(taxize)
# 
# species_list <- as.data.frame(table(ai_data_prj_area$Species))
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
# for (i in 26:nrow(not_sp)) {
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

# write.csv(species_list, "avian_flu_scripts/detecting_mammals.csv")

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

nrow(dt_pos[which(dt_pos$Species %in% unwanted_sp),])

# 73 mammal entries that we are removing. 

# unwanted_sp <- c("Unspecified Mammal", "Stone Marten", 
#                  "South American Coati (Nasua Nasua):Procyonidae-Carnivora",
#                  "Red Fox", "Polecat", "Polar Fox", "Nyctereutes Viverrinus (Japanese Racoon Dog)",
#                  "Mink", "Harbor Seal (Phoca Vitulina):Phocidae-Carnivora", 
#                  "Gray Seal (Halichoerus Grypus):Phocidae-Carnivora", "Fox",
#                  "European Pine Marten", "Civet", "Caspian Seal (Pusa Caspica)",
#                  "Halichoerus grypus", "Phoca vitulina", "Mustela putorius")

ai_pos_birds <- 
   dt_pos[-which(dt_pos$Species %in% unwanted_sp),]


table(ai_pos_birds$Serotype)
ai_pos_birds$serotype_HN <- ai_pos_birds$Serotype
ai_pos_birds$serotype_HN <-  gsub("HPAI","",as.character(ai_pos_birds$serotype_HN))
ai_pos_birds$serotype_HN <-  gsub("LPAI","",as.character(ai_pos_birds$serotype_HN))
ai_pos_birds$serotype_HN<-toupper(ai_pos_birds$serotype_HN) # so they match
ai_pos_birds$serotype_HN <- trimws(ai_pos_birds$serotype_HN) # tri, amy white space innthe entries to aid matching
table(ai_pos_birds$serotype_HN)

# manually change the one entry that has the serotype entered twice
ai_pos_birds[(which(ai_pos_birds$serotype_HN == "H7N2,H7N2")), "serotype_HN"] <- "H7N2"

# Difficult to split by LPAI and HPA as not specified directly. 
# FAO and WOAH only report HPAI so the ones listed in this data base should all be HPAI

table(ai_pos_birds$source, ai_pos_birds$serotype_HN)

serotype_source_df <- as.data.frame(table(ai_pos_birds$source, ai_pos_birds$serotype_HN))
serotype_source_df <- pivot_wider(serotype_source_df, names_from = Var1, values_from = "Freq")

sum(colSums(serotype_source_df[2:4]))

# if we filter by all the entries that don't have an entry in either FAO or WOAH then this
# may well leave only the LPAI? 

poss_lpai <- serotype_source_df[which(serotype_source_df$fao == 0 & serotype_source_df$woah == 0),]
sum(poss_lpai[,"bvbrc"]) # These are the ones that I am most confident are LPAI. But would need to exclude H5Nx as could be HPAI


# Add a year and month and year_month to the data

ai_pos_birds$year <- lubridate::year(ai_pos_birds$observation.date)
ai_pos_birds$month <- lubridate::month(ai_pos_birds$observation.date)
ai_pos_birds$month_year <-  format(as.Date(ai_pos_birds$observation.date), "%Y-%m")
ai_pos_birds$week <- format(ai_pos_birds$observation.date, "%Y Week %W")
ai_pos_birds$week_num <- lubridate::isoweek(ai_pos_birds$observation.date)

serotype_data <- as.data.frame(table(ai_pos_birds$serotype_HN))

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


# For quite a number of these there are only a few entries. 
# For the line plot by subtype, only use those that have >10 entries

serotype_data_ten_or_more <- serotype_data[which(serotype_data$Freq > 10),]
serotype_data_ten_or_more <- serotype_data_ten_or_more[which(serotype_data_ten_or_more$Var1 != ""),]

subtype_plot_data <- ai_pos_birds[which(ai_pos_birds$serotype_HN %in% serotype_data_ten_or_more$Var1),]
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

#ggsave("plots/subtypes/serotype_line_by_week.png")


## Start here
## If want to include year and month, will need to do another table to join with 
## so don't miss any months.

# subtype_counts_year <- subtype_plot_data %>%
#   group_by(year, serotype_HN) %>% count()
# 
# subtype_counts_year_all <- left_join()
# 
# ggplot(data = subtype_counts_year, aes(x = year, y = n, col = serotype_HN))+
#   geom_line()
# 
# # by month and year. 
# subtype_counts_monthyear <- subtype_plot_data %>%
#   group_by(month_year, serotype_HN) %>% count()
# subtype_counts_monthyear$log_count <- log(subtype_counts_monthyear$n)
# subtype_counts_monthyear$serotype_HN_factor <- as.factor(subtype_counts_monthyear$serotype_HN)
# 
# ggplot(data = subtype_counts_monthyear, 
#        aes(x = month_year, y = n, col = serotype_HN_factor, group = serotype_HN_factor)) +
#   geom_line()
# 
# 
# # by week 
# subtype_counts_weekyear <- subtype_plot_data %>%
#   group_by(week, serotype_HN) %>% count()
# 
# ggplot(data = subtype_counts_weekyear, 
#        aes(x = week, y = n, col = serotype_HN, group = serotype_HN)) +
#   geom_line(lwd = 0.8) + 
#   scale_color_brewer(palette="Paired")
# ggsave("plots/subtypes/serotype_line_by_week.png")
# 

# stacked bar plot coloured by serotype

subtype_counts_year$factor_HN <- as.factor(subtype_counts_year$serotype_HN)

ggplot(data = subtype_plot_data, aes(x = year, fill = serotype_HN))+
  geom_bar()


ggplot(data = subtype_plot_data, aes(x = month_year, fill = serotype_HN)) + 
  geom_bar()
ggsave("plots/subtypes/serotype_bar_by_month.png")

ggplot(data = subtype_plot_data, aes(x = week, fill = serotype_HN)) + 
  geom_bar()
ggsave("plots/subtypes/serotype_bar_by_week.png")



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



# ############################# MAKE RASTER  ############################
# 
# # Other manipulations needed to rasterize
# table(point_data$flu) # this is pos and neg data
# 
# points_sp <- sf::as_Spatial(point_data)
# head(points_sp)
# 
# points_vect <- terra::vect(points_sp)
# head(points_vect)
# plot(points_vect) # this is also still the global data set
# 
# #read in the raster we are using for the project
# euro_rast <- terra::rast("output/euro_rast.tif")
# euro_rast <- terra::rast("output/euro_rast_plus.tif")
# 
# euro_rast
# 
# ## use the vector within rasterize
# points_rast <- terra::rasterize(points_vect, euro_rast, field = "flu", fun = "max")
# # using max because if it has a positive in the cell (1) and a negative (0) we want
# # the cell to be classes as a negative
# plot(points_rast)
# points_rast
# head(points_rast)
# 
# plot(euro_shp)
# plot(points_rast, add = TRUE, col = "red")
# 
# table(terra::values(points_rast))
# 
# # Create the csv file for avian flu cases
# # make a points object using the centre of each pixel from the ref raster
# points_euro_rast <- terra::as.points(euro_rast)
# points_euro_rast
# 
# tictoc::tic()
# flu_res <- terra::extract(points_rast, points_euro_rast, method = "simple", xy = T)
# tictoc::toc()
# 
# table(flu_res$max)
# table(values(points_rast))
# euro_rast
# nrow(flu_res)
# points_rast
# 
# # why have we got fewer values in the results table rather than the raster? 
# # Looked at it with the 10km raster and wonder if some of them are in the sea/on land too small to be 
# # noted on our raster?? 
# 
# points_rast
# euro_rast
# table(values(points_rast))
# # we have the higher number in the points raster
# 
# points_euro_rast # the points are all within the boundaries of the euro-raster.
# # so that's not the issue
# 
# # so I think the issue is somewhere in extract or if there were NAs in the euroraster that
# # are not in the case data.
# 
# euro_nonNA_coords <- as.data.frame(terra::crds(euro_rast, na.rm = T))
# cases_nonNA_coords <- as.data.frame(terra::crds(points_rast, na.rm = T))
# sum(table(values(points_rast))) - sum(table(flu_res$max))
# diff_nas <- setdiff( cases_nonNA_coords, euro_nonNA_coords) # this is the same number
# # as are missing between the datasets
# 
# # let's turn these into spatial points and plot them 
# missing_points <- st_as_sf(x = diff_nas, 
#                            coords = c("x", "y"),
#                            crs = "3035")
# plot(euro_map)
# plot(missing_points, add = T, col = "red", pch = 18)
# 
# 
# 
# 
# 
# dev.off()
# plot(euro_map)
# plot(points_rast, xlim= c(4000000, 5000000), ylim = c(3000000, 4000000), 
#      col = "hot pink")#, background = "blue")
# plot(points_rast, col = "red")
# plot(points_rast, col = "red", add = T, axes = F)
# plot(euro_rast, col = "White", add = T, axes = F)
# #plot(euro_map, add = T, axes = F)
# 
# euro_rast
# points_rast
# euro_rast[1000:2000] # we can see that there are NAs
# points_rast[1:100]
# 
# # we want to know which points have NA in euro_rast that aren't NA in points rast. 
# # if it's a NA in the euro_rast but has a value in the points rast that might explain it
# 
# ################################################################
# ## Some trials to aid understanding/visualisation
# 
# ## Hard to see. Look at with larger raster
# 
# #read in the raster we are using for the project
# euro_rast_10k <- terra::rast("output/euro_rast_10k.tif")
# euro_rast_10k
# 
# points_rast_10k <- terra::rasterize(points_vect, euro_rast_10k,"flu", fun = max)
# plot(points_rast_10k)
# points_rast_10k
# head(points_rast_10k)
# table(values(points_rast_10k))
# # Easier to visualise that it probably is working OK. 
# # I think the 1km are just too small to see. 
# 
# dev.off()
# plot(euro_map)
# plot(points_rast_10k, col = "red")
# plot(euro_rast_10k, col = "White", add = T, axes = F)
# 
# 
# 
# # Create the csv file for avian flu cases
# # make a points object using the centre of each pixel from the ref raster
# points_euro_rast_10k <- terra::as.points(euro_rast_10k)
# points_euro_rast_10k
# 
# flu_res_10k <- terra::extract(points_rast_10k, points_euro_rast_10k,
#                               method = "simple", xy = T)
# 
# table(flu_res_10k$max)
# table(values(points_rast_10k))
# 
# dev.off()
# plot(points_rast_10k, col = "red")
# plot(euro_map, add = T, axes = F)
# 
# 
# 
# #Perhaps just look at UK? 
# GB_ext <- terra::ext(3300000, 3800000, 3100000, 4100000)
# GB_crop <- terra::crop(euro_map, GB_ext)
# plot(GB_crop)
# 
# GB_rast <- terra::crop(points_rast, GB_ext)
# plot(GB_crop)
# plot(GB_rast, add = T, axes = F)
# 
# 
# ## bit from fao data - useful for making mammal dataset later 
# 
# unwanted_sp <- c("Unspecified Mammal", "Stone Marten", 
#                  "South American Coati (Nasua Nasua):Procyonidae-Carnivora",
#                  "Red Fox", "Polecat", "Polar Fox", "Nyctereutes Viverrinus (Japanese Racoon Dog)",
#                  "Mink", "Harbor Seal (Phoca Vitulina):Phocidae-Carnivora", 
#                  "Gray Seal (Halichoerus Grypus):Phocidae-Carnivora", "Fox",
#                  "European Pine Marten", "Civet", "Caspian Seal (Pusa Caspica)")
# 
# fao_data <- fao_data[-which(fao_data$Species %in% unwanted_sp),]
