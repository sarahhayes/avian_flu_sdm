# In this script we investigate similiarities and differences between species-level seasons.

library(av)
require(ape)
library(BART)
library(DALEX)
library(data.table)
library(dplyr)
library(ebirdst)
library(fields)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(pdftools)
require(PhyloMeasures)
library(raster)
library(RColorBrewer)
library(readxl)
library(rnaturalearth)
library(sf)
library(stringr)
library(terra)
library(tidyr)

# Set season plot palette

pal <- c("#2271B2",
         "#F748A5",
         "#359B73",
         "#e69f00")

### First get list of birds which appear in Europe with eBird species codes, along
# with estimates of global population sizes. ####

# Get ebird species table
sp_df <- ebirdst_runs
sp_df$scientific_name <- sapply(sp_df$scientific_name,
                                tolower,
                                USE.NAMES = FALSE)

# Get codes for species in Europe
euro_bird_codes <- read.csv("ebird/codes_for_europe_clean.csv")

# euro_bird_codes <- read.csv("ebird/species_europe_2024_after_pkg_update.csv")
# euro_bird_codes <- euro_bird_codes[, c("common_name", "sci_code")]

no_euro_birds <- nrow(euro_bird_codes)

#Restrict to species in Europe list
sp_df <- sp_df[which(sp_df$species_code %in% euro_bird_codes$code), ]

p <- ggplot(sp_df) +
  geom_point(aes(x=breeding_start, y=breeding_end, colour="Breeding")) +
  geom_point(aes(x=nonbreeding_start, y=nonbreeding_end, colour="Nonbreeding")) +
  geom_point(aes(x=postbreeding_migration_start, y=postbreeding_migration_end, colour="Postbreeding")) +
  geom_point(aes(x=prebreeding_migration_start, y=prebreeding_migration_end, colour="Prereeding")) +
  xlab("Period starts") +
  ylab("Period ends") +
  labs(colour="Season")

nspecies <- nrow(sp_df)

date_seq <- seq(as.Date("2021/01/01"), by = "day", length.out = 365)

season_by_week <- matrix(nrow=nspecies, ncol=365)
for (i in 1:nspecies){
  starts <- c(sp_df$nonbreeding_start[i],
              sp_df$prebreeding_migration_start[i],
              sp_df$breeding_start[i],
              sp_df$postbreeding_migration_start[i])
  if (length(which(!is.na(starts)))>0){
    starts <- starts[which(!is.na(starts))]
    for (j in 1:365){
      if (min(starts)>date_seq[j]){
        season_by_week[i, j] <- which.max(starts)
      }
      else{
        preceding_dates <- starts[which(starts <= date_seq[j])]
        season_by_week[i, j] <- which(starts==max(preceding_dates))
      }
    }
  }
}
season_names <- c("Nonbreeding",
           "Prebreeding migration",
           "Breeding",
           "Postbreeding migration")
season_df <- season_by_week %>% data.frame()
colnames(season_df) <- date_seq
season_df <- season_df %>% pivot_longer(everything(),
                                                                           names_to = "Day",
                                                                           values_to = "Season")
season_df$Season <- season_names[season_df$Season]
season_table <- table(season_df)[, season_names] %>% data.frame()

# Check to make sure we get same number of species every day:
count_sp_by_day <- function(this_day){
  total_spec <- season_table %>%
    filter(Day == this_day) %>%
    select(Freq) %>%
    colSums()
  return(total_spec)
}
count_sp_by_day <- Vectorize(count_sp_by_day)
req_by_day <- data.frame(Day = unique(season_table$Day)) %>%
  mutate(Species_count = count_sp_by_day(Day))

# If the following has length 1, then we know we have consistent species
# numbers.
unique(req_by_day$Species_count)

h <- ggplot(season_table, aes(x=as.Date(Day), y=Freq, fill=Season)) +
  geom_bar(stat="identity", position="fill") +
  scale_x_date(date_breaks="months") +
  scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Date") +
  ylab("Proportion of species")

ggsave(paste("plots/species_seasons_by_week.png",
             sep=""),
       plot = h,
       width = 8.5,
       height = 4.5)

season_df_short <- pivot_wider(season_table, names_from = "Season", values_from = "Freq")
get_largest <- apply(season_df_short[,2:5], 1, max)
is_largest <- season_df_short[,2:5]==get_largest
breeding_start <- date_seq[is_largest[,1] %>% which() %>% min()]
nonbreeding_start <- date_seq[(is_largest[,3] %>% which() %>% max()) + 1]
postbreeding_start <- date_seq[is_largest[,3] %>% which() %>% min()]
prebreeding_start <- date_seq[is_largest[,4] %>% which() %>% min()]

cat("Based on plurality season limits are\n",
    nonbreeding_start %>% as.character(),
    "\n",
    prebreeding_start %>% as.character(),
    "\n",
    breeding_start %>% as.character(),
    "\n",
    postbreeding_start %>% as.character()
    )

is_maj <- season_df_short[,2:5]> 0.5 * rowSums(season_df_short[,2:5])
breeding_start <- date_seq[is_maj[,1] %>% which() %>% min()]
nonbreeding_start <- date_seq[(is_maj[,3] %>% which() %>% max()) + 1]
postbreeding_start <- date_seq[is_maj[,3] %>% which() %>% min()]
prebreeding_start <- date_seq[is_maj[,4] %>% which() %>% min()]

cat("Based on first majority date season limits are\n",
    nonbreeding_start %>% as.character(),
    "\n",
    prebreeding_start %>% as.character(),
    "\n",
    breeding_start %>% as.character(),
    "\n",
    postbreeding_start %>% as.character()
)


#### Now add time series plot ####

# Set ratio of pseudoabsences to positives (x:1)
pseud_ratio <- 1

# Do you need to initialise the environmental rasters for use in environmental thinning?
# If this is the first time this script is run, should be set to TRUE

INIT_ENV_VARS <- FALSE

# Set season plot palette

pal <- c("#2271B2",
         "#F748A5",
         "#359B73",
         "#e69f00")

#############
# Read data #
#############

# Read in extent of Europe mapped
base_map <- terra::rast("output/euro_rast_10k.tif")
euro_map_st <- st_read("output/euro_map.shp") # country borders

# Read in weighting layers for pseudoabsence sampling
# weight_layer <- terra::rast("variable_manipulation/variable_outputs/human_density.tif")
# layername <- "Human density"

weight_layer <- terra::rast("variable_manipulation/variable_outputs/ebird_records.tif")
layername <- "eBird record"

#weight_layer_log <- weight_layer %>% app(function(x) log(x+1))

# Define starts of ecological seasons
year_start <- as.Date("0000-01-01")
q2_start <- as.Date("0000-03-01") # pre-breeding migration
q3_start <- as.Date("0000-06-07") # breeding season
q4_start <- as.Date("0000-08-10") # post-breeding migration
q1_start <- as.Date("0000-11-30") # non-breeding season
year_end <- as.Date("0000-12-31")

# Read in positive flu site data 
pos_sites <- read.csv("hpai_pos_birds_nobvbrc.csv") %>% 
  rename(date = observation.date) %>%
  mutate(date = as.Date(date),
         Q = case_when((format(date, "%m%-%d") >= format(year_start, "%m%-%d")) & (format(date, "%m%-%d") < format(q2_start-1, "%m%-%d)")) ~ "Q1",
                       (format(date, "%m%-%d") >= format(q2_start, "%m%-%d")) & (format(date, "%m%-%d") < format(q3_start-1, "%m%-%d)")) ~ "Q2",
                       (format(date, "%m%-%d") >= format(q3_start, "%m%-%d")) & (format(date, "%m%-%d") < format(q4_start-1, "%m%-%d)")) ~ "Q3",
                       (format(date, "%m%-%d") >= format(q4_start, "%m%-%d")) & (format(date, "%m%-%d") < format(q1_start-1, "%m%-%d)")) ~ "Q4",
                       (format(date, "%m%-%d") >= format(q1_start, "%m%-%d")) & (format(date, "%m%-%d") <= format(year_end, "%m%-%d)")) ~ "Q1"))

pos_sites %>% pull(Q) %>% table
pos_sites %>% pull(serotype_HN) %>% table

# Data set A: 2.3.4.4b H5NX before H5N1  (includes H5N8, H5N6, retain ambiguous H5 or unlabelled HPAI [n = 80])
# Training set A: 2.3.4.4b H5NX in distinct 2016/2017 outbreak
# Test set A: 2.3.4.4b H5NX in distinct 2020/2021 outbreak
df_A <- pos_sites %>% 
  filter(date >= as.Date("2016-08-10") & date < as.Date("2021-08-10") & serotype_HN %in% c("H5N8", "H5N6", "H5", "")) %>%
  mutate(df = case_when(date < as.Date("2020-08-10") ~ "train_A",
                        date >= as.Date("2020-08-10") ~ "test_A"))

# Data set B: 2.3.4.4b H5N1 (retain ambiguous H5 or unlabelled HPAI [n = 64])
# Training set B: 2.3.4.4b H5N1 from Sep 21 - Mar 23 %>% 
# Test set B: H5N1 from Apr 23 - Mar 24 
df_B <- pos_sites %>% 
  filter(date >= as.Date("2021-08-10") & date < as.Date("2024-03-01") & serotype_HN %in% c("H5N1", "H5", "")) %>%
  mutate(df = case_when(date < as.Date("2023-03-01") ~ "train_B",
                        date >= as.Date("2023-03-01") ~ "test_B"))

# How many per quarter
pre_table <- bind_rows(df_A, df_B) %>% with(., table(df, Q))
pre_table

# Plot over time
# df_A %>% filter(grepl("train", df)) %>%
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#00BA38")
# 
# df_A %>% filter(grepl("test", df)) %>%
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#63db87")
# 
# df_B %>% filter(grepl("train", df)) %>%
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#619CFF")
# 
# df_B %>% filter(grepl("test", df)) %>%
#   ggplot(aes(x = date)) +
#   geom_histogram(fill = "#99bfff")

bind_rows(df_A, df_B) %>% 
  ggplot(aes(x = date, fill = df)) +
  geom_histogram(bins=100)

timeplot_seas <- pos_sites %>% 
  filter(date > as.Date("2016-08-10") & serotype_HN %in% c("H5N8", "H5N6", "H5N1", "H5", "")) %>%
  mutate(serotype_HN = case_when(
    serotype_HN == "H5N1" ~ "H5N1",
    serotype_HN == "H5N8" ~ "H5N8",
    serotype_HN == "H5N6" ~ "H5N6",
    TRUE ~ "H5NX"
  )) %>%
  ggplot(aes(x = date)) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2016-08-10"), xmax = as.Date("2016-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2017-08-10"), xmax = as.Date("2017-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2018-08-10"), xmax = as.Date("2018-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2019-08-10"), xmax = as.Date("2019-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2020-08-10"), xmax = as.Date("2020-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2021-08-10"), xmax = as.Date("2021-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2022-08-10"), xmax = as.Date("2022-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2023-08-10"), xmax = as.Date("2023-11-30"), ymin = -Inf, ymax = Inf), fill = pal[1], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2016-11-30"), xmax = as.Date("2017-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2017-11-30"), xmax = as.Date("2018-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2018-11-30"), xmax = as.Date("2019-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2019-11-30"), xmax = as.Date("2020-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2020-11-30"), xmax = as.Date("2021-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2021-11-30"), xmax = as.Date("2022-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2022-11-30"), xmax = as.Date("2023-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2023-11-30"), xmax = as.Date("2024-03-01"), ymin = -Inf, ymax = Inf), fill = pal[2], alpha = 0.2) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2017-03-01"), xmax = as.Date("2017-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2018-03-01"), xmax = as.Date("2018-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2019-03-01"), xmax = as.Date("2019-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2020-03-01"), xmax = as.Date("2020-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2021-03-01"), xmax = as.Date("2021-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2022-03-01"), xmax = as.Date("2022-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2023-03-01"), xmax = as.Date("2023-06-07"), ymin = -Inf, ymax = Inf), fill = pal[3], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2017-06-07"), xmax = as.Date("2017-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2018-06-07"), xmax = as.Date("2018-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2019-06-07"), xmax = as.Date("2019-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2020-06-07"), xmax = as.Date("2020-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2021-06-07"), xmax = as.Date("2021-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2022-06-07"), xmax = as.Date("2022-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_rect(data = slice(pos_sites,1), aes(xmin = as.Date("2023-06-07"), xmax = as.Date("2023-08-10"), ymin = -Inf, ymax = Inf), fill = pal[4], alpha = 0.2) + geom_vline(xintercept = as.Date("2020-08-10")) +
  geom_histogram(bins = round(as.numeric((max(pos_sites$date)-min(pos_sites$date))/30)), position = "stack", fill = "grey20") +
  geom_vline(xintercept = as.Date("2021-08-10")) +
  geom_vline(xintercept = as.Date("2023-03-01")) +
  geom_text(aes(x = as.Date("2016-08-10")+(as.Date("2020-08-10")-as.Date("2016-08-10"))/2, y = 330 ,label = "A. Training", hjust = 0.5)) +
  geom_text(aes(x = as.Date("2020-08-10")+(as.Date("2021-08-10")-as.Date("2020-08-10"))/2, y = 330 ,label = "A. Test", hjust = 0.5)) +
  geom_text(aes(x = as.Date("2021-08-10")+(as.Date("2023-03-01")-as.Date("2021-08-10"))/2, y = 330 ,label = "B. Training", hjust = 0.5)) +
  geom_text(aes(x = as.Date("2024-01-01"), y = 330 ,label = "B. Test", hjust = 0.5)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  ylab("Monthly reports") +
  xlab("Date") +
  theme_bw()

fig_combi <- (h)/(timeplot_seas) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(1,1)) &
  theme(legend.position = 'bottom')

ggsave(paste("plots/seasons_and_reports.png",
             sep=""),
       plot = fig_combi,
       width = 10.5,
       height = 8.5)
