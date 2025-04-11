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
  starts <- c(sp_df$breeding_start[i],
              sp_df$nonbreeding_start[i],
              sp_df$postbreeding_migration_start[i],
              sp_df$prebreeding_migration_start[i])
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
season_names <- c("Breeding",
           "Nonbreeding",
           "Postbreeding migration",
           "Prebreeding migration")
season_df <- season_by_week %>% data.frame()
colnames(season_df) <- date_seq
season_df <- season_df %>% pivot_longer(everything(),
                                                                           names_to = "Day",
                                                                           values_to = "Season")
season_df$Season <- season_names[season_df$Season]
season_table <- table(season_df) %>% data.frame()

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
