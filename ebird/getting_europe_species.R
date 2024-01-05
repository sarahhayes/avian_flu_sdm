# 04/01/2024
# Script to process the european bird species in a more reproducible format. 
# CSV that is read in comes from: https://science.ebird.org/en/status-and-trends/species?regionCode=eu
# These are the birds for which ebird had status and trends data that are present in Europe. 
# The list of species, after filtering to Europe, is cut and pasted into Notepad and then transferred to Excel. 

rm(list = ls())

library(tidyverse)
library(ebirdst)

raw_data <- read.csv("data/europe_bird_species.csv", header = F)
head(raw_data)

colnames(raw_data) <- "species"

# remove any repeating names
no_dups <- unique(raw_data, by = c("species"))

no_dups %>% filter(str_detect(species, "FORMES"))
no_dups %>% filter(str_detect(species, "dae"))
no_dups %>% filter(str_detect(species, "SPECIES"))


no_dups <- no_dups %>% filter(!str_detect(species, "SPECIES")) %>%
  filter(!str_detect(species, "FORMES")) %>%
  filter(!str_detect(species, "dae")) 

## Also want to remove the headings

heading_names <- c("Rheas", "Ducks, Geese, and Waterfowl", "Pheasants, Grouse, and Allies", 
                   "New World Quail", "Flamingos", "Grebes",
                   "Pigeons and Doves", "Sandgrouse",
                   "Cuckoos", "Nightjars and Allies", "Swifts", "Rails, Gallinules, and Coots", 
                   "Cranes", "Stilts and Avocets", "Oystercatchers", "Plovers and Lapwings",
                   "Jacanas", "Sandpipers and Allies", "Pratincoles and Coursers", "Skuas and Jaegers",
                   "Gulls, Terns, and Skimmers", "Loons", "Shearwaters and Petrels", "Storks",
                   "Frigatebirds", "Boobies and Gannets", "Cormorants and Shags", "Pelicans", 
                   "Herons, Egrets, and Bitterns", "Ibises and Spoonbills", "New World Vultures",
                   "Hawks, Eagles, and Kites", "Barn-Owls", "Owls", "Hoopoes", 
                   "Ground-Hornbills", "Kingfishers", "Bee-eaters", "Rollers", "Woodpeckers",
                   "Falcons and Caracaras", "Cockatoos", "Old World Parrots", "New World and African Parrots",
                   "Tyrant Flycatchers", "Vireos, Shrike-Babblers, and Erpornis", "Old World Orioles",
                   "Shrikes", "Crows, Jays, and Magpies", "Tits, Chickadees, and Titmice", 
                   "Penduline-Tits", "Larks", "Cisticolas and Allies", 
                   "Reed Warblers and Allies", "Grassbirds and Allies", "Swallows", "Bulbuls",
                   "Leaf Warblers", "Bush Warblers and Allies", "Sylviid Warblers and Allies",
                   "Parrotbills", "White-eyes, Yuhinas, and Allies", "Laughingthrushes and Allies",
                   "Kinglets", "Nuthatches", "Treecreepers", "Wrens", "Dippers",
                   "Starlings", "Mockingbirds and Thrashers", "Thrushes and Allies", "Old World Flycatchers",
                   "Waxwings", "Weavers and Allies", "Waxbills and Allies", "Indigobirds", 
                   "Accentors", "Old World Sparrows", "Wagtails and Pipits", "Finches, Euphonias, and Allies", 
                   "Longspurs and Snow Buntings", "Old World Buntings", "New World Sparrows",
                   "Troupials and Allies", "New World Warblers", "Cardinals and Allies",
                   "Tanagers and Allies")

no_dups <- no_dups %>% filter(!species %in% heading_names)

scientific_name <- no_dups %>% filter(row_number() %% 2 == 0) %>% ## Select even rows
  rename("sci_name" = "species")
scientific_name[,"sci_code"]<- ebirdst::get_species(scientific_name[,"sci_name"])

common_name <- no_dups %>% filter(row_number() %% 2 == 1) %>% ## Select odd rows
  rename("common_name" = "species") 
common_name[,"common_code"]<- ebirdst::get_species(common_name[,"common_name"])

species_europe <- cbind(common_name, scientific_name)

no_match <- species_europe[which(species_europe$common_code != species_europe$sci_code),]
  
## There are a few that are NA in the first column but give a code using the scientific name and vice versa

nrow(species_europe[is.na(species_europe$common_code),])
nrow(species_europe[is.na(species_europe$sci_code),])

# visual inspections shows not the same rows, so merge to get as many as possible.
# after package update this was no longer true

# species_europe[is.na(species_europe$common_code), "common_code"] <- 
#   species_europe[is.na(species_europe$common_code), "sci_code"]
# 
# species_europe[is.na(species_europe$sci_code), "sci_code"] <- 
#   species_europe[is.na(species_europe$sci_code), "common_code"]
# 
# nrow(species_europe[is.na(species_europe$common_code),])
# nrow(species_europe[is.na(species_europe$sci_code),])
# no_match <- species_europe[which(species_europe$common_code != species_europe$sci_code),] # checking haven't inadvertantly caused an issue!

#Save this file
# write.csv(species_europe, "ebird/species_europe_2024.csv")
# write.csv(species_europe, "ebird/species_europe_2024_after_pkg_update.csv")
 
# compare to the previous file

old_file <- read.csv("ebird/codes_for_europe_clean.csv")

old_codes <- old_file$code
new_codes <- scientific_name$sci_code

setdiff(old_codes, new_codes)
setdiff(new_codes, old_codes)


                   