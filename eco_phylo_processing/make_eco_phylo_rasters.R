# In this script we combine abundance estimates from eBird with species-specific
# ecological and phylogenetic factors to create rasters encoding the abundance
# of birds with different qualities.

PLOT <- FALSE

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

# First get list of birds which appear in Europewith eBird species codes, along
# with estimates of global population sizes.

# Get ebird species table
sp_df <- ebirdst_runs
sp_df$scientific_name <- sapply(sp_df$scientific_name,
                                tolower,
                                USE.NAMES = FALSE)

# Get codes for species in Europe
euro_bird_codes <- read.csv("ebird/codes_for_europe_clean.csv")

no_euro_birds <- nrow(euro_bird_codes)

#Restrict to species in Europe list
sp_df <- sp_df[which(sp_df$species_code %in% euro_bird_codes$code), ]

# Load in population sizes
pop_size_df <- read.csv("callaghan_pop_estimates.csv")
pop_size_df$Scientific.name <- sapply(pop_size_df$Scientific.name, tolower)

# Try to match up
cat(length(setdiff(sp_df$scientific_name, pop_size_df$Scientific.name)),
    "binomial names are missing from the population size list.")
cat(length(setdiff(sp_df$common_name, pop_size_df$Common.name)),
    "common names are missing from the population size list.")
sci_mismatch <- which(!(sp_df$scientific_name %in% pop_size_df$Scientific.name))
com_mismatch <- which(!(sp_df$common_name %in% pop_size_df$Common.name))
full_mismatch <- intersect(sci_mismatch, com_mismatch)
cat(length(full_mismatch),
  "species are missing binomial AND common names in the population size list.")
cat("Missing species common names are:",
    sp_df$common_name[full_mismatch],
    "scientific names are",
    sp_df$scientific_name[full_mismatch],
    sep = "\n")

# Also need to check matches are 1-1, i.e. when we have a match for both names
# they point us to the same place
eBird_to_pop_by_sci <- sapply(sp_df$scientific_name,
                              FUN = function(x){
                                which(pop_size_df$Scientific.name==x)})
eBird_to_pop_by_comm <- sapply(sp_df$common_name,
                              FUN = function(x){
                                which(pop_size_df$Common.name==x)})
consistent_matches <- lapply(
  1:no_euro_birds,
  FUN = function(i){
    identical(eBird_to_pop_by_sci$Value[i], eBird_to_pop_by_comm$Value[i]) |
      !(sp_df$scientific_name[i] %in% pop_size_df$Scientific.name) |
      !(sp_df$common_name[i] %in% pop_size_df$Common.name)
  }
)
cat("There are",
    length(which(consistent_matches==FALSE)),
    "common-scientific mismatches.")

# Bring in AVONET to see if we can fix the missing species using AVONET IDs
AVONET_df <- read_excel("data/AVONETSupplementarydataset1.xlsx",
                        sheet = "AVONET_Raw_Data")
AVONET_df <- AVONET_df[,
                       c("Avibase.ID",
                         "Species1_BirdLife",
                         "Species2_eBird",
                         "Species3_BirdTree")]
AVONET_df <- distinct(AVONET_df)
AVONET_df[, 2:4] <- sapply(AVONET_df[, 2:4], tolower)

missing_species_in_AVONET <- which(
  AVONET_df$Species2_eBird %in% sp_df$scientific_name[full_mismatch])

# Look at what AVONET has to say about missing species:
missing_species_df <- AVONET_df[missing_species_in_AVONET, ]
print(missing_species_df)

# See which alternative names show up in the population size data:
alts_from_BirdLife <- intersect(pop_size_df$Scientific.name,
                                missing_species_df$Species1_BirdLife)
alts_from_BirdTree <- intersect(pop_size_df$Scientific.name,
                                missing_species_df$Species3_BirdTree)

sub_strs <- sapply(union(alts_from_BirdLife, alts_from_BirdTree),
                   FUN = function(x){
                     paste(
                       missing_species_df$Species2_eBird[
                         which(missing_species_df$Species3_BirdTree==x)],
                       x,
                       sep = " -> ")
                   })
cat("The following substitutions are possible but need sense-checking:",
    sub_strs,
    sep = "\n"
    )

# Saxicola stejnegeri (Amur stonechat) is a vagrant species only rarely found in
# Europe, so we remove it from consideration
euro_bird_codes <- euro_bird_codes[
  !grepl("Amur Stonechat", euro_bird_codes$species), ]
sp_df <- sp_df[!grepl("saxicola stejnegeri", sp_df$scientific_name), ]

# Oenanthe melanoleuca (Eastern Black-eared Wheatear) is considered a subspecies
# of Oenanthe hispanica (Western Black-eared Wheatear). The two datasets
# conflict in their subspecies assignment so we remove both from consideration
euro_bird_codes <- euro_bird_codes[
  !grepl("Eastern Black-eared Wheatear", euro_bird_codes$species), ]
sp_df <- sp_df[!grepl("oenanthe melanoleuca", sp_df$scientific_name), ]
euro_bird_codes <- euro_bird_codes[
  !grepl("Western Black-eared Wheatear", euro_bird_codes$species), ]
sp_df <- sp_df[!grepl("oenanthe hispanica", sp_df$scientific_name), ]

# The classification of the entire Curruca genus is problematic, so we remove it
curruca_common_names <- sp_df$common_name[
  grepl("curruca",sp_df$scientific_name)]
for (name in curruca_common_names){
  euro_bird_codes <- euro_bird_codes[
    !grepl(name, euro_bird_codes$species), ]
}
sp_df <- sp_df[!grepl("curruca", sp_df$scientific_name), ]

# Remove both European green woodpecker and Iberian green woodpecker due to
# inconsistent subspecies assignment:
euro_bird_codes <- euro_bird_codes[
  !grepl("Iberian Green Woodpecker", euro_bird_codes$species), ]
sp_df <- sp_df[!grepl("picus sharpei", sp_df$scientific_name), ]
euro_bird_codes <- euro_bird_codes[
  !grepl("Eurasian Green Woodpecker", euro_bird_codes$species), ]
sp_df <- sp_df[!grepl("picus viridis", sp_df$scientific_name), ]

# The following species have no clear substitute:
true_missing <- missing_species_df$Species2_eBird[
  -which(
    (missing_species_df$Species3_BirdTree %in% pop_size_df$Scientific.name)|
        (missing_species_df$Species3_BirdTree %in% pop_size_df$Scientific.name))
    ]
cat("Can not identify matches for the following in AVONET:",
    true_missing,
    sep = "\n"
    )

# These two species appear to be missing with no obvious substitution so we
# remove them from the data
# Note that we need to use code in euro_bird_codes because the name has problem
# characters
euro_bird_codes <- euro_bird_codes[
  !grepl("krunut1", euro_bird_codes$code), ]
sp_df <- sp_df[!grepl("sitta krueperi", sp_df$scientific_name), ]
euro_bird_codes <- euro_bird_codes[
  !grepl("Northern Hawk Owl", euro_bird_codes$species), ]
sp_df <- sp_df[!grepl("surnia ulula", sp_df$scientific_name), ]

# Should now find that all species from eBird have population sizes:
cat(length(setdiff(sp_df$scientific_name, pop_size_df$Scientific.name)),
    "binomial names are missing from the population size list.")
cat(length(setdiff(sp_df$common_name, pop_size_df$Common.name)),
    "common names are missing from the population size list.")
sci_mismatch <- which(!(sp_df$scientific_name %in% pop_size_df$Scientific.name))
com_mismatch <- which(!(sp_df$common_name %in% pop_size_df$Common.name))
full_mismatch <- intersect(sci_mismatch, com_mismatch)
cat(length(full_mismatch),
    "species are missing binomial AND common names in the population size list.")

# Add pop sizes to species data:
no_euro_birds <- nrow(euro_bird_codes)
sp_df$pop_sizes <- 0
sp_df$pop_sizes <- sapply(1:no_euro_birds,
                          FUN = function(i){
                            if (i %in% sci_mismatch){
                              return(pop_size_df$Abundance.estimate[which(pop_size_df$Common.name==sp_df$common_name[i])])
                            }
                            else{
                              return(pop_size_df$Abundance.estimate[which(pop_size_df$Scientific.name==sp_df$scientific_name[i])])
                            }
                          })

################################################################################

# Control ambiguity by assigning only a single ID to each combination of names
# AVONET_df <- AVONET_df[!duplicated(AVONET_df[,2:4]), ]

# ID functions for matching species between datasets:
get_single_bird_id <- function(binomial_name, name_sources){
  
  # Return -100 if provided with na
  if (binomial_name=="na"){
    return("-100")
  }
  # Find IDs from each list and keep track of which species list they came from.
  ID_candidates <- c()
  ID_weightings <- c()
  ID_sources <- c()
  if (("BirdLife" %in% name_sources) & 
      (binomial_name %in% AVONET_df$Species1_BirdLife)){
    ID_candidates <- append(ID_candidates, AVONET_df$Avibase.ID[
      which(AVONET_df$Species1_BirdLife==binomial_name)])
    no_matches <- length(which(AVONET_df$Species1_BirdLife==binomial_name))
    for (i in 1:no_matches){
      ID_weightings <- append(ID_weightings, 1/no_matches)
      ID_sources <- append(ID_sources, "BirdLife")
    }
  }
  if (("eBird" %in% name_sources) & 
      (binomial_name %in% AVONET_df$Species2_eBird)){
    ID_candidates <- append(ID_candidates, AVONET_df$Avibase.ID[
      which(AVONET_df$Species2_eBird==binomial_name)])
    no_matches <- length(which(AVONET_df$Species2_eBird==binomial_name))
    for (i in 1:no_matches){
      ID_weightings <- append(ID_weightings, 1/no_matches)
      ID_sources <- append(ID_sources, "eBird")
    }
  }
  if (("BirdTree" %in% name_sources) & 
      (binomial_name %in% AVONET_df$Species3_BirdTree)){
    ID_candidates <- append(ID_candidates, AVONET_df$Avibase.ID[
      which(AVONET_df$Species3_BirdTree==binomial_name)])
    no_matches <- length(which(AVONET_df$Species3_BirdTree==binomial_name))
    for (i in 1:no_matches){
      ID_weightings <- append(ID_weightings, 1/no_matches)
      ID_sources <- append(ID_sources, "BirdTree")
    }
  }
  # How we work out the consensus depends on how many candidates are identified
  if (length(ID_candidates)==0){
    # If ID doesn't appear in any list, return -100
    return("-100")
  }
  else{
    return(unique(ID_candidates))
    # # Construct consensus distribution of candidate names
    # unique_candidates <- unique(ID_candidates)
    # if (length(unique_candidates)==1){
    #   # cat("Number of unique candidates is", length(unique_candidates), "\n")
    #   # cat(unique_candidates, "\n")
    #   return(unique_candidates)
    # }
    # candidate_dist <- c()
    # for (i in 1:length(unique_candidates)){
    #   # cat("i=", i, "\n")
    #   candidate_dist <- append(candidate_dist,
    #                            sum(ID_weightings[
    #                              which(ID_candidates==unique_candidates[i])]))
    # }
    # # cat(unique_candidates, "\n", ID_candidates, "\n", candidate_dist, "\n", ID_weightings, "\n")
    # if (length(which(candidate_dist==max(candidate_dist)))==1){
    #   return(unique_candidates[which(candidate_dist==max(candidate_dist))])
    # }
    # else{
    #   # cat("Multiple consensus options for species", binomial_name,"\n")
    #   return(sample(unique_candidates[which(candidate_dist==max(candidate_dist))], 1))
    # }
  }
}

get_bird_ids <- function(name_vect, name_sources){
  if (missing(name_sources)){
    name_sources <-c("BirdLife",
                     "eBird",
                     "BirdTree")
  }
  
  cols_needed <- 0
  # Check max number of options we will need
  if ("BirdLife" %in% name_sources){
    names_in_source <- AVONET_df$Species1_BirdLife[
      which(AVONET_df$Species1_BirdLife %in% name_vect)]
    if (max(table(names_in_source)) > cols_needed){
      cols_needed <- max(table(names_in_source))
    }
  }
  if ("eBird" %in% name_sources){
    names_in_source <- AVONET_df$Species2_eBird[
      which(AVONET_df$Species2_eBird %in% name_vect)]
    if (max(table(names_in_source)) > cols_needed){
      cols_needed <- max(table(names_in_source))
    }
  }
  if ("BirdTree" %in% name_sources){
    names_in_source <- AVONET_df$Species3_BirdTree[
      which(AVONET_df$Species3_BirdTree %in% name_vect)]
    if (max(table(names_in_source)) > cols_needed){
      cols_needed <- max(table(names_in_source))
    }
  }
  id_df <- data.frame(matrix(nrow = length(name_vect), ncol = (cols_needed+1)))
  id_df[, 1] <- name_vect
  for (n in 1:length(name_vect)){
    # cat("n=",n,"\n")
    id_vect <- get_single_bird_id(name_vect[n], name_sources)
    id_df[n, 2:(length(id_vect)+1)] <- id_vect
  }
  return(id_df)
}

# Find IDs for species in eBird
eBird_IDs <- get_bird_ids(sp_df$scientific_name)
# Check we got everything:
cat(length(which(eBird_IDs[, 2]=="-100")), "species were not ID'd.")
non_IDd_species <- sp_df$scientific_name[which(eBird_IDs[, 2]=="-100")]
cat("Unidentified species are:",
    non_IDd_species,
    sep = "\n")

# Porphyrio Poliocephalus doesn't have an obvious match in AVONET and is only
# found in small numbers at the very edge of the EPSG:3035 range so we just
# remove it from the list
euro_bird_codes <- euro_bird_codes[
  !grepl("Gray-headed Swamphen", euro_bird_codes$species), ]
sp_df <- sp_df[!grepl("porphyrio poliocephalus", sp_df$scientific_name), ]

# This is a simple alternative name:
sp_df$scientific_name[which(sp_df$scientific_name=="daptrius chimachima")] <- 
  "milvago chimachima"

# Try getting IDs again:
# Find IDs for species in eBird
eBird_IDs <- get_bird_ids(sp_df$scientific_name)
# Check we got everything:
cat(length(which(eBird_IDs[, 2]=="-100")), "species were not ID'd.")

# Add ID options to species dataframe
sp_df$Avibase_ID <- eBird_IDs[, 2:ncol(eBird_IDs)]

################################################################################
# Now introduce CLOVER database and identify all avian species known to be hosts
# of influenza A.

setwd("../clover")

# Load data
CLOVER_df <- read.csv(
  "clover/clover_1.0_allpathogens/CLOVER_1.0_Viruses_AssociationsFlatFile.csv")
# Restrict attention to samples from birds
CLOVER_df <- CLOVER_df[which(CLOVER_df$HostClass == "aves"), ]
# All pathogen species names are in lower case. Filtering for rows where this
# name contains "flu" should cover any variations on influenza.
CLOVER_df <- CLOVER_df[grepl("flu", CLOVER_df$Pathogen), ]

# Print list of pathogens in CLOVER_df - in practice this is now just influenza
# A.
unique(CLOVER_df$Pathogen)

# Now filter by a positive test result - we're only interested in host species.
CLOVER_df <- CLOVER_df[which((CLOVER_df$Detection_NotSpecified == "TRUE") |
                               (CLOVER_df$Detection_Serology == "TRUE") |
                               (CLOVER_df$Detection_Genetic == "TRUE") |
                               (CLOVER_df$Detection_Isolation == "TRUE")), ]
# In practice this last step doesn't do anything because the dataset is already
# restricted to species with known cases.

# Get list of all species where positive samples have been found
species_list <- unique(CLOVER_df$Host)
no_host_species <- length(species_list)

# Alternatively, we could count the number of times each species appears if we
# want to restrict to species with incidence above a certain window:

species_counts <- CLOVER_df %>% count(Host, sort = TRUE)

if (PLOT){
  # Plot species counts:
  p <- ggplot(data=CLOVER_df, aes(x=factor(Host))) + geom_bar(stat="count")
  p + coord_flip() + 
    scale_x_discrete(limits=species_counts$Host) + 
    xlab("Host species")
  
  # Try plotting only those species with at least 10 recorded cases:
  below_10_species <- species_counts$Host[which(species_counts$n<10)]
  p <- ggplot(data=CLOVER_df[-which(CLOVER_df$Host %in% below_10_species), ], aes(x=factor(Host))) + geom_bar(stat="count")
  p + coord_flip() + 
    scale_x_discrete(limits=species_counts$Host[-which(species_counts$Host %in% below_10_species)]) + 
    xlab("Host species")
}

# Do the same for higher taxonomic levels; this will be useful for identifying
# taxa of interest.
order_case_counts <- CLOVER_df %>% count(HostOrder, sort = TRUE)
family_case_counts <- CLOVER_df %>% count(HostFamily, sort = TRUE)
genus_case_counts <- CLOVER_df %>% count(HostGenus, sort = TRUE)

if (PLOT){
  # Plot order counts:
  p <- ggplot(data=CLOVER_df, aes(x=factor(HostOrder))) + 
    geom_bar(stat="count")
  p + coord_flip() + 
    scale_x_discrete(limits=order_case_counts$HostOrder) + 
    xlab("Host order") +
    ylab("Number of cases")
}

if (PLOT){
  # Plot family counts, but just do families with at least 5 host species:
  atleast5_families <- family_case_counts$HostFamily[which(family_case_counts$n>=5)]
  data_to_plot <- CLOVER_df[which(CLOVER_df$HostFamily %in% atleast5_families), ]
  plot_lims <- family_case_counts$HostFamily[which(family_case_counts$HostFamily %in% atleast5_families)]
  p <- ggplot(data=data_to_plot, aes(x=factor(HostFamily))) + 
    geom_bar(stat="count")
  p + coord_flip() + 
    scale_x_discrete(limits=plot_lims) + 
    xlab("Host family") +
    ylab("Number of cases")
}

if (PLOT){
  # Plot genus counts, but just do families with at least 10 host species:
  atleast10_genuses <- genus_case_counts$HostGenus[which(genus_case_counts$n>=10)]
  data_to_plot <- CLOVER_df[which(CLOVER_df$HostGenus %in% atleast10_genuses), ]
  plot_lims <- genus_case_counts$HostGenus[which(genus_case_counts$HostGenus %in% atleast10_genuses)]
  p <- ggplot(data=data_to_plot, aes(x=factor(HostGenus))) + 
    geom_bar(stat="count")
  p + coord_flip() + 
    scale_x_discrete(limits=plot_lims) + 
    xlab("Host genus") +
    ylab("Number of cases")
}

# Get species ID's of species in CLOVER
host_species_IDs <- get_bird_ids(species_list)

# See if any species were not successfully ID'd
cat(length(which(host_species_IDs=="-100")), "species were not ID'd.")
cat("Unidentified species are:",
    species_list[which(host_species_IDs=="-100")],
    sep = "\n")

# We manually replace these three species names with names found in AVONET.
# These are all unambiguous replacements - they are present in all three source
# datasets and associated with a single Avibase ID.
species_list[which(species_list=="coloeus monedula")] <- "corvus monedula"
species_list[which(species_list=="larus mongolicus")] <- "larus cachinnans"
species_list[which(species_list=="piaya minuta")] <- "coccycua minuta"

# Try again with ID's
host_species_IDs <- get_bird_ids(species_list)

# We now shouldn't see any failed ID's
cat(length(which(host_species_IDs=="-100")), "species were not ID'd.")

# Remove possible repeats:
host_species_IDs <- unique(host_species_IDs)
no_host_species <- length(host_species_IDs)

################################################################################
# Now introduce trait data

setwd('../avian_flu_sdm')

# Load in from file
EltonTraits_df <- 
  read.table("data/variables/bird_data/elton_traits/BirdFuncDat.txt",
             sep = '\t',
             quote="\"",
             header = TRUE)
# elton trait data has two rows of NA's at the bottom which we need to remove:
EltonTraits_df <- EltonTraits_df[1:9993, ]
EltonTraits_df$Scientific <- tolower(EltonTraits_df$Scientific)

# Search for IDs
EltonTraits_IDs <- get_bird_ids(EltonTraits_df$Scientific)

# Check we got everything:
cat(length(which(EltonTraits_IDs=="-100")), "species were not ID'd.")
non_IDd_species <- EltonTraits_df$Scientific[which(EltonTraits_IDs=="-100")]
cat("Unidentified species are:",
    non_IDd_species,
    sep = "\n")

# Rename to match AVONET:
EltonTraits_df$Scientific[
  which(EltonTraits_df$Scientific=="anthus longicaudatus")
] <- "anthus vaalensis"
EltonTraits_df$Scientific[
  which(EltonTraits_df$Scientific=="polioptila clementsi")
] <- "polioptila guianensis"
EltonTraits_df$Scientific[
  which(EltonTraits_df$Scientific=="hypositta perdita")
] <- "oxylabes madagascariensis"

# This next species just seems to be missing from AVONET so we remove it:
EltonTraits_df <- EltonTraits_df[
  which(EltonTraits_df$Scientific!="amaurospiza carrizalensis"), ]

EltonTraits_df$Scientific[
  which(EltonTraits_df$Scientific=="lophura hatinhensis")
] <- "lophura edwardsi"

# Try again:
EltonTraits_IDs <- get_bird_ids(EltonTraits_df$Scientific)
cat(length(which(EltonTraits_IDs=="-100")), "species were not ID'd.")
EltonTraits_df$Avibase_ID <- c(EltonTraits_IDs)

# Identify species that appear in CLOVER
host_IDs_in_elton <- unique(EltonTraits_IDs[which(EltonTraits_IDs %in% host_species_IDs)])

# Check if number of species left in trait database matches number from CLOVER:
length(host_IDs_in_elton)==no_host_species

# Create column indicating whether species is a host
host_indicator = c(EltonTraits_df$Avibase_ID %in% host_species_IDs)

# Append it to trait database
matched_data <- data.frame(EltonTraits_df, host_indicator)

# Now filter for species in Europe according to eBird:
# matched_data <- matched_data[which(EltonTraits_IDs %in% eBird_IDs), ]

# Remove fields we definitely won't want for fitting

# This line might be helpful if you want to print out all the variables in a
# clean format:
for (name in names(matched_data)){
  cat('"', name, '",\n', sep="")
}

fields_to_drop <- c("SpecID",
                    "BLFamilyEnglish",
                    "BLFamSequID",
                    "Taxo",
                    "English",
                    "Diet.Source",
                    "Diet.EnteredBy",
                    "ForStrat.Source",
                    "ForStrat.EnteredBy",
                    "BodyMass.Source",
                    "BodyMass.Comment",
                    "Record.Comment")
matched_data <- matched_data[, !(names(matched_data) %in% fields_to_drop)]

# Create a version filtered for certainty
certain_data <- matched_data
certain_data <- certain_data[
  which(certain_data$Diet.Certainty %in% c("A", "B")),
]
certain_data <- certain_data[
  which(certain_data$ForStrat.SpecLevel==1),
]
certain_data <- certain_data[
  which(certain_data$BodyMass.SpecLevel==1),
]

# Check proportion of species with a good level of certainty
cat(100 * nrow(certain_data) / nrow(matched_data),
    "% of species in Europe remain if we filter for certainty.\n")

# Check what this is in terms of number of species:
cat(sum(certain_data$host_indicator),
    "of",
    sum(certain_data$host_indicator),
    "host species remain if we filter for certainty.\n")

# Some of the ecological data is expressed as percentages which add up to 100.
# This makes them colinear and so to avoid this we will remove some of the
# percentage fields, specifically the ones that have the highest average value
# within their set of fields.

diet_data <- matched_data[, c("Diet.Inv",
                              "Diet.Vend",
                              "Diet.Vect",
                              "Diet.Vfish",
                              "Diet.Vunk",
                              "Diet.Scav",
                              "Diet.Nect",
                              "Diet.Fruit",
                              "Diet.Seed",
                              "Diet.PlantO")]
max_diet_comp <- names(diet_data)[which.max(colMeans(diet_data))]

forstrat_data <- matched_data[, c("ForStrat.watbelowsurf",
                                  "ForStrat.wataroundsurf",
                                  "ForStrat.aerial",
                                  "ForStrat.canopy",
                                  "ForStrat.ground",
                                  "ForStrat.understory",
                                  "ForStrat.midhigh")]
max_forstrat_comp <- names(forstrat_data)[which.max(colMeans(forstrat_data))]

################################################################################
# Incorporate additional ecological data from IUCN

IUCN_df <- rbind(
  read.csv("eco_phylo_processing/iucn-data-vol1/all_other_fields.csv"),
  read.csv("eco_phylo_processing/iucn-data-vol2/all_other_fields.csv"))
synonym_df <- rbind(
  read.csv("eco_phylo_processing/iucn-data-vol1/synonyms.csv"),
  read.csv("eco_phylo_processing/iucn-data-vol2/synonyms.csv"))
# IUCN_df <- IUCN_df[, c("scientificName",
#                        "Congregatory.value",
#                        "MovementPatterns.pattern")]
# There is lots of redundancy in this dataset because of the way we collated it:
IUCN_df <- distinct(IUCN_df)
IUCN_df$scientificName <- sapply(IUCN_df$scientificName, tolower)
synonym_df$scientificName <- sapply(synonym_df$scientificName, tolower)
IUCN_df$Avibase_ID <- get_bird_ids(IUCN_df$scientificName)
unmatched_names <- IUCN_df$scientificName[which(IUCN_df$Avibase_ID=="-100")]
cat("A total of", length(unmatched_names), "species where not ID'd.\n")

# Try using synonyms. First reduce synonyms to ones connected to species which
# were not ID'd.
synonym_df <- synonym_df[which(synonym_df$scientificName %in% unmatched_names), ]
# Convert synonyms to lower case and first two words only
synonym_df$name <- sapply(synonym_df$name, tolower)
synonym_df$name <- sapply(synonym_df$name, FUN = function(x) {str_extract(x, "[^ ]+ [^ ]+")})

# Find out which synonyms appear in AVONET and can be renamed
in_avonet <- sapply(synonym_df$name,
                    FUN = function(x) {(x %in% AVONET_df$Species1_BirdLife)|
                        (x %in% AVONET_df$Species2_eBird)|
                        (x %in% AVONET_df$Species3_BirdTree)})
species_in_avonet <- synonym_df$scientificName[which(in_avonet)]
cat("Of the",
    length(unmatched_names),
    "species which were not ID'd, a total of",
    length(unique(species_in_avonet)),
    "have a synonym which appears in AVONET.\n")

for (i in 1:length(which(in_avonet))){
  syn <- synonym_df$name[which(in_avonet)[i]]
  sp <- species_in_avonet[i]
  IUCN_df$scientificName[which(IUCN_df$scientificName==sp)] <- syn
}
IUCN_df$Avibase_ID <- get_bird_ids(IUCN_df$scientificName)
unmatched_names <- IUCN_df$scientificName[which(IUCN_df$Avibase_ID=="-100")]
cat("After relabelling with synonyms a total of",
    length(unmatched_names),
    "species where not ID'd.\n")

# Check if any of the remaining unmatched species are in the eBird list:
matched_IDs_not_in_IUCN <- setdiff(matched_data$Avibase_ID, IUCN_df$Avibase_ID)
cat("The following species of interest are not assigned ID's in the IUCN data:",
    matched_data$Scientific[
      which(matched_data$Avibase_ID %in% matched_IDs_not_in_IUCN)],
    sep="\n")

# Should find problem species are carduelis flammea/hornemanni and himantopus
# mexicanus
# We can fix the issue with carduelis flammea by renaming in the original IUCN
# data:
IUCN_df$scientificName[which(IUCN_df$scientificName=="acanthis flammea")] <-
  "carduelis flammea"

# himantopus mexicanus is a vagrant and should have very low numbers in Europe,
# so we remove it from the matched data:
matched_data <- matched_data[
  -which(matched_data$Scientific=="himantopus mexicanus"),
]

# Now try again
IUCN_df$Avibase_ID <- get_bird_ids(IUCN_df$scientificName)
matched_IDs_not_in_IUCN <- setdiff(matched_data$Avibase_ID, IUCN_df$Avibase_ID)
cat("There are now",
    length(matched_IDs_not_in_IUCN),
    "species of interest are not assigned ID's in the IUCN data.\n")

# Inspect coding of congregatory and migratory behaviours: 
cong_vals <- unique(IUCN_df$Congregatory.value)
move_vals <- unique(IUCN_df$MovementPatterns.pattern)
cat("Possible values of Congregatory.value are",
    cong_vals,
    sep = "\n")
cat("Possible values of MovementPatterns.pattern are",
    move_vals,
    sep = "\n")

# Set congregatory indicator to true if dispersive or year-round:
IUCN_df$is_congregatory <- (IUCN_df$Congregatory.value != "")

# Set migratory indicator to true if recorded as full migrant or not known:
IUCN_df$is_migratory <- (IUCN_df$MovementPatterns.pattern == "Full Migrant") |
  (IUCN_df$MovementPatterns.pattern == "Unknown")

# Attach migratory/congregatory indicators to matched_data
# Use any function so that if IUCN has multiple species for an ID we count it as
# congregatory/migatory if any of the matching species are.
# Set congregatory indicator to true if dispersive or year-round:
matched_data$is_congregatory <- sapply(matched_data$Avibase_ID,
                                       FUN = function(x){
                                         any(IUCN_df$Congregatory.value[which(IUCN_df$Avibase_ID==x)] != "")
                                       })
# Set migratory indicator to true if recorded as full:
matched_data$is_migratory <- sapply(matched_data$Avibase_ID,
                                    FUN = function(x){
                                      any(IUCN_df$MovementPatterns.pattern[
                                        which(IUCN_df$Avibase_ID==x)
                                      ] == "Full Migrant")
                                    })

################################################################################
# Bring in phylogenetic data

# If data already exists, load it in, otherwise calculate average phylogenetic
# distances across ABC samples.

{
  if (!file.exists("eco_phylo_processing/mean_phylo_distances.rds")){
    # Attempt to load in BirdTree data
    bird_tree <- read.tree("eco_phylo_processing/BirdzillaHackett10.tre")
    
    for (i in seq_along(bird_tree)){
      bird_tree[[i]]$tip.label <- tolower(bird_tree[[i]]$tip.label)
    }
    
    no_samples <- 100
    dmat_mean <- (1 / no_samples) * cophenetic.phylo(bird_tree[[1]])
    dmat_mean <- dmat_mean[order(rownames(dmat_mean)), order(colnames(dmat_mean))]
    for (i in 2:no_samples){
      start_time <- Sys.time()
      dmat <- cophenetic.phylo(bird_tree[[i]])
      dmat <- dmat[order(rownames(dmat)), order(colnames(dmat))]
      dmat_mean <- dmat_mean + (1 / no_samples) * dmat
      end_time <- Sys.time()
      cat("Matrix", i, "processed in", end_time - start_time, ".\n")
    }
    rownames(dmat_mean) <- gsub("_", " ", rownames(dmat_mean))
    colnames(dmat_mean) <- gsub("_", " ", colnames(dmat_mean))
    
    saveRDS(dmat_mean, "eco_phylo_processing/mean_phylo_distances.rds")
  }
  else{
    dmat_mean <- readRDS("eco_phylo_processing/mean_phylo_distances.rds")
  }
}

dmat_ids <- get_bird_ids(rownames(dmat_mean))

# See if any species were not successfully ID'd
cat("There are",
    length(setdiff(unique(EltonTraits_df$Scientific),
                   unique(rownames(dmat_mean)))),
    "species names in EltonTraits missing from the phylogeny and",
    length(setdiff(unique(rownames(dmat_mean)),
                   unique(EltonTraits_df$Scientific))),
    "species in the phylogeny missing from EltonTraits.\n"
)
# Should find that the phylogeny contains exactly the same species with the same
# names as EltonTraits, so we pass in the IDs from EltonTraits
for (i in 1:nrow(dmat_mean)){
  rownames(dmat_mean)[i] <- EltonTraits_df$Avibase_ID[which(EltonTraits_df$Scientific==rownames(dmat_mean)[i])[1]]
  colnames(dmat_mean)[i] <- EltonTraits_df$Avibase_ID[which(EltonTraits_df$Scientific==colnames(dmat_mean)[i])[1]]
}

dmat_host_cols <- dmat_mean[, which(dmat_ids %in% host_species_IDs)]
nearest_host_distance <- apply(dmat_host_cols, 1, min)
for (i in 1:nrow(matched_data)){
  species_in_dmat <- which(
    names(nearest_host_distance)==matched_data$Avibase_ID[i])
  if (length(species_in_dmat)==1){
    matched_data$nearest_host_distance[i] <- nearest_host_distance[species_in_dmat]
  }
  else if (length(species_in_dmat)==0){
    matched_data$nearest_host_distance[i] <- 1000
  }
  else{
    matched_data$nearest_host_distance[i] <- min(
      nearest_host_distance[species_in_dmat])
  }
}

# Do a quick check to make sure all host species have distance zero and all
# non-host species have nonzero distance
cat("Maximum distance for host species is",
    max(matched_data$nearest_host_distance[which(matched_data$host_indicator)]),
    ".\n")
cat("Minimum distance for non-host species is",
    min(matched_data$nearest_host_distance[which(!matched_data$host_indicator)]),
    ".\n")
# Check every species has had a valid distance assigned (i.e. not 1000))
cat("There are",
    length(which(matched_data$nearest_host_distance==1000)),
    "species without distances assigned.\n")

# Zeros will be too informative since they immediately tell us where the host
# species are, so we should really use the minimum non-zero value:
nearest_host_distance <- apply(dmat_host_cols, 1, FUN = function(x) {min(x)})
for (i in 1:nrow(matched_data)){
  species_in_dmat <- which(
    names(nearest_host_distance)==matched_data$Avibase_ID[i])
  if (length(species_in_dmat)==1){
    matched_data$nearest_host_distance[i] <- nearest_host_distance[species_in_dmat]
  }
  else if (length(species_in_dmat)==0){
    matched_data$nearest_host_distance[i] <- 1000
  }
  else{
    matched_data$nearest_host_distance[i] <- min(
      nearest_host_distance[species_in_dmat])
  }
}

cat("Mean distance to nearest confirmed host for confirmed hosts is",
    mean(matched_data$nearest_host_distance[which(matched_data$host_indicator)]),
    "and",
    mean(matched_data$nearest_host_distance[-which(matched_data$host_indicator)]),
    "for other species.\n")
cat("Median distance to nearest confirmed host for confirmed hosts is",
    median(matched_data$nearest_host_distance[which(matched_data$host_indicator)]),
    "and",
    median(matched_data$nearest_host_distance[-which(matched_data$host_indicator)]),
    "for other species.\n")

# Do a quick plot to see if there's an obvious difference in the distances for
# confirmed hosts:
if (PLOT){
  bin_size <- 10
  
  plot_ulim <- bin_size * ceiling(max(matched_data$nearest_host_distance) / bin_size)
  
  host_dhist <- hist(matched_data$nearest_host_distance[which(matched_data$host_indicator)],
                     breaks = seq(0, plot_ulim, by=bin_size),
                     plot = FALSE)
  nonhost_dhist <- hist(matched_data$nearest_host_distance[-which(matched_data$host_indicator)],
                        breaks = seq(0, plot_ulim, by=bin_size),
                        plot = FALSE)
  pal <- brewer.pal(6, "Dark2")
  plot(seq(0, plot_ulim-bin_size, by=bin_size),
       host_dhist$density,
       type = "l",
       lwd = 2,
       col = pal[1],
       xlab = "Distance to nearest\n confirmed host species",
       ylab = "Density")
  lines(seq(0, plot_ulim-bin_size, by=bin_size),
        nonhost_dhist$density,
        type = "l",
        lwd = 2,
        col = pal[2])
  legend("topright",
         inset=c(.025,0),
         legend = c("Confirmed hosts", "Other species"),
         col = pal,
         bty = "n",
         pch = 20,
         pt.cex = 2,
         cex = 0.8,
         xpd = TRUE)
}

################################################################################
# In this section we do some cleanup on genus/species names to make sure all our
# abundance data corresponds to something in matched_data

# Start by identifying all species that don't have either an unambiguous ID
# match of an unambiguous name match
problem_species <- which(sapply(1:nrow(sp_df),
                                FUN = function(i){
  (length(which(matched_data$Avibase_ID == sp_df$Avibase_ID[i]))!=1) &
  (length(which(matched_data$Scientific == sp_df$scientific_name[i]))!=1)}))

cat("The following are identified as problem species:",
    array(sp_df$scientific_name[problem_species]),
    sep = "\n"
    )

# Work through this list case-by-case
sp_df$scientific_name <- sub("melanitta deglandi",
                             "melanitta fusca",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("melanitta americana",
                             "melanitta nigra",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("himantopus mexicanus",
                             "himantopus himantopus",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("gallinago delicata",
                             "gallinago gallinago",
                             sp_df$scientific_name)
matched_data$Scientific <- sub("sterna nilotica",
                               "gelochelidon nilotica",
                               matched_data$Scientific)
matched_data$Scientific <- sub("sterna maxima",
                               "thalasseus maximus",
                               matched_data$Scientific)
sp_df$scientific_name <- sub("ardea intermedia",
                             "mesophoyx intermedia",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("circus hudsonius",
                             "circus cyaneus",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("picus sharpei",
                             "picus viridis",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("lanius borealis",
                             "lanius excubitor",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("cyanopica cooki",
                             "cyanopica cyanus",
                             sp_df$scientific_name)
matched_data$Scientific <- sub("parus lugubris",
                             "poecile lugubris",
                             matched_data$Scientific)
matched_data$Scientific <- sub("hirundo daurica",
                               "cecropis daurica",
                               matched_data$Scientific)
matched_data$Scientific <- sub("parus montanus",
                               "poecile montanus",
                               matched_data$Scientific)
sp_df$scientific_name <- sub("phylloscopus examinandus",
                             "phylloscopus borealis",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("saxicola rubicola",
                             "saxicola torquatus",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("saxicola maurus",
                             "saxicola torquatus",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("saxicola stejnegeri",
                             "saxicola torquatus",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("oenanthe melanoleuca",
                             "oenanthe hispanica",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("passer italiae",
                             "passer domesticus",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("motacilla tschutschensis",
                             "motacilla flava",
                             sp_df$scientific_name)
matched_data$Scientific <- sub("sturnus contra",
                               "gracupica contra",
                               matched_data$Scientific)
sp_df$scientific_name <- sub("curruca crassirostris",
                             "curruca hortensis",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("curruca subalpina",
                             "curruca cantillans",
                             sp_df$scientific_name)
sp_df$scientific_name <- sub("curruca balearica",
                             "curruca sarda",
                             sp_df$scientific_name)

# The taxonomy of Curruca proves to be quite problematic - 25 species were moved
# from Sylvia to Curruca in 2011 so we change the names in matched_data to
# reflect this

# Get species names only:
curruca_names <- sp_df$scientific_name %>%
                  str_subset(pattern = "curruca") %>%
                  sub(pattern = "curruca ", replacement = "")
names_to_replace <- lapply(curruca_names,
                           FUN = function(x){paste("sylvia ", x, sep = "")})
matched_data$Scientific[which(matched_data$Scientific %in% names_to_replace)] <-
  sub("sylvia",
      "curruca",
      matched_data$Scientific[
        which(matched_data$Scientific %in% names_to_replace)])

# Also need to change a few carduelis species to acanthis:
acanthis_names <- sp_df$scientific_name %>%
  str_subset(pattern = "acanthis") %>%
  sub(pattern = "acanthis ", replacement = "")
names_to_replace <- lapply(acanthis_names,
                           FUN = function(x){paste("carduelis ", x, sep = "")})
matched_data$Scientific[which(matched_data$Scientific %in% names_to_replace)] <-
  sub("carduelis",
      "acanthis",
      matched_data$Scientific[
        which(matched_data$Scientific %in% names_to_replace)])

# eBird is the only source to include acanthis cabaret, which was formerly
# considered a subspecies of acanthis flammea, so we treat counts of cabaret as
# counts of flammea
sp_df$scientific_name <- sub("acanthis cabaret",
             "acanthis flammea",
             sp_df$scientific_name)

# Another species adjustment, this time for anas zonorhyncha and anas
# poecilorhyncha
sp_df$scientific_name <- sub("anas zonorhyncha",
                             "anas poecilorhyncha",
                             sp_df$scientific_name)


# A genus name for American warblers was changed in early 2010a:
matched_data$Scientific <- sub("dendroica",
                               "setophaga",
                               matched_data$Scientific)

# See if any problem species remain:
problem_species <- which(sapply(1:nrow(sp_df),
                                FUN = function(i){
    (length(which(matched_data$Avibase_ID == sp_df$Avibase_ID[i]))!=1) &
      (length(which(matched_data$Scientific == sp_df$scientific_name[i]))!=1)}))
cat("There are now", length(problem_species), "problem species.\n")

################################################################################
# Now get abundance rasters for species and convolute with eco/phylo factors

# set the crs we want to use
crs <- "epsg:3035"

euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later

# Create a blank raster with appropriate projection and extent
blank_3035 <- rast(crs=crs, extent=euro_ext, res=9042.959)

# Create blank rasters for each of our quantities
# We will create rasters corresponding to the following traits:
# Congregatory behaviour
# Migatory behaviour
# Foraging around water surface
# Foraging >5cm below water surface
# Phylogenetic distance to a confirmed host

nlyrs <- 52

cong_rast <- rast(nlyrs=nlyrs,
                  crs=crs,
                  extent=euro_ext,
                  res=9042.959)
migr_rast <- rast(nlyrs=nlyrs,
                  crs=crs,
                  extent=euro_ext,
                  res=9042.959)
around_surf_rast <- rast(nlyrs=nlyrs,
                         crs=crs,
                         extent=euro_ext,
                         res=9042.959)
below_surf_rast <- rast(nlyrs=nlyrs,
                        crs=crs,
                        extent=euro_ext,
                        res=9042.959)

for (i in 1:nlyrs){
  cong_rast[, , i] <- 0
  migr_rast[, , i] <- 0
  around_surf_rast[, , i] <- 0
  below_surf_rast[, , i] <- 0
}

no_species <- nrow(sp_df)
sample_idx <- 1:no_species

# # Remove unmatched species from sp_df
# sp_df <- sp_df[which(sp_df$Avibase_ID %in% matched_data$Avibase_ID), ]

# set.seed(1)
# sample_idx <- sample(1:nrow(sp_df), no_species)

# Make sure timeout is set high enough so we can actually download the data
# Ten minutes per dataset should be enough
if (getOption('timeout') < 10 * 60 * no_species){
  options(timeout = 10 * 60 * no_species)
  cat("Setting timeout to",
      (1 / 60) * getOption('timeout'),
      "minutes.\n")
}

# Get number of species already downloaded by searching for codes in ebirdst
# data directory
starting_dls <- ebirdst_data_dir() %>% 
  paste("/2021", sep = "") %>% 
  list.files()
no_downloaded <- starting_dls %in% sp_df$species_code %>% 
  which() %>%
  length()
total_to_download <- no_species - no_downloaded
no_downloaded <- 0

not_downloaded <- setdiff((sp_df$species_code),
                          (data.frame(starting_dls)$starting_dls))
idx_to_download <- which(sp_df$species_code %in% not_downloaded)

# Loop over other first no_species species:
{
  loop.start <- Sys.time()
  mean_dl_time <- 0
  for (i in 1:no_species) {
    idx <- sample_idx[i]
    species_sel <- sp_df$species_code[idx]
    Avibase_ID <- sp_df$Avibase_ID[idx]
    species_factors <- matched_data[which(matched_data$Avibase_ID==Avibase_ID), ]
    if (nrow(species_factors)==0){
      cat("Species ID for",
          sp_df$scientific_name[idx],
          "is not in matched_data. Using scientific name to match.\n",
          sep = "\n")
      scientific_name <- sp_df$scientific_name[idx]
      species_factors <- matched_data[which(matched_data$Scientific==scientific_name), ]
    }
    if (nrow(species_factors)>1){
      cat("Same species ID for",
          species_factors$Scientific,
          "Using scientific name to match.\n",
          sep = "\n")
      scientific_name <- sp_df$scientific_name[idx]
      species_factors <- matched_data[which(matched_data$Scientific==scientific_name), ]
    }
    
    if (Avibase_ID %in% pop_size_ids){
      pop_size_scaler <- pop_size_df$Abundance.estimate[which()]
    }
    
    if (!(species_sel %in% starting_dls)){
      dl_start <- Sys.time()
      dl_flag <- TRUE
      attempt_count <- 0
      while (dl_flag){
        attempt_count <- attempt_count + 1
        cat("attempt = ", attempt_count, ".\n")
        path <- try(ebirdst_download(species = species_sel,
                                     pattern = "abundance_median_lr"))
        if (!inherits(path, "try-error")){
          dl_flag <- FALSE
        }
        if (attempt_count>50){
          print("Download failed")
          dl_flag <- FALSE
        }
      }
      no_downloaded <- no_downloaded + 1
      time.now <- Sys.time()
      elapsed <- as.numeric(difftime(time.now, dl_start, units = "mins"))
      mean_dl_time <- mean_dl_time + (1 / no_downloaded) * elapsed
      time_remaining <- (total_to_download - no_downloaded) * mean_dl_time
      cat(as.numeric(difftime(time.now, loop.start, units="mins")),
          " minutes elapsed since start, estimated ",
          time_remaining,
          " remaining.",
          no_downloaded,
          "of",
          total_to_download,
          "species downloaded\n")
    }
    else{
      path <- ebirdst_download(species = species_sel,
                               pattern = "abundance_median_lr")
    }
    this_rast <- load_raster(path = path,
                             product = "abundance",
                             period = "weekly",
                             resolution = "lr")
    this_rast <- project(x = this_rast, y = blank_3035, method = "near")
    # Get rid of NA's:
    this_rast <- replace(this_rast, is.na(this_rast), 0)
    
    # Fill in each field
    if (species_factors$is_congregatory){
      cong_rast <- cong_rast + this_rast
    }
    if (species_factors$is_migratory){
      migr_rast <- migr_rast + this_rast
    }
    around_surf_rast <- around_surf_rast +
      species_factors$ForStrat.wataroundsurf * this_rast
    below_surf_rast <- below_surf_rast +
      species_factors$ForStrat.watbelowsurf * this_rast
  }
}

av_capture_graphics(animate(cong_rast, n=1),
                    framerate = 5.2,
                    output = "cong_rast.mp4")
av_capture_graphics(animate(migr_rast, n=1),
                    framerate = 5.2,
                    output = "migr_rast.mp4")
av_capture_graphics(animate(around_surf_rast, n=1),
                    framerate = 5.2,
                    output = "around_surf_rast.mp4")
av_capture_graphics(animate(below_surf_rast, n=1),
                    framerate = 5.2,
                    output = "below_surf_rast.mp4")
