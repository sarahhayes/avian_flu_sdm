# In this script we combine abundance estimates from eBird with species-specific
# ecological and phylogenetic factors to create rasters encoding the abundance
# of birds with different qualities.

PLOT <- TRUE # Make plots/animations
HQ_ONLY <- TRUE # Determines whether to use only species with high-quality abundance data in eBird
QUARTERLY <- TRUE # Generate quarterly (weeks 1-13 etc) rasters, otherwise do weekly
SAVE_RAST <- FALSE # Save output rasters

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
library(pracma)
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
sp_df <- sp_df[, c(
  "species_code",
  "scientific_name",
  "common_name",
  "breeding_quality",
  "nonbreeding_quality",
  "postbreeding_migration_quality",
  "prebreeding_migration_quality",
  "resident_quality"
)]

# Get codes for species in Europe
euro_bird_codes <- read.csv("ebird/codes_for_europe_clean.csv")

no_euro_birds <- nrow(euro_bird_codes)

#Restrict to species in Europe list
sp_df <- sp_df[which(sp_df$species_code %in% euro_bird_codes$code), ]

# Load in population sizes
pop_size_df <- read.csv("data/eco_phylo_data/callaghan_pop_estimates.csv")
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
AVONET_df <- read_excel("data/eco_phylo_data/AVONETSupplementarydataset1.xlsx",
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
  synonyms <- c(binomial_name)
  if (("BirdLife" %in% name_sources) & 
      (binomial_name %in% AVONET_df$Species1_BirdLife)){
    ID_candidates <- append(ID_candidates, AVONET_df$Avibase.ID[
      which(AVONET_df$Species1_BirdLife==binomial_name)])
    synonyms <- append(synonyms, AVONET_df$Species2_eBird[
      which(AVONET_df$Species1_BirdLife==binomial_name)])
    synonyms <- append(synonyms, AVONET_df$Species3_BirdTree[
      which(AVONET_df$Species1_BirdLife==binomial_name)])
  }
  if (("eBird" %in% name_sources) & 
      (binomial_name %in% AVONET_df$Species2_eBird)){
    ID_candidates <- append(ID_candidates, AVONET_df$Avibase.ID[
      which(AVONET_df$Species2_eBird==binomial_name)])
    synonyms <- append(synonyms, AVONET_df$Species1_BirdLife[
      which(AVONET_df$Species2_eBird==binomial_name)])
    synonyms <- append(synonyms, AVONET_df$Species3_BirdTree[
      which(AVONET_df$Species2_eBird==binomial_name)])
  }
  if (("BirdTree" %in% name_sources) & 
      (binomial_name %in% AVONET_df$Species3_BirdTree)){
    ID_candidates <- append(ID_candidates, AVONET_df$Avibase.ID[
      which(AVONET_df$Species3_BirdTree==binomial_name)])
    synonyms <- append(synonyms, AVONET_df$Species1_BirdLife[
      which(AVONET_df$Species3_BirdTree==binomial_name)])
    synonyms <- append(synonyms, AVONET_df$Species2_eBird[
      which(AVONET_df$Species3_BirdTree==binomial_name)])
  }
  # How we work out the consensus depends on how many candidates are identified
  if (length(ID_candidates)==0){
    # If ID doesn't appear in any list, return -100
    ID_candidates <- c("-100")
  }
  # Remove NA's from synonyms
  synonyms <- synonyms[which(synonyms!="na")]
  return(list(ID = unique(ID_candidates), syn = unique(synonyms)))
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
  syn_df <- data.frame(matrix(nrow = length(name_vect), ncol = (cols_needed+1)))
  for (n in 1:length(name_vect)){
    # cat("n=",n,"\n")
    vect_list <- get_single_bird_id(name_vect[n], name_sources)
    id_df[n, 2:(length(vect_list$ID)+1)] <- vect_list$ID
    syn_df[n, 1:length(vect_list$syn)] <- vect_list$syn
  }
  return(list(id_df, syn_df))
}

# Find IDs for species in eBird
id_and_syn <- get_bird_ids(sp_df$scientific_name, name_sources = "eBird")
eBird_IDs <- id_and_syn[[1]]
eBird_syns <- id_and_syn[[2]]
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
id_and_syn <- get_bird_ids(sp_df$scientific_name, name_sources = "eBird")
eBird_IDs <- id_and_syn[[1]]
eBird_syns <- id_and_syn[[2]]
# Check we got everything:
cat(length(which(eBird_IDs[, 2]=="-100")), "species were not ID'd.")

# Add ID options to species dataframe
sp_df$Avibase_ID <- eBird_IDs[, 2:ncol(eBird_IDs)]
sp_df$synonyms <- eBird_syns

################################################################################
# Now introduce CLOVER database and identify all avian species known to be hosts
# of influenza A.

# Load data
CLOVER_df <- read.csv(
  "data/eco_phylo_data/CLOVER_1.0_Viruses_AssociationsFlatFile.csv")
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

species_counts <- CLOVER_df %>% dplyr::count(Host, sort = TRUE)

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
order_case_counts <- CLOVER_df %>% dplyr::count(HostOrder, sort = TRUE)
family_case_counts <- CLOVER_df %>% dplyr::count(HostFamily, sort = TRUE)
genus_case_counts <- CLOVER_df %>% dplyr::count(HostGenus, sort = TRUE)

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
host_species_IDs <- get_bird_ids(species_list)[[1]]

# See if any species were not successfully ID'd
cat(length(which(host_species_IDs[, 2]=="-100")), "species were not ID'd.")
cat("Unidentified species are:",
    species_list[which(host_species_IDs[, 2]=="-100")],
    sep = "\n")

# We manually replace these three species names with names found in AVONET.
# These are all unambiguous replacements - they are present in all three source
# datasets and associated with a single Avibase ID.
species_list[which(species_list=="coloeus monedula")] <- "corvus monedula"
species_list[which(species_list=="larus mongolicus")] <- "larus cachinnans"
species_list[which(species_list=="piaya minuta")] <- "coccycua minuta"

# Try again with ID's
host_species_IDs <- get_bird_ids(species_list)[[1]]

# We now shouldn't see any failed ID's
cat(length(which(host_species_IDs=="-100")), "species were not ID'd.")

# Remove possible repeats:
host_species_IDs <- unique(host_species_IDs)
no_host_species <- length(host_species_IDs)

host_indicator <- sapply(1:nrow(sp_df),
                         FUN = function(i){
                           (length(intersect(sp_df$synonyms[i, ], CLOVER_df$Host)) + 
                              length(intersect(sp_df$Avibase_ID[i, ], host_species_IDs))) > 0
                         })
sp_df$host_indicator <- host_indicator

################################################################################
# Now introduce trait data

# Load in from file
EltonTraits_df <- 
  read.table("data/variables/bird_data/elton_traits/BirdFuncDat.txt",
             sep = '\t',
             quote="\"",
             header = TRUE)
# elton trait data has two rows of NA's at the bottom which we need to remove:
EltonTraits_df <- EltonTraits_df[1:9993, ]
EltonTraits_df$Scientific <- tolower(EltonTraits_df$Scientific)

# Remove unnecessary fields:
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
EltonTraits_df <- EltonTraits_df[, !(names(EltonTraits_df) %in% fields_to_drop)]

# Some of the ecological data is expressed as percentages which add up to 100.
# This makes them colinear and so to avoid this we will remove some of the
# percentage fields, specifically the ones that have the highest average value
# within their set of fields.

diet_data <- EltonTraits_df[, c("Diet.Inv",
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

forstrat_data <- EltonTraits_df[, c("ForStrat.watbelowsurf",
                                    "ForStrat.wataroundsurf",
                                    "ForStrat.aerial",
                                    "ForStrat.canopy",
                                    "ForStrat.ground",
                                    "ForStrat.understory",
                                    "ForStrat.midhigh")]
max_forstrat_comp <- names(forstrat_data)[which.max(colMeans(forstrat_data))]

# Search for IDs
EltonTraits_IDs <- get_bird_ids(EltonTraits_df$Scientific)[[1]]

# Check we got everything:
cat(length(which(EltonTraits_IDs[, 2]=="-100")), "species were not ID'd.")
non_IDd_species <- EltonTraits_df$Scientific[which(EltonTraits_IDs[, 2]=="-100")]
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
EltonTraits_IDs <- get_bird_ids(EltonTraits_df$Scientific)[[1]]
cat(length(which(EltonTraits_IDs[, 2]=="-100")), "species were not ID'd.")

# Check that we can make at least one assignment for all species in eBird:
has_eco_data <- sapply(1:nrow(sp_df),
                       FUN = function(i){
                         (length(intersect(sp_df$synonyms[i, ], EltonTraits_df$Scientific))) > 0
                       })
cat("Elton traits can not be matched to",
    length(which(!has_eco_data)),
    "species in eBird using names only.")
# Should find there are no missing species!

# See which (if any) have ambiguous matches:
multiple_Elton_matches <- sapply(1:nrow(sp_df),
                                 FUN = function(i){
                                   (length(intersect(sp_df$synonyms[i, ], EltonTraits_df$Scientific))>1)
                                 })
cat("There are multiple matches for",
    length(which(multiple_Elton_matches)),
    "species in eBird.")
cat("Ambiguous species are:",
    sp_df$scientific_name[which(multiple_Elton_matches)],
    sep = "\n")

for (idx in which(multiple_Elton_matches)){
  cat("eBird name is ",
      sp_df$scientific_name[idx],
      ", Elton traits matches are ",
      sep = "")
  cat(EltonTraits_df$Scientific[which(EltonTraits_df$Scientific %in% sp_df$synonyms[idx, ])], sep = ", ")
  cat("\n")
}
# Should find for most problem species one of the synonyms is just the name
# from eBird. Only real problem species is Yellow Warbler, which we fix below:
EltonTraits_df$Scientific[
  which(EltonTraits_df$Scientific=="dendroica petechia")
] <- "setophaga petechia"

# Find position of each species from eBird in EltonTraits:
position_in_EltonTraits <- sapply(
  1:nrow(sp_df),
  FUN = function(i){
    if (sp_df$scientific_name[i] %in% EltonTraits_df$Scientific){
      return(which(EltonTraits_df$Scientific == sp_df$scientific_name[i]))
    }
    else{
      return(which(EltonTraits_df$Scientific %in% sp_df$synonyms[i, ]))
    }
  }
)

# Check to make sure name matching is correct:
names_in_Elton <- EltonTraits_df$Scientific[position_in_EltonTraits]
matched_locs <- sapply(1:nrow(sp_df),
                       FUN = function(i){
                         (names_in_Elton[i] %in% sp_df$synonyms[i, ])
                       })
cat("Name matching fails at", length(which(!matched_locs)), "locations.")

traits_by_eBird_species <- EltonTraits_df[position_in_EltonTraits, ]
sp_df$EltonTraits <- traits_by_eBird_species

################################################################################
# Bring in phylogenetic data

# If data already exists, load it in, otherwise calculate average phylogenetic
# distances across ABC samples.

{
  if (!file.exists("data/eco_phylo_data/mean_phylo_distances.rds")){
    # Attempt to load in BirdTree data
    bird_tree <- read.tree("data/eco_phylo_data/BirdzillaHackett10.tre")
    
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
    
    saveRDS(dmat_mean, "data/eco_phylo_data/mean_phylo_distances.rds")
  }
  else{
    dmat_mean <- readRDS("data/eco_phylo_data/mean_phylo_distances.rds")
  }
}

phylo_sp <- rownames(dmat_mean)
phylo_syns_and_IDs <- get_bird_ids(phylo_sp)
phylo_IDs <- phylo_syns_and_IDs[[1]]
phylo_syns <- phylo_syns_and_IDs[[2]]

# Identify host species in phylo data:
hosts_in_phylo <- sapply(1:length(phylo_sp),
                         FUN = function(i){
                           (length(intersect(phylo_syns[i, ], CLOVER_df$Host)) + 
                              length(intersect(phylo_IDs[i, ], host_species_IDs))) > 0
                         })

dmat_host_cols <- dmat_mean[, which(hosts_in_phylo)]
nearest_host_distance <- apply(dmat_host_cols, 1, min)
nearest_host_df <- data.frame(phylo_sp, nearest_host_distance)

# Now try to match up with eBird:

# Check that we can make at least one assignment for all species in eBird:
has_phylo_data <- sapply(1:nrow(sp_df),
                         FUN = function(i){
                           (length(intersect(sp_df$synonyms[i, ], nearest_host_df$phylo_sp))) > 0
                         })
cat("Phylogeny can not be matched to",
    length(which(!has_phylo_data)),
    "species in eBird using names only.")
# Should find there are no missing species!

# These species aren't in eBird but are needed for an EltonTraits-BirdTree
# matchup later:
rownames(dmat_mean)[
  which(rownames(dmat_mean)=="anthus longicaudatus")
] <- "anthus vaalensis"
rownames(dmat_mean)[
  which(rownames(dmat_mean)=="polioptila clementsi")
] <- "polioptila guianensis"
rownames(dmat_mean)[
  which(rownames(dmat_mean)=="hypositta perdita")
] <- "oxylabes madagascariensis"

# This next species just seems to be missing from AVONET so we remove it:
dmat_mean <- dmat_mean[
  -which(rownames(dmat_mean)=="amaurospiza carrizalensis"),
  -which(rownames(dmat_mean)=="amaurospiza carrizalensis") ]

rownames(dmat_mean)[
  which(rownames(dmat_mean)=="lophura hatinhensis")
] <- "lophura edwardsi"

# See which (if any) have ambiguous matches:
multiple_phylo_matches <- sapply(1:nrow(sp_df),
                                 FUN = function(i){
                                   (length(intersect(sp_df$synonyms[i, ], nearest_host_df$phylo_sp))>1)
                                 })
cat("There are multiple matches for",
    length(which(multiple_phylo_matches)),
    "species in eBird.")
cat("Ambiguous species are:",
    sp_df$scientific_name[which(multiple_phylo_matches)],
    sep = "\n")

# Should find it's the same problem species as for EltonTraits
identical(multiple_Elton_matches, multiple_phylo_matches)

nearest_host_df$phylo_sp[
  which(nearest_host_df$phylo_sp=="dendroica petechia")
] <- "setophaga petechia"

rownames(dmat_mean)[
  which(rownames(dmat_mean)=="dendroica petechia")
] <- "setophaga petechia"

# Find position of each species from eBird in phylo data:
position_in_phylo <- sapply(
  1:nrow(sp_df),
  FUN = function(i){
    if (sp_df$scientific_name[i] %in% nearest_host_df$phylo_sp){
      return(which(nearest_host_df$phylo_sp == sp_df$scientific_name[i]))
    }
    else{
      return(which(nearest_host_df$phylo_sp %in% sp_df$synonyms[i, ]))
    }
  }
)

phylo_dist_by_eBird_species <- nearest_host_df[position_in_phylo, ]
sp_df$nearest_host_distance <- phylo_dist_by_eBird_species$nearest_host_distance

# Set threshold to use for deciding a species is a risk factor
nhost_thresh <- 50.

################################################################################
# Filter for eBird data quality

sp_df$min_quality <- sapply(1:nrow(sp_df),
                            FUN = function(i){
                              min_qual <- sp_df[i,
                                                c("breeding_quality",
                                                  "nonbreeding_quality",
                                                  "postbreeding_migration_quality",
                                                  "prebreeding_migration_quality",
                                                  "resident_quality")] %>%
                                as.numeric %>%
                                min(na.rm = TRUE)
                              return(min_qual)
                            })

if (HQ_ONLY){
  sp_df <- sp_df[which(sp_df$min_quality > 1), ]
}

################################################################################
# Experiment: see if distance between species is always same up to order

orders_in_dmat <- sapply(rownames(dmat_mean),
                         FUN = function(x){
                           EltonTraits_df$IOCOrder[which(EltonTraits_df$Scientific==x)[1]]})
unique_ords <- unique(orders_in_dmat)
no_orders = n_distinct(orders_in_dmat)
order_idx_pairs <- meshgrid(1:no_orders, 1:no_orders)
order_idx_pairs$X <- c(order_idx_pairs$X)
order_idx_pairs$Y <- c(order_idx_pairs$Y)
order_combs <- data.frame(sp1=orders_in_dmat[order_idx_pairs$X],
                          sp2=orders_in_dmat[order_idx_pairs$Y])
order_dmat_min <- matrix(0, no_orders, no_orders)
rownames(order_dmat_min) <- unique_ords
colnames(order_dmat_min) <- unique_ords
order_dmat_max <- matrix(0, no_orders, no_orders)
rownames(order_dmat_max) <- unique_ords
colnames(order_dmat_max) <- unique_ords
for (i in 2:no_orders){
  cat("i=", unique_ords[i], "\n")
  for (j in 1:i-1){
    cat("j=", unique_ords[j], "\n")
    pairwise_dist_mat <- dmat_mean[which(orders_in_dmat==unique_ords[i]),
                           which(orders_in_dmat==unique_ords[j])]
    order_dmat_min[i, j] <- min(pairwise_dist_mat)
    order_dmat_min[j, i] <- min(pairwise_dist_mat)
    order_dmat_max[i, j] <- max(pairwise_dist_mat)
    order_dmat_max[j, i] <- max(pairwise_dist_mat)
  }
}

# Get difference between min and max distances by order:
order_dmat_diff <- order_dmat_max - order_dmat_min

# Should find there is one pair of orders where max and min don't match up:
problem_idx <- which(order_dmat_diff>1e-8, arr.ind=TRUE)
order_dmat_diff[problem_idx]

# Discrepancy is in distances restricted to hummingbirds and nightjars - these
# are known to be close and the discrepancy is fairly small, so we're safe
# assuming distances are uniform up to order level.
# So, to do max phylo span we just need to look at orders


# Plot demonstrating tree-hierarchical nature of distance values:
anser_dist_df <- data.frame(name=rownames(order_dmat_max),
                            distance=order_dmat_max["Anseriformes", ])

if (PLOT){
  # Plot order counts:
  p <- ggplot(data=anser_dist_df, aes(x=name, y=distance)) + 
    geom_bar(stat = "identity")
  p + coord_flip() + 
    scale_x_discrete(limits=rownames(anser_dist_df)) + 
    xlab("Order") +
    ylab("Distance")
}