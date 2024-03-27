# In this script we analyse species-level patterns regarding avian influenza
# host status. Rather than restricting to those species which appear in eBird,
# we attempt to use all the species which appear in EltonTraits.

# Some globals to decide if we want to do plots and full variable importance
# experiments:
PLOT <- TRUE
VARIMPS <- FALSE

library(av)
require(ape)
library(BART)
library(DALEX)
library(data.table)
library(dplyr)
library(ebirdst)
library(embarcadero)
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

# Bring in AVONET to see if we can fix the missing species using AVONET IDs
AVONET_df <- read_excel("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/eco_phylo_data/AVONETSupplementarydataset1.xlsx",
                        sheet = "AVONET_Raw_Data")
AVONET_df <- AVONET_df[,
                       c("Avibase.ID",
                         "Species1_BirdLife",
                         "Species2_eBird",
                         "Species3_BirdTree")]
AVONET_df <- distinct(AVONET_df)
AVONET_df[, 2:4] <- sapply(AVONET_df[, 2:4], tolower)

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

EltonTraits_df <- EltonTraits_df[,!names(EltonTraits_df) %in% c(max_diet_comp,
                                                            max_forstrat_comp)]

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
id_and_syn <- get_bird_ids(EltonTraits_df$Scientific)
EltonTraits_IDs <- id_and_syn[[1]]
EltonTraits_syns <- id_and_syn[[2]]
# Check we got everything:
cat(length(which(EltonTraits_IDs[, 2]=="-100")), "species were not ID'd.")

sp_df <- EltonTraits_df

# Add ID options to species dataframe
sp_df$Avibase_ID <- EltonTraits_IDs[, 2:ncol(EltonTraits_IDs)]
sp_df$synonyms <- EltonTraits_syns

################################################################################
# Now introduce CLOVER database and identify all avian species known to be hosts
# of influenza A.

# Load data
CLOVER_df <- read.csv(
  "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/eco_phylo_data/CLOVER_1.0_Viruses_AssociationsFlatFile.csv")
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
    ylab("Count")
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
    ylab("Count")
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
    ylab("Count")
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
                           (length(base::intersect(sp_df$synonyms[i, ], CLOVER_df$Host)) + 
                              length(base::intersect(sp_df$Avibase_ID[i, ], host_species_IDs))) > 0
                         })
sp_df$host_indicator <- host_indicator

################################################################################
# Incorporate additional ecological data from IUCN

IUCN_df <- rbind(
  read.csv("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/eco_phylo_data/iucn-data-vol1/all_other_fields.csv"),
  read.csv("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/eco_phylo_data/iucn-data-vol2/all_other_fields.csv"))
IUCN_df <- distinct(IUCN_df)
IUCN_df$scientificName <- sapply(IUCN_df$scientificName, tolower)

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

IUCN_df <- IUCN_df[, c("scientificName",
                       "is_congregatory",
                       "is_migratory")]

# Consolidate lines that differ only in fields we don't care about:
IUCN_df <- distinct(IUCN_df)

# Check if we can do IUCN to eBird matches by binomial name:
has_IUCN_data <- sapply(1:nrow(sp_df),
                        FUN = function(i){
                          (length(base::intersect(sp_df$synonyms[i, ], IUCN_df$scientificName))) > 0
                        })
cat("There are",
    length(which(!has_IUCN_data)),
    "species in eBird which can not be matched to the IUCN data.")
cat("Missing species are",
    sp_df$Scientific[which(!has_IUCN_data)],
    sep="\n")

# To save time we just ignore problem items
sp_df <- sp_df[has_IUCN_data, ]

multiple_IUCN_matches <- sapply(1:nrow(sp_df),
                                FUN = function(i){
                                  (length(base::intersect(sp_df$synonyms[i, ], IUCN_df$scientificName))) > 1
                                })
cat("There are",
    length(which(multiple_IUCN_matches)),
    "species in eBird with multiple matches in the IUCN data.")

# Reduce to cases where we can't let binomial name take priority:
problem_cases <- sapply(1:nrow(sp_df),
                        FUN = function(i){
                          (!(sp_df$Scientific[i] %in% IUCN_df$scientificName)) &
                            ((length(base::intersect(sp_df$synonyms[i, ], IUCN_df$scientificName))) > 1)
                        })
cat("There are",
    length(which(problem_cases)),
    "ambiguous cases which can not be resolved with default name.")
cat("Ambiguous species are:",
    sp_df$Scientific[which(problem_cases)],
    sep = "\n")

# Again, I'm just removing these rather than spending time correcting them
sp_df <- sp_df[-which(problem_cases), ]

# Find position of each species from eBird in IUCN data:
position_in_IUCN <- sapply(
  1:nrow(sp_df),
  FUN = function(i){
    if (sp_df$Scientific[i] %in% IUCN_df$scientificName){
      return(which(IUCN_df$scientificName == sp_df$Scientific[i]))
    }
    else{
      return(which(IUCN_df$scientificName %in% sp_df$synonyms[i, ]))
    }
  }
)

IUCN_by_species <- IUCN_df[position_in_IUCN, ]
sp_df <- cbind(sp_df, IUCN_by_species)

################################################################################
# Bring in phylogenetic data

# If data already exists, load it in, otherwise calculate average phylogenetic
# distances across ABC samples.

{
  if (!file.exists("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/eco_phylo_data/mean_phylo_distances.rds")){
    # Attempt to load in BirdTree data
    bird_tree <- read.tree("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/eco_phylo_data/BirdzillaHackett10.tre")
    
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
    
    saveRDS(dmat_mean, "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/eco_phylo_data/mean_phylo_distances.rds")
  }
  else{
    dmat_mean <- readRDS("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/eco_phylo_data/mean_phylo_distances.rds")
  }
}

phylo_sp <- rownames(dmat_mean)
phylo_syns_and_IDs <- get_bird_ids(phylo_sp)
phylo_IDs <- phylo_syns_and_IDs[[1]]
phylo_syns <- phylo_syns_and_IDs[[2]]

# Identify host species in phylo data:
hosts_in_phylo <- sapply(1:length(phylo_sp),
                         FUN = function(i){
                           (length(base::intersect(phylo_syns[i, ], CLOVER_df$Host)) + 
                              length(base::intersect(phylo_IDs[i, ], host_species_IDs))) > 0
                         })

dmat_host_cols <- dmat_mean[, which(hosts_in_phylo)]
nearest_host_distance <- apply(dmat_host_cols,
                               1,
                               FUN = function(x) {min(x[x>0])})
nearest_host_df <- data.frame(phylo_sp, nearest_host_distance)

# Now try to match up with EltonTraits:

# Check that we can make at least one assignment for all species in EltonTraits:
has_phylo_data <- sapply(1:nrow(sp_df),
                         FUN = function(i){
                           (length(base::intersect(sp_df$synonyms[i, ], nearest_host_df$phylo_sp))) > 0
                         })
cat("Phylogeny can not be matched to",
    length(which(!has_phylo_data)),
    "species using names only.")
# Should find there are no missing species!

# In fact, EltonTraits and BirdTree use the same taxonymy so we should be able
# to match by name:
has_phylo_data <- sapply(1:nrow(sp_df),
                         FUN = function(i){
                           (sp_df$Scientific[i] %in% nearest_host_df$phylo_sp)
                         })
cat("Phylogeny can not be matched to",
    length(which(!has_phylo_data)),
    "species using names only.")

# Find position of each species from EltonTraits in phylo data:
position_in_phylo <- sapply(
  1:nrow(sp_df),
  FUN = function(i){
      return(which(nearest_host_df$phylo_sp == sp_df$Scientific[i]))
  }
)

phylo_dist_by_species <- nearest_host_df[position_in_phylo, ]
sp_df$nearest_host_distance <- phylo_dist_by_species$nearest_host_distance

cat("Mean distance to nearest confirmed host for confirmed hosts is",
    mean(sp_df$nearest_host_distance[which(sp_df$host_indicator)]),
    "and",
    mean(sp_df$nearest_host_distance[-which(sp_df$host_indicator)]),
    "for other species.\n")
cat("Median distance to nearest confirmed host for confirmed hosts is",
    median(sp_df$nearest_host_distance[which(sp_df$host_indicator)]),
    "and",
    median(sp_df$nearest_host_distance[-which(sp_df$host_indicator)]),
    "for other species.\n")

# Do a quick plot to see if there's an obvious difference in the distances for
# confirmed hosts:
if (PLOT){
  bin_size <- 10
  
  plot_ulim <- bin_size * ceiling(max(sp_df$nearest_host_distance) / bin_size)
  
  host_dhist <- hist(sp_df$nearest_host_distance[which(sp_df$host_indicator)],
                     breaks = seq(0, plot_ulim, by=bin_size),
                     plot = FALSE)
  nonhost_dhist <- hist(sp_df$nearest_host_distance[-which(sp_df$host_indicator)],
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

# Repeat for foraging strategies
if (PLOT){
  bin_size <- 10.
  
  plot_ulim <- 100.
  
  host_dhist <- hist(sp_df$ForStrat.watbelowsurf[which(sp_df$host_indicator)],
                     breaks = seq(0, plot_ulim, by=bin_size),
                     plot = FALSE)
  nonhost_dhist <- hist(sp_df$ForStrat.watbelowsurf[-which(sp_df$host_indicator)],
                        breaks = seq(0, plot_ulim, by=bin_size),
                        plot = FALSE)
  pal <- brewer.pal(6, "Dark2")
  plot(seq(0, plot_ulim-bin_size, by=bin_size),
       nonhost_dhist$density,
       type = "l",
       lwd = 2,
       col = pal[2],
       xlab = "% of foraging time spent\n >5cm below water surface",
       ylab = "Density")
  lines(seq(0, plot_ulim-bin_size, by=bin_size),
        host_dhist$density,
        type = "l",
        lwd = 2,
        col = pal[1])
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
if (PLOT){
  bin_size <- 10.
  
  plot_ulim <- 100.
  
  host_dhist <- hist(sp_df$ForStrat.wataroundsurf[which(sp_df$host_indicator)],
                     breaks = seq(0, plot_ulim, by=bin_size),
                     plot = FALSE)
  nonhost_dhist <- hist(sp_df$ForStrat.wataroundsurf[-which(sp_df$host_indicator)],
                        breaks = seq(0, plot_ulim, by=bin_size),
                        plot = FALSE)
  pal <- brewer.pal(6, "Dark2")
  plot(seq(0, plot_ulim-bin_size, by=bin_size),
       nonhost_dhist$density,
       type = "l",
       lwd = 2,
       col = pal[2],
       xlab = "% of foraging time spent\n around below water surface",
       ylab = "Density")
  lines(seq(0, plot_ulim-bin_size, by=bin_size),
        host_dhist$density,
        type = "l",
        lwd = 2,
        col = pal[1])
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

if (PLOT){
  cong_by_host <- table(sp_df$is_congregatory,
                        sp_df$host_indicator,
                        dnn = c("Congregative", "Host")) %>%
    addmargins()
  mig_by_host <- table(sp_df$is_migratory,
                       sp_df$host_indicator,
                       dnn = c("Migratory", "Host")) %>%
    addmargins()
  
  png("cong_table.png")
  p <- tableGrob(cong_by_host,
                 rows = c("Not congregative", "Congregative", "Total"),
                 cols = c("Not confirmed", "Confirmed host", "Total"))
  grid.arrange(p)
  dev.off()
  
  png("mig_table.png")
  p <- tableGrob(mig_by_host,
                 rows = c("Not migratory", "Migratory", "Total"),
                 cols = c("Not confirmed", "Confirmed host", "Total"))
  grid.arrange(p)
  dev.off()
}

migr_by_family <- table(sp_df$BLFamilyLatin,
                        sp_df$is_migratory,
                        dnn = c("Family", "Migratory")) %>%
  addmargins() %>%
  as.data.frame.matrix()
migr_by_family$percent_migr <- 100 * (migr_by_family$`TRUE`) / (migr_by_family$Sum)
write.csv(migr_by_family, "data/migr_by_family.csv", row.names=TRUE)

################################################################################
# Work out proximity to key taxa

# Start by matching species in phylogeny to EltonTraits so we can add higher
# taxonomic levels to phylo distance matrix
missing_species <- which(!(rownames(dmat_mean) %in% EltonTraits_df$Scientific))
cat("There are",
    length(missing_species),
    "species in the phylogeny which aren't in EltonTraits.")

# Given the relatively small number of missing species we just remove them from
# the matrix.
dmat_mean <- dmat_mean[-missing_species, -missing_species]

orders_in_dmat <- lapply(rownames(dmat_mean),
                           FUN = function(x){
                             EltonTraits_df$IOCOrder[which(EltonTraits_df$Scientific==x)[1]]})
families_in_dmat <- lapply(rownames(dmat_mean),
                         FUN = function(x){
                           EltonTraits_df$BLFamilyLatin[which(EltonTraits_df$Scientific==x)[1]]})

# For genus, just take the first word of each species name.
genera_in_dmat <- lapply(rownames(dmat_mean),
                         FUN = function(x){
                           word(x, 1)})

# Corrections to genus based on AVONET:
genera_in_dmat[which(genera_in_dmat=="limicola")] <- "calidris"
genera_in_dmat[which(genera_in_dmat=="philomachus")] <- "calidris"
genera_in_dmat[which(genera_in_dmat=="eurynorhynchus")] <- "calidris"
genera_in_dmat[which(genera_in_dmat=="tryngites")] <- "calidris"
genera_in_dmat[which(genera_in_dmat=="aphriza")] <- "calidris"

genera_in_dmat[which(genera_in_dmat=="chen")] <- "anser"

# As an example, do distance to Laridae
laridae_cols <- dmat_mean[, which(families_in_dmat=="Laridae")]
laridae_distance <- apply(laridae_cols, 1, min)

# We then need to add this to sp_df, as we did with distance to host

################################################################################
# For interpretative purposes it is useful to have idiomatic versions of the
# Elton traits variable names
name_pairs <- data.frame(0)
name_pairs$Diet.Inv <- "Invertebrates"
name_pairs$Diet.Vend <- "Mammals and birds"
name_pairs$Diet.Vect <- "Amphibians and reptiles"
name_pairs$Diet.Vfish <- "Fish"
name_pairs$Diet.Vunk <- "Unknown vertebrates"
name_pairs$Diet.Scav <- "Scavenging"
name_pairs$Diet.Nect <- "Nectar"
name_pairs$Diet.Fruit <- "Fruit"
name_pairs$Diet.Seed <- "Seeds"
name_pairs$Diet.PlantO <- "Other plant matter"
name_pairs$ForStrat.watbelowsurf <- ">5cm below water surface"
name_pairs$ForStrat.wataroundsurf <- "Around water surface"
name_pairs$ForStrat.aerial <- "Aerial"
name_pairs$ForStrat.canopy <- "Canopy"
name_pairs$ForStrat.ground <- "Ground"
name_pairs$ForStrat.understory <- "Understory"
name_pairs$ForStrat.midhigh <- ">2m above ground, below canopy"
name_pairs$PelagicSpecialist <-  "Pelagic specialist"
name_pairs$Nocturnal <- "Nocturnal"
name_pairs$BodyMass.Value <- "Body mass"
name_pairs$nearest_host_distance <- "Phylo. distance to confirmed host"
name_pairs$anatidae_distance <- "Phylo. distance to Anatidae"
name_pairs$larus_distance <- "Phylo. distance to Larus"
name_pairs$is_congregatory <- "Congregatory"
name_pairs$is_migratory <- "Migratory"

expand_names <- function(name_list){
  expanded_names <- c()
  for (n in name_list){
    if (n %in% colnames(name_pairs)){
      expanded_names <- append(expanded_names, name_pairs[n])
    }
    else{
      expanded_names <- append(expanded_names, n)
    }
  }
  return(expanded_names)
}

################################################################################
# Now train BART

# We can half dataset size while keeping ~80% of the data by removing passerines
sp_df <- sp_df[-which(sp_df$IOCOrder=="Passeriformes"), ]

# Remove non-numerical data as well as stuff we don't want to fit
x <- sp_df[, -which(names(sp_df) %in% c("species_code",
                "scientific_name",
                "common_name",
                "breeding_quality",
                "nonbreeding_quality",
                "postbreeding_migration_quality",
                "prebreeding_migration_quality",
                "resident_quality",
                "Avibase_ID",
                "synonyms",
               "host_indicator",
               "PassNonPass",
               "IOCOrder",
               "BLFamilyLatin",
               "Scientific",
               "Diet.5Cat",
               "Diet.Certainty",
               "PelagicSpecialist",
               "ForStrat.SpecLevel",
               "Nocturnal",
               "BodyMass.Value",
               "BodyMass.SpecLevel",
               "scientificName",
               "nearest_host_distance",
               max_diet_comp,
               max_forstrat_comp))]
y <- as.numeric(sp_df$"host_indicator")

n_pts <- nrow(sp_df)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(x), n_training)
xtrain <- x[training, ]
ytrain <- y[training]
xtest <- x[-training, ]
ytest <- y[-training]

basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
varselect <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 full = TRUE,
                 quiet = TRUE)

varimp_raw <- varselect$varcount %>% 
  as.data.frame %>%
  mutate(./rowSums(.)) %>%
  mutate(./max(.))

varimp_summ <-  bind_rows(varimp_raw %>% summarise(across(everything(), mean)),
                               varimp_raw %>% summarise(across(everything(), sd))) %>%
  t() %>%
  data.frame(var = row.names(.)) %>%
  relocate(var) %>%
  rename("mean" = "X1", "sd" = "X2")

ordered_var <- varimp_summ$var[order(varimp_summ$mean, decreasing = TRUE)]
varimp_summ <- varimp_summ %>% mutate(var = fct_relevel(var, ordered_var),
                    upper = mean + sd, lower = mean - sd)

fig_varimp <- ggplot(varimp_summ, aes(x = var, y = mean, ymin = lower, ymax = upper)) + 
  geom_errorbar(width=0, position=position_dodge(0.6)) + 
  geom_point(position=position_dodge(0.6)) +
  geom_vline(xintercept=seq(1.5, nrow(varimp_summ)-0.5, 1), 
             lwd=.5, colour="grey75") + 
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_text(angle = 35, hjust = 1), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key.size=unit(0.8, "lines"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.background = element_rect(colour = "grey50"),
        axis.title.x=element_blank(),
        plot.margin = margin(10, 3, 3, 80),
        panel.grid.major.x = element_blank()) +
  ylab("Relative variable importance")

vars <- as.character(varimp_summ$var)
pd_summ <- vector("list", length(vars))

for(i in 1:length(vars)){
  
  # Adapted from embarcadero::partial
  fullvarname <- varselect$fit$data@x %>% as.data.frame %>% dplyr::select(matches(vars[i])) %>% names
  if (length(fullvarname) == 1) {
    raw <- varselect$fit$data@x[, fullvarname]
    lev <- list(seq(min(raw), max(raw), ((max(raw) - min(raw))/15)))
    pd <- pdbart(varselect, xind = fullvarname, levs = lev, pl = FALSE)
    
    pd_summ[[i]] <- data.frame(var = vars[i],
                                    x = unlist(lev), 
                                    y = pd$fd[[1]] %>% apply(., 2, median) %>% pnorm,
                                    upper = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.025) %>% pnorm,
                                    lower = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.975) %>% pnorm
                             )
  }
}

fig_pd_chosen <- pd_summ %>%
  bind_rows %>%
  ggplot(aes(x = x, y = y, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = 0.08, colour = NA) +
  geom_line(lwd = 0.8, alpha = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  ylab("Probability") +
  theme_bw() +
  theme(axis.title.x=element_blank()) +
  facet_wrap(~ var, scales = "free_x")
