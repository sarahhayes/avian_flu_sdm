# In this script we analyse species-level patterns regarding avian influenza
# host status.

library(BART)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(RColorBrewer)
library(readxl)

# First use the AVONET data to create a mapping between names and alphanumeric
# identifiers
AVONET_df <- read_excel("data/AVONETSupplementarydataset1.xlsx",
                        sheet = "AVONET_Raw_Data")
AVONET_df <- AVONET_df[,
                        c("Avibase.ID",
                          "Species1_BirdLife",
                          "Species2_eBird",
                          "Species3_BirdTree")]
AVONET_df <- distinct(AVONET_df)
AVONET_df[, 2:4] <- sapply(AVONET_df[, 2:4], tolower)

# Assign ID based on consensus
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
    # Construct consensus distribution of candidate names
    unique_candidates <- unique(ID_candidates)
    if (length(unique_candidates)==1){
      # cat("Number of unique candidates is", length(unique_candidates), "\n")
      # cat(unique_candidates, "\n")
      return(unique_candidates)
    }
    candidate_dist <- c()
    for (i in 1:length(unique_candidates)){
      # cat("i=", i, "\n")
      candidate_dist <- append(candidate_dist,
                               sum(ID_weightings[
                                 which(ID_candidates==unique_candidates[i])]))
    }
    # cat(unique_candidates, "\n", ID_candidates, "\n", candidate_dist, "\n", ID_weightings, "\n")
    if (length(which(candidate_dist==max(candidate_dist)))==1){
      return(unique_candidates[which(candidate_dist==max(candidate_dist))])
    }
    else{
      # cat("Multiple consensus options for species", binomial_name,"\n")
      return(sample(unique_candidates[which(candidate_dist==max(candidate_dist))], 1))
    }
  }
}

get_bird_ids <- function(name_vect, name_sources){
  if (missing(name_sources)){
    name_sources <-c("BirdLife",
                     "eBird",
                     "BirdTree")
  }
  id_vect <- vector(length = length(name_vect))
  for (n in 1:length(name_vect)){
    # cat("n=",n,"\n")
    id_vect[n] <- get_single_bird_id(name_vect[n], name_sources)
  }
  return(id_vect)
}

################################################################################
# Now introduce CLOVER database and identify all avian species known to be hosts
# of influenza A.

setwd("../clover")

# Load data
CLOVER_df <- read.csv("clover/clover_1.0_allpathogens/CLOVER_1.0_Viruses_AssociationsFlatFile.csv")
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
# Now introduce eBird data

setwd('../avian_flu_sdm')

# Read in codes
eBird_names <- read.csv("eBird/eBird_species_europe_copy.csv", header = F)

# Remove number and URL rows:
remove_rows <- sort(c(seq(from = 1, to = nrow(eBird_names), by = 4),
                      seq(from = 4, to = nrow(eBird_names), by = 4)))
eBird_names <- eBird_names[-remove_rows, ]

# The data now has a line with the common name followed by a line with the the
# common name and the binomial name for each species. The following loop removes
# the common name from the second row for each species.
for (i in seq(2, length(eBird_names), 2)){
  eBird_names[i] <- substr(eBird_names[i],
                           nchar(eBird_names[i-1]) + 2,
                           1000)
}
eBird_names <- eBird_names[seq(2, length(eBird_names), 2)]
eBird_names <- sapply(eBird_names, tolower)

# Now filter CLOVER data for names in eBird
eBird_names_in_clover <- species_list[which(species_list %in% eBird_names)]

# Check it's symmetric in terms of which way round we do it - it should be, but
# if it wasn't that might tell us something interesting!
clover_names_in_eBird <- eBird_names[which(eBird_names %in% species_list)]
setequal(eBird_names_in_clover, clover_names_in_eBird)

# Not all the names from CLOVER appear in the eBird data - this could be because
# these species aren't found in Europe, but it could also be because the two
# datasets are using inconsistent binomial names. There should only be 50 names
# from CLOVER missing in the eBird data, so we'll inspect visually
clover_names_not_in_eBird <- species_list[-which(species_list %in% eBird_names)]
for (c in clover_names_in_eBird){
  cat(c, "\n")
}

# We now want the ID's for the species in eBird. 
eBird_IDs <- get_bird_ids(eBird_names, c("eBird"))

# See if any species were not successfully ID'd
cat(length(which(eBird_IDs=="-100")), "species were not ID'd.")
non_IDd_species <- eBird_names[which(eBird_IDs=="-100")]
cat("Unidentified species are:",
    non_IDd_species,
    sep = "\n")
for (spe in non_IDd_species){
  if (spe %in% AVONET_df$Species1_BirdLife){
    cat(spe,
        "is in BirdLife and cross-referenced with",
        AVONET_df$Species2_eBird[which(AVONET_df$Species1_BirdLife==spe)],
        "in eBird.\n")
  }
  else{
    cat(spe, "does not appear in BirdLife.\n")
  }
  if (spe %in% AVONET_df$Species2_eBird){
    cat(spe,
        "is in eBird.\n")
  }
  else{
    cat(spe, "does not appear in eBird\n")
  }
  if (spe %in% AVONET_df$Species3_BirdTree){
    cat(spe,
        "is in BirdTree and cross-referenced with",
        AVONET_df$Species2_eBird[which(AVONET_df$Species3_BirdTree==spe)],
        "in eBird.\n")
  }
  else{
    cat(spe, "does not appear in BirdTree.\n")
  }
}

# Should find a matchup between Cuculus Optatus and Cuculus Solitarius, but
# if you look this up online it appears to be spurious.
# We match up Cuculus Optatus with Cuculus Saturatus because C. Optatus was
# formerly classed as a subspecies of C. Saturatus.
eBird_names[which(eBird_names=="cuculus optatus")] <- "cuculus saturatus"

# Porphyrio Poliocephalus doesn't have an obvious match in AVONET and is only
# found in small numbers at the very edge of the EPSG:3035 range so we just
# remove it from the list
eBird_names <- eBird_names[which(eBird_names!="porphyrio poliocephalus")]

# This is a simple alternative name:
eBird_names[which(eBird_names=="daptrius chimachima")] <- "milvago chimachima"

# Now try doing ID matchup again:
eBird_IDs <- get_bird_ids(eBird_names, c("eBird"))

# Check we got everything:
cat(length(which(eBird_IDs=="-100")), "species were not ID'd.")

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

# Identify species that appear in CLOVER
host_IDs_in_elton <- unique(EltonTraits_IDs[which(EltonTraits_IDs %in% host_species_IDs)])

# Check if number of species left in trait database matches number from CLOVER:
length(host_IDs_in_elton)==no_host_species

# Create column indicating whether species is a host
host_indicator = c(EltonTraits_df$Scientific %in% species_list)

# Append it to trait database
matched_data <- data.frame(EltonTraits_df, host_indicator)

# Now filter for species in Europe according to eBird:
matched_data <- matched_data[which(EltonTraits_IDs %in% eBird_IDs), ]

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

# For interpretative purposes it is useful to have idiomatic versions of the
# elton traits variable names
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

train <- sample(1:nrow(matched_data), nrow(matched_data)/2)

# Remove non-numerical data as well as stuff we don't want to fit
x <- matched_data[,!(names(matched_data) %in% c("PassNonPass",
                                                "IOCOrder",
                                                "BLFamilyLatin",
                                                "Scientific",
                                                "Diet.5Cat",
                                                "Diet.Certainty",
                                                "ForStrat.SpecLevel",
                                                "BodyMass.SpecLevel",
                                                "host_indicator",
                                                max_diet_comp,
                                                max_forstrat_comp))]
y <- matched_data[, "host_indicator"]
xtrain <- x[train, ]
ytrain <- y[train]
xtest <- x[-train, ]
ytest <- y[-train]

bartfit <- lbart(xtrain,
                 ytrain,
                 x.test = xtest,
                 ntree = 1000)

# Use a majority judgement to decide whether to accept
yhat.train <- plogis(bartfit$yhat.train - bartfit$binaryOffset)
yhat.maj <- colSums(yhat.train) >= nrow(yhat.train)/2

# Calculate error rates
false_neg_rate <- length(which((!yhat.maj)&ytrain))/length(which(ytrain))
false_pos_rate <- length(which(yhat.maj&!ytrain))/length(which(!ytrain))
misclass_rate <- (length(which((!yhat.maj)&ytrain)) +
  length(which(yhat.maj&!ytrain))) / length(ytrain)
cat("Training false negative rate is",
    false_neg_rate,
    ".\n Training false positive rate is",
    false_pos_rate,
    ".\n Training misclassification rate is",
    misclass_rate,
    ".\n")

# And for test data:
yhat.test <- plogis(bartfit$yhat.test - bartfit$binaryOffset)
yhat.maj <- colSums(yhat.test) >= nrow(yhat.test)/2

# Calculate error rates
false_neg_rate <- length(which((!yhat.maj)&ytest))/length(which(ytest))
false_pos_rate <- length(which(yhat.maj&!ytest))/length(which(!ytest))
misclass_rate <- (length(which((!yhat.maj)&ytest)) +
                    length(which(yhat.maj&!ytest))) / length(ytest)
cat("Test false negative rate is",
    false_neg_rate,
    ".\n Test false positive rate is",
    false_pos_rate,
    ".\n Test misclassification rate is",
    misclass_rate,
    ".\n")

# How many times to variables appear?
ord <- order(bartfit$varcount.mean, decreasing = T)
varcount_df <- data.frame((1/1000) * bartfit$varcount.mean)

par(mar= c(14, 5, 1, 8))
varimp_ord <- data.frame(varcount_df[ord, ])
rownames(varimp_ord) <- rownames(varcount_df)[ord]
short_varnames <- expand_names(rownames(varimp_ord))
diet_vars <- which(rownames(varimp_ord) %like% "Diet.")
forstrat_vars <- which(rownames(varimp_ord) %like% "ForStrat.")
other_vars <- 1:18
other_vars <- other_vars[-which(other_vars %in% diet_vars)]
other_vars <- other_vars[-which(other_vars %in% forstrat_vars)]
bar_cols <- vector(length = 18)
pal <- brewer.pal(3, "Dark2")
bar_cols[diet_vars] <- pal[1]
bar_cols[forstrat_vars] <- pal[2]
bar_cols[other_vars] <- pal[3]
varimp_bars <- barplot(varimp_ord[, 1],
                       border = F,
                       las = 2,
                       names.arg = short_varnames,
                       col = bar_cols,
                       ylab = "Mean appearances\n per tree")
legend("topright",
       inset=c(-0.225,0),
       legend = c("Dietary (% of calories)", "Foraging (% strategy)", "Other"),
       col = pal,
       bty = "n",
       pch = 20,
       pt.cex = 2,
       cex = 0.8,
       xpd = TRUE)

# Do a more thorough variable importance
tree_nos <- c(10,
              20,
              50,
              100,
              150,
              200)

for (i in 1:(length(tree_nos)-1)){
  this_bartfit <- lbart(xtrain,
                        ytrain,
                        x.test = xtest,
                        ntree = tree_nos[i])
  if (i==1){
    varimp_df <- data.frame(this_bartfit$varcount.mean/sum(this_bartfit$varcount.mean))
  }
  else{
    varimp_df <- data.frame(varimp_df, this_bartfit$varcount.mean/sum(this_bartfit$varcount.mean))
  }
}
varimp_df <- data.frame(varimp_df, varcount_df)

pal <- brewer.pal(6, "Dark2")

par(mar= c(10, 8, 1, 5))
matplot(varimp_df[ord, ],
        type = "l",
        xlab = "",
        ylab = "Normalised mean\n appearances per tree",
        xaxt = "n",
        col=pal,
        lty = 1)
for (i in 1:length(tree_nos)){
  points(x = 1:length(ord),
         y = varimp_df[ord, i],
         col = pal[i],
         pch = 16)
}
axis(1, at = 1:18, labels = names(xtrain)[ord], las = 2)
legend("topright",
       inset = c(-0.0,0),
       legend = tree_nos,
       title = "Number of trees",
       col = pal,
       lty = 1,
       xpd = TRUE)

################################################################################
# Now try some PCA.

matched_numeric <- data.frame(x, as.numeric(y))
names(matched_numeric)[21] <- 'host_indicator'

pc_all <- prcomp(matched_numeric,
             center = TRUE,
             scale. = TRUE)
summary(pc_all)

layout(matrix(1:100, ncol = 10), respect = TRUE)

for (i in 1:5){
  for (j in 1:5){
    if (i > j){
      autoplot(pc_all,
               data = matched_numeric,
               colour = 'host_indicator',
               x=i,
               y=j)
    }
  }
}

host_species_traits <- data.frame(x[which(y==TRUE), ])

pc_host <- prcomp(host_species_traits,
             center = TRUE,
             scale. = TRUE)
summary(pc_host)
autoplot(pc_host, data = host_species_traits)
