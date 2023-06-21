# Load in CLOVER database and identify all avian species known to be hosts of
# influenza A.

library(BART)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(gridExtra)

setwd("Github/clover")

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

################################################################################
# Now introduce trait data

setwd('..')
setwd('avian_flu_sdm')
# Load in from file
elton_trait_df <- 
  read.table("data/variables/bird_data/elton_traits/BirdFuncDat.txt",
             sep = '\t',
             quote="\"",
             header = TRUE)
# elton trait data has two rows of NA's at the bottom which we need to remove:
elton_trait_df <- elton_trait_df[1:9993, ]

# Set scientific names to lower case to match with CLOVER data
elton_trait_df$Scientific <- tolower(elton_trait_df$Scientific)
# Identify species that appear in CLOVER
hosts_in_elton <- elton_trait_df[which(elton_trait_df$Scientific %in% species_list), 'Scientific']

# Check if number of species left in trait database matches number from CLOVER:
length(hosts_in_elton)==no_host_species

# The numbers don't match because elton_traits uses outdated binomial names
# and doesn't account for some recently classified species. We make a look-up
# table to match the names used by the two datasets.

# Important caveat: I generated the list of alternative names to match with the
# ones from CLOVER using ChatGPT. I've cleaned it up enough to make sure all the
# names are present in elton_trait but I'm not totally certain that all the
# matches are correct. For instance, ChatGPT tried to replace the binomial name
# for the Northern waterthrush with that of the Kentucky warbler.

CLOVER_names <- c(
  "anas carolinensis",
  "anser caerulescens",
  "ardea alba",
  "ardenna pacifica",
  "bubo scandiacus",
  "chroicocephalus brunnicephalus",
  "chroicocephalus novaehollandiae",
  "chroicocephalus ridibundus",
  "coloeus monedula",
  "geothlypis tolmiei",
  "ichthyaetus ichthyaetus",
  "ichthyaetus melanocephalus",
  "larus mongolicus",
  "leucophaeus atricilla",
  "mareca americana",
  "mareca penelope",
  "mareca strepera",
  "onychoprion fuscatus",
  "parkesia noveboracensis",
  "piaya minuta",
  "sibirionetta formosa",
  "spatula clypeata",
  "spatula discors"
)

elton_traits_names <- c(
  "anas crecca",
  "chen caerulescens",
  "casmerodius albus",
  "puffinus pacificus",
  "bubo scandiaca",
  "larus brunnicephalus",
  "larus novaehollandiae",
  "larus ridibundus",
  "corvus monedula",
  "oporornis tolmiei",
  "larus ichthyaetus",
  "larus melanocephalus",
  "larus cachinnans",
  "larus atricilla",
  "anas americana",
  "anas penelope",
  "anas strepera",
  "sterna fuscata",
  "seiurus noveboracensis",
  "crotophaga sulcirostris",
  "anas formosa",
  "anas clypeata",
  "anas discors"
)

name_lookup <- data.frame(CLOVER = CLOVER_names,
                          elton_traits = elton_traits_names)

# Replace names in species_list:
for (name_idx in 1:nrow(name_lookup)){
  species_list[which(species_list==name_lookup$CLOVER[name_idx])] <- 
    name_lookup$elton_traits[name_idx]
}
# elton_traits treats the Green-winged and Eurasian teals as the same species,
# meaning the species_list now has a repeated element. We reduce to its unique
# elements, noting that this won't affect our analysis since based on the
# elton_traits data we'd just have two identical species.
species_list <- unique(species_list)
no_host_species <- length(species_list)

# Identify species that appear in CLOVER
hosts_in_elton <- elton_trait_df[which(elton_trait_df$Scientific %in% species_list), 'Scientific']

# Check if number of species left in trait database matches number from CLOVER:
length(hosts_in_elton)==no_host_species
# We now have matching numbers!

# Create column indicating whether species is a host
host_indicator = c(elton_trait_df$Scientific %in% species_list)

# Append it to trait database
matched_data <- data.frame(elton_trait_df, host_indicator)

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

# Should find around 80% of species have a good level of certainty
proportion_certain <- nrow(certain_data) / nrow(matched_data)

# Should also find all but one confirmed host species stay in the data when
# filtering for certainty:
sum(certain_data$host_indicator)

# Identify this species and see how many samples there are from it
hosts_with_certainty <- certain_data[which(certain_data$host_indicator),
                                     "Scientific"]
uncertain_host <- species_list[-which(species_list %in% hosts_with_certainty)]
species_counts[which(species_counts$Host==uncertain_host), "n"]
# Should find there are 8 samples from this species - so not totally negligible
# but not a massive chunk of the data either!

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

bartfit <- lbart(xtrain, ytrain, x.test = xtest)

# Use a majority judgement to decide whether to accept
yhat.train <- plogis(bartfit$yhat.train - bartfit$binaryOffset)
yhat.maj <- colSums(yhat.train) >= nrow(yhat.train)/2

# Calculate error rates
false_neg_rate <- length(which((!yhat.maj)&ytrain))/length(which(ytrain))
false_pos_rate <- length(which(yhat.maj&!ytrain))/length(which(!ytrain))
misclass_rate <- (length(which((!yhat.maj)&ytrain)) +
  length(which(yhat.maj&!ytrain))) / length(ytrain)

# And for test data:
yhat.test <- plogis(bartfit$yhat.test - bartfit$binaryOffset)
yhat.maj <- colSums(yhat.test) >= nrow(yhat.test)/2

# Calculate error rates
false_neg_rate <- length(which((!yhat.maj)&ytest))/length(which(ytest))
false_pos_rate <- length(which(yhat.maj&!ytest))/length(which(!ytest))
misclass_rate <- (length(which((!yhat.maj)&ytest)) +
                    length(which(yhat.maj&!ytest))) / length(ytest)

# How manyu times to variables appear?
ord <- order(bartfit$varcount.mean, decreasing = T)
varcount_df <- data.frame(bartfit$varcount.mean[ord])
colnames(varcount_df) <- c("MeanCountsPerTree")
plot(varcount_df$MeanCountsPerTree,
     xlab = "",
     ylab = "Mean appearances per tree",
     xaxt = "n")
axis(1, at = 1:18, labels = names(xtrain), las = 2)

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
  this_ord <- order(this_bartfit$varcount.mean, decreasing = T)
  if (i==1){
    varimp_df <- data.frame(this_bartfit$varcount.mean[this_ord])
  }
  else{
    varimp_df <- data.frame(varimp_df, this_bartfit$varcount.mean[this_ord])
  }
}
varimp_df <- data.frame(varimp_df, varcount_df)

# Putting this in as a bit of a patch because I accidentally made the data frame
# the wrong way round!
varimp_df_t <- transpose(varimp_df)
names(varimp_df_t) <- rownames(varimp_df)
rownames(varimp_df_t) <- names(varimp_df)

matplot(varimp_df,
        type = "l",
        xlab = "",
        ylab = "Mean appearances per tree",
        xaxt = "n")
axis(1, at = 1:18, labels = names(xtrain), las = 2)
legend("right", legend = tree_nos)
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