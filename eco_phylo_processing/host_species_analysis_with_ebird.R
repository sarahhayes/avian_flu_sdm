# Load in CLOVER database and identify all avian species known to be hosts of
# influenza A.

library(BART)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(readxl)

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
# Now introduce eBird data

setwd('..')
setwd('avian_flu_sdm')

# Read in codes
ebird_names <- read.csv("ebird/ebird_species_europe_copy.csv", header = F)

# Remove number and URL rows:
remove_rows <- sort(c(seq(from = 1, to = nrow(ebird_names), by = 4),
                      seq(from = 4, to = nrow(ebird_names), by = 4)))
ebird_names <- ebird_names[-remove_rows, ]

# The data now has a line with the common name followed by a line with the the
# common name and the binomial name for each species. The following loop removes
# the common name from the second row for each species.
for (i in seq(2, length(ebird_names), 2)){
  ebird_names[i] <- substr(ebird_names[i],
                           nchar(ebird_names[i-1]) + 2,
                           1000)
}
ebird_names <- ebird_names[seq(2, length(ebird_names), 2)]
ebird_names <- sapply(ebird_names, tolower)

# Now filter CLOVER data for names in ebird
ebird_names_in_clover <- species_list[which(species_list %in% ebird_names)]

# Check it's symmetric in terms of which way round we do it - it should be, but
# if it wasn't that might tell us something interesting!
clover_names_in_ebird <- ebird_names[which(ebird_names %in% species_list)]
setequal(ebird_names_in_clover, clover_names_in_ebird)

# Not all the names from CLOVER appear in the ebird data - this could be because
# these species aren't found in Europe, but it could also be because the two
# datasets are using inconsistent binomial names. There should only be 50 names
# from CLOVER missing in the ebird data, so we'll inspect visually
clover_names_not_in_ebird <- species_list[-which(species_list %in% ebird_names)]
for (c in clover_names_in_ebird){
  cat(c, "\n")
}

################################################################################
# Now introduce trait data

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

# Try rechecking CLOVER-eBird matchups to see if we do anything better with new
# binomial names
for (name_idx in 1:nrow(name_lookup)){
  ebird_names[which(ebird_names==name_lookup$CLOVER[name_idx])] <- 
    name_lookup$elton_traits[name_idx]
}
# Now filter CLOVER data for names in ebird
ebird_names_in_clover <- species_list[which(species_list %in% ebird_names)]

# Check it's symmetric in terms of which way round we do it - it should be, but
# if it wasn't that might tell us something interesting!
clover_names_in_ebird <- ebird_names[which(ebird_names %in% species_list)]
setequal(ebird_names_in_clover, clover_names_in_ebird)

# Identify species that appear in CLOVER
hosts_in_elton <- elton_trait_df[which(elton_trait_df$Scientific %in% species_list), 'Scientific']

# Check if number of species left in trait database matches number from CLOVER:
length(hosts_in_elton)==no_host_species
# We now have matching numbers!

# Create column indicating whether species is a host
host_indicator = c(elton_trait_df$Scientific %in% species_list)

# Append it to trait database
matched_data <- data.frame(elton_trait_df, host_indicator)

# Now filter for species in Europe according to eBird:
matched_data <- matched_data[which(matched_data$Scientific %in% ebird_names), ]

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
varcount_df <- data.frame((1/1000) * bartfit$varcount.mean)
par(mar= c(10, 8, 1, 5))
matplot(varcount_df[ord, ],
        type = "l",
        xlab = "",
        ylab = "Normalised mean\n appearances per tree",
        xaxt = "n",
        lty = 1)
points(x = 1:length(ord),
         y = varcount_df[ord, ],
         pch = 16)
axis(1, at = 1:18, labels = rownames(varcount_df)[ord], las = 2)

par(mar= c(14, 5, 1, 8))
varimp_ord <- data.frame(varcount_df[ord, ])
rownames(varimp_ord) <- rownames(varcount_df)[ord]
diet_vars <- which(rownames(varimp_ord) %like% "Diet.")
forstrat_vars <- which(rownames(varimp_ord) %like% "ForStrat.")
other_vars <- 1:18
other_vars <- other_vars[-which(other_vars %in% diet_vars)]
other_vars <- other_vars[-which(other_vars %in% forstrat_vars)]
short_varnames <- c("Around water surface",
                    "Seeds",
                    ">5cm below water surface",
                    "Nocturnal",
                    ">2m above ground, below canopy",
                    "Pelagic specialist",
                    "Other plant matter",
                    "Mammals and birds",
                    "Fish",
                    "Canopy",
                    "Understory",
                    "Fruit",
                    "Nectar",
                    "Amphibians and reptiles",
                    "Scavenging",
                    "Body mass",
                    "Aerial",
                    "Unknown vertebrates"
                    )
bar_cols <- vector(length = 18)
bar_cols[diet_vars] <- rgb(0.3,0.1,0.4,0.6)
bar_cols[forstrat_vars] <- rgb(0.3,0.5,0.4,0.6)
bar_cols[other_vars] <- rgb(0.3,0.9,0.4,0.6)
varimp_bars <- barplot(varimp_ord[, 1],
                       border = F,
                       las = 2,
                       names.arg = short_varnames,
                       col = bar_cols,
                       ylab = "Mean appearances\n per tree")
legend("topright",
       inset=c(-0.5,0),
       legend = c("Dietary (% of calories)", "Foraging (% strategy)", "Other"),
       col = c(rgb(0.3,0.1,0.4,0.6), rgb(0.3,0.5,0.4,0.6), rgb(0.3,0.9,0.4,0.6)),
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

par(mar= c(10, 8, 1, 5))
matplot(varimp_df[ord, ],
        type = "l",
        xlab = "",
        ylab = "Normalised mean\n appearances per tree",
        xaxt = "n",
        col=1:6,
        lty = 1)
for (i in 1:length(tree_nos)){
  points(x = 1:length(ord),
         y = varimp_df[ord, i],
         col = i,
         pch = 16)
}
axis(1, at = 1:18, labels = names(xtrain)[ord], las = 2)
legend("topright",
       inset=c(-0.3,0),
       legend = tree_nos,
       col = 1:6,
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
