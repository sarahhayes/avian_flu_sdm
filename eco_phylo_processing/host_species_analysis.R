# Load in CLOVER database and identify all avian species known to be hosts of
# influenza A.

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
no_species <- length(species_list)

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

# Load in from file
elton_trait_df <- 
  read.table("data/variables/bird_data/elton_traits/BirdFuncDat.txt",
             sep = '\t',
             quote="\"",
             header = TRUE)
# Set scientific names to lower case to match with CLOVER data
elton_trait_df$Scientific <- tolower(elton_trait_df$Scientific)
# Filter by birds that appear in CLOVER
elton_trait_df <- elton_trait_df[which(elton_trait_df$Scientific %in% species_list), ]

# Check if number of species left in trait database matches number from CLOVER:
nrow(elton_trait_df)==no_species

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
no_species <- length(species_list)

# Now reload elton_traits data but filter for changed species names:
elton_trait_df <- 
  read.table("data/variables/bird_data/elton_traits/BirdFuncDat.txt",
             sep = '\t',
             quote="\"",
             header = TRUE)
# Set scientific names to lower case to match with CLOVER data
elton_trait_df$Scientific <- tolower(elton_trait_df$Scientific)
# Filter by birds that appear in CLOVER
elton_trait_df <- elton_trait_df[which(elton_trait_df$Scientific %in% species_list), ]

# Check if number of species left in trait database matches number from CLOVER:
nrow(elton_trait_df)==no_species
# We now have matching numbers!