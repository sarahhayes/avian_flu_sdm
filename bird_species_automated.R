# trying to automate the process of identifying the taxons

rm(list = ls())

library(tidyverse)
#install.packages("taxize")
library(taxize)

## bv-brc data

pos_data <- read.csv("data/flu_data/BVBRC_wild_captivewild_positive.csv")
columns_needed <- c("Collection.Date", "Collection.Year", "Collection.Country", 
                    "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                    "Host.Species", "Host.Group", "Host.Natural.State", "Host.Health")

pos_data <- select(pos_data, all_of(columns_needed))

## if want to just select the wild
pos_data <- pos_data %>%
  filter(Host.Natural.State == "Wild")

bvbrc_species <-as.data.frame(table(pos_data$Host.Species))
colnames(bvbrc_species)<- c("species", "count")

bvbrc_species$species <- as.character(bvbrc_species$species)

#example <- bvbrc_species[1:10, "species"]
#exampled_df <- tax_name(example, get = 'family', db = 'ncbi')

bvbrc_species[which(bvbrc_species$species == "Cygnus sp."),"species"] <- "Cygnus"

bvbrc_species[which(bvbrc_species$species == "Gallinago spp"),"species"] <- "Gallinago"

bvbrc_species[which(bvbrc_species$species == "Bucephala sp."),"species"] <- "Bucephala"

bvbrc_species[which(bvbrc_species$species == "Calidris sp."),"species"] <- "Calidris"

bvbrc_species[which(bvbrc_species$species == "Larus sp."),"species"] <- "Larus"

bvbrc_species[which(bvbrc_species$species == "Larus sp"),"species"] <- "Larus"




# reports bad request intermittently. 
# Connection issue? Can just restart from where it had got to. 

# first do for order
for (i in 1:nrow(bvbrc_species)) {
  bvbrc_species[i,"order"] <- tax_name(bvbrc_species[i,"species"], 
                                        get = "order",
                                        db = "ncbi")$order 
  print(i)
}

# Manually add that charadrii are charadriiformes


bvbrc_species[which(bvbrc_species$species == "Charadrii"),"order"] <- "Charadriiformes"

bvbrc_orders <- bvbrc_species %>%
  group_by(order) %>%
  summarise(bvbrc_order_count = sum(count))


## now do for families
for (i in 1:nrow(bvbrc_species)) {
  bvbrc_species[i,"family"] <- tax_name(bvbrc_species[i,"species"], 
                                         get = "family",
                                         db = "ncbi")$family 
  print(i)
}

bvbrc_families <- bvbrc_species %>%
  group_by(family) %>%
  summarise(bvbrc_family_count = sum(count))


# now genus
for (j in 1:nrow(bvbrc_species)) {
  bvbrc_species[j,"genus"] <- tax_name(bvbrc_species[j,"species"], 
                                        get = "genus",
                                        db = "ncbi")$genus
  print(j)
}

bvbrc_genus <- bvbrc_species %>%
  group_by(genus) %>%
  summarise(bvbrc_genus_count = sum(count))


## repeat for the negative data as useful to have to get % positive

neg_data <- read.csv("data/flu_data/BVBRC_wild_captivewild_negative.csv")

neg_data <- select(neg_data, all_of(columns_needed))
# use below if just want the wild and not captive wild
neg_data <- neg_data %>%
  filter(Host.Natural.State == "Wild")


bvbrc_species_neg <-as.data.frame(table(neg_data$Host.Species))
colnames(bvbrc_species_neg)<- c("species", "count_neg")

bvbrc_species_neg$species <- as.character(bvbrc_species_neg$species)

# Really want this by species so just use the info from those that have the full species in
# both positive and negative. 

bvbrc_species_pos <- bvbrc_species
colnames(bvbrc_species_pos) <- c("species", "count_pos")

bvbrc_pos_neg <- merge(bvbrc_species_neg, bvbrc_species_pos) # default for merge is to just include 
# columns that have entries for both. Probably most suitable here

bvbrc_pos_neg[,"total_tested"] <- bvbrc_pos_neg[,"count_pos"] + bvbrc_pos_neg[,"count_neg"]
bvbrc_pos_neg[,"percent_pos"] <- (bvbrc_pos_neg[,"count_pos"]/bvbrc_pos_neg[,"total_tested"])*100

# quite a few with high % because so few were tested. 
# Filter to those where at least 10 tested

bvbrc_pos_neg_fifty <- bvbrc_pos_neg %>%
  filter(total_tested >= 50)


###################################################################
### empres-i data

empres_europe_pos <- read.csv("data/flu_data/empresi_wild_confirmed_europe.csv")

# we don't need to humans affected of human deaths columns. 
empres_europe_pos <- dplyr::select(empres_europe_pos, !c(Humans.affected, Human.deaths))

# The species entries are much less standardised in these data
emp_species <-  as.data.frame(table(empres_europe_pos$Species))
colnames(emp_species) <- c("species", "counts")


## empres i data are entered slightly chaotically
## Try removing everything from ( onwards

emp_species <- emp_species %>% mutate(sp_edit = str_remove(species, ' \\(.*')) %>%
  mutate(sp_edit =  str_remove(sp_edit, '\\:.*')) %>%
  mutate(sp_edit = str_remove(sp_edit, '\\?.*'))

# let's now combine species where possible

emp_species <- emp_species %>%
  group_by(sp_edit) %>%
  summarise(counts = sum(counts)) 

#colnames(emp_species) <- c("species", "counts")

# first order
for (i in 1:nrow(emp_species)) {
  emp_species[i,"order"] <- tax_name(emp_species[i,"sp_edit"], 
                                        get = "order",
                                        db = "ncbi")$order 
  print(i)
}

emp_orders <- emp_species %>%
  group_by(order) %>%
  summarise(emp_order_count = sum(counts))

## now families

for (i in 1:nrow(emp_species)) {
  emp_species[i,"family"] <- tax_name(emp_species[i,"sp_edit"], 
                                      get = "family",
                                      db = "ncbi")$family 
  print(i)
}

emp_families <- emp_species %>%
  group_by(family) %>%
  summarise(emp_family_count = sum(counts))

## and finally genus
for (i in 67:nrow(emp_species)) {
  emp_species[i,"genus"] <- tax_name(emp_species[i,"sp_edit"], 
                                      get = "genus",
                                      db = "ncbi")$genus 
  print(i)
}

emp_genus <- emp_species %>%
  group_by(genus) %>%
  summarise(emp_genus_count = sum(counts))


## combine the databases results


order_counts <- full_join(bvbrc_orders, emp_orders)
family_counts <- full_join(bvbrc_families, emp_families)
genus_counts <- full_join(bvbrc_genus, emp_genus)


#write.csv(emp_species, "output/empres_full_data.csv")
#write.csv(bvbrc_species, "output/bvbrc_full_data.csv")
#write.csv(order_counts, "output/order_counts.csv")
#write.csv(family_counts, "output/family_counts.csv")
#write.csv(genus_counts, "output/genus_counts.csv")
