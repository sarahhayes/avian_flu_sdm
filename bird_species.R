## Considering which species we might want to include in the bird distribution data.
## One aspect to take into account is which species have positive results 
## for avian influenza reported. 

## Will look at both data sets but remember that EMPRES-i only contains HPAI and thus
## bv-brc MAY be more useful? 

rm(list = ls())

library(tidyverse)

pos_data <- read.csv("data/flu_data/BVBRC_wild_captivewild_positive.csv")
columns_needed <- c("Collection.Date", "Collection.Year", "Collection.Country", 
                    "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                    "Host.Species", "Host.Group", "Host.Natural.State", "Host.Health")

pos_data <- select(pos_data, all_of(columns_needed))


bvbrc_species <-as.data.frame(table(pos_data$Host.Species))
colnames(bvbrc_species)<- c("species", "count")

range(bvbrc_species$count)
table(bvbrc_species$count)
hist(bvbrc_species$count, breaks = seq(0,17500,10), xlim = c(0,1000), col = 'orange',
     main = "", xlab = "Species count")

# For these data, splitting the name into two columns will probably help us. 

bvbrc_grp_sp <- bvbrc_species %>%
  separate(species, c("genus","species_detail"), sep = " ") %>%
  dplyr::select(genus, count) %>%
  group_by(genus) %>%
  summarise(total_count = sum(count))

range(bvbrc_grp_sp$total_count)
table(bvbrc_grp_sp$total_count)
hist(bvbrc_grp_sp$total_count, breaks = seq(0,26000,10), xlim = c(0,1000), col = 'turquoise',
     main = "Counts by genus (x-axis truncated at 1000)", xlab = "Genus count")

# these aren't actually all genus - some are order i.e. there is a charadriiformes count 
# what I really need to know is which of these aren't anseriformes or charadriiformes

#write.csv(bvbrc_grp_sp, "output/bvbrc_species_groups.csv", row.names = F)

# I have manually added the order to the table and saved it as a different name
# Now read that back in

bvbrc_order <- read.csv("output/bvbrc_species_groups_order_added_manually.csv")

order_counts <- bvbrc_order %>%
  dplyr::select(total_count, order) %>%
  group_by(order) %>%
  summarise(order_count = sum(total_count))

#write.csv(order_counts, "output/order_counts.csv", row.names = F)

################################################################################
## Now let's look at Empres-i

empres_europe_pos <- read.csv("data/flu_data/empresi_wild_confirmed_europe.csv")

# we don't need to humans affected of human deaths columns. 
empres_europe_pos <- dplyr::select(empres_europe_pos, !c(Humans.affected, Human.deaths))

# The species entries are much less standardised in these data
emp_species <-  as.data.frame(table(empres_europe_pos$Species))
colnames(emp_species) <- c("species", "counts")

range(emp_species$counts)
hist(emp_species$counts, breaks = seq(0,600,10), xlab = "Count of a species", 
     main = "Counts of species from Empres-i", col = "light blue")

# There are 220 unique entries for the species. The histogram shows that a lot of these 
# have fewer than 10 entries. Let's look at those with 10 or more entries only

emp_species_ten <- emp_species %>%
  filter(counts>= 10)

# write.csv(emp_species_ten, "output/empres_species_counts.csv", row.names = F)

empres_order <- read.csv("output/empres_species_counts_order_manually_added.csv")

empres_order_counts <- empres_order %>%
  dplyr::select(counts, order) %>%
  group_by(order) %>%
  summarise(order_count = sum(counts))

#write.csv(empres_order_counts, "output/empres_order_counts.csv", row.names = F)

