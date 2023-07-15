# trying to automate the process of identifying the taxons

species_table$species <- as.character(species_table$species)

## automatically remove any labels such as 'sp' or 'spp' in species
species_table$species <- str_remove(species_table$species, pattern = " sp.")

# reports bad request intermittently. 
# Connection issue? Can just restart from where it had got to. 

library(taxize)

# first add order to table
for (i in 1:nrow(species_table)) {
  species_table[i,"order"] <- tax_name(species_table[i,"species"], 
                                       get = "order",
                                       db = "ncbi")$order 
  print(i)
}


## now add families
for (i in 15:nrow(species_table)) {
  species_table[i,"family"] <- tax_name(species_table[i,"species"], 
                                        get = "family",
                                        db = "ncbi")$family 
  print(i)
}

# finally add genus
for (j in 40:nrow(species_table)) {
  species_table[j,"genus"] <- tax_name(species_table[j,"species"], 
                                       get = "genus",
                                       db = "ncbi")$genus
  print(j)
}

## For bvbrc data there are relatively few that haven't been assigned.
## We can add these manually

species_table[which(species_table$species == "Anas clypeata"), c("order", "family", "genus")] <-
  c("Anseriformes", "Anatidae", "Spatula")

species_table[which(species_table$species == "Anser fabalis serrirostris"), c("order", "family", "genus")] <-
  c("Anseriformes", "Anatidae", "Anser")

species_table[which(species_table$species == "Cygnus columbianus bewickii"), c("order", "family", "genus")] <-
  c("Anseriformes", "Anatidae", "Cygnus")


## repeat for the negative data as useful to have to get % positive

bvbrc_species_neg$species <- as.character(bvbrc_species_neg$species)

# Really want this by species so just use the info from those that have the full species in
# both positive and negative. 

bvbrc_species_pos <- bvbrc_species
colnames(bvbrc_species_pos) <- c("species", "count_pos")

bvbrc_pos_neg <- merge(bvbrc_species_neg, bvbrc_species_pos) # default for merge is to just include 
# columns that have entries for both. Probably most suitable here

bvbrc_pos_neg <- full_join(bvbrc_species_neg, bvbrc_species_pos)


bvbrc_pos_neg[,"total_tested"] <- bvbrc_pos_neg[,"count_pos"] + bvbrc_pos_neg[,"count_neg"]
bvbrc_pos_neg[,"percent_pos"] <- (bvbrc_pos_neg[,"count_pos"]/bvbrc_pos_neg[,"total_tested"])*100

# quite a few with high % because so few were tested. 
# Filter to those where at least 10 tested

bvbrc_pos_neg_fifty <- bvbrc_pos_neg %>%
  filter(total_tested >= 50)

# the list below isn't exhaustive as yet
# bvbrc_pos_neg[which(bvbrc_pos_neg$species %in% c("Anas sp.", "Anas sp", "Anas spp.")),"species"] <- "Anas"
# bvbrc_pos_neg[which(bvbrc_pos_neg$species == "Cygnus sp."),"species"] <- "Cygnus"
# bvbrc_pos_neg[which(bvbrc_pos_neg$species == "Gallinago spp"),"species"] <- "Gallinago"
# bvbrc_pos_neg[which(bvbrc_pos_neg$species == "Anas sp."),"species"] <- "Anas"
# bvbrc_pos_neg[which(bvbrc_pos_neg$species == "Columba sp."),"species"] <- "Columba"
# bvbrc_pos_neg[which(bvbrc_pos_neg$species == "Larus sp."),"species"] <- "Larus"
# bvbrc_pos_neg[which(bvbrc_pos_neg$species == "Larus sp"),"species"] <- "Larus"
# 

## now do for families
for (i in 1077:nrow(bvbrc_pos_neg)) {
  bvbrc_pos_neg[i,"family"] <- tax_name(bvbrc_pos_neg[i,"species"], 
                                        get = "family",
                                        db = "ncbi")$family 
  print(i)
}

bvbrc_families_pos_neg <- bvbrc_pos_neg %>%
  group_by(family) %>%
  summarise(bvbrc_family_count = sum(count))

unique(bvbrc_pos_neg$family)

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


emp_species <- emp_species %>%
  group_by(sp_edit) %>%
  summarise(counts = sum(counts)) 

#colnames(emp_species) <- c("species", "counts")

# I want to know how many species across the different data frames
# However bvbrc as scientific and empres are common

emp_species_df <- as.data.frame(emp_species)

sp_uids <- get_uid(emp_species_df[,"sp_edit"])   
uids.found <- as.uid(sp_uids[!is.na(sp_uids)])

sci.names <- comm2sci(uids.found, db = 'ncbi')

sci_names_emp <- unlist(sci.names)
sci_names_bvbrc <- bvbrc_pos_neg$species

sci_names_comb <- unique(c(sci_names_bvbrc, sci_names_emp))
sci_names_comb_df <- as.data.frame(sci_names_comb)
# let's now combine databases where possible

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
