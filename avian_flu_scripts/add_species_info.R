## Lines 1-50 are for use within create_flu_csv scripts if updates to raw data have occurred

# use Taxize to see how many can be identified automatically
# If no updates to data can skip this part and just read in the csv below.

library(taxize)

#for Europe
#species_list <- as.data.frame(table(ai_data_prj_area$Species))
#for asia and americas
species_list <- as.data.frame(table(americas_asia_data$Species))

colnames(species_list) <- c("species", "freq")
species_list$species <- as.character(species_list$species)

for (i in 410:nrow(species_list)) {
  species_list[i,"class"] <- tax_name(species_list[i,"species"],
                                      get = "class",
                                      db = "ncbi")$class
  print(i)
}

not_sp <- species_list[which(is.na(species_list$class)),]
# manual inspection suggests some of these may not be found due to things like having 'incognita' added to name

not_sp$species_edit <- not_sp$species
not_sp$species_edit <- sub("\\(.*", "", not_sp$species_edit)
not_sp$species_edit <- sub("\\:.*", "", not_sp$species_edit)

for (i in 71:nrow(not_sp)) {
  not_sp[i,"class"] <- tax_name(not_sp[i,"species_edit"],
                                get = "class",
                                db = "ncbi")$class
  print(i)
}


# still quite a few without a class which may have to manually look through

not_sp <- rename(not_sp, class_2 = class)
species_list <- left_join(species_list, not_sp)
species_list[which(is.na(species_list$class)), "class"] <- species_list[which(is.na(species_list$class)), "class_2"]
species_list <- dplyr::select(species_list, c("species", "freq", "class"))

# Peacock has been mislabelled as an insect so change to Aves

species_list[which(species_list$species == "Peacock"),"class"] <- "Aves"


# select the ones labelled "Mammalia" for removal
# save the species list so don't have to generate every time'

#write.csv(species_list, "avian_flu_scripts/detecting_mammals.csv", row.names = F)
#write.csv(species_list, "avian_flu_scripts/detecting_mammals_asia_americas.csv", row.names = F)

#####################################################################################

# trying to automate the process of identifying the taxons

species_list$species <- as.character(species_list$species)

## automatically remove any labels such as 'sp' or 'spp' in species
species_list$species <- str_remove(species_list$species, pattern = " sp.")

library(taxize)

# remove the brackets as this helped earlier
species_list$species <- sub("\\(.*", "", species_list$species)
species_list$species <- sub("\\:.*", "", species_list$species)
species_list$species <- str_squish(species_list$species)

#for (i in 1:20) {
for (i in 407:nrow(species_list)) {
  species_list[i,"order"] <- tax_name(species_list[i,"species"], 
                                       get = "order",
                                       db = "ncbi")$order 
  species_list[i,"family"] <- tax_name(species_list[i,"species"], 
                                       get = "family",
                                       db = "ncbi")$family
  species_list[i,"genus"] <- tax_name(species_list[i,"species"], 
                                      get = "genus",
                                      db = "ncbi")$genu
  print(i)
}


no_family <- species_list[which(is.na(species_list$family)), "species"]
no_family

species_list[which(species_list$species == "Aixnsa"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Aix")

species_list[which(species_list$species == "Anas clypeata"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Spatula")

species_list[which(species_list$species == "Anser fabalis serrirostris"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Anser")

species_list[which(species_list$species == "Arctic Skua"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "Stercorariidae", "Stercorarius")

species_list[which(species_list$species == "Barn Owl"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Strigiformes", "Tytonidae", "Tyto")

species_list[which(species_list$species == "Bearded Vulture"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Accipitriformes", "Accipitridae", "Gypaetus")

species_list[which(species_list$species == "Black Backed Gull"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "Laridea", "Larus")

species_list[which(species_list$species == "Black Guillemot"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "Alcidea", "Cepphus")

species_list[which(species_list$species == "Brent Goose"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidea", "Branta")

species_list[which(species_list$species == "Buzzard"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Accipitriformes", "Accipitridae", "")

species_list[which(species_list$species == "Common Wood-Pigeon"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Columbiformes", "Columbidae", "Columba")

species_list[which(species_list$species == "Crane"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Gruiiformes", "Gruoidae", "")

species_list[which(species_list$species == "Crow"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Passeriformes", "Corvidae", "")

species_list[which(species_list$species == "Curlew"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "Scolopacidae", "Numenius")

species_list[which(species_list$species == "Cygnus columbianus bewickii"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Cygnus")

species_list[which(species_list$species == "Egret"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Pelecaniformes", "Ardeidae", "")

species_list[which(species_list$species == "Eurasian Buzzard"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Accipitriformes", "Accipitridae", "Buteo")

species_list[which(species_list$species == "Eurasian Jackdaw"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Passeriformes", "Corvidae", "Coloeus")

species_list[which(species_list$species == "European Turtle Dove"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Columbiformes", "Columbidae", "Sreptopelia")

species_list[which(species_list$species == "Falcon"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Falconiformes", "Falconidae", "Falco")

species_list[which(species_list$species == "Goosander"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Mergus")

species_list[which(species_list$species == "Greater White-Fronted Goose"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Anser")

species_list[which(species_list$species == "Grebe"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Podicipediformes", "Podicipedidae", "")

species_list[which(species_list$species == "Hawk"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Accipitriformes", "Accipitridae", "")

species_list[which(species_list$species == "Heron"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Pelecaniformes", "Ardeidae", "")

species_list[which(species_list$species == "Kite"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Accipitriformes", "Accipitridae", "")

species_list[which(species_list$species == "Larus brunicephalus"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "Laridea", "Chroicocephalus")

species_list[which(species_list$species == "Magpie"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Passeriformes", "Corvidae", "")

species_list[which(species_list$species == "Northern Crested Caracara"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Falconiformes", "Falconidae", "Caracara")

species_list[which(species_list$species == "Owl"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Strigiformes", "", "")

species_list[which(species_list$species == "Oyster Catcher"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "Haematopodidae", "Haematopus")

species_list[which(species_list$species == "Peregrin Falcon"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Falconiformes", "Falconidae", "Falco")

species_list[which(species_list$species == "Phalacrocorax pygmaeus"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Suliformes", "Phalacrocoracidae", "Microcarbo")

species_list[which(species_list$species == "Pigeon"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Columbiformes", "Columbidae", "Columba")

species_list[which(species_list$species == "Rose Pelican"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Pelecaniformes", "Pelecanidae", "Pelecanus")

species_list[which(species_list$species == "Silver Teal"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Spatula")

species_list[which(species_list$species == "Stork"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Ciconiiformes", "Ciconiidae", "")

species_list[which(species_list$species == "Sulids"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Suliformes", "Sulidae", "")

species_list[which(species_list$species == "Swallow"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Passeriformes", "Hirundinidae", "")

species_list[which(species_list$species == "Swan"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Cygnus")

species_list[which(species_list$species == "Swift Tern"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "Laridea", "Thalasseus")

species_list[which(species_list$species == "Teal"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Anas")

species_list[which(species_list$species == "Tern"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "Laridea", "")

species_list[which(species_list$species == "Thrush"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Passeriformes", "Turdidae", "")

species_list[which(species_list$species == "Variable Hawk"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Accipitriformes", "Accipitridae", "")

species_list[which(species_list$species == "Vulture"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Accipitriformes", "Accipitridae", "")

species_list[which(species_list$species == "Wader"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "", "")

species_list[which(species_list$species == "Water Hen"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Gruiformes", "Rallidae", "Gallinula")

species_list[which(species_list$species == "Whooper Swan"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Cygnus")

species_list[which(species_list$species == "Wigeons"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "Mareca")

species_list[which(species_list$species == "Wild Duck"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Anseriformes", "Anatidae", "")



#write.csv(species_list, "avian_flu_scripts/species_with_family.csv")

##############################
## May not need the below. 


## For bvbrc data there are relatively few that haven't been assigned.
## We can add these manually

species_list[which(species_list$species == "Anas clypeata"), c("order", "family", "genus")] <-
  c("Anseriformes", "Anatidae", "Spatula")

species_list[which(species_list$species == "Anser fabalis serrirostris"), c("order", "family", "genus")] <-
  c("Anseriformes", "Anatidae", "Anser")

species_list[which(species_list$species == "Cygnus columbianus bewickii"), c("order", "family", "genus")] <-
  c("Anseriformes", "Anatidae", "Cygnus")



### old code - may no longer need

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
