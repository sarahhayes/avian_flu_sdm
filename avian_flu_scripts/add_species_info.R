## Lines 1-50 are for use within create_flu_csv scripts if updates to raw data have occurred

# use Taxize to see how many can be identified automatically
# If no updates to data can skip this part and just read in the csv below.

library(taxize)

#for Europe
species_list <- as.data.frame(table(ai_data_prj_area$Species))
#for asia and americas
#species_list <- as.data.frame(table(americas_asia_data$Species))

colnames(species_list) <- c("species", "freq")
species_list$species <- as.character(species_list$species)

for (i in 294:nrow(species_list)) {
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

for (i in 121:nrow(not_sp)) {
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
for (i in 410:nrow(species_list)) {
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

species_list[which(species_list$species == "Cape Gannet Shy Albatross"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Procellariiformes", "Diomedeidae", "")

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

species_list[which(species_list$species == "Partridge"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Galliformes", "Phasianidae", "")

species_list[which(species_list$species == "Peregrin Falcon"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Falconiformes", "Falconidae", "Falco")

species_list[which(species_list$species == "Phalacrocorax pygmaeus"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Suliformes", "Phalacrocoracidae", "Microcarbo")

species_list[which(species_list$species == "Pigeon"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Columbiformes", "Columbidae", "Columba")

species_list[which(species_list$species == "Purple Sandpiper"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Charadriiformes", "Scolopacidae", "Calidris")

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

species_list[which(species_list$species == "Umbrella Cockatoo"), 
             c("class", "order", "family", "genus")] <- 
  c("Aves", "Psittaciformes", "Cacatuidae", "Cacatua")

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

# there are two lines with black guillemot that haven't been combined. Seven birds in total. 
# Combine these and add manually. 

species_list[which(species_list$species == "Black Guillemot"),]

species_list <- species_list %>% dplyr::filter(species != "Black Guillemot")
species_list <- rbind(species_list, c("Black Guillemot", 7, "Aves", "Charadriiformes", "Alcidea", "Cepphus"))

#write.csv(species_list, "avian_flu_scripts/species_with_family.csv")


