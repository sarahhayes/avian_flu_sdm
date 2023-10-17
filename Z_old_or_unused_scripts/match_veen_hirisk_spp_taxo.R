library(tidyverse)
library(magrittr)

# High risk spp from Veen et al tables 2.1 - 2.3
hirisk <- structure(list(Scientific = c("Anas acuta", "Anas clypeata", "Anas crecca", 
                              "Anas penelope", "Anas platyrhynchos", "Anas querquedula", "Anser albifrons", 
                              "Anser anser", "Anser brachyrhynchus", "Anser erythropus", "Anser fabalis", 
                              "Ardea alba", "Ardea cinerea", "Ardea purpurea", "Ardeola ralloides", 
                              "Aythya ferina", "Aythya fuligula", "Branta bernicla", "Branta canadensis", 
                              "Branta leucopsis", "Branta ruficollis", "Bubulcus ibis", "Ciconia ciconia", 
                              "Columba oenas", "Columba palumbus", "Corvus frugilegus", "Corvus monedula", 
                              "Cygnus columbianus", "Cygnus olor", "Egretta garzetta", "Fringilla coelebs", 
                              "Fringilla montifringilla", "Fulica atra", "Fulica cristata", 
                              "Hirundo rustica", "Larus canus", "Larus ridibundus", "Limosa limosa", 
                              "Marmaronetta angustirostris", "Netta rufina", "Nycticorax nycticorax", 
                              "Passer domesticus", "Passer hispaniolensis", "Pelecanus crispus", 
                              "Pelecanus onocrotalus", "Phalacrocorax carbo", "Phalacrocorax pygmeus", 
                              "Philomachus pugnax", "Platalea leucorodia", "Plegadis falcinellus", 
                              "Pluvialis apricaria", "Podiceps cristatus", "Riparia riparia", 
                              "Streptopelia decaocto", "Sturnus unicolor", "Sturnus vulgaris", 
                              "Turdus iliacus", "Turdus pilaris", "Vanellus vanellus")), class = "data.frame", row.names = c(NA, 
                                                                                                                             -59L)) %>%
  mutate(Scientific = tolower(Scientific))

# Read EltonTraits for taxonomy
EltonTraits_df <- 
  read.table("data/variables/bird_data/elton_traits/BirdFuncDat.txt",
             sep = '\t',
             quote="\"",
             header = TRUE)
# elton trait data has two rows of NA's at the bottom which we need to remove:
EltonTraits_df <- EltonTraits_df[1:9993, ]
EltonTraits_df$Scientific <- tolower(EltonTraits_df$Scientific)
EltonTraits_df %<>% tidyr::separate(Scientific, c("Genus", "Spp"), " ", remove = FALSE)

hirisk %>% left_join(EltonTraits_df %>% select(IOCOrder, BLFamilyLatin, Genus, Spp, Scientific)) %>%
  mutate(BLFamilyLatin = ifelse(Scientific == "ardea alba", "Ardeidae", BLFamilyLatin), # fill in mismatch spp manually
         IOCOrder = ifelse(Scientific == "ardea alba", "Pelecaniformes", IOCOrder)) %>%
  arrange(IOCOrder, BLFamilyLatin, Genus, Spp) %>%
  write.csv("hirisk.csv")

