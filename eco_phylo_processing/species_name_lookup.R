# In this script we load in the collection of species names from AVONET and
# explore the ways they conflict

library(readxl)

# Load the page from AVONET data which includes the multiple name formats
AVONET_df <- read_excel("data/AVONETSupplementarydataset1.xlsx",
                        sheet = "AVONET_Raw_Data")

# Extract names under each naming scheme
BirdLife_names <- AVONET_df$Species1_BirdLife
eBird_names <- AVONET_df$Species2_eBird
BirdTree_names <- AVONET_df$Species3_BirdTree

# Identify all species where the three data sources are consistent
nonconflict_rows <- which(
  (BirdLife_names==eBird_names)&(BirdLife_names==BirdTree_names))

# We're only interested in inconsistent species so we remove the consistent ones
# from the dataframe:
AVONET_df <- AVONET_df[-nonconflict_rows,
                       c("Avibase.ID",
                         "Species1_BirdLife",
                         "Species2_eBird",
                         "Species3_BirdTree")]

# Reduce to unique combinations
AVONET_df <- distinct(AVONET_df)
cat("Total number of conflicts is", nrow(AVONET_df))
