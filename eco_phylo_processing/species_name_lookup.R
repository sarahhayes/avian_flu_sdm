# In this script we load in the collection of species names from AVONET and
# explore the ways they conflict

library(readxl)

# Load the page from AVONET data which includes the multiple name formats
AVONET_df <- read_excel("data/AVONETSupplementarydataset1.xlsx",
                        sheet = "AVONET_Raw_Data")
AVONET_df <- AVONET_df[,
                       c("Avibase.ID",
                         "Species1_BirdLife",
                         "Species2_eBird",
                         "Species3_BirdTree")]
AVONET_df <- distinct(AVONET_df)
AVONET_df[, 2:4] <- sapply(AVONET_df[, 2:4], tolower)

# Extract names under each naming scheme
BirdLife_names <- AVONET_df$Species1_BirdLife
eBird_names <- AVONET_df$Species2_eBird
BirdTree_names <- AVONET_df$Species3_BirdTree

# Identify all species where the three data sources are consistent
nonconflict_rows <- which(
  (BirdLife_names==eBird_names)&(BirdLife_names==BirdTree_names))

# We're only interested in inconsistent species so we remove the consistent ones
# from the dataframe:
AVONET_conflicts <- AVONET_df[-nonconflict_rows, ]

# Reduce to unique combinations
AVONET_conflicts <- distinct(AVONET_conflicts)
cat("Total number of conflicts is", nrow(AVONET_conflicts))

# Assign ID based on consensus
get_single_bird_id <- function(binomial_name){
  # Return -100 if provided with na
  if (binomial_name=="na"){
    return("-100")
  }
  # Find IDs from each list and keep track of which species list they came from.
  ID_candidates <- c()
  ID_weightings <- c()
  ID_sources <- c()
  if (binomial_name %in% AVONET_df$Species1_BirdLife){
    ID_candidates <- append(ID_candidates, AVONET_df$Avibase.ID[
      which(AVONET_df$Species1_BirdLife==binomial_name)])
    no_matches <- length(which(AVONET_df$Species1_BirdLife==binomial_name))
    for (i in 1:no_matches){
      ID_weightings <- append(ID_weightings, 1/no_matches)
      ID_sources <- append(ID_sources, "BirdLife")
    }
  }
  if (binomial_name %in% AVONET_df$Species2_eBird){
    ID_candidates <- append(ID_candidates, AVONET_df$Avibase.ID[
      which(AVONET_df$Species2_eBird==binomial_name)])
    no_matches <- length(which(AVONET_df$Species2_eBird==binomial_name))
    for (i in 1:no_matches){
      ID_weightings <- append(ID_weightings, 1/no_matches)
      ID_sources <- append(ID_sources, "eBird")
    }
  }
  if (binomial_name %in% AVONET_df$Species3_BirdTree){
    ID_candidates <- append(ID_candidates, AVONET_df$Avibase.ID[
      which(AVONET_df$Species3_BirdTree==binomial_name)])
    no_matches <- length(which(AVONET_df$Species3_BirdTree==binomial_name))
    for (i in 1:no_matches){
      ID_weightings <- append(ID_weightings, 1/no_matches)
      ID_sources <- append(ID_sources, "BirdTree")
    }
  }
  # How we work out the consensus depends on how many candidates are identified
  if (length(ID_candidates)==0){
    # If ID doesn't appear in any list, return -100
    return("-100")
  }
  else{
    # Construct consensus distribution of candidate names
    unique_candidates <- unique(ID_candidates)
    if (length(unique_candidates)==1){
      # cat("Number of unique candidates is", length(unique_candidates), "\n")
      # cat(unique_candidates, "\n")
      return(unique_candidates)
    }
    candidate_dist <- c()
    for (i in 1:length(unique_candidates)){
      # cat("i=", i, "\n")
      candidate_dist <- append(candidate_dist,
                               sum(ID_weightings[
                                 which(ID_candidates==unique_candidates[i])]))
    }
    # cat(unique_candidates, "\n", ID_candidates, "\n", candidate_dist, "\n", ID_weightings, "\n")
    if (length(which(candidate_dist==max(candidate_dist)))==1){
      return(unique_candidates[which(candidate_dist==max(candidate_dist))])
    }
    else{
      # cat("Multiple consensus options for species", binomial_name,"\n")
      return(sample(unique_candidates[which(candidate_dist==max(candidate_dist))], 1))
    }
  }
}

get_bird_ids <- function(name_vect){
  id_vect <- vector(length = length(name_vect))
  for (n in 1:length(name_vect)){
    # cat("n=",n,"\n")
    id_vect[n] <- get_single_bird_id(name_vect[n])
  }
  return(id_vect)
}

# Need to test this is doing the right thing:
test_id_list1 <- get_bird_ids(AVONET_df$Species1_BirdLife)
test_id_list2 <- get_bird_ids(AVONET_df$Species2_eBird)
test_id_list3 <- get_bird_ids(AVONET_df$Species3_BirdTree)

list1_mistmatch <- which(test_id_list1!=AVONET_df$Avibase.ID)
mismatch_species1 <- AVONET_df$Species1_BirdLife[list1_mistmatch]
list2_mistmatch <- which(test_id_list2!=AVONET_df$Avibase.ID)
mismatch_species2 <- AVONET_df$Species2_eBird[list2_mistmatch]
list3_mistmatch <- which(test_id_list3!=AVONET_df$Avibase.ID)
mismatch_species3 <- AVONET_df$Species3_BirdTree[list3_mistmatch]

# Compare differing number of mismatches to differing numbers of species:
cat("BirdLife has",
    length(unique(AVONET_df$Species1_BirdLife))-
      length(unique(AVONET_df$Species2_eBird)),
    "more species than eBird, and eBird produces",
    length(mismatch_species2) - length(mismatch_species1),
    "more mismatches.")
cat("BirdLife has",
    length(unique(AVONET_df$Species1_BirdLife))-
      length(unique(AVONET_df$Species3_BirdTree)),
    "more species than BirdTree, and BirdTree produces",
    length(mismatch_species3) - length(mismatch_species1),
    "more mismatches.")
# These values are close together, suggesting the number of mismatches is
# closely related to the degree of splitting/lumping.
# Try comparing directly with number of ID's instead:
cat("BirdLife has",
    length(unique(AVONET_df$Avibase.ID))-
      length(unique(AVONET_df$Species1_BirdLife)),
    "fewer species than AVIBASE, and produces",
    length(mismatch_species1),
    "mismatches.")
cat("eBird has",
    length(unique(AVONET_df$Avibase.ID))-
      length(unique(AVONET_df$Species2_eBird)),
    "fewer species than AVIBASE, and produces",
    length(mismatch_species2),
    "mismatches.")
cat("BirdTree has",
    length(unique(AVONET_df$Avibase.ID))-
      length(unique(AVONET_df$Species3_BirdTree)),
    "fewer species than AVIBASE, and produces",
    length(mismatch_species3),
    "mismatches.")
# How does it compare to number of distinct combinations of ID and species name?
cat("There are",
    nrow(distinct(AVONET_df[, c("Avibase.ID", "Species1_BirdLife")])),
    "distinct combinations of Avibase ID and BirdLife names, and",
    nrow(AVONET_df)-length(mismatch_species1),
    "successful matches.")
cat("There are",
    nrow(distinct(AVONET_df[, c("Avibase.ID", "Species2_eBird")])),
    "distinct combinations of Avibase ID and BirdLife names, and",
    nrow(AVONET_df)-length(mismatch_species2),
    "successful matches.")
cat("There are",
    nrow(distinct(AVONET_df[, c("Avibase.ID", "Species3_BirdTree")])),
    "distinct combinations of Avibase ID and BirdLife names, and",
    nrow(AVONET_df)-length(mismatch_species3),
    "successful matches.")

length(unique(test_id_list1))
length(unique(test_id_list2))
length(unique(test_id_list3))
