# In this script we load in 1000 MCMC samples of trees for Aves downloaded from
# https://data.vertlife.org/?basetree=birdtree&start_folder=Stage2/ and use ape
# to get a single consensus tree.

require(ape)
require(dplyr)
require(TreeTools)
bird_tree <- read.tree("eco_phylo_processing/BirdzillaHackett10.tre")

# Do a benchmarking analysis to get an estimate on how long the consensus tree
# construction from all 1000 trees will take:
subset_sizes <- c(200, 100, 50, 10)
consensus_times <- c()
for (s in subset_sizes) {
  this_tree <- bird_tree[1:s]
  start.time <- Sys.time()
  consensus_tree <- consensus(this_tree, rooted = TRUE)
  end.time <- Sys.time()
  elapsed <- as.numeric(difftime(end.time, start.time, units="secs"))
  cat("Built consensus tree from",
      s,
      "subtrees in",
      elapsed,
      "seconds.\n")
  consensus_times <- append(consensus_times, elapsed)
}

benchmark_df <- data.frame(subset_sizes, consensus_times)
bm_fit <- lm(consensus_times~subset_sizes, benchmark_df)

plot(subset_sizes,
     consensus_times,
     xlab = "Number of trees",
     ylab = "Execution time in seconds")
abline(bm_fit)

cat("Expected execution time for 1000 trees is",
    bm_fit$coefficients[1] + bm_fit$coefficients[2]*1000,
    "seconds.\n")

# Now get full consensus tree and save to file
{
  start.time <- Sys.time()
  consensus_tree <- consensus(bird_tree,
                              p=0.5,
                              rooted = TRUE)
  end.time <- Sys.time()
  elapsed <- as.numeric(difftime(end.time, start.time, units="secs"))
  cat("Built consensus tree from 1000 subtrees in",
      elapsed,
      "seconds.\n")
}
write.nexus(consensus_tree,
            file = "data/phylogeny/consensus_bird_tree.nex",
            translate = FALSE)
consensus_copy <- read.nexus("data/phylogeny/consensus_bird_tree.nex")

# Check it writes and reads properly by working out distance matrices for both:
dmat <- cophenetic.phylo(consensus_tree)

# Check for NaN's in the distance matrix
rows_to_inspect <- 1:50
rownames(dmat)[which(is.nan(rowSums(dmat[rows_to_inspect, ])))]
# Should find all the first fifty rows have a NaN somewhere!

# Try identifying the species where we have NaN's:
nan_sp_list <- list()
for (sp in 1:nrow(dmat)){
  nan_sp <- colnames(dmat)[which(is.nan(dmat[sp, ]))]
  nan_sp_list <- append(nan_sp_list, list(nan_sp))
}

# See if it's the same entries for all species:
unique_nan_sp_list <- unique(nan_sp_list)
unique_nan_sp_list %>% lengths %>% table
# We should see that there's a list of 43 "bad species" that consistently show
# up, plus 37 of length 9992 (all species except the one in question) plus three
# of length 9991 (all species except two). This suggests there are 37 completely
# isolated species and three pairs of species that are only connected to each
# other.
sp_pair_lists <- unique_nan_sp_list[which(lengths(unique_nan_sp_list)==9991)]
sp_names <- rownames(dmat)
for (i in 1:3){
  sp_pair <- setdiff(sp_names, unlist(sp_pair_lists[i]))
  cat(sp_pair[1], "and", sp_pair[2], "are paired.\n")
}

# Go back to bird tree, and check if there are any unrooted or non-binary trees:
for (i in seq_along(bird_tree)){
  if (!is.rooted(bird_tree[i])){
    cat("Tree", i, "is unrooted.\n")
  }
  if (!is.binary(bird_tree[i])){
    cat("Tree", i, "is not binary\n")
  }
}

# We have forced the consensus tree to be rooted, but it may not be binary:
{
  if (is.binary(consensus_tree)){
    print("Consensus tree is binary.")
  }
  else{
    print("Consensus tree is not binary.")
  }
}
