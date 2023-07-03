# In this script we load in 1000 MCMC samples of trees for Aves downloaded from
# https://data.vertlife.org/?basetree=birdtree&start_folder=Stage2/ and use ape
# to get a single consensus tree.

require(ape)
bird_tree <- read.tree("eco_phylo_processing/BirdzillaHackett10.tre")

# Do a benchmarking analysis to get an estimate on how long the consensus tree
# construction from all 1000 trees will take:
subset_sizes <- c(200, 100, 50, 10)
consensus_times <- c()
for (s in subset_sizes) {
  this_tree <- bird_tree[1:s]
  start.time <- Sys.time()
  consensus_tree <- consensus(this_tree)
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
  consensus_tree <- consensus(bird_tree)
  end.time <- Sys.time()
  elapsed <- as.numeric(difftime(end.time, start.time, units="secs"))
  cat("Built consensus tree from 1000 subtrees in",
      elapsed,
      "seconds.\n")
}
write.nexus(consensus_tree,
            file = "data/phylogeny/consensus_bird_tree",
            translate = FALSE)
