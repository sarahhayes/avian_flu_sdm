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
            file = "data/phylogeny/consensus_bird_tree.nex",
            translate = FALSE)
consensus_copy <- read.nexus("data/phylogeny/consensus_bird_tree.nex")

# Check it writes and reads properly by working out distance matrices for both:
dmat <- cophenetic.phylo(consensus_tree)
dmat_copy <- cophenetic.phylo(consensus_copy)
dmat_order <- dmat[, order(colnames(dmat))]
dmat_order <- dmat_order[order(rownames(dmat)), ]
dmat_copy_order <- dmat_copy[, order(colnames(dmat_copy))]
dmat_copy_order <- dmat_copy_order[order(rownames(dmat_copy)), ]
cat("Max difference in sorted distance matrices is",
    max(dmat_order - dmat_copy_order),
    ".\n")
# The maximum difference is NaN, so there must be NaN's in one or both of the
# distance matrices. What species do these correspond to?
rows_to_inspect <- 1:50
rownames(dmat_order)[which(is.nan(rowSums(dmat_order[rows_to_inspect, ])))]
# Should find all the first fifty rows have a NaN somewhere!

# Try identifying the species where we have NaN's:
nan_sp_list <- list()
for (sp in 1:5){
  nan_sp <- colnames(dmat_order)[which(is.nan(dmat_order[sp, ]))]
  cat("NaN distance between",
      rownames(dmat_order)[sp],
      "and",
      nan_sp,
      ".\n",
      sep = "\n")
  nan_sp_list <- append(nan_sp_list, list(nan_sp))
}

# Appears to be the same species consistently, so let's keep going:
for (sp in 6:nrow(dmat_order)){
  nan_sp <- colnames(dmat_order)[which(is.nan(dmat_order[sp, ]))]
  nan_sp_list <- append(nan_sp_list, list(nan_sp))
}

# See if it's the same entries for all species:
unique_nan_sp_list <- unique(nan_sp_list)

# We should see that there's a list of 43 "bad species" that consistently show
# up, plus 40 other lists consisting of all the species minus one or two.
# To me this suggests that we have 39 species that don't link to everything plus
# 4 that are in two 2-species clusters.

sp_names <- rownames(dmat)
n_samples <- 10
sp_sample <- sample(sp_names, 2 * n_samples)
dsamples <- c()
dsamples_copy <- c()
for (i in 1:n_samples){
  sp1_idx <- which(rownames(dmat)==sp_sample[i])
  sp2_idx <- which(rownames(dmat)==sp_sample[n_samples + i])
  dsamples <- append(dsamples, dmat[sp1_idx, sp2_idx])
  
  sp1_idx <- which(rownames(dmat_copy)==sp_sample[i])
  sp2_idx <- which(rownames(dmat_copy)==sp_sample[n_samples + i])
  dsamples_copy <- append(dsamples_copy, dmat_copy[sp1_idx, sp2_idx])
}
identical(dsamples, dsamples_copy)
