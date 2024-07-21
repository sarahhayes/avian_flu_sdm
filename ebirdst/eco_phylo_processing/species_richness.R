# Create blank rasters for each of our quantities
# We will create rasters corresponding to the following traits:
# Congregatory behaviour
# Migatory behaviour
# Foraging around water surface
# Foraging >5cm below water surface
# Phylogenetic distance to a confirmed host

if (QUARTERLY){
  nlyrs <- 4
  qrtr_bds <- c(0, 13, 26, 39, 52)
}else{
  nlyrs <- 52
}

pop_rast <- rast(nlyrs=nlyrs,
                  crs=crs,
                  extent=euro_ext,
                  res=res(blank_3035))

species_rast <- rast(nlyrs=nlyrs,
                 crs=crs,
                 extent=euro_ext,
                 res=res(blank_3035))


for (i in 1:nlyrs){
  pop_rast[[i]] <- 0
  species_rast <-0
}




# Loop over other first no_species species:
{
  loop.start <- Sys.time()
  mean_dl_time <- 0
  mean_load_time <- 0
  mean_process_time <- 0
  no_processed <- 0
  for (i in 1:no_species) {
#   for (i in 1:10) {
#    i = 1
  idx <- sample_idx[i]
    species_sel <- sp_df$species_code[idx]
    species_factors <- sp_df[idx, ]
    
    if (!(species_sel %in% starting_dls)){
      dl_start <- Sys.time()
      dl_flag <- TRUE
      attempt_count <- 0
      while (dl_flag){
        attempt_count <- attempt_count + 1
        cat("attempt = ", attempt_count, ".\n")
        path <- try(ebirdst_download(species = species_sel,
                                     pattern = "_lr_"))
        if (!inherits(path, "try-error")){
          dl_flag <- FALSE
        }
        if (attempt_count>50){
          print("Download failed")
          dl_flag <- FALSE
        }
      }
      no_downloaded <- no_downloaded + 1
      time.now <- Sys.time()
      elapsed <- as.numeric(difftime(time.now, dl_start, units = "mins"))
      mean_dl_time <- (1 / no_downloaded) * 
        ((no_downloaded - 1) * mean_dl_time + elapsed)
    }else{
      load_start <- Sys.time()
      path <- ebirdst_download(species = species_sel,
                               pattern = "_lr_")
      no_loaded <- no_loaded + 1
      time.now <- Sys.time()
      elapsed <- as.numeric(difftime(time.now, load_start, units = "mins"))
      mean_load_time <- (1 / no_loaded) * 
        ((no_loaded - 1) * mean_load_time + elapsed)
    }
    
    process_start <- Sys.time()
    
    this_rast <- load_raster(path = path,
                             product = "percent-population",
                             period = "weekly",
                             resolution = "lr")
    this_rast <- project(x = this_rast, y = blank_3035, method = "near")
    
    # Get rid of NA's:
    this_rast <- replace(this_rast, is.na(this_rast), 0)
    
    if (QUARTERLY){
      this_rast <- lapply(1:nlyrs,
                          FUN = function(i){
                            app(this_rast[[qrtr_bds[i]:qrtr_bds[i+1]]], mean)}
      ) %>%
        rast
      set.names(this_rast, c("Qrt1", "Qrt2", "Qrt3", "Qrt4"))
    }
    
    # Fill in each field
    pop_rast <- species_factors$pop_sizes * this_rast
    
      species_rast <- species_rast + (pop_rast >= 1) 
  
    no_processed <- no_processed + 1
    
    time.now <- Sys.time()
    elapsed <- as.numeric(difftime(time.now, process_start, units = "mins"))
    mean_process_time <- (1 / no_processed) * 
      ((no_processed - 1) * mean_process_time + elapsed)
    time_remaining <- (total_to_download - no_downloaded) * mean_dl_time +
      (total_to_load - no_loaded) * mean_load_time +
      (no_species - i) * mean_process_time
    cat(as.numeric(difftime(time.now, loop.start, units="mins")),
        " minutes elapsed since start, estimated ",
        time_remaining,
        " remaining.",
        i,
        "of",
        no_species,
        "species processed.\n")
   }
}  

dev.off()
plot(species_rast)


if (SAVE_RAST){
  writeRaster(species_rast, "output/species_richness_rast.tif", overwrite=TRUE)
}
