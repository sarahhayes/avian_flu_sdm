 # In this script we train a BART model to classify sites for avian flu
# presence/absence based on environmental and species abundance factors.
# This script includes cross-validation of BART parameters.

SAVE_FITS <- FALSE
SAVE_PLOTS <- FALSE
BUILD_COVS <- FALSE

# Set path to folder containing data, and where output will be stored
PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

################################################################################
# Snippet adapted from
# https://stackoverflow.com/questions/14334970/convert-latitude-and-longitude-coordinates-to-country-name-in-r
library(sp)
library(rworldmap)

# The single argument to this function, points, is a data.frame in which:
#   - column 1 contains the longitude in degrees
#   - column 2 contains the latitude in degrees
coords2country = function(sample_data)
{  
  # Convert from EPSG to lat long decimal:
  coordinates(sample_data) <- c("X", "Y")
  proj4string(sample_data) <- CRS("epsg:3035")
  sample_data <- spTransform(sample_data, CRS("+init=epsg:4326")) %>% as.data.frame()
  
  # Extract coordinates as dataframe with appropriate labels:
  points <- data.frame(lon = sample_data$coords.x1, lat = sample_data$coords.x2)
  countriesSP <- getMap(resolution='low')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # convert our list of points to a SpatialPoints object
  
  # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  #setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  # return the ADMIN names of each country
  indices$ADMIN  
  #indices$ISO3 # returns the ISO3 code 
  #indices$continent   # returns the continent (6 continent model)
  #indices$REGION   # returns the continent (7 continent model)
}

################################################################################

# EMBARCADERO FUNCTIONS ADJUSTED TO ALLOW CUSTOM BART PARAMETERS:

bart.flex <- function(x.data, y.data, ri.data = NULL,
                      y.name = NULL, ri.name = NULL,
                      n.trees = 200,
                      k = 2.0, power = 2.0, base = 0.95) {
  
  if(is.null(ri.data)) {
    train <- cbind(y.data, x.data) 
    if(!is.null(y.name)) {colnames(train)[1] <- y.name}
    train <- na.omit(train)
    model <- bart(y.train = train[,1], 
                  x.train = train[,2:ncol(train)], 
                  ntree = n.trees, keeptrees=TRUE)
  } else { 
    train <- cbind(y.data, x.data, ri.data) 
    if(!is.null(y.name)) {colnames(train)[1] <- y.name}
    if(!is.null(ri.name)) {colnames(train)[ncol(train)] <- ri.name}
    f <- as.formula(paste(paste(colnames(train)[1],paste(colnames(train)[2:(ncol(train)-1)], 
                                                         collapse=' + '), sep = ' ~ '), 
                          colnames(train)[ncol(train)], sep=' - '))
    
    train <- na.omit(train) 
    model <- rbart_vi(f, 
                      group.by=train[,ncol(train)],
                      data=train,
                      n.samples=1000,
                      n.burn=100,
                      n.chains=1,
                      n.threads=1,
                      n.trees = n.trees,
                      keepTrees = TRUE,
                      k = k,
                      power = power,
                      base = base) 
  }
  return(model)
}

varimp.diag <- function(x.data, y.data, ri.data=NULL, iter=50, quiet=FALSE,
                        k = 2.0, power = 2.0, base = 0.95) {
  
  nvars <- ncol(x.data)
  varnums <- c(1:nvars)
  varlist <- colnames(x.data)
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY 4 THIS CODE :) 
  
  ###############
  
  # auto-drops 
  
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               ri.data = ri.data,
                               n.trees = 200,
                               k = k,
                               power = power,
                               base = base))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  
  if(length(dropnames)==0) {} else{
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames,collapse = ' '), ' \n')
  }
  
  x.data %>% dplyr::select(-dropnames) -> x.data  
  
  ###############
  
  for (n.trees in c(10, 20, 50, 100, 150, 200)) {
    
    cat(paste('\n', n.trees, 'tree models:', iter, 'iterations\n'))
    if(!quiet){pb <- txtProgressBar(min = 0, max = iter, style = 3)}
    
    for(index in 1:iter) {
      quietly(model.j <- bart.flex(x.data = x.data[,varnums], y.data = y.data, 
                                   ri.data = ri.data,
                                   n.trees = n.trees,
                                   k = k,
                                   power = power,
                                   base = base))
      
      vi.j <- varimp(model.j)
      if(index==1) {
        vi.j.df <- vi.j
      } else {
        vi.j.df[,index+1] <- vi.j[,2]
      }
      if(!quiet){setTxtProgressBar(pb, index)}
    }
    vi.j <- data.frame(vi.j.df[,1],
                       rowMeans(vi.j.df[,-1]))
    
    if(n.trees==10) { vi <- vi.j } else {  vi <- cbind(vi,vi.j[,2])  }
  }
  
  colnames(vi) <- c('variable','10','20','50','100','150','200')
  vi <- reshape::melt(vi, "variable")
  colnames(vi) <- c('variable','trees','imp')
  
  vi %>% group_by(variable) %>% summarise(max = max(imp)) -> vi.fac
  vi.fac <- vi.fac[order(-vi.fac$max),]
  
  vi$names <- factor(vi$variable, levels=vi.fac$variable)
  
  g1 <- ggplot2::ggplot(vi, aes(y=imp, x=names, group=trees)) +
    geom_line(aes(color=trees)) + geom_point(size=3) + theme_classic() +
    ylab("Relative contribution\n") + xlab("\nVariables dropped") +
    ggpubr::rotate_x_text(angle = 35) + 
    theme(axis.text = element_text(size=10),
          axis.title = element_text(size=14,face="bold")); print(g1)
  
}

variable.step <- function(x.data, y.data, ri.data=NULL, n.trees=10, iter=50, quiet=FALSE,
                          k = 2.0, power = 2.0, base = 0.95) {
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY
  
  comp <- complete.cases(x.data)
  
  if(length(comp) < (nrow(x.data))) {
    message("Some rows with NA's have been automatically dropped. \n")
  }
  x.data <- x.data[comp,]
  y.data <- y.data[comp]
  
  ###############
  
  # auto-drops 
  
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               ri.data = ri.data,
                               n.trees = 200,
                               k = k,
                               power = power,
                               base = base))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  
  if(length(dropnames) > 0) {
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames,collapse = ' '), ' \n')
  }
  
  x.data %>% dplyr::select(-any_of(dropnames)) -> x.data  
  
  ###############
  
  nvars <- ncol(x.data)
  varnums <- c(1:nvars)
  varlist.orig <- varlist <- colnames(x.data)
  
  rmses <- data.frame(Variable.number=c(),RMSE=c())
  dropped.varlist <- c()
  
  for(var.j in c(nvars:3)) {
    
    print(noquote(paste("Number of variables included:",var.j)))
    print(noquote("Dropped:"))
    print(if(length(dropped.varlist)==0) {noquote("")} else {noquote(dropped.varlist)})
    
    rmse.list <- c()
    
    if(!quiet){pb <- txtProgressBar(min = 0, max = iter, style = 3)}
    for(index in 1:iter) {
      quietly(model.j <- bart.flex(x.data = x.data[,varnums], y.data = y.data, 
                                   ri.data = ri.data,
                                   n.trees = n.trees,
                                   k = k,
                                   power = power,
                                   base = base))
      
      quietly(vi.j <- varimp(model.j))
      if(index==1) {
        vi.j.df <- vi.j
      } else {
        vi.j.df[,index+1] <- vi.j[,2]
      }
      
      pred.p <- colMeans(pnorm(model.j$yhat.train))[y.data==1]
      pred.a <- colMeans(pnorm(model.j$yhat.train))[y.data==0]
      #e <- evaluate(p=pred.p,
      #              a=pred.a)
      #aucs <- rbind(aucs,c(var.j,e@auc)); colnames(aucs) <- c('Vars','AUC')
      
      pred.c <- c(pred.p, pred.a)
      true.c <- c(rep(1,length(pred.p)), rep(0,length(pred.a)))
      rmsej.i <- Metrics::rmse(true.c,pred.c)
      rmse.list <- c(rmse.list,rmsej.i)
      if(!quiet){setTxtProgressBar(pb, index)}
    }
    
    vi.j <- data.frame(vi.j.df[,1],
                       rowMeans(vi.j.df[,-1]))
    vi.j <- vi.j[order(vi.j[,2]),]
    
    drop.var <- vi.j[1,1]
    dropped.varlist <- c(dropped.varlist,as.character(drop.var))
    
    rmsej <- mean(rmse.list)
    
    rmses <- rbind(rmses,c(nvars-var.j,rmsej)); colnames(rmses) <- c('VarsDropped','RMSE')
    
    varnums <- varnums[!(varnums==which(varlist.orig==drop.var))]
    varlist <- varlist.orig[varnums]
    print(noquote("---------------------------------------"))
  }
  
  g1 <- ggplot2::ggplot(rmses, aes(y=RMSE, x=VarsDropped)) +
    geom_line(color="black") + geom_point(size=3) + theme_bw() +
    ylab("RMSE of model\n") + xlab("\nVariables dropped") +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14,face="bold")) +
    scale_x_discrete(limits=c(0:(nrow(rmses)))); print(g1)
  
  print(noquote("---------------------------------------"))
  print(noquote("Final recommended variable list"))
  varlist.final <- varlist.orig[!(varlist.orig %in% dropped.varlist[0:(which(rmses$RMSE==min(rmses$RMSE))-1)])]
  print(noquote(varlist.final))
  invisible(varlist.final)
}

bart.step <- function(x.data, y.data, ri.data=NULL,
                      iter.step=100, tree.step=10,
                      iter.plot=100,
                      k = 2.0, power = 2.0, base = 0.95,
                      full=FALSE,
                      quiet=FALSE) {
  
  ###############
  
  # auto-drops 
  
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }  # THANKS HADLEY
  
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               ri.data = ri.data,
                               n.trees = 200,
                               k = k,
                               power = power,
                               base = base))
  
  if(class(model.0)=='rbart') {
    fitobj <- model.0$fit[[1]]
  }
  if(class(model.0)=='bart') {
    fitobj <- model.0$fit
  }
  
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
  
  if(length(dropnames)==0) {} else{
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames,collapse = ' '), ' \n')
  }
  
  x.data %>% dplyr::select(-dropnames) -> x.data  
  
  ###############
  
  quiet2 <- quiet
  if(full==TRUE){varimp.diag(x.data, y.data, ri.data, iter=iter.plot, quiet=quiet2)}
  vs <- variable.step(x.data,
                      y.data,
                      ri.data,
                      n.trees=tree.step,
                      iter=iter.step,
                      quiet=quiet2)
  
  invisible(best.model <- bart.flex(x.data = x.data[,vs], y.data = y.data, 
                                    ri.data = ri.data, n.trees=200,
                                    k = k,
                                    power = power,
                                    base = base))
  if(full==TRUE){varimp(best.model, plots=TRUE)}
  if(full==TRUE) {p <- summary(best.model, plots=TRUE)
  print(p)} else 
  {p <- summary(best.model, plots=FALSE)
  print(p)}
  invisible(best.model)
}

if (BUILD_COVS){
  # set the crs we want to use
  crs <- "epsg:3035"
  
  euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later
  
  eco_paths <- list.files(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Eco-Rasters", sep = ""),
                          pattern = "*.tif",
                          full.names = TRUE)
  eco_lyrnames <- eco_paths %>%
    sub(pattern = paste(PATH_TO_DATA, "AI_S2_SDM_storage/Eco-Rasters/", sep = ""),
        replacement = "") %>%
    sub(pattern = "_rast.tif",
        replacement = "")
  
  eco_layers <- rast(lapply(eco_paths, rast))
  
  # Separate by quarters:
  q1_eco_layers <- eco_layers[[seq(1, nlyr(eco_layers), 4)]]
  q2_eco_layers <- eco_layers[[seq(2, nlyr(eco_layers), 4)]]
  q3_eco_layers <- eco_layers[[seq(3, nlyr(eco_layers), 4)]]
  q4_eco_layers <- eco_layers[[seq(4, nlyr(eco_layers), 4)]]
  set.names(q1_eco_layers, eco_lyrnames)
  set.names(q2_eco_layers, eco_lyrnames)
  set.names(q3_eco_layers, eco_lyrnames)
  set.names(q4_eco_layers, eco_lyrnames)
  
  # Now do environmental:
  env_paths <- list.files(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental rasters", sep = ""),
                          pattern = "*.tif",
                          full.names = TRUE)
  env_paths <- env_paths[-grep("landcover_output_full", env_paths)]
  env_paths <- env_paths[-grep("chicken_density_2015", env_paths)]
  env_paths <- env_paths[-grep("duck_density_2015", env_paths)]
  env_paths <- env_paths[-grep("mean_tmax", env_paths)]
  env_paths <- env_paths[-grep("mean_tmin", env_paths)]
  env_lyrnames <- env_paths %>%
    sub(pattern = paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental rasters/", sep = ""),
        replacement = "") %>%
    sub(pattern = ".tif",
        replacement = "")
  
  env_layers <- rast(lapply(env_paths, rast))
  names(env_layers) <- env_lyrnames
  
  all_excludes <- grep("quart", env_lyrnames, value = TRUE)
  q1_excludes <- grep("first", all_excludes, value = TRUE, invert = TRUE)
  q2_excludes <- grep("second", all_excludes, value = TRUE, invert = TRUE)
  q3_excludes <- grep("third", all_excludes, value = TRUE, invert = TRUE)
  q4_excludes <- grep("fourth", all_excludes, value = TRUE, invert = TRUE)
  
  q1_env_layers <- subset(env_layers, setdiff(env_lyrnames, q1_excludes))
  q2_env_layers <- subset(env_layers, setdiff(env_lyrnames, q2_excludes))
  q3_env_layers <- subset(env_layers, setdiff(env_lyrnames, q3_excludes))
  q4_env_layers <- subset(env_layers, setdiff(env_lyrnames, q4_excludes))
  
  cov_coast <- read.csv(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental variable csvs/dist_to_coast_output.csv", sep = "")) %>%
    rename(x = X, y = Y) %>%
    select(-X.1) %>% 
    relocate(x,y) %>%
    rast(type = "xyz", crs = "epsg:3035")
  cov_water <- read.csv(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental variable csvs/dist_to_water_output.csv", sep = "")) %>%
    select(-ID) %>% 
    relocate(x,y) %>%
    rast(type = "xyz", crs = "epsg:3035")
  cov_alt <- read.csv(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental variable csvs/elevation_outputs.csv", sep = "")) %>%
    select(-ID, -X) %>% 
    relocate(x,y) %>%
    rast(type = "xyz", crs = "epsg:3035")
  
  landcover_rast <- rast(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental rasters/landcover_output_full.tif", sep = ""))
  n_landtypes <- 17
  landcover_layers <- lapply(1:n_landtypes,
                             FUN = function(i){landcover_rast == i}) %>%
    rast()
  names(landcover_layers) <- lapply(1:n_landtypes,
                                    FUN = function(i){paste("lc_", i, sep = "")})
  lc_lookup <- pairlist("lc_1" = "Water_bodies",
                        "lc_2" = "Evergreen_Needleleaf_Forests",
                        "lc_3" = "Evergreen_Broadleaf_Forests",
                        "lc_4" = "Deciduous_Needleleaf_Forests",
                        "lc_5" = "Deciduous_Broadleaf_Forests",
                        "lc_6" = "Mixed_Forests",
                        "lc_7" = "Closed_Shrublands",
                        "lc_8" = "Open_Shrublands",
                        "lc_9" = "Woody_Savannas",
                        "lc_10" = "Savannas",
                        "lc_11" = "Grasslands",
                        "lc_12" = "Permanent_Wetlands",
                        "lc_13" = "Croplands",
                        "lc_14" = "Urban_and_Built-up_Lands",
                        "lc_15" = "Cropland/Natural_Vegetation_Mosaics",
                        "lc_16" = "Non-Vegetated_Lands",
                        "lc_17" = "Unclassified")
  names(landcover_layers) <- lapply(names(landcover_layers),
                                    FUN = function(n){lc_lookup[[n]]})
  # Remove layers that are either all True or all False:
  hom_layers <- sapply(1:n_landtypes,
                       FUN = function(i){diff(minmax(landcover_layers[[i]]))==0}) %>%
    which()
  if (length(hom_layers)>0){
    landcover_layers <- landcover_layers[[-hom_layers]]
  }
  
  q1_covs <- c(resample(q1_eco_layers, env_layers, method = "near"),
                         q1_env_layers,
                         landcover_layers,
                         cov_coast,
                         cov_water,
                         cov_alt)  
  writeRaster(q1_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs.tif", sep = ""), overwrite = TRUE)
  rm(q1_covs)
  gc()
  q2_covs <- c(resample(q2_eco_layers, env_layers, method = "near"),
                         q2_env_layers,
                         landcover_layers,
                         cov_coast,
                         cov_water,
                         cov_alt)
  writeRaster(q2_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs.tif", sep = ""), overwrite = TRUE)
  rm(q2_covs)
  gc()
  q3_covs <- c(resample(q3_eco_layers, env_layers, method = "near"),
                        q3_env_layers,
                        landcover_layers,
                        cov_coast,
                        cov_water,
                        cov_alt)
  writeRaster(q3_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs.tif", sep = ""), overwrite = TRUE)
  rm(q3_covs)
  gc()
  q4_covs <- c(resample(q4_eco_layers, env_layers, method = "near"),
                         q4_env_layers,
                         landcover_layers,
                         cov_coast,
                         cov_water,
                         cov_alt)
  writeRaster(q4_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs.tif", sep = ""), overwrite = TRUE)
  rm(q4_covs)
  gc()

}

################################################################################
# Do the Q1 analysis

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

# Load training data
training_coords <- readRDS("training_sets/training_coords_Q1.RDS")
sample_countries <- coords2country(training_coords) %>% as.data.frame()
cov_df <- data.frame(raster::extract(covstack, training_coords))

n_pts <- nrow(training_coords)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(cov_df), n_training)
xtrain <- cov_df[training, ]
ytrain <- training_coords$pos[training]
xtest <- cov_df[-training, ]
ytest <- training_coords$pos[-training]

if (SAVE_FITS){
  save(xtrain, ytrain, xtest, ytest, file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/q1_train_test_data.rds", sep = ""))
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_Q1.rds", sep = ""))
}

retuned_model <- retune(basic_model)

k_retune = retuned_model$fit$model@node.hyperprior@k
power_retune = retuned_model$fit$model@tree.prior@power
base_retune = retuned_model$fit$model@tree.prior@base

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = k_retune,
                 power = power_retune,
                 base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_sdm_Q1.rds", sep = ""))
}

covstack_lores <- aggregate(covstack, fact = 10)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
                      )
if (SAVE_FITS){
  save(pred_layer,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_prediction_Q1.rds", sep = ""))
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q1_map.png")
}
# Plot map
plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Mean prediction, Q1')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q1_lbound.png")
}
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q1')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q1_ubound.png")
}
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q1')
if (SAVE_PLOTS){
  dev.off()
}

################################################################################
# Do the Q2 analysis

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")


# Load training data
training_coords <- readRDS("training_sets/training_coords_Q2.RDS")
cov_df <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

n_pts <- nrow(training_coords)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(cov_df), n_training)
xtrain <- cov_df[training, ]
ytrain <- training_coords$pos[training]
xtest <- cov_df[-training, ]
ytest <- training_coords$pos[-training]

if (SAVE_FITS){
  save(xtrain, ytrain, xtest, ytest, file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/q2_train_test_data.rds", sep = ""))
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_Q2.rds", sep = ""))
}

retuned_model <- retune(basic_model)

k_retune = retuned_model$fit$model@node.hyperprior@k
power_retune = retuned_model$fit$model@tree.prior@power
base_retune = retuned_model$fit$model@tree.prior@base

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = k_retune,
                 power = power_retune,
                 base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_sdm_Q2.rds", sep = ""))
}

covstack_lores <- aggregate(covstack, fact = 10)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
)
if (SAVE_FITS){
  save(pred_layer,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_prediction_Q2.rds", sep = ""))
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q2_map.png")
}
# Plot map
plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Mean prediction, Q2')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q2_lbound.png")
}
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q2')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q2_ubound.png")
}
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q2')
if (SAVE_PLOTS){
  dev.off()
}

################################################################################
# Do the Q3 analysis

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")


# Load training data
training_coords <- readRDS("training_sets/training_coords_Q3.RDS")
cov_df <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

n_pts <- nrow(training_coords)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(cov_df), n_training)
xtrain <- cov_df[training, ]
ytrain <- training_coords$pos[training]
xtest <- cov_df[-training, ]
ytest <- training_coords$pos[-training]

if (SAVE_FITS){
  save(xtrain, ytrain, xtest, ytest, file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/q3_train_test_data.rds", sep = ""))
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_Q3.rds", sep = ""))
}

retuned_model <- retune(basic_model)

k_retune = retuned_model$fit$model@node.hyperprior@k
power_retune = retuned_model$fit$model@tree.prior@power
base_retune = retuned_model$fit$model@tree.prior@base

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = k_retune,
                 power = power_retune,
                 base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_sdm_Q3.rds", sep = ""))
}

covstack_lores <- aggregate(covstack, fact = 10)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
)
if (SAVE_FITS){
  save(pred_layer,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_prediction_Q3.rds", sep = ""))
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q3_map.png")
}
# Plot map
plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Mean prediction, Q3')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q3_lbound.png")
}
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q3')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q3_ubound.png")
}
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q3')
if (SAVE_PLOTS){
  dev.off()
}

################################################################################
# Do the Q4 analysis

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")


# Load training data
training_coords <- readRDS("training_sets/training_coords_Q4.RDS")
cov_df <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

n_pts <- nrow(training_coords)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(cov_df), n_training)
xtrain <- cov_df[training, ]
ytrain <- training_coords$pos[training]
xtest <- cov_df[-training, ]
ytest <- training_coords$pos[-training]

if (SAVE_FITS){
  save(xtrain, ytrain, xtest, ytest, file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/q4_train_test_data.rds", sep = ""))
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_Q4.rds", sep = ""))
}

retuned_model <- retune(basic_model)

k_retune = retuned_model$fit$model@node.hyperprior@k
power_retune = retuned_model$fit$model@tree.prior@power
base_retune = retuned_model$fit$model@tree.prior@base

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = k_retune,
                 power = power_retune,
                 base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_sdm_Q4.rds", sep = ""))
}

covstack_lores <- aggregate(covstack, fact = 10)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
)
if (SAVE_FITS){
  save(pred_layer,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_prediction_Q4.rds", sep = ""))
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q4_map.png")
}
# Plot map
plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Mean prediction, Q4')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q4_lbound.png")
}
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q4')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q4_ubound.png")
}
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q4')
if (SAVE_PLOTS){
  dev.off()
}