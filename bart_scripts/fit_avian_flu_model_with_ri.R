 # In this script we train a BART model to classify sites for avian flu
# presence/absence based on environmental and species abundance factors.
# This script includes cross-validation of BART parameters.

SAVE_FITS <- FALSE
SAVE_PLOTS <- FALSE
BUILD_COVS <- FALSE
PLOT_COUNTRY_VALIDATION <- TRUE

# Set path to folder containing data, and where output will be stored
PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(dbarts)
library(rasterVis)
library(embarcadero)
library(raster)
library(sp)
library(terra)
library(rworldmap)
set.seed(12345)

################################################################################
# Try making raster of country idx's, adapted from
# https://gis.stackexchange.com/questions/44139/convert-a-spatialpolygonsdataframe-to-raster-using-rasterize-function


# Create a blank raster with appropriate projection and extent
blank_3035 <- terra::rast("output/euro_rast_10k.tif")
euro_ext <- ext(blank_3035)
crs <- "epsg:3035"

countriesSP <- getMap(resolution='low') %>% spTransform(CRS(crs))
countriesSP$Grd_ranks <- rank(countriesSP$ADMIN)
country_lookup <- data.frame(country=countriesSP$ADMIN,
                             val=rank(countriesSP$Grd_ranks))

r <- rast(nlyrs=1,
          crs=crs,
          extent=extent(countriesSP),
          res=res(blank_3035)) %>% raster()

country_rast <- rasterize(countriesSP, r, field="Grd_ranks", fun="first") %>%
                terra::rast() %>%
                crop(euro_ext)
# country_fact <- as.factor(country_rast)
# tar<-levels(country_fact)[[1]]
# tar[["Country"]]<-country_lookup$country[as.numeric(freq(country_fact)$value)]
# levels(country_fact)<-tar
# levelplot(country_fact)

# ################################################################################
# # Snippet adapted from
# # https://stackoverflow.com/questions/14334970/convert-latitude-and-longitude-coordinates-to-country-name-in-r
# 
# 
# # The single argument to this function, points, is a data.frame in which:
# #   - column 1 contains the longitude in degrees
# #   - column 2 contains the latitude in degrees
# coords2country = function(sample_data)
# {  
#   # Convert from EPSG to lat long decimal:
#   coordinates(sample_data) <- c("X", "Y")
#   proj4string(sample_data) <- CRS("epsg:3035")
#   sample_data <- spTransform(sample_data, CRS("+init=epsg:4326")) %>% as.data.frame()
#   
#   # Extract coordinates as dataframe with appropriate labels:
#   points <- data.frame(lon = sample_data$coords.x1, lat = sample_data$coords.x2)
#   countriesSP <- getMap(resolution='high')
#   #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
#   
#   # convert our list of points to a SpatialPoints object
#   
#   # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
#   
#   #setting CRS directly to that from rworldmap
#   pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))
#   
#   
#   # use 'over' to get indices of the Polygons object containing each point 
#   indices = over(pointsSP, countriesSP)
#   
#   # return the ADMIN names of each country
#   indices$ADMIN  
#   #indices$ISO3 # returns the ISO3 code 
#   #indices$continent   # returns the continent (6 continent model)
#   #indices$REGION   # returns the continent (7 continent model)
# }

################################################################################

# EMBARCADERO FUNCTIONS ADJUSTED TO ALLOW CUSTOM BART PARAMETERS:
# 
# bart.flex <- function(x.data, y.data, ri.data = NULL,
#                       y.name = NULL, ri.name = NULL,
#                       n.trees = 200,
#                       k = 2.0, power = 2.0, base = 0.95) {
#   
#   if(is.null(ri.data)) {
#     train <- cbind(y.data, x.data) 
#     if(!is.null(y.name)) {colnames(train)[1] <- y.name}
#     train <- na.omit(train)
#     model <- bart(y.train = train[,1], 
#                   x.train = train[,2:ncol(train)], 
#                   ntree = n.trees, keeptrees=TRUE)
#   } else { 
#     train <- cbind(y.data, x.data, ri.data) 
#     if(!is.null(y.name)) {colnames(train)[1] <- y.name}
#     if(!is.null(ri.name)) {colnames(train)[ncol(train)] <- ri.name}
#     f <- as.formula(paste(paste(colnames(train)[1],paste(colnames(train)[2:(ncol(train)-1)], 
#                                                          collapse=' + '), sep = ' ~ '), 
#                           colnames(train)[ncol(train)], sep=' - '))
#     
#     train <- na.omit(train) 
#     model <- rbart_vi(f, 
#                       group.by=train[,ncol(train)],
#                       data=train,
#                       n.samples=1000,
#                       n.burn=100,
#                       n.chains=1,
#                       n.threads=1,
#                       n.trees = n.trees,
#                       keepTrees = TRUE,
#                       k = k,
#                       power = power,
#                       base = base) 
#   }
#   return(model)
# }
# 
# varimp.diag <- function(x.data, y.data, ri.data=NULL, iter=50, quiet=FALSE,
#                         k = 2.0, power = 2.0, base = 0.95) {
#   
#   nvars <- ncol(x.data)
#   varnums <- c(1:nvars)
#   varlist <- colnames(x.data)
#   
#   quietly <- function(x) {
#     sink(tempfile())
#     on.exit(sink())
#     invisible(force(x))
#   }  # THANKS HADLEY 4 THIS CODE :) 
#   
#   ###############
#   
#   # auto-drops 
#   
#   quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
#                                ri.data = ri.data,
#                                n.trees = 200,
#                                k = k,
#                                power = power,
#                                base = base))
#   
#   if(class(model.0)=='rbart') {
#     fitobj <- model.0$fit[[1]]
#   }
#   if(class(model.0)=='bart') {
#     fitobj <- model.0$fit
#   }
#   dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
#   
#   if(length(dropnames)==0) {} else{
#     message("Some of your variables have been automatically dropped by dbarts.")
#     message("(This could be because they're characters, homogenous, etc.)")
#     message("It is strongly recommended that you remove these from the raw data:")
#     message(paste(dropnames,collapse = ' '), ' \n')
#   }
#   
#   x.data %>% dplyr::select(-dropnames) -> x.data  
#   
#   ###############
#   
#   for (n.trees in c(10, 20, 50, 100, 150, 200)) {
#     
#     cat(paste('\n', n.trees, 'tree models:', iter, 'iterations\n'))
#     if(!quiet){pb <- txtProgressBar(min = 0, max = iter, style = 3)}
#     
#     for(index in 1:iter) {
#       quietly(model.j <- bart.flex(x.data = x.data[,varnums], y.data = y.data, 
#                                    ri.data = ri.data,
#                                    n.trees = n.trees,
#                                    k = k,
#                                    power = power,
#                                    base = base))
#       
#       vi.j <- varimp(model.j)
#       if(index==1) {
#         vi.j.df <- vi.j
#       } else {
#         vi.j.df[,index+1] <- vi.j[,2]
#       }
#       if(!quiet){setTxtProgressBar(pb, index)}
#     }
#     vi.j <- data.frame(vi.j.df[,1],
#                        rowMeans(vi.j.df[,-1]))
#     
#     if(n.trees==10) { vi <- vi.j } else {  vi <- cbind(vi,vi.j[,2])  }
#   }
#   
#   colnames(vi) <- c('variable','10','20','50','100','150','200')
#   vi <- reshape::melt(vi, "variable")
#   colnames(vi) <- c('variable','trees','imp')
#   
#   vi %>% group_by(variable) %>% summarise(max = max(imp)) -> vi.fac
#   vi.fac <- vi.fac[order(-vi.fac$max),]
#   
#   vi$names <- factor(vi$variable, levels=vi.fac$variable)
#   
#   g1 <- ggplot2::ggplot(vi, aes(y=imp, x=names, group=trees)) +
#     geom_line(aes(color=trees)) + geom_point(size=3) + theme_classic() +
#     ylab("Relative contribution\n") + xlab("\nVariables dropped") +
#     ggpubr::rotate_x_text(angle = 35) + 
#     theme(axis.text = element_text(size=10),
#           axis.title = element_text(size=14,face="bold")); print(g1)
#   
# }
# 
# variable.step <- function(x.data, y.data, ri.data=NULL, n.trees=10, iter=50, quiet=FALSE,
#                           k = 2.0, power = 2.0, base = 0.95) {
#   
#   quietly <- function(x) {
#     sink(tempfile())
#     on.exit(sink())
#     invisible(force(x))
#   }  # THANKS HADLEY
#   
#   comp <- complete.cases(x.data)
#   
#   if(length(comp) < (nrow(x.data))) {
#     message("Some rows with NA's have been automatically dropped. \n")
#   }
#   x.data <- x.data[comp,]
#   y.data <- y.data[comp]
#   
#   ###############
#   
#   # auto-drops 
#   
#   quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
#                                ri.data = ri.data,
#                                n.trees = 200,
#                                k = k,
#                                power = power,
#                                base = base))
#   
#   if(class(model.0)=='rbart') {
#     fitobj <- model.0$fit[[1]]
#   }
#   if(class(model.0)=='bart') {
#     fitobj <- model.0$fit
#   }
#   
#   dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
#   
#   if(length(dropnames) > 0) {
#     message("Some of your variables have been automatically dropped by dbarts.")
#     message("(This could be because they're characters, homogenous, etc.)")
#     message("It is strongly recommended that you remove these from the raw data:")
#     message(paste(dropnames,collapse = ' '), ' \n')
#   }
#   
#   x.data %>% dplyr::select(-any_of(dropnames)) -> x.data  
#   
#   ###############
#   
#   nvars <- ncol(x.data)
#   varnums <- c(1:nvars)
#   varlist.orig <- varlist <- colnames(x.data)
#   
#   rmses <- data.frame(Variable.number=c(),RMSE=c())
#   dropped.varlist <- c()
#   
#   for(var.j in c(nvars:3)) {
#     
#     print(noquote(paste("Number of variables included:",var.j)))
#     print(noquote("Dropped:"))
#     print(if(length(dropped.varlist)==0) {noquote("")} else {noquote(dropped.varlist)})
#     
#     rmse.list <- c()
#     
#     if(!quiet){pb <- txtProgressBar(min = 0, max = iter, style = 3)}
#     for(index in 1:iter) {
#       quietly(model.j <- bart.flex(x.data = x.data[,varnums], y.data = y.data, 
#                                    ri.data = ri.data,
#                                    n.trees = n.trees,
#                                    k = k,
#                                    power = power,
#                                    base = base))
#       
#       quietly(vi.j <- varimp(model.j))
#       if(index==1) {
#         vi.j.df <- vi.j
#       } else {
#         vi.j.df[,index+1] <- vi.j[,2]
#       }
#       
#       pred.p <- colMeans(pnorm(model.j$yhat.train))[y.data==1]
#       pred.a <- colMeans(pnorm(model.j$yhat.train))[y.data==0]
#       #e <- evaluate(p=pred.p,
#       #              a=pred.a)
#       #aucs <- rbind(aucs,c(var.j,e@auc)); colnames(aucs) <- c('Vars','AUC')
#       
#       pred.c <- c(pred.p, pred.a)
#       true.c <- c(rep(1,length(pred.p)), rep(0,length(pred.a)))
#       rmsej.i <- Metrics::rmse(true.c,pred.c)
#       rmse.list <- c(rmse.list,rmsej.i)
#       if(!quiet){setTxtProgressBar(pb, index)}
#     }
#     
#     vi.j <- data.frame(vi.j.df[,1],
#                        rowMeans(vi.j.df[,-1]))
#     vi.j <- vi.j[order(vi.j[,2]),]
#     
#     drop.var <- vi.j[1,1]
#     dropped.varlist <- c(dropped.varlist,as.character(drop.var))
#     
#     rmsej <- mean(rmse.list)
#     
#     rmses <- rbind(rmses,c(nvars-var.j,rmsej)); colnames(rmses) <- c('VarsDropped','RMSE')
#     
#     varnums <- varnums[!(varnums==which(varlist.orig==drop.var))]
#     varlist <- varlist.orig[varnums]
#     print(noquote("---------------------------------------"))
#   }
#   
#   g1 <- ggplot2::ggplot(rmses, aes(y=RMSE, x=VarsDropped)) +
#     geom_line(color="black") + geom_point(size=3) + theme_bw() +
#     ylab("RMSE of model\n") + xlab("\nVariables dropped") +
#     theme(axis.text = element_text(size=12),
#           axis.title = element_text(size=14,face="bold")) +
#     scale_x_discrete(limits=c(0:(nrow(rmses)))); print(g1)
#   
#   print(noquote("---------------------------------------"))
#   print(noquote("Final recommended variable list"))
#   varlist.final <- varlist.orig[!(varlist.orig %in% dropped.varlist[0:(which(rmses$RMSE==min(rmses$RMSE))-1)])]
#   print(noquote(varlist.final))
#   invisible(varlist.final)
# }
# 
# bart.step <- function(x.data, y.data, ri.data=NULL,
#                       iter.step=100, tree.step=10,
#                       iter.plot=100,
#                       k = 2.0, power = 2.0, base = 0.95,
#                       full=FALSE,
#                       quiet=FALSE) {
#   
#   ###############
#   
#   # auto-drops 
#   
#   quietly <- function(x) {
#     sink(tempfile())
#     on.exit(sink())
#     invisible(force(x))
#   }  # THANKS HADLEY
#   
#   quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
#                                ri.data = ri.data,
#                                n.trees = 200,
#                                k = k,
#                                power = power,
#                                base = base))
#   
#   if(class(model.0)=='rbart') {
#     fitobj <- model.0$fit[[1]]
#   }
#   if(class(model.0)=='bart') {
#     fitobj <- model.0$fit
#   }
#   
#   dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x,"drop"))==FALSE)))]
#   
#   if(length(dropnames)==0) {} else{
#     message("Some of your variables have been automatically dropped by dbarts.")
#     message("(This could be because they're characters, homogenous, etc.)")
#     message("It is strongly recommended that you remove these from the raw data:")
#     message(paste(dropnames,collapse = ' '), ' \n')
#   }
#   
#   x.data %>% dplyr::select(-dropnames) -> x.data  
#   
#   ###############
#   
#   quiet2 <- quiet
#   if(full==TRUE){varimp.diag(x.data, y.data, ri.data, iter=iter.plot, quiet=quiet2)}
#   vs <- variable.step(x.data,
#                       y.data,
#                       ri.data,
#                       n.trees=tree.step,
#                       iter=iter.step,
#                       quiet=quiet2)
#   
#   invisible(best.model <- bart.flex(x.data = x.data[,vs], y.data = y.data, 
#                                     ri.data = ri.data, n.trees=200,
#                                     k = k,
#                                     power = power,
#                                     base = base))
#   if(full==TRUE){varimp(best.model, plots=TRUE)}
#   if(full==TRUE) {p <- summary(best.model, plots=TRUE)
#   print(p)} else 
#   {p <- summary(best.model, plots=FALSE)
#   print(p)}
#   invisible(best.model)
# }

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
training_coords <- readRDS("training_sets/training_coords_A_Q1.RDS")
# training_coords$country <- coords2country(training_coords)
# 
# # Identify and plot samples where a country could not be assigned - should find
# # it's all on-water samples.
# 
country_df <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))
mismatch_samples <- which(is.na(rowSums(country_df)))
mismatch_pts <- terra::vect(training_coords[mismatch_samples, ], geom=c("X", "Y"),
                            crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
crs <- "epsg:3035"
euro_ext <- extent(covstack[[1]])

# # change projection and extent. 
# # using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs)
euro_map_crop <- terra::crop(euro_map, euro_ext)
plot(euro_map_crop,
     col = "white",
     background = "azure2",
     axes = FALSE,
     buffer = FALSE,
     mar = c(0, 0, 0, 0))
plot(mismatch_pts, add = T, col = "red", pch = 16, cex = 1.)
# 
# # Remove NA's from data - although we should think about how to fix this
training_coords <- training_coords[-mismatch_samples,]
country_df <- country_df[-mismatch_samples,]

cov_df <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

if (PLOT_COUNTRY_VALIDATION){
  # Quick digression to validate coords2country function:
  
  plot(euro_map_crop,
       col = "white",
       background = "azure2",
       axes = FALSE,
       buffer = FALSE,
       xmin = euro_ext@xmin,
       mar = c(0, 0, 0, 0))
  pts_to_plot <- sample(1:nrow(training_coords), 25)
  for (i in 1:25){
    pts_pos <- terra::vect(training_coords[pts_to_plot[i], ], geom=c("X", "Y"),
                           crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
    plot(pts_pos, add = T, col = "red", pch = 16, cex = .3)
    this_country <- country_lookup$country[which(country_lookup$val==country_df$layer[pts_to_plot[i]])]
    text(pts_pos, labels=this_country)
    cat("This point is in",
        this_country,
        ".\n")
    Sys.sleep(1)
  }
}

# FROM https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
sampled = apply(X = training_coords[, 1:2], MARGIN = 1, FUN = function(xy) r@data@values[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])

# # Breaking construction of cov_df into steps speeds up execution time - this
# # appears to be a quirk of how raster extraction and/or data frame construction
# # is done.
# {
#   n_blocks <- floor((1/250) * nrow(training_coords))
#   time_now <- Sys.time()
#   cov_df <- data.frame(raster::extract(covstack, training_coords[1:250, 1:2]))
#   for (i in 2:n_blocks){
#                     new_block <- data.frame(
#                     raster::extract(
#                       covstack,
#                       training_coords[(250*(i-1)+1:250*i), 1:2]
#                     )
#                   )
#                   cat("New block has",
#                       length(which(is.na(rowSums(new_block)))),
#                       "rows with NAs",
#                       "\n")
#                   print(head(new_block[which(is.na(rowSums(new_block))), ]))
#                   cov_df <- rbind(cov_df, new_block)
#   cat("nrow(cov_df) =",nrow(cov_df),"\n")
#   }
#   cov_df <- rbind(cov_df,
#                   data.frame(
#                     raster::extract(
#                       covstack,
#                       training_coords[(250*n_blocks+1):nrow(training_coords), 1:2]
#                       )
#                     )
#                   )
#   cov_build_time <- as.numeric(difftime(Sys.time(), time_now, units="mins"))
# }

# Try to plot out coordinates where we can't extract covariates:
bad_rows <- which(is.na(rowSums(cov_df)))
bad_pts <- terra::vect(training_coords[bad_rows, ], geom=c("X", "Y"),
                       crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
plot(euro_map_crop,
     col = "white",
     background = "azure2",
     axes = FALSE,
     buffer = FALSE,
     mar = c(0, 0, 0, 0))
plot(bad_pts, add = T, col = "red", pch = 16, cex = .3)

# Remove bad points
cov_df <- cov_df[-bad_rows, ]
training_coords <- training_coords[-bad_rows, ]
# Just redraw countries:
country_df <- data.frame(raster::extract(country_rast, training_coords[, 1:2]))

n_pts <- nrow(training_coords)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(cov_df), n_training)
xtrain <- cov_df[training, ]
ytrain <- training_coords$pos[training]
xtest <- cov_df[-training, ]
ytest <- training_coords$pos[-training]
countrytrain <- country_df$layer[training]
countrytest <- country_df$layer[-training]

combined_data <- cov_df
combined_data$pos <- training_coords$pos
combined_data$country <- sapply(1:nrow(combined_data),
                                FUN=function(i){
                                  country_lookup$country[which(country_lookup$val==country_df$layer[i])]})

df_train <- combined_data[training, ]
df_test <- combined_data[-training, ]

if (SAVE_FITS){
  save(xtrain, ytrain, xtest, ytest, file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/q1_train_test_data_ri.rds", sep = ""))
}

# Initialise model
basic_model <- rbart_vi(pos ~ . - country,
                    df_train,
                    group.by = country,
                    group.by.test = country,
                    test = df_test,
                    keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

# Eyeball RI terms:
ranef_means <- data.frame(name = names(basic_model$ranef.mean),
                          rem = unname(basic_model$ranef.mean))
ggplot(ranef_means) + geom_col(aes(rem, name)) + theme(axis.title=element_blank())

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_Q1.rds", sep = ""))
}

# retuned_model <- retune(basic_model)
# 
# k_retune = retuned_model$fit$model@node.hyperprior@k
# power_retune = retuned_model$fit$model@tree.prior@power
# base_retune = retuned_model$fit$model@tree.prior@base

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = df_train$country,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_sdm_Q1.rds", sep = ""))
}

# Eyeball RI terms:
vs_ranef_means <- data.frame(name = names(sdm$ranef.mean),
                          rem = unname(sdm$ranef.mean))
p <- ggplot(vs_ranef_means) + geom_col(aes(rem, name)) + theme(axis.title=element_blank())
p

covstack <- aggregate(covstack, fact = 10)
gc()

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack,
                      ri.data = country_rast,
                      ri.name = "country_rast",
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