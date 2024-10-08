# code to install all the packages listed below
# install.packages(c("dismo", "tidyterra", "colorRamps", "usdm", "ecospat", 
# "corrplot", "maxnet", "ENMeval", "predicts", "parallel", "randomForest", "biomod2"))

## ---- Load libraries ---------------------------------------------------------
library(geodata)
library(dplyr)
library(sf)
library(terra)
library(rnaturalearth)
library(ggplot2)
library(dismo)
library(tidyterra)
library(colorRamps)
library(usdm) # variable corr testing
library(ecospat) # variable corr testing, boyce index
library(corrplot) # variable corr testing
library(maxnet)
library(ENMeval)
library(predicts)
library(parallel) # parallel processing
library(randomForest)
library(PRROC) # for AUC-PR


## ---- Download & prepare species data-----------------------------------------
#data(bradypus) # brown-throated sloth (Bradypus variegatus)
bradypus <- sp_occurrence(genus='Bradypus', species='variegatus',
                                   download=TRUE, geom=TRUE, 
                                   removeZeros = TRUE) %>%
  # select columns of interest
  select(species, country, lon, lat, year, coordinateUncertaintyInMeters) %>%
  # remove high uncertainty records
  filter(coordinateUncertaintyInMeters < 10000) %>%
  # remove records before 2000
  filter(year > 2000) %>%
  # remove duplicates
  distinct(lon, lat, .keep_all = TRUE) %>%
  # convert to sf object
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# ggplot of records over South American countries
bvCountries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf", 
                                           country=unique(bradypus$country))
ggplot() +
  geom_sf(data = bvCountries) +
  geom_sf(data = bradypus, aes(color = species)) +
  theme_minimal()


## ---- Get env data-------------------------------------------------------
# download the bioclim variables
bioRasts <- worldclim_global(var="bio", res=5, path=getwd())

# rename the layers for convenience
names(bioRasts) <- paste0("bio", 1:19)

# select a subset of the bioclim variables
bioRasts <- bioRasts[[c("bio4", "bio5", "bio6", "bio15", "bio16", "bio17")]]

# we want an additional raster with biome information that is 
# provided with the dismo package
# file to path to where the rasters are saved
path <- file.path(system.file(package="dismo"), 'ex')
# list the files
list.files(path) # want want only the rasters
rastFiles <- list.files(path, pattern='biome.grd$', full.names=TRUE )
rastFiles

# read the biome raster
biome <- rast(rastFiles) 
biome
# plot the biome data
plot(biome)

# crop bioRasts to the extent of the biome data
bioRasts <- crop(bioRasts, biome)

# Let's stack the raster files together, but they have different resolutions!
try(bioRasts <- c(bioRasts, biome))

# we need to get the rasters to the same resolution
# first, let's check the resolution of the rasters
res(bioRasts)
res(biome)
res(biome)/res(bioRasts)

# let's also check their extent and origin - these all match
ext(bioRasts)
ext(biome)
origin(bioRasts)
origin(biome)

# resample the rasters to the same resolution
?terra::resample

# let's disaggregate the biome data to the same resolution as the bioclim data
biome <- disagg(biome, fact=res(biome)/res(bioRasts), method='near')
plot(biome)

# now we can stack the rasters...dang, we get a WARNING about the CRSs
bioRasts <- c(bioRasts, biome)
bioRasts

# remove the biome layer from bioRasts
bioRasts <- subset(bioRasts, "biome", negate=TRUE)

# need to project the rasters to the same CRS
biome <- project(biome, crs(bioRasts), method='near')

#...and now we can stack the rasters
bioRasts <- c(bioRasts, biome)
plot(bioRasts)

# lets plot the data
ggplot() + 
  geom_spatraster(data = bioRasts, aes(fill = bio5)) +
  scale_fill_gradientn(colors = rgb.tables(1000), na.value = "transparent", 
                       name = "bio5") + 
  geom_sf(data = bradypus) +  
  geom_sf(data = bvCountries, fill=NA) +
  theme_minimal() 


## ---- Extract env data-------------------------------------------------------
# check dimensions of bradypus
dim(bradypus)

# Next, we want to extract the environmental values at each
# of the occurrence locations.
envPres <- terra::extract(bioRasts, bradypus, cells=TRUE, xy=TRUE) %>%
  # remove spatial duplicates
  distinct(cell, .keep_all = TRUE)
head(envPres)
dim(envPres)


## ---- Create background (pseudo-absence) data --------------------------------
# setting a random seed to always create the same
# random set of points for this example
set.seed(896)

# create 10,000 random background points 
envBg <- spatSample(bioRasts, size=10000, method='random', na.rm=TRUE, 
                     cells=T, xy=TRUE)

# create a vector of 1's (presence) and 0's (background / pseudo-absence)
bvOcc <- c(rep(1, nrow(envPres)), rep(0, nrow(envBg)))

# select cell, and x-y columns
xyPres <- select(envPres, x, y)
xyBg <- select(envBg, x, y)

# select env columns
envPres.x <- select(envPres, names(bioRasts))
envBg.x <- select(envBg, names(bioRasts))  

# now we bind everything together into a data.frame
sdmData <- data.frame(cbind(rbind(xyPres, xyBg), # coords
                            bvOcc, # presence/background
                            rbind(envPres.x, envBg.x))) # env data

# biome is a factor, so define it that way
sdmData[,'biome'] = as.factor(sdmData[,'biome'])

# and have a look at the data
head(sdmData)
summary(sdmData)


## ---- Assess correlation structure -------------------------------------------
varCor <- cor(na.omit(sdmData[,-c(1:3,10)]))
corrplot::corrplot(varCor)
ecospat.cor.plot(na.omit(sdmData[,-c(1:3,10)]))

# clustering with dendrogram
allDistNew <- abs(as.dist(cor(sdmData[,-c(1:3,10)])))
allClusNew <- hclust(1 - allDistNew)
plot(allClusNew, hang=-1)

# Variance Inflation Factor
# report VIF by variable
usdm::vif(sdmData[,-c(1:3,10)])
# variable selection using VIF
usdm::vifstep(sdmData[,-c(1:3,10)])
# variable selection using correlation and VIF
usdm::vifcor(sdmData[,-c(1:3,10)], th=0.8) 

# remove problematic variable(s)
remVars <- vifstep(sdmData[,-c(1:3,10)])@excluded
sdmData <- sdmData[,-which(names(sdmData) %in% remVars)]


## ---- Create training / test datasets ----------------------------------------
set.seed(896)

# k-fold cross-validation partitioning
bv.kfold <- ENMeval::get.randomkfold(occs=xyPres,
                       bg=xyPres,
                       k=4) # 4 random folds, allows for a 75/25 split

# look at the structure, etc
str(bv.kfold)
table(bv.kfold$occs.grp)

# plot the data
evalplot.grps(pts=xyPres, pts.grp=bv.kfold$occs.grp, envs=stack(bioRasts))

# # spatial block partition: checkerboard
# bv.spBlock <- get.checkerboard1(occs=xyPres, envs=stack(bioRasts), 
#                          bg=xyBg, 
#                          aggregation.factor=100)
# evalplot.grps(pts = xyBg, pts.grp = cb1$bg.grp, 
#               envs = stack(bioRasts))
# evalplot.grps(pts = xyPres, pts.grp = cb1$occs.grp, 
#               envs = stack(bioRasts))

# make life a little easier
selTrainTest <- as.numeric(unlist(bv.kfold))

# create a training dataset
sdmData.train <- subset(sdmData, selTrainTest != 1)
dim(sdmData.train)
# create a testing dataset
sdmData.test <- subset(sdmData, selTrainTest <= 1)
dim(sdmData.test)

# make sf versions for plotting
sdmData.train.sf <- st_as_sf(sdmData.train, coords = c("x", "y"), crs=4326)
sdmData.test.sf <- st_as_sf(sdmData.test, coords = c("x", "y"), crs=4326)


## ---- Mahalanobis Distance ---------------------------------------------------
?mahal

# remove the categorical variable
bioRasts.noBiome <- subset(bioRasts, "biome", negate=TRUE)

# fit the model
mm <- mahal(stack(bioRasts.noBiome), # raster stack
            sdmData.train[sdmData.train$bvOcc==1,c("x","y")]) #presence-only data

# prediction is slow...
# predict the distribution
pm <- raster::predict(stack(bioRasts.noBiome), # raster stack
              mm, # model
              progress='text')

# raw output is 1-distance, see ?mahal
plot(pm) 

# let's convert the distance values to a p-value
# Mahal distances (D^2) are Chi-square distributed
mm.prob <- app(rast(pm), function(x, k=nlyr(bioRasts.noBiome)){
  x <- 1-x
  x <- x^2
  p_value <- pchisq(x, df = k, lower.tail = FALSE)
  return(p_value)
})

# plot the p-values
ggplot() + 
  geom_spatraster(data = mm.prob, aes(fill = lyr.1)) +
  scale_fill_gradientn(colors = rgb.tables(1000), na.value = "transparent", 
                       name = "Prob") + 
  geom_sf(data = subset(sdmData.train.sf, bvOcc==1)) +  
  geom_sf(data = bvCountries, fill=NA) +
  theme_minimal() 


## ---- Evaluate Mahalanobis Model ---------------------------------------------
# evaluate the model (presences + background)
?predicts::pa_evaluate
predAtTrain <- terra::extract(mm.prob, sdmData.train[sdmData.train$bvOcc==1,c("x","y")])$lyr.1
predAtTest <- terra::extract(mm.prob, sdmData.test[sdmData.test$bvOcc==1,c("x","y")])$lyr.1
predAtBg <- terra::extract(mm.prob, sdmData.train[sdmData.train$bvOcc==0,c("x","y")])$lyr.1

# evaluate using training data
evalTrain <- predicts::pa_evaluate(p=predAtTrain, # presences
              a=predAtBg) # background / absences
evalTrain

# evaluate using test data
evalTest <- predicts::pa_evaluate(p=predAtTest, # presences
              a=predAtBg) # background / absences
evalTest

# let's check model quality using the Boyce Index
boyceTrain <- ecospat::ecospat.boyce(mm.prob, obs=sdmData.train[sdmData.train$bvOcc==1,c("x","y")])
boyceTrain

boyceTest <- ecospat::ecospat.boyce(mm.prob, obs=sdmData.test[sdmData.test$bvOcc==1,c("x","y")])
boyceTest

# AUC-PR
prTrain <- PRROC::pr.curve(scores.class0 = predAtTrain, 
                    scores.class1 = predAtBg, 
                    curve = TRUE)
prTrain

prTest <- PRROC::pr.curve(scores.class0 = predAtTest, 
                   scores.class1 = predAtBg, 
                   curve = TRUE)
prTest

# To plot the PR curve
plot(prTrain)
plot(prTest)


## ---- Fit MaxEnt -------------------------------------------------------------
MaxEnt()
?MaxEnt

# folder to save the output
dir.create("/Users/mfitzpatrick/maxentOut")
filePath <- "/Users/mfitzpatrick/maxentOut"
#dir.create("/home/mfitzpatrick/maxentOut")
#filePath <- "/home/mfitzpatrick/maxentOut"

mx <- MaxEnt(x=bioRasts, # env data as a raster stack
             p=sdmData.train[sdmData.train$bvOcc==1,c("x","y")], # presence data
             a=sdmData.train[sdmData.train$bvOcc==0,c("x","y")], # background data
             factors='biome', # biome is categorical
             path=filePath) # where to save all the output
plot(mx)
mx

# predict back to geography
mxPred <- predict(mx, bioRasts)
plot(mxPred, col=rgb.tables(1000))

# # Predict in 'raw' format and save raster of prediction
# mxPred <- predict(mx, stack(bioRasts), args=c("outputformat=raw"), overwrite=TRUE,
#                   filename=paste0(filePath, '/maxent_predictionRAW.tif'))
# 
# # let's check model quality using the Boyce Index
# mxPred <- rast(paste0(filePath, '/maxent_predictionRAW.tif'))

# evaluate the model
predAtTrain <- terra::extract(mxPred, sdmData.train[sdmData.train$bvOcc==1,c("x","y")])$maxent
predAtTest <- terra::extract(mxPred, sdmData.test[sdmData.test$bvOcc==1,c("x","y")])$maxent
predAtBg <- terra::extract(mxPred, sdmData.train[sdmData.train$bvOcc==0,c("x","y")])$maxent

# evaluate using training data
evalTrain <- pa_evaluate(p=predAtTrain, # presences
                         a=predAtBg) # background / absences
evalTrain

# evaluate using test data
evalTest <- pa_evaluate(p=predAtTest, # presences
                        a=predAtBg) # background / absences
evalTest

# let's check model quality using the Boyce Index
boyceTrain <- ecospat.boyce(mxPred, obs=sdmData.train[sdmData.train$bvOcc==1,c("x","y")])
boyceTrain

boyceTest <- ecospat.boyce(mxPred, obs=sdmData.test[sdmData.train$bvOcc==1,c("x","y")])
boyceTest

# AUC-PR
prTrain <- pr.curve(scores.class0 = predAtTrain, 
                    scores.class1 = predAtBg, 
                    curve = TRUE)
prTrain

prTest <- pr.curve(scores.class0 = predAtTest, 
                   scores.class1 = predAtBg, 
                   curve = TRUE)
prTest

# To plot the PR curve
plot(prTrain)
plot(prTest)

# using the maxnet package instead
#?maxnet
# mx <- maxnet(p=sdmData.train$bvOcc, # pres-bg vector
#              data=sdmData.train[,names(bioRasts)])
# 
# # function to predict to a spatRast stack
# predict.maxnet.spatRasts <- function(model, rasters, nCores=NULL, 
#                                      type="cloglog", clamp=T){
#   # convert spatRast stack to a data.frame
#   rastTable <- na.omit(as.data.frame(rasters))
#   cells <- as.numeric(row.names(rastTable))
#   
#   # make empty mask raster to hold prediction
#   mask <- rasters[[1]]
#   mask[] <- NA
#   names(mask) <- paste0(type, '_prediction')
#   
#   # predict the model
#   pred <- predict(model, newdata=rastTable, clamp=clamp, type)
#   
#   # assign the prediction to the mask
#   mask[cells] <- pred
#   return(mask)
# }
# 
# # predict the model
# pred <- predict.maxnet.spatRasts(mx, bioRasts, type="cloglog")
# 
# # plot the prediction
# plot(pred, main='MaxEnt prediction')


## ---- Random Forest ----------------------------------------------------------
# the model formula
model <- bvOcc ~ bio4 + bio5 + bio15 + bio16 + bio17 + biome

# fit the RF model
?randomForest # lots of arguments, but let's keep it simple
rf.reg <- randomForest(model, data=sdmData.train, importance=T)

# note the warning - the function assumed regression instead of 
# classification
rf.reg
importance(rf.reg) # variable importance summary

# plot response curves
partialPlot(rf.reg, sdmData.train, bio5)
partialPlot(rf.reg, sdmData.train, biome)

# with pa as a factor to perform classification
model <- factor(bvOcc) ~ bio4 + bio5 + bio15 + bio16 + bio17 + biome
rf.class <- randomForest(model, data=sdmData.train, importance=T)
rf.class
importance(rf.class)

# look at variable importance
par(mfrow=c(1,2))
varImpPlot(rf.reg)
varImpPlot(rf.class)
dev.off()

# plot the prediction
par(mfrow=c(1,2))
pHS <- predict(stack(bioRasts), rf.reg, type="response")
plot(pHS, main='RF: Regression')

pPA <- predict(stack(bioRasts), rf.class, type="response")
plot(pPA, main='RF: Classification')
dev.off()

# evaluate the model
predAtTrain <- terra::extract(pHS, sdmData.train[sdmData.train$bvOcc==1,c("x","y")])
predAtTest <- terra::extract(pHS, sdmData.test[sdmData.test$bvOcc==1,c("x","y")])
predAtBg <- terra::extract(pHS, sdmData.train[sdmData.train$bvOcc==0,c("x","y")])

# evaluate using training data
evalTrain <- pa_evaluate(p=predAtTrain, # presences
                         a=predAtBg) # background / absences
evalTrain

# evaluate using test data
evalTest <- pa_evaluate(p=predAtTest, # presences
                        a=predAtBg) # background / absences
evalTest

# let's check model quality using the Boyce Index
boyceTrain <- ecospat.boyce(rast(pHS), obs=sdmData.train[sdmData.train$bvOcc==1,c("x","y")])
boyceTrain

boyceTest <- ecospat.boyce(rast(pHS), obs=sdmData.test[sdmData.train$bvOcc==1,c("x","y")])
boyceTest

# AUC-PR
prTrain <- pr.curve(scores.class0 = predAtTrain, 
                    scores.class1 = predAtBg, 
                    curve = TRUE)
prTrain

prTest <- pr.curve(scores.class0 = predAtTest, 
                   scores.class1 = predAtBg, 
                   curve = TRUE)
prTest

# To plot the PR curve
plot(prTrain)
plot(prTest)


## ---- Random Forest, weighted ------------------------------------------------
# equal weighting of presences and absences
sum(sdmData.train$bvOcc==1)
sum(sdmData.train$bvOcc==0)

weights <- c(rep(0.1, sum(sdmData.train$bvOcc)),
             rep(1/sum(sdmData.train$bvOcc==0),sum(sdmData.train$bvOcc==0)))

# fit, predict, plot the RF model
model <- bvOcc ~ bio4 + bio5 + bio15 + bio16 + bio17 + biome
rf.reg.w <- randomForest(model, data=sdmData.train, importance=T, weights=weights)
pHS.w <- predict(stack(bioRasts), rf.reg.w, type="response")
plot(pHS.w, main='RF: Regression, weighted')

# plot response curves
partialPlot(rf.reg.w, sdmData.train, bio5)
partialPlot(rf.reg.w, sdmData.train, biome)

# evaluate the model
predAtTrain <- terra::extract(pHS.w, sdmData.train[sdmData.train$bvOcc==1,c("x","y")])
predAtTest <- terra::extract(pHS.w, sdmData.test[sdmData.test$bvOcc==1,c("x","y")])
predAtBg <- terra::extract(pHS.w, sdmData.train[sdmData.train$bvOcc==0,c("x","y")])

# evaluate using training data
evalTrain <- pa_evaluate(p=predAtTrain, # presences
                         a=predAtBg) # background / absences
evalTrain

# evaluate using test data
evalTest <- pa_evaluate(p=predAtTest, # presences
                        a=predAtBg) # background / absences
evalTest

# let's check model quality using the Boyce Index
boyceTrain <- ecospat.boyce(rast(pHS.w), obs=sdmData.train[sdmData.train$bvOcc==1,c("x","y")])
boyceTrain

boyceTest <- ecospat.boyce(rast(pHS.w), obs=sdmData.test[sdmData.train$bvOcc==1,c("x","y")])
boyceTest

# AUC-PR
prTrain <- pr.curve(scores.class0 = predAtTrain, 
                    scores.class1 = predAtBg, 
                    curve = TRUE)
prTrain

prTest <- pr.curve(scores.class0 = predAtTest, 
                   scores.class1 = predAtBg, 
                   curve = TRUE)
prTest

# To plot the PR curve
plot(prTrain)
plot(prTest)

