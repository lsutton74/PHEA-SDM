## Species Distribution Model for the Philippine Eagle
## Author: Luke J. Sutton
## The Peregrine Fund
## lsutton@peregrinefund.org
## 30-09-2022

## Preprint DOI: https://www.biorxiv.org/content/10.1101/2021.11.29.470363v2

# Citation: Sutton, L.J., Ibañez, J.C., Salvador, D.I., Taraya, R.L., Opiso, G.S., Senarillos, T.L.P. & McClure, C.J.W. (2022). 
# Priority conservation areas and a global population estimate for the Critically Endangered Philippine Eagle derived from modelled range metrics using remote sensing habitat characteristics. 
# bioRxiv. DOI: https://doi.org/10.1101/2021.11.29.470363 [preprint]


## R packages not on CRAN
devtools::install_github("narayanibarve/ENMGadgets")
remotes::install_github('adamlilith/legendary', dependencies=TRUE)
remotes::install_github('adamlilith/omnibus', dependencies=TRUE)
remotes::install_github('adamlilith/statisfactory', dependencies=TRUE)
remotes::install_github("adamlilith/enmSdm")

library(raster) 
library(sp) 
library(rgdal)
library(rgeos)
library(gtools)

library(glmnet) 
library(maxnet)
library(ENMeval)
library(enmSdm)
library(dismo) 
library(ENMGadgets) 

library(prettymapr) 
library(maptools) 
library(maps) 
library(mapdata)
library(rangeBuilder)
library(MetBrewer)

library(omnibus)
library(statisfactory)
library(legendary)
library(scales)
library(plotrix)

setwd("C:/PHEA-SDM")
rasterOptions(progress="text")
memory.limit()

# load covariate raster layers for Mindanao
list.env <- mixedsort(list.files(paste(getwd(), "/data/MIN", sep = ""), 
                                 full.names = TRUE, 
                                 pattern = ".tif"))
list.env 

# stack raster files
env <- stack(list.env) 
env
names(env)
plot(env)

# load range extent shapefile for Mindanao
ext <- readOGR(dsn="C:/PHEA-SDM/data/GADM", layer="Mindanao")
summary(ext)

# randomly sample n=100 occurrences from first raster
library(virtualspecies)

vs <- sampleOccurrences(x=env[[1]], 
                        n=100,
                        type = "presence only",
                        extract.probability = FALSE,
                        sampling.area = ext,
                        detection.probability = 1,
                        correct.by.suitability = FALSE,
                        error.probability = 0,
                        bias = "no.bias",
                        bias.strength = 50,
                        bias.area = NULL,
                        weights = NULL,
                        sample.prevalence = NULL,
                        replacement = TRUE,
                        plot = TRUE)
str(vs)

locs <- vs$sample.points[,1:2]
head(locs)

names <- c("lon", "lat")
names(locs) <- names
head(locs)

plot(ext, col="lightgrey", border="grey60")
points(locs$lon, 
       locs$lat, 
       col="steelblue", 
       pch=21, 
       cex=0.7)

# Select best fit model
# All models AICc <2 and consider model with lowest OR10

# Block cross-validation
# Model transfer in space

evalb <- ENMevaluate(occ=locs, 
                     env=env, 
                     RMvalues=seq(1, 10, 0.5), # grid search to 10
                     fc=c("L", "LQ"),          # linear and quadratic functions
                     algorithm="maxnet",
                     method='block',
                     n.bg=1000,                # 1:10 ratio presences:background points
                     rasterPreds=TRUE,
                     clamp=TRUE, 
                     parallel = TRUE, 
                     numCores = 4)
evalb
aicmods <- which(evalb@results$AICc == min(na.omit(evalb@results$AICc)))
evalb@results[aicmods,]
evalb@results
write.csv(evalb@results, "./results/block.csv")

# fit penalized logistic regression model
library(glmnet)
library(maxnet)
library(enmSdm)
library(ENMeval)

# withold a 20% sample for testing 
fold <- kfold(x=locs, k=5)
test <- locs[fold == 1, ]
train <- locs[fold != 1, ]

# create training environment
TrainEnv <- extract(x=env, y=train)
head(TrainEnv)
nrow(TrainEnv)

# generate random background/absence points 
bg <- randomPoints(mask=env, 
                   n=1000,  # 1:10 ratio
                   ext=ext,
                   lonlatCorrection=TRUE)
head(bg)
nrow(bg)

# extract values from model environment for bg points
absvals <- extract(x=env, y=bg)
head(absvals)

# presence/absence column 
presabs <- c(rep(1, nrow(TrainEnv)), rep(0, nrow(absvals)))
head(presabs)

# combine presence/absence data with model environment values
sdmdata <- data.frame(cbind(presabs, rbind(TrainEnv, absvals)))
head(sdmdata)

# check for NA values
row.has.na <- apply(sdmdata, 1, function(x){any(is.na(x))})
sum(row.has.na)

# remove NA values
sdmdata <- sdmdata[!row.has.na,]
head(sdmdata)

# subset of dataset without the presence and absence values
env.data <- sdmdata[ ,-1]
head(env.data)

# input model parameters from ENMevaluate function to maxnet function arguments

# run maxnet function to fit the SDM using glmnet
mn <- maxnet(p = sdmdata$presabs, 
             data = env.data,
             f = maxnet.formula(p = sdmdata$presabs, 
                                data = env.data, 
                                classes = "lq"), # feature classes or functions
             regmult = 3) # level of coefficient penalty

mn$call
names(mn)
str(mn)
summary(mn)

mn$alpha
mn$entropy

# extract beta coefficients
round(mn$betas,3)
write.csv(round(mn$betas,3), "./results/beta-coefficients.csv")

# save model as RData
saveRDS(mn, file = './results/maxnet.RData')

# load saved model
mn <- readRDS('./results/maxnet.RData')
summary(mn)
round(mn$betas,3)

# covariate response curves
plot(x=mn,
     common.scale=TRUE,
     type="cloglog",
     ylab="Predicted value",
     mod=mn)

par(mfrow=c(1,1))

# plot coefficient paths
plot.glmnet(x = mn, xvar = "norm", label = TRUE)
plot.glmnet(x = mn, xvar = "lambda", label = TRUE)
plot.glmnet(x = mn, xvar = "dev", label = TRUE)

# predict
library(ENMeval)

pred.mn <- maxnet.predictRaster(mod=mn,
                                env=env, 
                                type="cloglog", 
                                clamp=TRUE) 
pred.mn

# save raster
writeRaster(pred.mn, 
            "./results/Mindanao_Continuous.tif", 
            overwrite=TRUE)

# re-load raster if needed
pred.mn <- raster("./results/Mindanao_Continuous.tif")
pred.mn

# colour ramp palette
MetBrewer::colorblind_palettes

pal <- rev(met.brewer("Johnson", n=100, type="continuous"))

# or...
#pal <- rev(met.brewer("Tiepolo", n=100, type="continuous"))

# plot map
plot(pred.mn,
     maxpixels=500000,
     col=pal,
     interpolate=FALSE,
     axes=F,
     zlim=c(0.0,1.0), 
     legend.only=F,
     bty="n",
     box=F)

points(locs$lon, 
       locs$lat, 
       pch=21, 
       col="white",
       bg="black",
       cex=0.6) 

# Model Evaluation

# partial ROC
library(ENMGadgets)

write.csv(test, "./pROC/testdata.csv")
test <- read.csv("./pROC/testdata.csv")

# save raster to pROC file
writeRaster(pred.mn, "./pROC/predmn.tif", overwrite=TRUE)

proc <- PartialROC(PresenceFile="./pROC/testdata.csv", 
                   PredictionFile="./pROC/predmn.tif", 
                   OmissionVal=0.1, 
                   RandomPercent=50, 
                   NoOfIteration=1000,
                   OutputFile="./pROC/TestRoc.csv") 

head(proc)
round(mean(proc$AUC_ratio),3)  
round(sd(proc$AUC_ratio),3)    
round(range(proc$AUC_ratio),3)

# Continuous Boyce Index (CBI)
predPres <- extract(x = pred.mn, 
                    y = cbind(locs$lon, locs$lat),
                    na.rm=TRUE)

predPres <- na.omit(predPres)
head(predPres)
min(predPres)

bg <- as.data.frame(bg)

predBg <- extract(x = pred.mn, y = cbind(bg$x, bg$y))
predBg <- na.omit(predBg)
head(predBg)

# get proportion of presences or background sites in each bin
presDistrib <- hist(predPres, 
                    plot=FALSE, 
                    breaks=seq(0, 1, length.out=11))$counts

bgDistrib <- hist(predBg, 
                  plot=FALSE, 
                  breaks=seq(0, 1, length.out=11))$counts

# convert to proportion of sites
presDistrib <- presDistrib / sum(presDistrib)
bgDistrib <- bgDistrib / sum(bgDistrib)

# P/E plot
pe <- presDistrib / bgDistrib

plot(pe,
     xlab='Prediction class (1 to 10)',
     ylab='Predicted-to-Expected Ratio', 
     main='P/E Plot',
     col='red', 
     pch=16)

cbi <- contBoyce(pres=predPres, 
                 bg=predBg,
                 numBins = 101, 
                 binWidth = 0.1,
                 autoWindow = TRUE, 
                 method = "spearman", 
                 dropZeros = TRUE,
                 graph = FALSE,
                 na.rm=TRUE)
round(cbi,3)


# CBI Random 5-fold cross-validation
head(locs)
occnames <- c("x", "y")
names(locs) <- occnames
head(locs)

envdata <- raster::extract(env, locs)
head(envdata)

dv <- cbind(locs, envdata)
head(dv)
tail(dv)

# calculate k-folds for presences and background sites
head(dv)
nrow(dv)
kPres <- kfold(x = dv, k=5)

bg2 <- cbind(bg, absvals)
head(bg2)
nrow(bg2)

kBg <- kfold(x = bg2, k=5)

head(kPres)
head(kBg)

# map
plot(ext, main='K-fold #1')

points(bg2$x,#[kBg==1],
       bg2$y,#[kBg==1],
       col='lightgrey',
       cex=0.7,
       pch=20)

points(dv$x, dv$y)

points(dv$x[kPres==1],
       dv$y[kPres==1],
       bg='red',
       pch=21)

legend('topright',
       legend=c('Training presence', 'Test presence'),
       pch=c(1, 16),
       col=c('black', 'red'),
       bg='white',
       cex=0.8)


# for storing CBI
cbiRandom <- rep(NA, 5)
names(env)

predictors <- c("B1Red",
                "B2NIR",
                "B7SWIR",
                "EvergreenForest",
                "HumanFootprintIndex",
                "LeafAreaIndex")


# cycle through each k-fold
for (i in 1:5) {
  
  omnibus::say('K-fold ', i, ':', post=0)
  
  # make training data frame with predictors and vector of 1/0 for
  # presence/background... using only points not in this k-fold
  envData <- rbind(
    dv[kPres!=i, predictors],
    bg2[kBg!=i, predictors]
  )
  
  presBg <- c(rep(1, sum(kPres!=i)), rep(0, sum(kBg!=i))
  
  trainData <- cbind(presBg, envData)
  
  # check for NA values
  row.has.na <- apply(trainData, 1, function(x){any(is.na(x))})
  sum(row.has.na)
  
  # remove NA values
  trainData.filtered <- trainData[!row.has.na,]
  trainData <- trainData.filtered
  
  # tuned model
  model <- trainMaxNet(
    data=trainData,
    regMult=3,
    classes='lq',
    verbose=FALSE
  )
  
  # predict to presences and background sites
  predPres <- raster::predict(model, 
                              newdata=dv[kPres==i, ], 
                              type='cloglog')
  
  predBg <- raster::predict(model, 
                            newdata=bg2[kBg==i, ], 
                            type='cloglog')
  
  thisCbi <- contBoyce(pres=predPres, 
                       bg=predBg,
                       numBins = 101, 
                       binWidth = 0.1,
                       autoWindow = TRUE, 
                       method = "spearman", 
                       dropZeros = TRUE,
                       graph = FALSE,
                       na.rm = TRUE)
  
  omnibus::say('CBI = ', round(thisCbi, 3))
  
  cbiRandom[i] <- thisCbi
  
}

omnibus::say('Mean CBI:', round(mean(cbiRandom), 3))


# Project to Eastern Visayas and Luzon

# load Eastern Visayas and Luzon shapefiles
setwd("C:/PHEA-SDM/data/GADM")

ev <- readOGR(dsn="C:/PHEA-SDM/data/GADM", layer="EasternVisayas")
plot(ev, col="lightgrey")

lz <- readOGR(dsn="C:/PHEA-SDM/data/GADM", layer="Luzon")
plot(lz, col="lightgrey")

# Philippines shapefile
shp <- raster::getData(country='PHL', level=0)
plot(shp, col="lightgrey", border="grey60")

# load covariate rasters for Philippines
setwd("C:/PHEA-SDM")

list.phl <- mixedsort(list.files(paste(getwd(), "./data/PHL/", sep = ""), 
                                 full.names = TRUE, 
                                 pattern = ".tif"))
list.phl

# stack raster files
phlenv <- stack(list.phl)
phlenv
names(phlenv)
plot(phlenv)

phl_names <- c("B1Red",
               "B2NIR",
               "B7SWIR",
               "EvergreenForest",
               "HumanFootprintIndex",
               "LeafAreaIndex")

names(phlenv) <- phl_names 
plot(phlenv)

# define projection system of data
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs ")
proj4string(phlenv) <- crs.geo 
phlenv

# crop to EV shapefile
evcrop <- raster::crop(phlenv, ev)
evcrop
plot(evcrop)

# mask
evmask <- raster::mask(x= evcrop, mask = ev)
evmask
plot(evmask)

# stack
ev.env <- stack(evmask)
ev.env

# project/predict
ev.mn <- maxnet.predictRaster(mod=mn,
                              env=ev.env, 
                              type="cloglog", 
                              clamp=TRUE) 
ev.mn

par(mfrow=c(1,1))
plot(ev.mn, col=pal)

# save model raster
writeRaster(ev.mn, 
            "./results/EasternVisayas_Continuous.tif",
            overwrite=TRUE)

#ev.mn <- raster("./results/EasternVisayas_Continuous.tif")
#ev.mn

# Project to Luzon

# crop to shapefile
lzcrop <- crop(phlenv, lz)
lzcrop
plot(lzcrop)

# mask
lzmask <- raster::mask(x= lzcrop, mask = lz)
lzmask
plot(lzmask)

# stack
lz.env <- stack(lzmask)
lz.env

# project/predict
lz.mn <- maxnet.predictRaster(mod=mn,
                              env=lz.env, 
                              type="cloglog", 
                              clamp=TRUE) 
lz.mn

par(mfrow=c(1,1))
plot(lz.mn, col=pal)

# save raster
writeRaster(lz.mn, 
            "./results/Luzon_Continuous.tif",
            overwrite=TRUE)

#lz.mn <- raster("./results/Luzon_Continuous.tif")
#plot(lz.mn, col=pal)

# plot continuous range map
r1 <- raster("./results/Mindanao_Continuous.tif")
r2 <- raster("./results/EasternVisayas_Continuous.tif")
r3 <- raster("./results/Luzon_Continuous.tif")

c <- merge(r1, r2, r3)
plot(c, col=pal)

writeRaster(c, "./results/PHEA_Continuous.tif", overwrite=TRUE)

## Reclassify continuous map to binary threshold

# evaluation statistics
eval <- dismo::evaluate(p=as.vector(predPres), 
                        a=as.vector(predBg), 
                        tr=seq(0, 1, by=0.01))
eval
str(eval)

threshold(eval)

# Threshold that maximizes the sum of sensitivity and specificity
# maxTSS
maxSeSp <- eval@t[which.max(eval@TPR + eval@TNR)]
maxSeSp

# define matrix
maxt <- c(0,0.490, 0,0.490, 1, 1)
rclmat <- matrix(maxt, ncol=3, byrow=TRUE)
rclmat

# reclassify continuous to binary
rc <- reclassify(c, rclmat)
rc

# define projection system of data
proj4string(rc) <- crs.geo 
rc

plot(rc, 
     xlab="Longitude", 
     ylab="Latitude", las=1,
     legend=F,
     main = "maxTSS")

writeRaster(rc, "./results/PHEA_maxTSS.tif", overwrite=T)


# Calculate IUCN range metrics
library(redlistr)

# project to Transverse cylindrical equal area - IUCN standard
crs.proj <- "+proj=tcea +lon_0=122.2558465 +datum=WGS84 +units=m +no_defs"

rcp <- projectRaster(rc, crs=crs.proj)
rcp

# calculate AOO
aoo <- getAOO(input.data = rcp, 
              grid.size = 2000,
              min.percent.rule = T, 
              percent = 1)
aoo

# calculate model AOH and run Gap Analysis 

# Raster to polygon conversion
# maxTSS threshold
rc <- reclassify(c, rclmat)
rc

# 8-neighbor rule
patchID8 <- clump(x = rc, 
                  directions=8,
                  gaps=TRUE)
plot(patchID8)

# number of patches identified
cellStats(patchID8, max)

# convert raster layer to shapefile
poly <- rasterToPolygons(x = patchID8, 
                         n = 4,
                         na.rm = TRUE,
                         digits = 12,
                         dissolve = TRUE)
poly
#plot(poly)

# simplify raw polygon
sim <- gSimplify(spgeom = poly, 
                 tol = 0.01, 
                 topologyPreserve=FALSE)
sim
#plot(sim)

# smooth polygon
library(smoothr)

sm <- smooth(x = sim, method = "chaikin")
sm
plot(sm)

aoh <- as(sm, "SpatialPolygonsDataFrame")
aoh

# write shapefile object
writeOGR(obj=aoh, 
         dsn="C:/PHEA-SDM/results", 
         layer="PHEA_AOH", 
         overwrite_layer=T,
         driver="ESRI Shapefile")

# load new shapefile if needed
aoh <- readOGR(dsn="C:/PHEA-SDM/results", layer="PHEA_AOH")
summary(aoh)

crs(aoh)  
crs(aoh) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(aoh)

plot(shp, 
     col = "lightgrey",
     border="grey70", 
     lwd=.5)

plot(aoh, 
     add=TRUE, 
     #xlab="Longitude", 
     #ylab="Latitude", las=1,
     axes=F,
     col= "darkkhaki",
     border = "darkkhaki",
     bty="n",
     box=F,
     main = "")

# plot MCP for max EOO
ch <- gConvexHull(aoh)
ch
plot(ch, add=T, lty=3)

ch <- as(ch, "SpatialPolygonsDataFrame")

# write shapefile object
writeOGR(obj=ch, 
         dsn="C:/PHEA-SDM/results", 
         layer="PHEA_EOO_max", 
         overwrite_layer=T,
         driver="ESRI Shapefile")

# mask to land area within MCP
int <- raster::intersect(shp, ch)
plot(int)
plot(ch, add=T, lty = 4, border = "red")

# write shapefile object
writeOGR(obj=int, 
         dsn="C:/PHEA-SDM/results", 
         layer="PHEA_EOO_min",
         overwrite_layer=T,
         driver="ESRI Shapefile")

# load range extent
ext2 <- readOGR(dsn="C:/PHEA-SDM/data/GADM", layer="PHEA_extent")
ext2

# project to Transverse cylindrical equal area
crs.proj <- "+proj=tcea +lon_0=122.2558465 +datum=WGS84 +units=m +no_defs"
aoh2 <- spTransform(aoh, crs.proj)
shp2 <- spTransform(shp, crs.proj)

# calculate AOH in km2
aoh.area <- round(raster::area(aoh2) / 1000000, 1)
aoh.area               # km2 for each polygon
round(sum(aoh.area),0) # total km2


# Territorial area in km2 based on median Home Range Estimates
med <- 73  # 95% KDE
min <- 64  # r-LoCoH
max <- 90  # 99% MCP

# Median number of pairs
T.hat <- round((sum(aoh.area)/med)*2,1)/2
T.hat

# Minimum number of pairs
T.hat.min <- round((sum(aoh.area)/max)*2,1)/2
T.hat.min

# Maximum number of pairs
T.hat.max <- round((sum(aoh.area)/min)*2,1)/2
T.hat.max

# max EOO
ch <- gConvexHull(aoh2)
ch.area <- round(raster::area(ch) / 1000000, 1)
round(sum(ch.area),0) # total km2

# min EOO
int <- raster::intersect(shp2, ch)
int.area <- round(raster::area(int) / 1000000, 1)
round(sum(int.area),0) # total km2


## GAP ANALYSIS

# set WDPA target representation based on AOH
t <- max(0.1, min(1, -0.375*log10(sum(aoh.area))+2.126))
round(t,2)

# load WDPA shapefile
wdpa <- readOGR(dsn="C:/PHEA-SDM/data/PA", layer="PHEA_WDPA")
wdpa

# intersect WDPA with AOH
wdpa.aoh <- raster::intersect(x = aoh, y = wdpa)
wdpa.aoh

plot(aoh, 
     axes=F,
     col= "darkkhaki",
     border = "darkkhaki",
     bty="n",
     box=F,
     main = "")

plot(wdpa.aoh, 
     add=T, 
     border="red")

wdpa.aoh2 <- spTransform(wdpa.aoh, crs.proj)

# calculate suitability in WDPA network area in km2
poly.area <- round(area(wdpa.aoh2) / 1000000, 1)

# km2 for each polygon
poly.area  

# total km2
sum(poly.area) 

# % of WDPA in AOH area
round((sum(poly.area)/sum(aoh.area)*100),1)

# load KBA shapefile
kba <- readOGR(dsn="C:/PHEA-SDM/data/PA", layer="PHEA_KBA")
crs(kba)

# intersect KBA with AOH
kba.aoh <- raster::intersect(x = aoh, y = kba)
kba.aoh

plot(aoh, col="darkkhaki", border="darkkhaki")

plot(kba.aoh, 
     add=T, 
     border="grey70")

kba.aoh2 <- spTransform(kba.aoh, crs.proj)

# calculate suitability in KBA network area in km2
poly.area <- round(area(kba.aoh2) / 1000000, 1)

# km2 for each polygon
poly.area  

# total km2
sum(poly.area) 

# % of KBA in AOH area
round((sum(poly.area)/sum(aoh.area)*100),1)
