#packages <- c("dismo", "gbm", "rgeos", "ENMeval", "spThin", "maps", "maptools", "raster", "rgbif", "adehabitatHR", "sp", "tidyverse")
#lapply(packages, install.packages)
#sapply(packages, library, character.only = T)

install.packages("dismo")
install.packages("gbm")
install.packages("rgeos")
install.packages("ENMeval")
install.packages("spThin")
install.packages("maps")
install.packages("maptools")
install.packages("raster")
install.packages("rgbif")
install.packages("adehabitatHR")
install.packages("sp")
install.packages("tidyverse")
install.packages("sqldf")
install.packages("roxygen2")
install.packages("testthat")

library(dismo)
library(gbm)
library(rgeos)
library(ENMeval)
library(spThin)
library(maps)
library(maptools)
library(raster)
library(rgbif)
library(adehabitatHR)
library(sp)
library(tidyverse)
library(sqldf)
library(roxygen2)
library(testthat)

setwd("~/Desktop/CastellanosOSOS2018/ECOLOGICAL_NICHE_MODELING_IN_R")

#####INTRODUCTION TO ECOLOGICAL NICHE MODELS#####
###Step 1: Gather coordinate data

SPEC <- occ_search(scientificName = "Dipodomys compactus", fields = c('name', 'stateProvince', 'country', 'county', 'locality', 'decimalLongitude', 'decimalLatitude', 'geodeticDatum', 'datasetName', 'institutionCode', 'year', 'issues'), limit = 1000)

#GBIF is a great place to start your search for coordinate information (simply because of all the places it grabs data from)
#a few fields that I find most helpful for determining reasonable coordinates to use are the ones used above
#stateProvince, country, county, and locality are always good to have for georeferencing purposes and to determine that coordinates are placed where they should be
#issues is a good field that you can whittle records down by

SPEC
#There are over 600 records of Dipodomys compactus on GBIF
#Not all of them are likely great, though

plot(SPEC$data[, c("decimalLongitude", "decimalLatitude")], pch = 19, col = "red")
maps::map(add = T)
#the maps::map function is helpful for providing quick boundaries
levels(as.factor(SPEC$data$stateProvince))
#mathwoman.gif
#What is the range of D. compactus, again?
#Oftentimes, a placeholder of 0, 0 is used for records, which can be a pain to deal with
#Still, doesn't explain the Oklahoma or California records. Always check the accepted/acceptable range of your study species.

#`dplyr::filter` is a powerful tool that acts like the `base::subset` function but becomes very useful if you start creating pipelines with the magrittr package
OCC <- filter(SPEC$data, !is.na(decimalLongitude))
#For the purposes of our workshop today, we will get rid of all the records without coordinate information. In an actual study, you would set these to the side and georeference them (if possible)

DUPL <- duplicated(OCC[ , c("decimalLongitude", "decimalLatitude")]) 
#base::duplicated determines which rows have duplicated information in the decimalLongitude and decimalLatitude columns, which is quite a lot

sum(DUPL) 
#the base::duplicated function returns a logical vector of the number of rows of the data.frame we entered. As such, you can use the base::sum function to find out how many duplicated coordinates we have

head(OCC[DUPL, ])
#take a look at the duplicate rows just to make sure everything worked

OCC <- OCC[!DUPL, ] 
#this takes out the duplicated rows (remember, ! is the logical NOT operator; you can use it to take the opposite of something)

dim(OCC)
#Well...
#Those 600 records certainly went away fast

#But you can still do so much more to get so much less info

gbif_issues()
summary(as.factor(OCC$issues))
#The best way I've found to filter out specific issues (due to how they are placed in the issue column) is to run a `grep` function and use the `|` symbol to specify everything you want to get rid of.

grep("cdiv|cdrepf|cdreps|txmatnon|cucdmis|zerocd", OCC$issues)
#Because it returns a vector of row numbers that match the pattern, we can use - to negate these numbers and grab everything else.

dim(OCC[-grep("cdiv|cdrepf|cdreps|txmatnon|cucdmis|zerocd", OCC$issues), ])

#There are still a bunch of things you can do to clean up your data (check for outliers, make sure the georeferencing was done correctly, etc.).

#HOWEVER, let's all start from the same point
PRES <- read.csv("Dipo_presence_OSOS.csv", header = T)

plot(PRES, pch = 19, col = "red")
maps::map(add = T)

#Now that we have our presence data, what else do we need

###Step 2: gather your raster data
PRED <- raster::getData("worldclim", var = "bio", res = 2.5)
#We have used the raster::getData function to grab GDAL data for administrative boundaries, but it can also be used to grab the oft used Bioclim dataset of 19 bioclimatic variables

PRED
#can see that its class is a RasterStack along with its resolution, dimensions, extent, and the names of the variables
#neat.jpg

###Step 3: gather your background data
###The background points give the algorithm information about the environmental space available in your study region and thus should reflect a reasonable amount of your study region
###If you increase the size of your study region, you can increase your AUC (increase the possibility of "correctly" identifying areas as absent of the species)
##Current reasonable recommendations are to choose the current extent of your points and one half to one degree (but see Anderson, R.P and A. Raza. 2010. J Biogeogr. The effect of the extent of the study region on GIS models of species geographic distributions and estimates of niche evolution: Preliminary tests with montane rodents (genus Nephelomys) in Venezuela.)

MCP <- mcp(SpatialPoints(PRES, proj4string = CRS("+proj=longlat +datum=WGS84")), percent = 100) #creates a polygon that encapsulates all of your chosen points

#Take a look at your polygon and points
plot(MCP) + points(PRES) 
maps::map(add = T)

extent(MCP)
#raster::extent gives the geographical extent of the chosen spatial object (in this case, our polygon)

attributes(MCP)
#take a look at what makes up the SpatialPolygonsDataFrame that is MCP. In particular, we are looking at $bbox
#a SpatialPolygonsDataFrame is considered an S4 class in R, which means it is indexed using @

EXT <- extent(c(MCP@bbox[, 1] - 1, MCP@bbox[, 2] + 1)[c(1, 3, 2, 4)]) 
#and now we have the extent of our study!

PNTS <- randomPoints(PRED, n = 10000, ext = EXT, extf = 1)
#the dismo::randomPoints function is useful for generating random background points within an area
#it uses a raster mask restricted by an extent object to generate n number of points within our study region

#A note about background points: The number chosen is usually recommended to be a large number (e.g. 10,000), especially for MaxEnt. However, the number of background points that works best seems to differ with algorithm choice
#For more information, see Barbet-Massin, M. et al. 2012. Methods in Ecology and Evolution. Selection psuedo-absences for species distribution models: How, where, and how many?

#When using natural history collection specimens, it is often recommended to apply the same bias in selecting background points as is found in collection patterns
#But how do we determine this bias?
#You can see collection bias by obtaining records of similar species (for this example, we will grab records of all heteromyid and cricetid rodents) 
#We will use these records to create a background mask and then bias the selection of background points to within this mask


#DO NOT RUN the below lines of code
#The rgbif::occ_search functions here are grabbing all cricetid and heteromyid records with coordinate information from Texas and Tamaulipas. This will run for a while.
#####
CRIC <- occ_search(taxonKey = 3240723, stateProvince = c("Texas", "Tamaulipas"), hasCoordinate = T, fields = c("name", "decimalLongitude", "decimalLatitude"), limit = 200000)
HMYD <- occ_search(taxonKey = 5504, stateProvince = c("Texas", "Tamaulipas"), hasCoordinate = T, fields = c("name", "decimalLongitude", "decimalLatitude"), limit = 200000)
#when querying for multiple options (such as here for both Texas and Tamaulipas), the result is in a list, requiring indexing of the different elements of the list to grab the $data portion
MASK <- data.frame(rbind(CRIC[[1]]$data, CRIC[[2]]$data, HMYD[[1]]$data, HMYD[[2]]$data))
coordinates(MASK) <- ~decimalLongitude + decimalLatitude
#I often find it easier to use the sp::coordinates function to turn a data.frame into a SpatialPointsDataFrame when I am trying to crop it to a certain geographic extent
MASK <- crop(MASK, EXT)
#the raster::crop function is handy because it allows us to restrict the points to the extent of our study area (the EXT object)
#####


###DO RUN this below code
#since searching for 50,000 records takes a slight amount of time, just read in the results from the above few lines of code and turn them into a SpatialPointsDataFrame
MASK <- read.csv("enm_mask_OSOS.csv", header = T, row.names = 1)
coordinates(MASK) <- ~decimalLongitude + decimalLatitude
 
plot(MASK, pch = 19, col = "blue", cex = 0.4)
maps::map(add = T)
maps::map(add = T, "county", "texas")
#the maps::map function unfortunately doesn't have Mexican municipalities, so we can only show country boundaries and Texas county boundaries
#It is slightly messy, but it is easy and doesn't take absurd amounts of time to plot

CIRC <- circles(MASK, d = 50000, lonlat = T)
#the dismo::circles function will take the points from the MASK object and create circles of a 50km radius around them

plot(CIRC)
maps::map(add = T)
#You can see that overlapping polygons are dissolved into one another

CIRC@polygons
#The result of the function is one SpatialPolygon

RAST <- rasterToPolygons(crop(PRED, EXT))
#This line uses both the raster::rasterToPolygons and raster::crop functions to create a smaller RasterStack the same size as our study extent and then turn it into a SpatialPolygonsDataFrame
#"But why did you do this? We have the circle polygon and can just grab random data points from within it."
#Well you can see that a not so insignificant portion of the polygon is over the Gulf. We don't want points landing in the ocean.
#Also the border from say a polygon from GDAL isn't the same border as the Bioclim raster data, so we can still get points with NA info even if they technically aren't in the Gulf

POL <- gIntersection(CIRC@polygons, RAST)
#the rgeos::gIntersection grabs the intersecting bits from the circle polygon (CIRC@polygons) and our raster polygon (RAST)

plot(POL)
#Check out our new background mask!

PNTS <- spsample(POL, 10000, type = "random", iter = 25)
#here our sp::spsample function will sample 10,000 points randomly within our POL object. The iter argument just gives it the number of times to try grabbing random points before declaring it a failure

plot(PNTS, pch = 19, cex = 0.2, col = "blue")
maps::map(add = T)
#Look at that spread of random background points!

ABSV <- as.data.frame(PNTS)
colnames(ABSV) <- c("lon", "lat")
#Now, let's turn it into a data.frame and give it more useful column names than just x and y

#So now we have coordinate data and raster data. Where do we go from here?

###Step 4: building a model
#The BIOCLIM algorithm for environmental niche modeling isn't particularly used anymore (other algorithms often outperform it except in a few special cases)
#However, it is easy to use and doesn't require you to go through some sort of R version of a SAW movie's torture chamber to run it (looking at you, MaxEnt)

BCLM <- bioclim(PRED, PRES)
#all that the dismo::bioclim function needs is a raster of predictor variables and the coordinates. The function will extract the needed information for each of the records (probably using raster::extract, but I haven't checked)

plot(BCLM)
#Cool! It's a plot of something.
#Actually it appears to be a plot of values for each of the records for certain bioclimatic variables (default appears to be the first two) with a box around a certain percentage of the points

plot(BCLM, a = 1, b = 2, p = 0.85) 
#You can seemingly adjust the axes and the size of the box using the arguments a, b, and p

BCP <- predict(PRED, BCLM, ext = EXT, progress = '')
#Now things get slightly weird here. We are using the raster::predict function to create a raster prediction of the model object (BCLM) restricted to the extent (EXT). Simple enough. However, there is also the dismo::predict function that does the same thing but with the model object coming before the raster object in the function
#The raster and dismo packages have a few authors in common, which likely led to this. I believe that the raster version of the predict function is more applicable to a wider variety of model objects, so I normally use that version

BCP
#Look! It's a RasterLayer. Let's plot it
plot(BCP)
points(PRES, pch = 19, col = "blue", cex = 0.5)
maps::map(add = T)
#Plot the raster prediction, add our points for context, and add country boundaries
#Congrats, you have a visualization of an ENM that you can use (but that is likely worthless because it's BIOCLIM and we've done the minimum because this is a 2 hour workshop)

###Step 5: evaluating your model
ME <- evaluate(PRES, ABSV, BCLM, PRED)
#the dismo::evaluate function takes presence coordinates, background coordinates, a model object and a RasterStack of predictor variables to generate an evaluation of the model

ME
#You can see that it gives the number of presence and absence/background points used and the AUC of the model

plot(ME, "ROC")
#AUC (area under the curve) is probably the most used evaluation statistic for ecological niche models. It ranges from 0 to 1 with a score of 0.5 meaning the model is essentially random and a score of 1 is a perfect prediction of a species distribution
#It is also a highly contentious statistic, especially in recent literature. Unfortunately given the prevalence of presence-only modeling techniques, it seems to be necessary for now despite its flaws

threshold(ME, "no_omission")
#Thresholds are something used in ENMs to display the predictions in a binary way (presence vs. absence). If you type in ?threshold and check under the Values section, you will get an idea of the options available to set as a threshold. For more insight into the use of thresholds in ENMs, see Norris, D. 2014. Tropical Conservation Science. Model Thresholds are more important than presence location type: Understanding the distribution of lowland tapir (Tapirus terrrestris) in a continuous Atlantic forest of southeast Brazil

plot(BCP > threshold(ME, "no_omission"))
maps::map(add = T)
#Here I used the minimum presence threshold ("no_omission") which sets the threshold at the rate of the lowest presence point used in the model

attributes(ME)
#this just sent up a lot of horrifying numbers
#one of the most important things to look at here is the confusion table (ME@confusion) which shows the range of true positives, false positives, false "negatives", true "negatives", for a given threshold (ME@t)
#I used quotes around "negatives" because in a presence-only framework there are no absences, only background points

#k-fold cross validation is often used to evaulate the models

group <- list(kfold(PRES, 5), kfold(ABSV, 5))
#dismo::kfold will randomly, equally assign each record with a number from 1 to k (given here as 5)
#I put these in a list because I really like lists

DATA <- lapply(seq(1, 5, 1), function(x) list(PRES[which(group[[1]] != x), ], PRES[which(group[[1]] == x), ], ABSV[which(group[[2]] != x), ], ABSV[which(group[[2]] == x), ]))
#this is obnoxious, but is a fairly simple way of assigning presence and background points to k groups of training and test datasets
#For each number from 1 to k, the base::lapply function will create a list with 4 elements in it: 1) the training presence dataset with 4/5 of the records in it, 2) the test presence dataset with 1/5 of the records in it, 3) the training background dataset with 4/5 of the background points, and 4) the test background dataset with the remaining 1/5 background points

STATS <- rep(0, 5)


for(i in 1:5) {
M <- bioclim(PRED, DATA[[i]][[1]])
#creates the model using the training data
ME <- evaluate(DATA[[i]][[2]], DATA[[i]][[4]], M, PRED)
#evaluates the model using the test data
STATS[i] <- ME@auc
#grabs the AUC for the partition and places it in the i position of the empty STATS vector
}


STATS
mean(STATS)
#And now we have the AUC results of each partition of our k-fold cross validation

###Additional options: Filtering presence points
#Oftentimes, occurence data will be biased in some way towards locations that are easier to access or nearer to population centers. In these events, you can get large clusters of occurences that the model can overfit to, thus reducing the usefulness of the model
#See Boria, R.A. et al. 2014. Ecological Modelling. Spatial filtering to reduce sampling bias can improve the performance of ecological niche models. 
#We will be using an environmental filter (described in Varela, S.A. et al. 2014. Ecography. Environmental filters reduce the effects of sampling bias and improve predictions of ecological niche models.)

envSample<- function (coord, filters, res, do.plot=TRUE){
  
  n<- length (filters)
  pot_points<- list ()
  for (i in 1:n){
    k<- filters [[i]] [!is.na(filters[[i]])]
    ext1<- range (k)
    ext1 [1]<- ext1[1]- 1
    x<- seq(ext1[1],ext1[2], by=res[[i]])
    pot_points[[i]]<- x
  }
  pot_p<- expand.grid(pot_points)
  
  ends<- NULL
  for (i in 1:n){
    fin<- pot_p [,i] + res[[i]]
    ends<- cbind (ends, fin)
  }
  
  pot_pp<- data.frame (pot_p, ends)
  pot_pp<- data.frame (pot_pp, groupID=c(1:nrow (pot_pp)))
  rows<- length (filters[[1]])
  filter<- data.frame(matrix(unlist(filters), nrow=rows))
  real_p<- data.frame (coord, filter)
  
  names_real<- c("lon", "lat")
  names_pot_st<- NULL
  names_pot_end<- NULL
  sql1<- NULL
  for (i in 1:n){
    names_real<- c(names_real, paste ("filter", i, sep=""))
    names_pot_st<- c(names_pot_st, paste ("start_f", i, sep=""))
    names_pot_end<- c(names_pot_end, paste ("end_f", i, sep=""))
    sql1<- paste (sql1, paste ("real_p.filter", i, sep=""), sep=", ")   
  }
  
  names (real_p)<- names_real
  names (pot_pp)<- c(names_pot_st, names_pot_end, "groupID")
  
  conditions<- paste ("(real_p.filter", 1, "<= pot_pp.end_f", 1,") and (real_p.filter", 1, "> pot_pp.start_f", 1, ")", sep="")
  for (i in 2:n){
    conditions<- paste (conditions, 
                        paste ("(real_p.filter", i, "<= pot_pp.end_f", i,") and (real_p.filter", i, "> pot_pp.start_f", i, ")", sep=""), 
                        sep="and")
  }
  
  selection_NA<- sqldf(paste ("select real_p.lon, real_p.lat, pot_pp.groupID",   
                              sql1, "from pot_pp left join real_p on", conditions, sep=" "))
  
  selection<- selection_NA [complete.cases(selection_NA),]
  
  final_points<- selection[!duplicated(selection$groupID), ]
  coord_filter<- data.frame (final_points$lon, final_points$lat) 
  names (coord_filter)<- c("lon", "lat")
  
  if (do.plot==TRUE){
    par (mfrow=c(1,2), mar=c(4,4,0,0.5))
    plot (filters[[1]], filters[[2]], pch=19, 
          col="grey50", xlab="Filter 1", ylab="Filter 2")
    points (final_points$filter1, final_points$filter2, 
            pch=19, col="#88000090")
    plot (coord, pch=19, col="grey50")
    maps::map(add=T)
    points (coord_filter, pch=19, col="#88000090")
    
  }
  coord_filter
}

#If you read closely, you'll see that it requires coordinates, environmental variables to filter by, and a resolution to filter by. 
DATA <- cbind.data.frame(PRES, raster::extract(PRED, PRES))
head(DATA)

EFILT <- envSample(DATA[, 1:2], filters = list(DATA$bio1, DATA$bio12), res = list(diff(range(DATA$bio1))/10, diff(range(DATA$bio12))/10))

#Now create a model and prediction using the spatially filtered points
ECLM <- bioclim(PRED, EFILT)
ECP <- predict(PRED, ECLM, ext = EXT, progress = '')
plot(ECP)
maps::map(add = T)

###Additional options: Reducing correlation of bioclimatic variables
#Many studies will check to see if the bioclimatic variables they are using are highly correlated (For an example see Cooper, D.M. et al. 2016. Diversity Distrib. Predicted Pleistocene-Holocene range shifts of the tiger (Panthera tigris))
#Some have found that this reduces over-parameterization in the models

COR <- layerStats(crop(PRED, EXT), "pearson", na.rm = T)
#raster::layerStats can be used to compute correlation and covariance for Raster objects. Here we are using it to compute the Pearson correlation coefficient for each of the bioclimatic variables
#the na.rm = T argument removes the cells of the RasterStack object with NA values (often the ocean) so that they don't interfere with the calculations
#the result is a list with the correlation coefficients and the mean values for each variable

symnum(COR[[1]])
#If you wish to visualize this, run stats::symnum on the first element of the list
#Sometimes, you'll have multiple highly correlated variables and have to choose one. You can do this via jacknife tests or by looking at a PCA

PREDI <- crop(PRED, EXT)
PCA <- prcomp(na.omit(getValues(PREDI)))
summary(PCA)
plot(PCA$x[, 1:2])
PREDI[[1]][!is.na(PREDI[[1]])] <- PCA$x[, 1]
PREDI[[2]][!is.na(PREDI[[2]])] <- PCA$x[, 2]
PCAR <- stack(PREDI[[1:2]])
plot(PCAR[[1]])
plot(PCAR[[2]])

###Additional options: Schoener's D for comparison
#sometimes you have predictions from multiple algorithms or from the use of multiple arguments within the same algorithm and want to see how the prediction has changed
MAPS <- stack(BCP, ECP)
calc.niche.overlap(MAPS)
#ENMeval::calc.niche.overlap takes a RasterStack of the predictive maps created for each model and calculates the Schoener's D statistic of niche similarity
#A result of 0 means that the two maps have no similarity and a result of 1 means that they are exactly the same

###Additional options: MaxEnt
#Maxent is a slight pain to setup to run in R
#it requires the rJava package to run which will let you know in the R console window if you do absolutely anything on your computer
#First you will have to actually download Maxent on your computer
#https://biodiversityinformatics.amnh.org/open_source/maxent
#From this download, you will grab the maxent.jar file and place it in the folder that comes up when you run the below line
system.file("java", package = "dismo")
#Chances are you will possibly run into more errors involving pathways
#As always, googling the error is the best way out of the problem

#When everything is finally setup, you can run maxent easily as with dismo::bioclim
install.packages("rJava")
library(rJava)

MXNT <- maxent(PRED, PRES)
MXP <- predict(PRED, MXNT, ext = EXT, progress = '')
maps::map(add = T)

#Good introductory reading on Maxent can be found in:
#Merow, C. et al. 2013. Ecography. A practical guide to MaxEnt for modeling species' distributions: what it does, and why inputs and settings matter.
#Elith, J. et al. 2011. Diversity Distrib. A statistical explanation of MaxEnt for ecologists

###Additional options: Boosted Regression Trees
#You can use boosted regression trees using the gbm package
#There is a nice vignette for the dismo package that goes over the use of boosted regression trees in creating ENMs
#Additionally, for a good source of information, see Elith, J. et al. 2008. Journal of Animal Ecology. A working guide to boosted regression trees.