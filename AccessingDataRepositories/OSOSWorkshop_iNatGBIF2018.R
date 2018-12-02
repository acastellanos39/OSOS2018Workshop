#packages <- c("rgbif", "raster", "maps", "rgdal", "viridis", "picante", "rinat", "sp", "broom", "maptools", "ape", "tidyverse", "magrittr")
#lapply(packages, install.packages)
#sapply(packages, library, character.only = T)

install.packages("rgbif")
install.packages("raster")
install.packages("maps")
install.packages("rgdal")
install.packages("viridis")
install.packages("picante")
install.packages("rinat")
install.packages("sp")
install.packages("broom")
install.packages("maptools")
install.packages("ape")
install.packages("tidyverse")
install.packages("magrittr")

library(rgbif)
library(raster)
library(maps)
library(rgdal)
library(viridis)
library(picante)
library(rinat)
library(sp)
library(broom)
library(maptools)
library(ape)
library(tidyverse)
library(magrittr)


setwd("~/Desktop/CastellanosOSOS2018/ACCESSING DATA REPOSITORIES ") #change to whatever you desire

#####INTRODUCTION TO iNATURALIST
###package::rinat
###TAXON-LEVEL QUERIES
OCC <- get_inat_obs(taxon_name = "Lontra canadensis", maxresults = 25) 
#to query for specific species, use the rinat::get_inat_obs function

#Can query  by taxon name, taxon id (from iNat), or from a general query
#the max number of results is set to 100 by default. Here I set it to 25 just to show a quick example

##Take a look at the data
class(OCC)
head(OCC)
summary(OCC) 
#OR
str(OCC) #gives structure of object, can look messy at time but does give the class of each column (which can prove helpful)


OCC <- get_inat_obs(taxon_name = "Lontra canadensis", year = 2018, maxresults = 1000, quality = "research", geo = T)
#queries can be further filtered by year (can only give one year, not a range, but there is a way around this for enterprising individuals...), month(s), status of georeferenced coordinates, and quality 
#Here I set quality to only those observations deemed research quality (determined when the wider iNaturalist community agrees on a species level ID). If you are using iNat to determine where a species is, you should probably only look at these results.

#try plotting the coordinates of 2018 research grade otter observations
plot(OCC[, c("longitude", "latitude")], pch = 19, col = "red", cex = 0.5) 
maps::map(add = T) 
#the maps::map function can prove helpful if you need quick boundaries. The database argument is set to "world" by default, but you can choose a few US options by changing it to "state" or "county"
maps::map(add = T, "state")
#*gasp* US state boundaries


OCC <- get_inat_obs(taxon_name = "Lontra canadensis", geo = T, bounds = c(25, -106, 37, -94), maxresults = 1000)
#can apply a bounding box to filter out observations (S-LAT, W-LON, N-LAT, E-LON)
#I created a bounding box that essentially covers Texas

plot(OCC[, c("longitude", "latitude")], pch = 19, col = "red")
maps::map(add = T, database = "county", "texas") 
#can see a general plot of otter observations in/around Texas


###GROUP-LEVEL QUERIES
#Besides querying for a specific species or taxonomic group, one can query for a specific iNaturalist group project (e.g., Herps of Texas or Mammals of Texas). 


#DON'T run the below code, it can take a while 
MAMM <- get_inat_obs_project(grpid = "mammals-of-texas", type = "observations")

#DO RUN the below line to read a .csv file with the Mammals of Texas Project information
MAMM <- read.csv("iNatMammals.csv", header = T)
summary(MAMM) 
#Look at all those squirrels!
#check Iconic.taxon.name for a laugh
#check summary for User.login, Quality.grade, Scientific.name

plot(MAMM[, c("Longitude", "Latitude")], pch = 19, col = "red") 
#well that seems to look like Texas


#Let's take a look at trends by each county
#iNat doesn't have explicit county data for each observation, so let's create it
#First, take out all those records without coordinates to make things easier on us. 
MAMM <- filter(MAMM, !is.na(Longitude))
#^the above line checks for which rows have NA values (this the is.na function) for longitude and filters them out to only select those rows that are the opposite of it. The ! symbol acts as a negator for logical values (the result of is.na), so it grabs those with coordinate data

CNTY <- raster::getData(country = "USA", level = 2) #the raster::getData function will grab the R equivalent of a shapefile of US county data (level 0 = country, level 1 = stateProvince, level 2 = county/municipality)

class(CNTY) 
#it is a SpatialPolygonsDataFrame

head(CNTY) 
#with an associated data.frame (you can access the data.frame only by using CNTY@data if you wish)
#Take a look at the NAME_1 and NAME_2 columns

CNTY <- CNTY[CNTY$NAME_1 %in% "Texas", ]
#This subsets the polygons to just Texas *ComeAndTakeItFlagWaving.gif*
#When subsetting a data.frame, you have a few options
#You can use the base::subset function
#subset(CNTY, CNTY$NAME_1 == "Texas")
#OR you can subset within square brackets
#CNTY[CNTY$NAME_1 == "Texas", ] (or how I did it above)

CNTY$NAME_2 
#Oh, wow! A vector of Texas county names

coordinates(MAMM) <- ~Longitude + Latitude 
#turns the MAMM data.frame into a SpatialPointsDataFrame by specifying the coordinate columns. I found that this is easier than using sp::SpatialPointsDataFrame

crs(MAMM) <- crs(CNTY) 
MAMM <- spTransform(MAMM, crs(CNTY))
#all coordinates from iNat are stored in unprojected lat/lon coordinates using the WGS84 datum, so you can grab the coordinate reference system from the county data

OVR <- over(MAMM, CNTY) 
#the raster::over function is useful for determining what specific polygons (y argument, shown here as CNTY) a point in the x argument (here shown as MAMM) fall into 
dim(OVR) 
head(OVR)
#this means you now have county information for each iNat record by simply writing out OVR$NAME_2

MAMM@data <- cbind(MAMM@data, county = OVR$NAME_2) #let's add that county information to our SpatialPointsDataFrame

OBS <- sapply(CNTY$NAME_2, function(x) nrow(filter(MAMM@data, county == x))) 
#grabs the number of observations by county
#Hoo boy, time for an apply function explanation
#for the above line, we are grabbing the vector of Texas county names and indexing the iNat mammal data by each county (this is what everything after function(x) is concerned with). `sapply` (results in a vector) and `lapply` (results in a list) have the following syntax: sapply(x, function(x) whatever function you can think up).
#so this is going county by county and counting up the number of observations

SR <- sapply(CNTY$NAME_2, function(x) MAMM@data %>% filter(county == x) %$% Scientific.name %>% unique %>% length) 
#grabs the number of species found in each county
#this is similar to the code written to count observations, but instead further indexes by only grabbing the Scientific.name column and then grabs only the number of unique names in this column
#Soooo, you may be wondering what kind of mystical arts I'm attempting with `%$% and `%>%`. I'm lazy, so I constructed a pipeline. Pipelines in R are created in using the handy `%>%` batch of symbols from the magrittr package. Essentially what it does is take the result from the left and use it as the first argument in whatever is after `%>%` (`%$%` is a special usage that only grabs a specific column or index from the left)
#So essentially, we are doing `length(unique(filter(MAMM\@data, county = x)$Scientific.name))` but in an easier to read way (you can actually understand it from left to right)
#Ta-dah! Species richness

USER <- sapply(CNTY$NAME_2, function(x) MAMM@data %>% filter(county == x) %$% User.login %>% unique %>% length)
#one more thing...
#Lets document how many unique users there are recording observations in each county

CNTY$inat.obs <- OBS 
CNTY$inat.sr <- SR
CNTY$inat.user <- USER
#place our observation, species richness, and user number vectors as columns within our SpatialPolygonsDataFrame of Texas counties

#How about we visualize what we just wrought?
#ggplot2 is a fickle companion, so we will use the broom::tidy function to break down our observations and create a sort of choropleth map
OBSR <- tidy(CNTY, region = "inat.obs")
#broom::tidy will turn our SpatialPointsDataFrame into a regular data.frame and ID each row according to the region specified. This will allow the fill argument of the ggplot2::aes function to fill each county appropriately

class(OBSR$id) #character

OBSR$id <- as.numeric(OBSR$id)
#we don't want these numbers to be characters when filling in each polygon, so change them to numeric

OBSR[which(OBSR$id == 0), "id"] <- NA
#those polygons with no observations will be set to NA so as to not confuse them with counties that only have a handful of observations

SRR <- tidy(CNTY, region = "inat.sr")
#do the same thing with species richness
SRR$id <- as.numeric(SRR$id)
SRR[which(SRR$id == 0), "id"] <- NA

#and the same thing with the user information
USR <- tidy(CNTY, region = "inat.user")
USR$id <- as.numeric(USR$id)
USR[which(USR$id == 0), "id"] <- NA

INAT <- list(OBSR, SRR, USR)
#Lists are great classes to store things in. Here we have the data.frames featuring county observation, species richness, and user numbers and can access each one by using [[]] 

ggplot(data = INAT[[2]], aes(x = long, y = lat, fill = id, group = group)) + geom_polygon(color = "black", size = 0.2) + scale_fill_viridis(name = "iNaturalist Species Richness", option = "D", na.value = NA)

filter(MAMM@data, county == "Brazos")
#Now that we have county data associated with each record, you can subset the records by a specific county (such as Brazos Co.)
MAMM@data %>% filter(county == "Brazos") %$% Scientific.name %>% unique
#You can use the base::unique function to determine which species are found in Brazos Co.

#####INTRODUCTION TO GBIF
###package::rgbif

CHAET <- occ_search(scientificName = "Chaetodipus hispidus", limit = 10)
class(CHAET) 
#hmmm, it is a gbif class. Can't really review the information we want (presumably you want info regarding the record(s))

attributes(CHAET)
CHAET$data 
#WHAT. THE. HELL. is a tibble?
#trimmed down data frame
#relax, breathe and change it to something you are more comfortable with if you desire. You are in control

CHAET <- data.frame(CHAET$data)
head(CHAET)
#Well...
#That is a lot of data. Do we really need that much data? Do you REALLY want the vernacularName, license, and publishingCountry for example?

summary(as.factor(CHAET$county))
#the summary of a vector of factors will give you a handy breakdown of what the vector is comprised of. Neat!

CHAET <- data.frame(occ_search(scientificName = "Chaetodipus hispidus", fields = c('name', 'infraspecificEpithet', 'stateProvince', 'country', 'county', 'locality', 'decimalLongitude', 'decimalLatitude', 'geodeticDatum', 'datasetName', 'institutionCode', 'recordNumber', 'catalogNumber', 'year', 'preparations'), limit = 10)$data)
#here we not only listed the set of fields that we wanted, but also only grabbed the $data portion and turned it into a data.frame in one line

#But what if you wanted more than 200,000 records (the limit imposed by the draconian creators of the package)? *thinkingfaceemoji*
CNTRY <- c("US", "MX")
CHAET <- lapply(1:2, function(x) as.data.frame(occ_search(scientificName = "Chaetodipus hispidus", country = CNTRY[x], fields = c('name', 'stateProvince', 'country', 'county', 'locality', 'decimalLongitude', 'decimalLatitude', 'geodeticDatum', 'datasetName', 'institutionCode', 'recordNumber', 'catalogNumber', 'year', 'preparations'), limit = 15)$data))
#this lapply function will call the rgbif::occ_search function twice, each with the country specified as each element of the CNTRY vector

length(CHAET)
CHAET 
#This results in a list composed of Chaetodipus hispidus records from the US and from Mexico

#"But I don't want a list"
CHAET <- do.call(rbind, CHAET)
dim(CHAET)
#the base::do.call function will apply the function specified in the first part of the argument to each portion of the list, thus binding the rows of each list element together
#There are other methods to use a function like base::rbind or base::cbind on each element of a list (e.g., plyr::rbind.fill) but base::do.call is accessible to everyone and is simple to use

HMYD <- occ_search(taxonKey = 5504, limit = 20)
#searching using the taxonKey (always found in the url of the GBIF page of the group you are looking at) can be very helpful as well

HMYD
#There are over 200,000 heteromyid records on GBIF 
#Look at all those cool species of Dipodomys, Perognathus, and Chaetodipus!

#What if you weren't interested in a particular species or family, but wanted to know what was in a particular institution? rgbif has you covered there as well

###DO NOT RUN the below code
BRTC <- occ_search(taxonKey = 359, basisOfRecord = "PRESERVED_SPECIMEN", institutionCode = "TCWC", fields = c('name', 'stateProvince', 'country', 'county', 'locality', 'decimalLongitude', 'decimalLatitude', 'geodeticDatum', 'datasetName', 'institutionCode', 'recordNumber', 'catalogNumber', 'year', 'preparations'), limit = 200000)$data
###
#The above code searches for taxonKey 359 (Mammalia) of PRESERVED_SPECIMENS (so not videos, pictures, or hearsay of a species) at the institutionCode formerly known as the TCWC (our Biodiversity Research and Teaching Collections is still considered the TCWC by GBIF)

###DO RUN the below code instead
BRTC <- read.csv("BRTCMammals.csv", header = T)
dim(BRTC)
summary(BRTC)
#Look at all those Peromyscus
#These data aren't up to date, but do include much of the BRTC Mammal collection

summary(BRTC$country) 
length(unique(BRTC$country))
#The BRTC Mammal collection houses specimens from over 70 countries

plot(BRTC[, c("decimalLongitude", "decimalLatitude")])
maps::map(add = T)
#specimens with coordinate information are very US biased

TEX <- filter(BRTC, stateProvince == "Texas")
#index the BRTC records to only pick those from Texas 
#you can use dplyr::filter or base::subset or %in% to do this

summary(TEX)
#Hmmm that doesn't seem quite right. Why is Aguascalientes and Alabama showing up in stateProvince

TEX <- droplevels(filter(BRTC, stateProvince == "Texas"))
#use base::droplevels to get rid of unneeded levels when using the base::summary function

dim(TEX) 
#over 18,000 specimens from Texas
summary(TEX)

plot(TEX[, c("decimalLongitude", "decimalLatitude")])
#Looks like some specimens weren't georeferenced correctly. Let's ignore those.

#However, we have more county records (only 311 NA values) than coordinate records (7,953 NA values)
#Since georeferencing is not the goal of this module, we will visualize records and species richness at the county level

SR <- sapply(CNTY$NAME_2, function(x) TEX %>% filter(county == x) %$% name %>% unique %>% length)
#You've seen this before (hopefully you have been paying some modicum of attention). The base::sapply function takes each county name, subsets the BRTC Texas data with it and counts the unique species names found

OBS <- sapply(CNTY$NAME_2, function(x) TEX %>% filter(county == x) %>% nrow)
#same as above, but uses base::nrow to count the number of records per county

CNTY$brtc.sr <- SR
CNTY$brtc.obs <- OBS

TSR <- tidy(CNTY, region = "brtc.sr")
TSR$id <- as.numeric(TSR$id)
TSR[which(TSR$id == 0), "id"] <- NA

TOBS <- tidy(CNTY, region = "brtc.obs")
TOBS$id <- as.numeric(TOBS$id)
TOBS[which(TOBS$id == 0), "id"] <- NA

DIV <- list(TSR, TOBS)
#Again, the above several lines are creating a data.frame of the species richness and record values for each county in a way that ggplot2 can easily use to fill

ggplot(data = DIV[[1]], aes(x = long, y = lat, fill = id, group = group)) + geom_polygon(color = "black", size = 0.2) + scale_fill_viridis(name = "BRTC Mammal\nSpecies Richness", option = "D", na.value = NA)
#You can use \n to create a return in a title in ggplot2

#But not all species richness is created the same
#5 Peromyscus species are different than a community that has marsupials, rodents, and cervids, for example
#But how do we show this...
#Great question, follow below for more!

TREE <- read.nexus("ELE_1307_sm_SA1.nex")
#the ape::read.nexus function will allow us to read in a nexus file that contains a phylogenetic tree
#In this case, it is a tree of 5,020 mammal species from S.A. Fritz et al. 2009. Ecology Letters. Geographical variation in predictors of mammalian extinction risk: Big is bad, but only in the tropics

class(TREE) 
summary(TREE)
attributes(TREE)
#this multiPhylo class is essentially a list that contains three different trees with different dates
#For our purposes, we will use TREE$mammalST_MSW05_bestDates

summary(TREE$mammalST_MSW05_bestDates)
#shows number of tips, nodes, and branch length distribution

attributes(TREE$mammalST_MSW05_bestDates)
head(TREE$mammalST_MSW05_bestDates$tip.label)
#if you notice, all tip label names for the tree have a "_" within them 

head(TEX$name)
#not so with our species names from GBIF
#could be a problem but it is a simple fix, so let's fix this!

TEX$name <- sapply(strsplit(as.character(TEX$name), " "), paste, collapse = "_")
#here we use the base::sapply function to split all names by a space and replace it with a _ using the base::paste function
#base::strsplit only works on characters, so I had to use base::as.character to change the class of the TEX$name vector

#Now comes the tricky part. Unfortunately the picante::pd function to determine phylogenetic diversity requires a community matrix, so we have to create one

SAMP <- matrix(0, nrow = length(CNTY$NAME_2), ncol = length(unique(TEX$name)))
#create a matrix of 0s with each row being a county name and each column being a species name

row.names(SAMP) <- CNTY$NAME_2
colnames(SAMP) <- unique(TEX$name)
#populate the empty row and column names with the names of the counties and species

for(i in 1:nrow(SAMP)) {
	SAMP[i, unique(TEX[TEX$county %in% CNTY$NAME_2[i], "name"])] <- 1
}
#then we create a for loop that goes to each county row and places a 1 in each species column that is in that county (because the species names are the column names, we can use "name" to reference both)

SAMP <- SAMP[, colnames(SAMP)[!colnames(SAMP) %in% setdiff(colnames(SAMP), TREE$mammalST_MSW05_bestDates$tip.label)]]
#the picante::pd function does not play well with species being found in the matrix that are not found in the tree
#This also looks horrifying but is really quite simple
#We find the difference between species names from the BRTC Mammal Collection and the tree tip labels using the base::setdiff function and then exclude the names that are not found in the tree

PD <- pd(SAMP, TREE$mammalST_MSW05_bestDates)
head(PD)
#the picante::pd function returns a data.frame comprised of vectors of phylogenetic diversity (PD) and species richness (SR) for each county

CNTY$pd <- log(PD$PD)

TPD <- tidy(CNTY, region = "pd")
TPD$id <- as.numeric(TPD$id)
TPD[which(TPD$id == 0), "id"] <- NA

DIV <- list(TSR, TOBS, TPD)
#add to our list

ggplot(data = DIV[[3]], aes(x = long, y = lat, fill = id, group = group)) + geom_polygon(color = "black", size = 0.2) + scale_fill_viridis(name = "BRTC Phylogenetic Diversity", option = "D", na.value = NA)
#Plot!