library(reshape2)
library(readr)
library(dplyr)
library(data.table)
library(readxl)
library(rgeos)
library(sp)
library(maptools)
library(rgdal)
library(sf)
library(FNN)
library(raster)
library(future)
library(tcltk2)

readPoints <- function(csvfile){
  NRD <- read.csv(csvFile)
  NRD <- na.omit(melt(NRD, id = c("Pixel_id", "X", "Y")))
  PixelCoords <- NRD[c("Pixel_id", "X", "Y")]
  return(PixelCoords)
}

readgridPoints<- function(Modelgrd.Path,Model.Shape){
  Modelgrd %<-% readOGR(Modelgrd.Path,Model.Shape)
  gridCentroids <-gCentroid(Modelgrd,byid=TRUE)
  print(proj4string(Modelgrd))
  if (!compareCRS(HARNSP17ft,proj4string(Modelgrd))) {
    gridCentroids <-spTransform(gridCentroids,HARNSP17ft)
  }
  #-------------------------------------------------
  # Convert model cell points to a dataframe
  #-------------------------------------------------
  asSFGC <- st_as_sf(gridCentroids)
  ModelGridCoords <- do.call(rbind,st_geometry(asSFGC))
  ModelGridCoords <-cbind(Modelgrd$Row, Modelgrd$Column_,ModelGridCoords)
  
  return(ModelGridCoords)
}

#-------------------------------------------------
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
#-------------------------------------------------
HARNSP17ft  = CRS("+init=epsg:2881")
HARNUTM17Nm  = CRS("+init=epsg:3747")
latlongs = CRS("+proj=longlat +datum=WGS84")

basePath <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS/"
basePath <-  tk_choose.dir(default = basePath, caption = "Select input directory for biasNRD files")

#-------------------------------------------------
# Read one year for Pixels
#-------------------------------------------------
csvFile <- paste0(basePath, paste0("biasNRD2003.csv"))
csvFile <- choose.files(default=paste(basePath,'*.csv',sep='/'))

plan(multiprocess)

cat(paste('Reading',csvFile,'\n'))
f<- future({readPoints(csvFile)})

#-------------------------------------------------
# Read Model grid shapefile as polys 
# and convert to points and convert to data frame
#-------------------------------------------------
Modelgrd.Path <- "//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/CFWI/Data/From_SW_SJ"
Model.Shape <-"ECFTX_GRID_V3"
setwd(Modelgrd.Path)

cat(paste('Reading Model shapefile',Model.Shape,'\n'))
g<-future({readgridPoints(Modelgrd.Path,Model.Shape)})

#-------------------------------------------------
# Wait for values from futures
#-------------------------------------------------
cat(paste('Waiting for background processing to complete','\n'))
PixelCoords <-future::value(f)
ModelGridCoords <-future::value(g)

cat(paste('Background processing is complete','\n'))

#-------------------------------------------------
# Find closest points
#-------------------------------------------------
head(ModelGridCoords)
NearNeighbor = get.knnx(PixelCoords[,2:3],ModelGridCoords[,3:4],1)
ClosestPixels <- NearNeighbor[["nn.index"]]
ClosestDistance <-NearNeighbor[["nn.dist"]]
ModelGridCoords <-cbind(PixelCoords[ClosestPixels,1],ClosestDistance,ModelGridCoords)

#-------------------------------------------------
# Output data
#-------------------------------------------------
ModelGCdf <- as.data.frame(ModelGridCoords)
head (ModelGCdf)
names(ModelGCdf)<- c("PixelID","distance","row","col","X","Y")
ModelCellfile<-paste0(Modelgrd.Path,'/Testing.ModelCells.csv')
write_csv(as.data.frame(ModelGCdf),ModelCellfile)
gc()

#-------------------------------------------------
# FNN::get.knnx simple example
#
# x1 is like Pixels
# x2 is like modelgrid
#-------------------------------------------------
#   x1 = cbind(runif(5),runif(5))
#   x2 = cbind(runif(10),runif(10))
#   nn = get.knnx(x1,x2,1)
#-------------------------------------------------
# returns index and distance of x1 closest to x2
#-------------------------------------------------
#   nn[["nn.index"]]
#   nn[["nn.dist"]]
