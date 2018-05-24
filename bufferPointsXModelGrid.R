library(data.table)
library(rgeos)
library(sp)
library(rgdal)
library(sf)
library(FNN)
library(raster)
library(future)
library(tcltk2)

readPoints <- function(csvfile){
  #-------------------------------------------------
  #  Read NexRad Pixels coordinates into
  #  data frame
  #-------------------------------------------------
  NRD <- utils::read.csv(csvFile)
  NRD <- stats::na.omit(data.table::melt(NRD, id = c("Pixel_id", "X", "Y")))
  PixelCoords <- NRD[c("Pixel_id", "X", "Y")]
  return(PixelCoords)
}

readgridPoints<- function(Modelgrd.Path,Model.Shape){
  #-------------------------------------------------
  #  -Read arcGIS shape file of model mesh
  #  -Convert mesh polygons to points 
  #     referencing the center of the model cells
  #  -Change spatial reference if necessary
  #  -Convert point coordinates into a data frame
  #-------------------------------------------------
  Modelgrd %<-% rgdal::readOGR(Modelgrd.Path,Model.Shape)
  gridCentroids <- rgeos::gCentroid(Modelgrd,byid=TRUE)
  print(sp::proj4string(Modelgrd))
  if (!raster::compareCRS(HARNSP17ft,sp::proj4string(Modelgrd))) {
    gridCentroids <- sp::spTransform(gridCentroids,HARNSP17ft)
  }
  
  # Convert model cell points to a dataframe
  asSFGC <- sf::st_as_sf(gridCentroids)
  ModelGridCoords <- do.call(base::rbind,sf::st_geometry(asSFGC))
  ModelGridCoords <-base::cbind(Modelgrd$Row, Modelgrd$Column_,ModelGridCoords)
  
  return(ModelGridCoords)
}

#-------------------------------------------------
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
#-------------------------------------------------
HARNSP17ft  <- CRS("+init=epsg:2881")
HARNUTM17Nm  <- CRS("+init=epsg:3747")
latlongs <- CRS("+proj=longlat +datum=WGS84")

#-------------------------------------------------
# Provide default basepath with option to change
#-------------------------------------------------
basePath <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS/"
basePath <-  tk_choose.dir(default = basePath, 
                           caption = "Select input directory for biasNRD files")

#-------------------------------------------------
#  Prepare environment to support multiprocessing
#-------------------------------------------------
future::plan(multiprocess)

#-------------------------------------------------
#  Select input and Read one year for Pixels
#  using a background process
#-------------------------------------------------

#csvFile <- paste0(basePath, paste0("biasNRD2003.csv"))
csvFile <- utils::choose.files(default=paste(basePath,'*.csv',sep='/'))
cat(paste('"+" indicates Progress Reading',csvFile,'\n'))
f <- future({readPoints(csvFile)})

#-------------------------------------------------
# Read Model grid shapefile as polys 
# and convert to points and convert to data frame
# using a background process  
# - readgridPoints() runs async to readPoints()
#-------------------------------------------------
Modelgrd.Path <- 
  "//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/CFWI/Data/From_SW_SJ"
Model.Shape <-"ECFTX_GRID_V3"
setwd(Modelgrd.Path)

cat(paste('":" indicates Progress Reading Model shapefile',Model.Shape,'\n'))
g <- future({readgridPoints(Modelgrd.Path,Model.Shape)})

#-------------------------------------------------
# Wait for values from futures
#-------------------------------------------------
cat(paste('Waiting for background processing to complete','\n'))

while (!resolved(g)){
  if (!resolved(f)){
  cat("+")
  }
  cat(":")
}
cat("\n")

PixelCoords <-future::value(f)
ModelGridCoords <-future::value(g)

cat(paste('Background processing is complete','\n'))

#-------------------------------------------------
# Find closest points using nearest neighbor
#-------------------------------------------------
NearNeighbor <- FNN::get.knnx(PixelCoords[,2:3],ModelGridCoords[,3:4],1)
ClosestPixels <- NearNeighbor[["nn.index"]]
ClosestDistance <- NearNeighbor[["nn.dist"]]
ModelGridCoords <- cbind(PixelCoords[ClosestPixels,1],ClosestDistance,ModelGridCoords)

#-------------------------------------------------
# Output data
#-------------------------------------------------
ModelGCdf <- as.data.frame(ModelGridCoords)
head(ModelGCdf)
names(ModelGCdf) <- c("PixelID","distance","row","col","X","Y")
ModelCellfile <- paste(Modelgrd.Path,'ModelCells.csv',sep='/')
utils::write.csv(as.data.frame(ModelGCdf),ModelCellfile,row.names =FALSE)
cat(paste('Crosswalk/Lookup table for NexRad Pixel to Model ROw and Column','\n'))
cat(paste('saved to: ',ModelCellfile,'\n'))

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
