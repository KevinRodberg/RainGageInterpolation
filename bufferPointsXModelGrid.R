
library(reshape2)
library(readr)
library(dplyr)
library(data.table)
library(readxl)
library(rgeos)
library(sp)
library(maptools)
library(sf)
library(FNN)


basePAth <-  "//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/"

#-------------------------------------------------
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
#-------------------------------------------------
HARNSP17ft  = CRS("+init=epsg:2881")

#-------------------------------------------------
# Read one year for Pixels
#-------------------------------------------------
LWC_NRD1999 <- read.csv(paste0(basePAth, "LWC_NRD1999.csv"))
NRD1999 <- na.omit(melt(LWC_NRD1999, id = c("Pixel_id", "X", "Y")))
PixelCoords <- NRD1999[c("Pixel_id", "X", "Y")]
PixelPoints<-SpatialPointsDataFrame(coords = PixelCoords[, c("X", "Y")],
                                    data = PixelCoords,proj4string = HARNSP17ft)

#-------------------------------------------------
# Read Model grid shapefile as polys 
# and convert to points
#-------------------------------------------------
LWCSIMgrd.Path <- "//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/Model_shapefiles"
LWCSIM.Shape <- "LWCSIM_Model_Grid.shp"
LWCSIMgrd.Path <- "//ad.sfwmd.gov/dfsRoot/data/wsd/MOD/LWCSASIAS/Model/Felipe"
LWCSIM.Shape <-"LWCSIM_IBOUND1_Layer01_Cells.shp"
setwd(LWCSIMgrd.Path)
LWCSIMgrd <- readShapePoly(LWCSIM.Shape, proj4string = HARNSP17ft)
gridCentroids <-gCentroid(LWCSIMgrd,byid=TRUE)

#-------------------------------------------------
# Convert Pixel shapefile to dataframe
#-------------------------------------------------
asSFPP <- st_as_sf(PixelPoints)
PixelCoordinates <- do.call(rbind,st_geometry(asSFPP))

#-------------------------------------------------
# Convert model cell points to a dataframe
#-------------------------------------------------
asSFGC <- st_as_sf(gridCentroids)
ModelGridCoords <- do.call(rbind,st_geometry(asSFGC))
ModelGridCoords <-cbind(LWCSIMgrd$ROW, LWCSIMgrd$COL,ModelGridCoords)

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
# returns index and sistance of x1 closest to x2
#-------------------------------------------------
#   nn[["nn.index"]]
#   nn[["nn.dist"]]


#-------------------------------------------------
# Find closest points
#-------------------------------------------------
NearNeighbor = get.knnx(PixelCoordinates,ModelGridCoords[,3:4],1)
ClosestPixels <- NearNeighbor[["nn.index"]]
ClosestDistance <-NearNeighbor[["nn.dist"]]
ModelGridCoords <-cbind(PixelCoords[ClosestPixels,1],ClosestDistance,ModelGridCoords)

#-------------------------------------------------
# Output data
#-------------------------------------------------
ModelGCdf <- as.data.frame(ModelGridCoords)
names(ModelGCdf)<- c("PixelID","distance","row","col","X","Y")
ModelCellfile<-paste0(LWCSIMgrd.Path,'ModelCells.csv')
write_csv(as.data.frame(ModelGCdf),ModelCellfile)

