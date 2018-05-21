#-------------------------------------------------
#
# Script:     ECFTPrepDataUsingBiasMP.R
#   based upon      PrepDataUsingBiasMP.R
#
# Programmer: Kevin A. Rodberg
# Date:       05/21/2018
#
# Implemented with Multiprocessing functionality
# Execution time approx 10 minutes. with 4 processors
# 
# Current configuration has been set up and tested for 2003.
#
#->>> Not there are some missing pixels from input NexRad data.   <<<<<<<<<<<
#
#   These haven't been accounted for.
#
# The program:
# Calculates annual bias multipliers from RainGage vs NexRad Pixels
# Bias factors at RainGages are interpolated using Ordinary Kriging
# producing an annual 'bias' raster.  Daily NexRad Pixels are rasterized
# and multiplid by the 'bias' raster to produce 'bias adjusted' Daily NexRad 
# rasters, which are converted back to pixels by extracting raster values
# from the rasters.  Daily pixels values by row are combined as columns 
# and exported to csv with an Annual total column added.
# Raster plots of Annual data are also produced.
#-------------------------------------------------

list.of.packages <-  c("reshape2","readr","dplyr", "data.table","readxl",
                       "rgeos","sp","dismo","lattice","rasterVis","maptools",
                       "raster","fields","automap","gstat","future","listenv")

new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

if (length(new.packages))
  install.packages(new.packages)

library(reshape2)
#library(readr)
library(dplyr)
library(data.table)
#library(readxl)
library(rgeos)
library(sp)
#library(dismo)
library(lattice)
library(rasterVis)
library(maptools)
library(raster)
#library(fields)
library(automap)
library(latticeExtra)
#library(gstat)
library(future)
library(listenv)

myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))

basePath <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS/"

NRDpixels <- 
    read.csv("//whqhpc01p/hpcc_shared/krodberg/NexRadTS/Rain/PixelsByWMD2016.csv")
names(NRDpixels)

#-------------------------------------------------
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
#-------------------------------------------------
HARNSP17ft  = CRS("+init=epsg:2881")
HARNUTM17Nm  = CRS("+init=epsg:3747")
latlongs = CRS("+proj=longlat +datum=WGS84")
PixelCoords <- NRDpixels[c("PixelID", "Longitude", "Latitude")]
names(PixelCoords) <-c("Pixel_id", "X", "Y")
coordinates(PixelCoords) <-~X+Y
proj4string(PixelCoords) <-latlongs
plot(PixelCoords)

#-------------------------------------------------
# Add additional melted NRD data sets 
# to calculate bias based on more than one year
#-------------------------------------------------
#NRD <- do.call("rbind", list(NRD2016))

#-------------------------------------------------
# Set up county boundry shapefile for overlay 
# on raster maps
#-------------------------------------------------
gClip <- function(shp, bb) {
  if (class(bb) == "matrix")
    b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else
    b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

WMDbnd.Path <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS"
WMDbnd.Shape <- "CntyBnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape, proj4string = HARNSP17ft)

#-------------------------------------------------
# Calculate raster extents
#-------------------------------------------------
data.proj <-spTransform(PixelCoords, HARNSP17ft)
#xmin = floor(min(PixelCoords[c('X')]))
#xmax = ceiling(max(PixelCoords[c('X')]))
#ymin = floor(min(PixelCoords[c('Y')]))
#ymax = ceiling(max(PixelCoords[c('Y')]))
xmin = floor(xmin(data.proj))
xmax = ceiling(xmax(data.proj))
ymin = floor(ymin(data.proj))
ymax = ceiling(ymax(data.proj))

#-------------------------------------------------
# calculate number of rows and columns
# 6561.679 = feet in 2 Kilometer pixel spacing
#-------------------------------------------------
rasRows <- round(((ymax - ymin) / 6561.679),0)
rasCols <- round(((xmax - xmin) / 6561.679),0)

#-------------------------------------------------
# define raster and map extents using NRD pixel data extents
#-------------------------------------------------
ras <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,
               ymn=ymin,ymx=ymax,crs=HARNSP17ft)
rasExt <- extent(ras)
clpBnds2 <- gClip(WMDbnd, rasExt)
#-------------------------------------------------
# FUNCTION: dayBiasFn
# Multiplies Daily NRD rasters by bias  
#   calculating adjusted daily NexRad Raster
# and returns 
#-------------------------------------------------
dayBiasFn <- function(DailyNRD,biasRas){
  DailyNRD.pnts <-SpatialPointsDataFrame(coords = DailyNRD[, c("Longitude", "Latitude")],
                                         data = DailyNRD,proj4string = latlongs)
  DailyNRD.pnts <- spTransform(DailyNRD.pnts,HARNSP17ft)
  NRDras <-rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = mean) * biasRas
  NRDBiasPnts <- extract(NRDras,DailyNRD.pnts,fun=mean,df=TRUE)
  NRDBiasPnts$Pixel_id <- DailyNRD.pnts$PixelID
#  NRDBiasPnts$X<-DailyNRD.pnts$X
#  NRDBiasPnts$Y<-DailyNRD.pnts$Y
  return(NRDBiasPnts[c(3,2)])
}

#-------------------------------------------------
# FUNCTION: biasByYear
#   [Works well with future function for multiprocessing]
# Processes RainVsGage data by year
# Creating CSV files with NRDstat, biasNRD and 
# updates NRD with annual totals
#-------------------------------------------------
yearStr = '2003'
biasByYear <-function(yearStr,RvsG){
  ECFTX_NRDbyYr <- read.csv(paste0(basePath, "nr", yearStr, ".csv"))
  justPixels <-ECFTX_NRDbyYr[, -which(names(ECFTX_NRDbyYr) 
                                %in% c("SEQNUM","Row","Column_","ROWCOL"))]
  uniquePixels <- justPixels[!duplicated(justPixels[1:2]),]
  names(uniquePixels)
  dateList<-names(uniquePixels)[-c(1,2)]


  #-------------------------------------------------
  #  Define SignificantDays to
  #  Filter out Raing Gage data with fewer than 90%
  #-------------------------------------------------
#  TotDays = endDay - startDay + 1
#  SignificantDays = as.numeric(.9 * TotDays)
  
#  RvsG[c("biasNRD")] <- NA
#  RainStats <- NULL
#  stationList = unique(RvsG$RainGage)
  RGpath <- paste0("//ad.sfwmd.gov/dfsRoot/data/wsd/GIS/GISP_2012/DistrictAreaProj/",
                   "CFWI/Data/RainGageVsNexRad/"  )
  RGfile <-"AdjFact_2003.csv"
  RGfilepath <-paste0(RGpath,RGfile)
  #  RGdata <- read.csv(paste0(basePath, "nr", yearStr, ".csv"))
  RGdata <- read.csv(RGfilepath)
  names(RGdata)
  stationList = unique(RGdata$RGAGE_ID)
  
  NRDwCoords<-inner_join(NRDpixels,uniquePixels, by= c("PixelID"= "Pixel"))

  names(NRDwCoords)
  NRDbyYr <- melt(NRDwCoords, id = c("X", "PixelID", "Latitude", "Longitude", "WMD","District"))
  NRDbyYr$X <-NULL
  NRDbyYr$District <-NULL
  
  #-------------------------------------------------
  # read and organize daily NexRad data
  #-------------------------------------------------
  for (d in dateList[1]) {
    DailyNRD <- NRDbyYr[NRDbyYr$variable == d, ]
    names(DailyNRD)
    DailyNRD.pnts <-SpatialPointsDataFrame(coords = DailyNRD[, c("Longitude", "Latitude")],
                                           data = DailyNRD, 
                                           proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    DailyNRD.pnts <- spTransform(DailyNRD.pnts,HARNSP17ft)
    NRDras <-rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = mean) 
    plot(NRDras)
    NRDBiasPnts <- data.frame(extract(NRDras,DailyNRD.pnts))
    DailyNRD.grid <- as(NRDras, "SpatialGridDataFrame")
  }
  NRDPixel<- data.frame(DailyNRD.pnts$PixelID,DailyNRD.pnts$Longitude,DailyNRD.pnts$Latitude)
  
  #-------------------------------------------------
  # Interpolate Bias from RainGages to NexRad pixels
  # Make NRD data correction using Bias
  #-------------------------------------------------
#  biasData <- na.omit(combineData[, c("XCOORD", "YCOORD", "bias")])
  biasData <- na.omit(RGdata[, c("X", "Y", "ADJ_FACT")])
  
  coordinates(biasData) =  ~ X + Y
  proj4string(biasData) = HARNUTM17Nm
  
  biasData$XCOORD <- coordinates(biasData)[, 1]
  biasData$YCOORD <- coordinates(biasData)[, 2]
  biasData <- spTransform(biasData,HARNSP17ft)
  names(biasData)
  rainGage.pnts <-
    SpatialPointsDataFrame(coords = RGdata[, c("X", "Y")],
                           data = RGdata,proj4string = HARNUTM17Nm)
  rainGage.pnts <-spTransform(rainGage.pnts,HARNSP17ft)
  NRDras <- rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = mean)
  #-------------------------------------------------
  #  Theisen Polygon and raster code:
  # theisPoly <- voronoi(rainGage.pnts)
  # TheisRas <- rasterize(theisPoly, NRDras, theisPoly$bias, fun = mean)
  
  myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))
  
  #-------------------------------------------------
  #  autoKrige implemented for Ordinary kriging
  #-------------------------------------------------
  surf <- autoKrige(formula=ADJ_FACT ~ 1, input_data=biasData, new_data = DailyNRD.grid)
  biasRas <- raster(surf$krige_output)
  plot(biasRas)
  #-------------------------------------------------
  #  IDW raster code:
  # gs <- gstat(formula=bias~1, locations=biasData)
  # idw <- interpolate(NRDras, gs)

  NRDbiasPixels <-NRDPixel
  names(NRDbiasPixels)<-c( "Pixel_id" ,"X" ,"Y")
  #-------------------------------------------------
  #  Process NexRad for each day calling "dayBiasFn"
  #  NRD with bias results are merged into a single table 
  #-------------------------------------------------
  for (d in dateList) {
    #   cat(paste("."))
    table = dayBiasFn(NRDbyYr[NRDbyYr$variable == d, ],biasRas)
    names(table) <- c( "Pixel_id" ,d  )
    NRDbiasPixels <- merge(NRDbiasPixels, table, by =c("Pixel_id"))
  }
  
  #-------------------------------------------------
  # Add final column for annual total NRD with Bias correction
  # and export to csv
  #-------------------------------------------------
  NRDbiasPixels$Annual<-rowSums(NRDbiasPixels[,-c(1,2,3)],na.rm=TRUE)
  csvFile <- paste0(basePath, paste0("biasNRD",yearStr,".csv"))
  fwrite(NRDbiasPixels, csvFile) 
  
  #-------------------------------------------------
  # Create Raster for bias corrected Annual NexRAD rain
  #-------------------------------------------------
  Annuals <-NRDbiasPixels[,c(1,2,3,ncol(NRDbiasPixels))]
  names(Annuals)
  xy <- Annuals[,c(2,3)]
  Annualspdf <-SpatialPointsDataFrame(coords=xy,data=Annuals,proj4string=HARNSP17ft)
  AnnRas <- rasterize(Annualspdf,ras,Annualspdf$Annual,fun=mean)

  #-------------------------------------------------
  # Add final column for annual total NRD
  # and export to csv
  #-------------------------------------------------

  NRDwCoords$Annual<-rowSums(NRDwCoords[,-c(1,2,3,4,5,6)],na.rm=TRUE)
  csvFile <- paste0(basePath, paste0("NRD",yearStr,".csv"))
  fwrite(NRDwCoords, csvFile) 
  
  #-------------------------------------------------
  # Create Raster for Annual NEXRad rain
  #-------------------------------------------------
  AnnNRDs<-NRDwCoords[,c(2,3,4,ncol(NRDwCoords))]
  xy <- AnnNRDs[,c(3,2)]
  AnnNRDsspdf <-SpatialPointsDataFrame(coords=xy,data=AnnNRDs,proj4string=latlongs)
  AnnNRDsspdf <-spTransform(AnnNRDsspdf,HARNSP17ft)
  NRDAnnRas <- rasterize(AnnNRDsspdf,ras,AnnNRDsspdf$Annual,fun=mean)
  plot(NRDAnnRas)
  #-------------------------------------------------
  #  Return a list of 4 spatial objects 
  #     (3 rasters & 1 points)
  #-------------------------------------------------
  biasStuff <-list("B_ras"=biasRas,
                   "R_pnts"=rainGage.pnts,
                   "annualRas"=AnnRas,
                   "annNRDRas"=NRDAnnRas )
  return(biasStuff)
}

#-------------------------------------------------
# Set up for multiprocessing function calls
#-------------------------------------------------
plan(multiprocess)
processed= listenv(NULL)
yrList=list()

yearStr <- as.character(2003)
#-------------------------------------------------
# Define range of years to process
#-------------------------------------------------
#processYears <- seq(1999, 2016)
processYears <- seq(2003, 2003)
x=0

for (yr in processYears) {
  yearStr <- as.character(yr)
  x=x+1
  #---------------------------
  # Compare NexRad vs RainGage
  #---------------------------
  startDay = as.Date(paste0("01/01/",yearStr), "%m/%d/%Y")
  endDay = as.Date(paste0("12/31/",yearStr), "%m/%d/%Y")
#  RvsG <- as.data.frame(NRDvsGage[NRDvsGage$Rainfall < 20 &
#                                             as.Date(NRDvsGage$DailyDate) >= startDay &
#                                             as.Date(NRDvsGage$DailyDate) <= endDay, ])
  cat (paste(yearStr,"\n"))
  #-------------------------------------------------
  # Call FUNCTION "biasByYear" with futures multiprocessing
  # wrapper function
  #-------------------------------------------------
  processed[[x]] <- future({biasByYear(yearStr,RvsG)})
}

mpList<-list()
rList<-list()
stackList <-list()

#-------------------------------------------------
# value function waits for results to become available
# for each process
#-------------------------------------------------
for (i in seq(1:x)){
  mpList <-value(processed[[i]])
  rList[[i]]<-mpList
}

#-------------------------------------------------
# unlist results returned from FUNCTION "biasByYear"
#-------------------------------------------------
stackList = unlist(lapply(rList,"[[",1))
pointList = unlist(lapply(rList,"[[",2))
AnnRainList = unlist(lapply(rList,"[[",3))
AnnNRDList = unlist(lapply(rList,"[[",4))

rasStack <-stack()
rasStack <-stack(stackList)
AnnrasStack <-stack()
AnnrasStack <-stack(AnnRainList)
AnnNRDStack <-stack()
AnnNRDStack <-stack(AnnNRDList)

yearList <-as.character(processYears)
names(rasStack)<- yearList
names(AnnrasStack)<- yearList
names(AnnNRDStack)<- yearList

#-------------------------------------------------
# Create plot files for each raster type by year
#-------------------------------------------------
for (i in 1:nlayers(rasStack)){
  filename=paste(basePath,"/rasterPlots/bias",names(rasStack)[i],".png",sep="")
  myplot=(  levelplot(rasStack[[i]], par.settings = myTheme, main=names(rasStack)[i], 
                      at=seq(.8,1.3,length=18),layout=c(1,1),contour=FALSE) + 
              layer(sp.polygons(clpBnds2)) +
              layer(sp.points(pointList[[i]], col = "red")))
  labels <- layer(sp.text(coordinates(pointList[[i]]),
                  txt=pointList[[i]]$RGAGE_ID,pos=1,cex=.5 ))
  
  print (myplot+labels)
  dev.copy(png, filename, width=2400,height=2400,units="px",res=300)
  dev.off()
  #   Plot Annaul Adjusted NRD
  filename=paste(basePath,"/rasterPlots/AnnNRDwBias",names(AnnrasStack)[i],".png",sep="")
  myplot=(  levelplot(AnnrasStack[[i]], par.settings = myTheme, main=names(AnnrasStack)[i], 
                      at=seq(30,90,length=18),layout=c(1,1),contour=FALSE) + 
              layer(sp.polygons(clpBnds2)) +
              layer(sp.points(pointList[[i]], col = "red")))
  labels <- layer(sp.text(coordinates(pointList[[i]]),
                  txt=pointList[[i]]$RGAGE_ID,pos=1,cex=.5 ))
  
  print (myplot+labels)
  dev.copy(png, filename, width=2400,height=2400,units="px",res=300)
  dev.off()
  # Plot Annual NRD
  filename=paste(basePath,"/rasterPlots/AnnNRD",names(AnnNRDStack)[i],".png",sep="")
  myplot=(  levelplot(AnnNRDStack[[i]], par.settings = myTheme, main=names(AnnNRDStack)[i], 
                      at=seq(30,90,length=18),layout=c(1,1),contour=FALSE) + 
              layer(sp.polygons(clpBnds2)) +
              layer(sp.points(pointList[[i]], col = "red")))
  labels <- layer(sp.text(coordinates(pointList[[i]]),
                  txt=pointList[[i]]$RGAGE_ID,pos=1,cex=.5 ))
  
  print (myplot+labels)
  dev.copy(png, filename, width=2400,height=2400,units="px",res=300)
  dev.off()
}