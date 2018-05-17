#-------------------------------------------------
#
# Script:     PrepDataUsingBiasMP.R
#
# Programmer: Kevin A. Rodberg
# Date:       05/11/2018
#
# Implemented with Multiprocessing functionality
# Execution time approx 10 minutes. with 4 processors
# 
# Calculates annual bias multipliers from RainGage vs NexRad Pixels
# Bias factors at RainGages are interpolated using Ordinary Kriging
# producing an annual 'bias' raster.  Daily NexRad Pixels are rasterized
# and multiplid by the 'bias' raster to produce 'bias adjusted' Daily NexRad 
# rasters, which are converted back to pixels by extracting raster values
# from the rasters.  Daily pixels values by row are combined as columns 
# and exported to csv with an Annual total column added.
# Raster plots of Annual data are also produced.
#-------------------------------------------------

list.of.packages <-  c("reshape2","readr","dplyr", "data.table","readxl","rgeos","sp","dismo","lattice",
                       "rasterVis","maptools","raster","fields","automap","gstat","future","listenv")

new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

if (length(new.packages))
  install.packages(new.packages)

library(reshape2)
library(readr)
library(dplyr)
library(data.table)
library(readxl)
library(rgeos)
library(sp)
library(dismo)
library(lattice)
library(rasterVis)
library(maptools)
library(raster)
library(fields)
library(automap)
library(latticeExtra)
library(gstat)
library(future)
library(listenv)

myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))

basePAth <-  "//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/"
LWCNRDvsGage99to16 <-  read_excel(paste0(basePAth, "Compared_/LWC_NRDvsGageLY1999-2016.xlsx"))

LWC_NRD2016 <- read.csv(paste0(basePAth, "LWC_NRD2016.csv"))
NRD2016 <- na.omit(melt(LWC_NRD2016, id = c("Pixel_id", "X", "Y")))
PixelCoords <- NRD2016[c("Pixel_id", "X", "Y")]

#-------------------------------------------------
# Add additional melted NRD data sets 
# to calculate bias based on more than one year
#-------------------------------------------------
NRD <- do.call("rbind", list(NRD2016))

#-------------------------------------------------
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
#-------------------------------------------------
HARNSP17ft  = CRS("+init=epsg:2881")

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
# Read Model grid shapefile
#-------------------------------------------------
#LWCSIMgrd.Path <- "//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/Model_shapefiles"
#LWCSIM.Shape <- "LWCSIM_Model_Grid.shp"
#setwd(LWCSIMgrd.Path)
#LWCSIMgrd <- readShapePoly(LWCSIM.Shape, proj4string = HARNSP17ft)
#-------------------------------------------------
# Calculate raster extents
#-------------------------------------------------
xmin = floor(min(PixelCoords[c('X')]))
xmax = ceiling(max(PixelCoords[c('X')]))
ymin = floor(min(PixelCoords[c('Y')]))
ymax = ceiling(max(PixelCoords[c('Y')]))

#-------------------------------------------------
# calculate number of rows and columns
# 6561.679 = feet in 2 Kilometer pixel spacing
#-------------------------------------------------
rasRows <- (ymax - ymin) / 6561.679
rasCols <- (xmax - xmin) / 6561.679

#-------------------------------------------------
# define raster and map extents using NRD pixel data extents
#-------------------------------------------------
ras <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=HARNSP17ft)
rasExt <- extent(ras)
clpBnds2 <- gClip(WMDbnd, rasExt)

#-------------------------------------------------
# FUNCTION: dayBiasFn
# Multiplies Daily NRD rasters by bias  
#   calculating adjusted daily NexRad Raster
# and returns 
#-------------------------------------------------
dayBiasFn <- function(DailyNRD,biasRas){
  DailyNRD.pnts <-SpatialPointsDataFrame(coords = DailyNRD[, c("X", "Y")],
                                         data = DailyNRD,proj4string = HARNSP17ft)
  NRDras <-rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = max) * biasRas
  NRDBiasPnts <- extract(NRDras,DailyNRD.pnts,fun=max,df=TRUE)
  NRDBiasPnts$Pixel_id <- DailyNRD.pnts$Pixel_id
  NRDBiasPnts$X<-DailyNRD.pnts$X
  NRDBiasPnts$Y<-DailyNRD.pnts$Y
  return(NRDBiasPnts[c(3,4,5,2)])
}

#-------------------------------------------------
# FUNCTION: biasByYear
#   [Works well with future function for multiprocessing]
# Processes RainVsGage data by year
# Creating CSV files with NRDstat, biasNRD and 
# updates NRD with annual totals
#-------------------------------------------------
biasByYear <-function(yearStr,RvsG){
  LWC_NRDbyYr <- read.csv(paste0(basePAth, "LWC_NRD", yearStr, ".csv"))
  dateList<-names(LWC_NRDbyYr)[-c(1,2,3)]
  NRDbyYr <- melt(LWC_NRDbyYr, id = c("Pixel_id", "X", "Y"))

  #-------------------------------------------------
  #  Define SignificantDays to
  #  Filter out Raing Gage data with fewer than 90%
  #-------------------------------------------------
  TotDays = endDay - startDay + 1
  SignificantDays = as.numeric(.9 * TotDays)
  
  RvsG[c("biasNRD")] <- NA
  RainStats <- NULL
  stationList = unique(RvsG$RainGage)
  
  #-------------------------------------------------
  #  Calculate bias at each rainGage vs NRD Pixel
  #  where NRD data > 0  --->>>>  may need to be .25
  #  Ignore bias factors < .8 and  > 1.2
  #-------------------------------------------------
  for (gage in stationList) {
    if ((nrow(na.omit(RvsG[RvsG$RainGage == gage & RvsG$NRD >= 0.0, -c(14)]))) > 0) {
      nobsSigRain <-nrow(na.omit(RvsG[RvsG$RainGage == gage & RvsG$NRD >= 0.0, -c(14)]))
      nobsNRD <- nrow(RvsG[RvsG$RainGage == gage & !is.na(RvsG$NRD), ])
      nobsRain <-nrow(RvsG[RvsG$RainGage == gage & !is.na(RvsG$Rainfall), ])
      meanNRD <-mean(as.numeric(RvsG[RvsG$RainGage == gage & !is.na(RvsG$Rainfall), 11]),na.rm=TRUE)
      meanRain <-mean(as.numeric(RvsG[RvsG$RainGage == gage & !is.na(RvsG$Rainfall), 10]),na.rm=TRUE)
      bias = meanRain / meanNRD
      if (bias > .8 & bias < 1.2 &  nobsRain > SignificantDays) {
        RvsG[RvsG$RainGage == gage &!is.na(RvsG$NRD), 14] <- 
          na.omit((RvsG[RvsG$RainGage == gage  & !is.na(RvsG$NRD), 11]) * bias)
        
        ValsList = data.frame(gage,as.numeric(bias),as.numeric(nobsSigRain),nobsNRD,nobsRain)
        RainStats <- rbind.data.frame(RainStats, ValsList)
      }
    }
  }
  
  
  #-------------------------------------------------
  # Calculate Rain Gage vs NexRad pixel Stats by year
  #-------------------------------------------------
  names(RainStats) <-c("RainGage", "bias", "SigRainObs", "NRDobs", "RainObs")
  NRDSum <-aggregate(RvsG$NRD,list(RvsG$RainGage, RvsG$XCOORD, RvsG$YCOORD),FUN = sum)
  names(NRDSum) <- c("RainGage", "XCOORD", "YCOORD", "NRD")
  
  GageSum <- aggregate(RvsG$Rainfall,list(RvsG$RainGage, RvsG$XCOORD, RvsG$YCOORD),FUN = sum)
  names(GageSum) <- c("RainGage", "XCOORD", "YCOORD", "Rainfall")
  BiasSum <- aggregate(RvsG$biasNRD,list(RvsG$RainGage, RvsG$XCOORD, RvsG$YCOORD),FUN = sum)
  names(BiasSum) <- c("RainGage", "XCOORD", "YCOORD", "biasNRD")
  
  compareSums <- merge(merge(BiasSum, NRDSum, all = TRUE), GageSum, all =TRUE)
  combineData <- plyr::join(compareSums,RainStats,by = 'RainGage',type = 'left',match = 'all')
  combineData <- na.omit(combineData)
  
  #-------------------------------------------------
  # Export Stats to csv files by year
  #-------------------------------------------------
  csvFile <- paste0(basePAth, paste0("newNRDstats",yearStr,".csv"))
  fwrite(na.omit(combineData), csvFile)
  
  #-------------------------------------------------
  # read and organize daily NexRad data
  #-------------------------------------------------
  for (d in dateList[1]) {
    DailyNRD <- NRDbyYr[NRDbyYr$variable == d, ]
    DailyNRD.pnts <-SpatialPointsDataFrame(coords = DailyNRD[, c("X", "Y")],
                                           data = DailyNRD,proj4string = HARNSP17ft)
    
    NRDras <-rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = max) 
    NRDBiasPnts <- data.frame(extract(NRDras,DailyNRD.pnts))
    DailyNRD.grid <- as(NRDras, "SpatialGridDataFrame")
  }
  NRDPixel<- data.frame(DailyNRD.pnts$Pixel_id,DailyNRD.pnts$X,DailyNRD.pnts$Y)
  
  #-------------------------------------------------
  # Interpolate Bias from RainGages to NexRad pixels
  # Make NRD data correction using Bias
  #-------------------------------------------------
  biasData <- na.omit(combineData[, c("XCOORD", "YCOORD", "bias")])
  
  coordinates(biasData) =  ~ XCOORD + YCOORD
  proj4string(biasData) = HARNSP17ft
  
  biasData$X <- coordinates(biasData)[, 1]
  biasData$Y <- coordinates(biasData)[, 2]
  
  rainGage.pnts <-
    SpatialPointsDataFrame(coords = combineData[, c("XCOORD", "YCOORD")],
                           data = combineData,proj4string = HARNSP17ft)
  
  NRDras <- rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = max)
  
  #-------------------------------------------------
  #  Theisen Polygon and raster code:
  # theisPoly <- voronoi(rainGage.pnts)
  # TheisRas <- rasterize(theisPoly, NRDras, theisPoly$bias, fun = mean)
  
  myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))
  
  #-------------------------------------------------
  #  autoKrige implemented for Ordinary kriging
  #-------------------------------------------------
  surf <- autoKrige(formula=bias ~ 1, input_data=biasData, new_data = DailyNRD.grid)
  biasRas <- raster(surf$krige_output)

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
    names(table) <- c( "Pixel_id" ,"X" ,"Y",d  )
    NRDbiasPixels <- merge(NRDbiasPixels, table, by =c("Pixel_id","X","Y"))
  }
  
  #-------------------------------------------------
  # Add final column for annual total NRD with Bias correction
  # and export to csv
  #-------------------------------------------------
  NRDbiasPixels$Annual<-rowSums(NRDbiasPixels[,-c(1,2,3)],na.rm=TRUE)
  csvFile <- paste0(basePAth, paste0("biasNRD",yearStr,".csv"))
  fwrite(NRDbiasPixels, csvFile) 
  
  #-------------------------------------------------
  # Create Raster for bias corrected Annual NexRAD rain
  #-------------------------------------------------
  Annuals <-NRDbiasPixels[,c(1,2,3,ncol(NRDbiasPixels))]
  xy <- Annuals[,c(2,3)]
  Annualspdf <-SpatialPointsDataFrame(coords=xy,data=Annuals,proj4string=HARNSP17ft)
  AnnRas <- rasterize(Annualspdf,ras,Annualspdf$Annual,fun=mean)

  #-------------------------------------------------
  # Add final column for annual total NRD
  # and export to csv
  #-------------------------------------------------
  LWC_NRDbyYr$Annual<-rowSums(LWC_NRDbyYr[,-c(1,2,3)],na.rm=TRUE)
  csvFile <- paste0(basePAth, paste0("NRD",yearStr,".csv"))
  fwrite(LWC_NRDbyYr, csvFile) 
  
  #-------------------------------------------------
  # Create Raster for Annual NEXRad rain
  #-------------------------------------------------
  AnnNRDs<-LWC_NRDbyYr[,c(1,2,3,ncol(LWC_NRDbyYr))]
  xy <- AnnNRDs[,c(2,3)]
  AnnNRDsspdf <-SpatialPointsDataFrame(coords=xy,data=AnnNRDs,proj4string=HARNSP17ft)
  NRDAnnRas <- rasterize(AnnNRDsspdf,ras,AnnNRDsspdf$Annual,fun=mean)
  
  #-------------------------------------------------
  #  Return a list of 4 spatial objects 
  #     (3 rasters & 1 points)
  #-------------------------------------------------
  biasStuff <-list("B_ras"=biasRas,"R_pnts"=rainGage.pnts,"annualRas"=AnnRas,"annNRDRas"=NRDAnnRas )
  return(biasStuff)
}

#-------------------------------------------------
# Set up for multiprocessing function calls
#-------------------------------------------------
plan(multiprocess)
processed= listenv(NULL)
yrList=list()

yearStr <- as.character(1999)
#-------------------------------------------------
# Define range of years to process
#-------------------------------------------------
processYears <- seq(1999, 2016)
x=0

for (yr in processYears) {
  yearStr <- as.character(yr)
  x=x+1
  #---------------------------
  # Compare NexRad vs RainGage
  #---------------------------
  startDay = as.Date(paste0("01/01/",yearStr), "%m/%d/%Y")
  endDay = as.Date(paste0("12/31/",yearStr), "%m/%d/%Y")
  RvsG <- as.data.frame(LWCNRDvsGage99to16[LWCNRDvsGage99to16$Rainfall < 20 &
                                             as.Date(LWCNRDvsGage99to16$DailyDate) >= startDay &
                                             as.Date(LWCNRDvsGage99to16$DailyDate) <= endDay, ])
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
  filename=paste(basePAth,"/rasterPlots/bias",names(rasStack)[i],".png",sep="")
  myplot=(  levelplot(rasStack[[i]], par.settings = myTheme, main=names(rasStack)[i], 
                      at=seq(.8,1.2,length=18),layout=c(1,1),contour=FALSE) + 
              layer(sp.polygons(clpBnds2)) +
              layer(sp.points(pointList[[i]], col = "red")))
  labels <- layer(sp.text(coordinates(pointList[[i]]),txt=pointList[[i]]$RainGage,pos=1,cex=.5 ))
  
  print (myplot+labels)
  dev.copy(png, filename, width=2400,height=2400,units="px",res=300)
  dev.off()
  #   Plot Annaul Adjusted NRD
  filename=paste(basePAth,"/rasterPlots/AnnNRDwBias",names(AnnrasStack)[i],".png",sep="")
  myplot=(  levelplot(AnnrasStack[[i]], par.settings = myTheme, main=names(AnnrasStack)[i], 
                      at=seq(30,90,length=18),layout=c(1,1),contour=FALSE) + 
              layer(sp.polygons(clpBnds2)) +
              layer(sp.points(pointList[[i]], col = "red")))
  labels <- layer(sp.text(coordinates(pointList[[i]]),txt=pointList[[i]]$RainGage,pos=1,cex=.5 ))
  
  print (myplot+labels)
  dev.copy(png, filename, width=2400,height=2400,units="px",res=300)
  dev.off()
  # Plot Annual NRD
  filename=paste(basePAth,"/rasterPlots/AnnNRD",names(AnnNRDStack)[i],".png",sep="")
  myplot=(  levelplot(AnnNRDStack[[i]], par.settings = myTheme, main=names(AnnNRDStack)[i], 
                      at=seq(30,90,length=18),layout=c(1,1),contour=FALSE) + 
              layer(sp.polygons(clpBnds2)) +
              layer(sp.points(pointList[[i]], col = "red")))
  labels <- layer(sp.text(coordinates(pointList[[i]]),txt=pointList[[i]]$RainGage,pos=1,cex=.5 ))
  
  print (myplot+labels)
  dev.copy(png, filename, width=2400,height=2400,units="px",res=300)
  dev.off()
}