list.of.packages <-  c("reshape2","readr","dplyr", "data.table","readxl","rgeos","sp","dismo","lattice",
                       "rasterVis","maptools","raster","fields","automap","gstat","future")

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
detach("package:ggplot2", unload=TRUE)

basePAth <-  "//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/"
LWCNRDvsGage99to16 <-  read_excel(paste0(basePAth, "Compared_/LWC_NRDvsGageLY1999-2016.xlsx"))

LWC_NRD2016 <- read.csv(paste0(basePAth, "LWC_NRD2016.csv"))
NRD2016 <- na.omit(melt(LWC_NRD2016, id = c("Pixel_id", "X", "Y")))
PixelCoords <- NRD2016[c("Pixel_id", "X", "Y")]

NRD <- do.call("rbind", list(NRD2016))

# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
HARNSP17ft  = CRS("+init=epsg:2881")
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

# Calculate raster extents
xmin = floor(min(PixelCoords[c('X')]))
xmax = ceiling(max(PixelCoords[c('X')]))
ymin = floor(min(PixelCoords[c('Y')]))
ymax = ceiling(max(PixelCoords[c('Y')]))

# calculate number of rows and columns
# 6561.679 = feet in 2 Kilometer pixel spacing
rasRows <- (ymax - ymin) / 6561.679
rasCols <- (xmax - xmin) / 6561.679

# define raster template using NRD pixel data set extents
ras <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=HARNSP17ft)
rasExt <- extent(ras)
clpBnds2 <- gClip(WMDbnd, rasExt)
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
yr = 2003

for (yr in seq(2001, 2003)) {
#for (yr in seq(1999, 2016)) {
  year <- as.character(yr)
  #---------------------------
  # Compare NexRad vs RainGage
  #---------------------------
  startDay = as.Date(paste0("01/01/",year), "%m/%d/%Y")
  endDay = as.Date(paste0("12/31/",year), "%m/%d/%Y")
  RvsG <- as.data.frame(LWCNRDvsGage99to16[LWCNRDvsGage99to16$Rainfall < 20 &
                                               as.Date(LWCNRDvsGage99to16$DailyDate) >= startDay &
                                               as.Date(LWCNRDvsGage99to16$DailyDate) <= endDay, ])
  LWC_NRDbyYr <- read.csv(paste0(basePAth, "LWC_NRD", year, ".csv"))
  dateList<-names(LWC_NRDbyYr)[-c(1,2,3)]
 # NRDbyYr <- na.omit(melt(LWC_NRDbyYr, id = c("Pixel_id", "X", "Y")))
  NRDbyYr <- melt(LWC_NRDbyYr, id = c("Pixel_id", "X", "Y"))
  
  TotDays = endDay - startDay + 1
  SignificantDays = as.numeric(.9 * TotDays)
  
  RvsG[c("biasNRD")] <- NA
  diffDays =    RainStats <- NULL
  stationList = unique(RvsG$RainGage)
  
  for (gage in stationList) {
    if ((nrow(na.omit(RvsG[RvsG$RainGage == gage & RvsG$NRD >= 0.0, -c(14)]))) > 0) {
      nobsSigRain <-nrow(na.omit(RvsG[RvsG$RainGage == gage & RvsG$NRD >= 0.0, -c(14)]))
      nobsNRD <- nrow(RvsG[RvsG$RainGage == gage & !is.na(RvsG$NRD), ])
      nobsRain <-nrow(RvsG[RvsG$RainGage == gage & !is.na(RvsG$Rainfall), ])
      meanNRD <-mean(as.numeric(RvsG[RvsG$RainGage == gage & !is.na(RvsG$Rainfall), 11]),na.rm=TRUE)
      meanRain <-mean(as.numeric(RvsG[RvsG$RainGage == gage & !is.na(RvsG$Rainfall), 10]),na.rm=TRUE)
      bias = meanRain / meanNRD
      if (bias > .8 & bias < 1.2 &  nobsRain > SignificantDays) {
        cat(paste(gage, ",", bias, ",", nrow(na.omit(RvsG[RvsG$RainGage == gage , -c(14)]))), "\n")
        
        RvsG[RvsG$RainGage == gage &!is.na(RvsG$NRD), 14] <- 
          na.omit((RvsG[RvsG$RainGage == gage  & !is.na(RvsG$NRD), 11]) * bias)
        
        ValsList = data.frame(gage,as.numeric(bias),as.numeric(nobsSigRain),nobsNRD,nobsRain)
        RainStats <- rbind.data.frame(RainStats, ValsList)
      }
    }
  }
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
  csvFile <- paste0(basePAth, paste0("newNRDstats",year,".csv"))
  fwrite(na.omit(combineData), csvFile)
  
  #------------------------------
  # read and organize daily NexRad data
  #------------------------------
  rasStack <- stack()
  for (d in dateList[1]) {
    DailyNRD <- NRDbyYr[NRDbyYr$variable == d, ]
    DailyNRD.pnts <-SpatialPointsDataFrame(coords = DailyNRD[, c("X", "Y")],
                                          data = DailyNRD,proj4string = HARNSP17ft)
    
    NRDras <-rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = max) 
    NRDBiasPnts <- data.frame(extract(NRDras,DailyNRD.pnts))
    DailyNRD.pix <-SpatialPixels(SpatialPoints(coordinates(NRDras)[!is.na(values(NRDras)), ]))
    DailyNRD.grid <- as(NRDras, "SpatialGridDataFrame")
  }
  NRDPixel<- data.frame(DailyNRD.pnts$Pixel_id,DailyNRD.pnts$X,DailyNRD.pnts$Y)

  #------------------------------------------------------
  # Interpolate Bias from RainGages to NexRad pixels
  # Make NRD data correction using Bias
  #------------------------------------------------------
  
  biasData <- na.omit(combineData[, c("XCOORD", "YCOORD", "bias")])

  coordinates(biasData) =  ~ XCOORD + YCOORD
  proj4string(biasData) = HARNSP17ft
  
  biasData$X <- coordinates(biasData)[, 1]
  biasData$Y <- coordinates(biasData)[, 2]
  
  rainGage.pnts <-
    SpatialPointsDataFrame(coords = combineData[, c("XCOORD", "YCOORD")],
                           data = combineData,proj4string = HARNSP17ft)
  
  theisPoly <- voronoi(rainGage.pnts)
  NRDras <- rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = max)
  TheisRas <- rasterize(theisPoly, NRDras, theisPoly$bias, fun = mean)

  #Semivariogram
  #vgOK <- autofitVariogram(bias ~ 1, biasData)
  #predict
  #predOK <-
  # krige(bias ~ 1, biasData, NRD2007.grid, model = vgOK$var_model)
  #plot(predOK)
  #biasRas <- raster(predOK)
  
  myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))
  
  surf <- autoKrige(formula=bias ~ 1, input_data=biasData, new_data = DailyNRD.grid)
  biasRas <- raster(surf$krige_output)
  #plot(biasRas)
  
  gs <- gstat(formula=bias~1, locations=biasData)
  idw <- interpolate(NRDras, gs)
  # plot(idw)
  
  filename=paste(basePAth,"/rasterPlots/_bias",year,".png",sep="")
#  myplot=(  levelplot(surfRas, par.settings = myTheme, contour=FALSE,main=year,at= seq(.8,1.2,by=.01)) + 
  myplot=(  levelplot(biasRas, par.settings = myTheme, contour=FALSE,main=year) + 
                          layer(sp.polygons(clpBnds2)) +
              layer(sp.points(rainGage.pnts, col = "red")))  
  print (myplot)
  
  dev.copy(png, filename, width=2400,height=2400,units="px",res=300)
  dev.off()
  

  rasStack <- stack()
  NRDbiasPixels <-NRDPixel
  names(NRDbiasPixels)<-c( "Pixel_id" ,"X" ,"Y")
  for (d in dateList) {
    cat(paste("."))
    table = dayBiasFn(NRDbyYr[NRDbyYr$variable == d, ],biasRas)
    names(table) <- c( "Pixel_id" ,"X" ,"Y",d  )
    NRDbiasPixels <- merge(NRDbiasPixels, table, by =c("Pixel_id","X","Y"))

  #  rasStack <- stack(rasStack, NRDras)
  }
  NRDbiasPixels$Annual<-rowSums(NRDbiasPixels[,-c(1,2,3)])
  csvFile <- paste0(basePAth, paste0("biasNRD",year,".csv"))
  fwrite(NRDbiasPixels, csvFile)
}
