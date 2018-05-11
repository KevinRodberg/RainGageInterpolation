library(RODBC)
library(reshape2)
library(raster)
library(dplyr)
library(rasterVis)
library(data.table)
library(sp)
library(rgdal)
library(maps)
library(maptools)
library(animation)
library(readr)
library(ggplot2)
geomSeries <- function(base, max) {
  (base ^ (0:floor(log(max, base))) - 1) / 100
}

clipToExtent <- function(sp, extent) {
  require(rgeos)
  keep <-
    gContains(extent, sp, byid = TRUE)
    gOverlaps(extent, sp, byid = TRUE)
  stopifnot(ncol(keep) == 1)
  sp[drop(keep), ]
}
gClip <- function(shp, bb) {
  if (class(bb) == "matrix") {
    b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
    } else{
      b_poly <- as(extent(bb), "SpatialPolygons")
      print (shp)
      print (b_poly)
      gIntersection(shp, b_poly, byid = T)
    }
}
DistrictBnds <- function() {
  #  color=rgb(0,0,0,alpha=.01)
  #  plot(clpBnds2,col=color )
  plot(clpBnds2, bg = "transparent", add = TRUE)
}

# Calculate raster extents
UTM17m = CRS("+init=epsg:26917") # UTM Zone 17 Meters
llwgs84  = CRS("+init=epsg:4326") # lat-long projection
HARNSP17ft  = CRS("+init=epsg:2881") # NAD_1983_HARN_StatePlane_Florida_East_FIPS_0901_Feet

# Calculate raster extents
xmin = floor(min(WMD.rain[c('X')]))
xmax = ceiling(max(WMD.rain[c('X')]))
ymin = floor(min(WMD.rain[c('Y')]))
ymax = ceiling(max(WMD.rain[c('Y')]))

# calculate number of rows and columns
# 6561.679 = feet in 2 Kilometer pixel spacing
rasRows <- (ymax - ymin) / 6561.679
rasCols <- (xmax - xmin) / 6561.679

lwcRAS <-
  raster(
    nrow = rasRows,
    ncol = rasCols,
    xmn = xmin,
    xmx = xmax,
    ymn = ymin,
    ymx = ymax,
    crs = NA
  )


WMDbnd.Path <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS/Rain"
WMDbnd.Shape <- "WMD_bnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape, proj4string = HARNSP17ft)
WMDbndUTM <- spTransform(WMDbnd, UTM17m)


#keep <- gContains( extent, sp,byid=TRUE ) | gOverlaps( extent, sp,byid=TRUE )
rasPoly <- as(extent(lwcRAS), 'SpatialPolygons')
proj4string(rasPoly) <- HARNSP17ft
#keep <- gContains( rasPoly, WMDbndUTM,byid=TRUE ) | gOverlaps( rasPoly, WMDbndUTM,byid=TRUE )

clpdBnds <- clipToExtent(WMDbnd, rasPoly)
clpBnds2 <- gClip(WMDbnd, rasPoly)

#Reads the rain gauge data query result from DBHydro
LWC_Rain_G_data_2003to2018 <-
  read_csv(
    "//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/LWC_Rain_G_data_2003to2018.csv",
    col_types = cols(
      DAILY_DATE = col_date(format = "%m/%d/%Y"),
      VALUE = col_double(),
      XCOORD = col_double(),
      YCOORD = col_double()
    ),
    skip = 0
  )

yearlist <-c(1997, 2003, 2005, 2006, 2007 ,2009, 2010 ,2011, 2013, 2014, 2015, 2016)
#for (year in c(2004,2008,2012,2016)) {
#for (year in yearlist) {
for (year in seq(1997,2016)) {
year = 1999
#print(paste("The year is", year))
  gageByyear <-
    subset(LWC_Rain_G_data_2003to2018, format(as.Date(DAILY_DATE), "%Y") == year, VALUE < 20)
  WideGage <-
    dcast(gageByyear,
          STATION + DATA_TYPE + XCOORD + YCOORD + DBKEY ~ DAILY_DATE,
          value.var = "VALUE")

# READS the NEXRDAD Files fro the NRD_Query.R file
  
  file <- paste("//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/LWC_NRD",
                year,".csv",sep = "")
  #
  #   Add an extra "d" when doing LEAP YEAR
  #
  Colstypes = "idddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd"
  
  WMD.rain <- read_csv(file, col_types = Colstypes)
  names(WMD.rain)
  # Clear raster Stack
  rasStack <- stack()
  if (is.element(year,  c(2000, 2004, 2008, 2012, 2016))) {
    dailyColNames = names(WMD.rain[, -c(1, 2, 3)])
  }  else {
    dailyColNames = names(WMD.rain[, -c(1, 2, 3, length(WMD.rain))])
  }
  #dailyColNames = names(nrd[,-c(1,2,3,4,5,6)])
  
  # create points from each column and convert
  # them to rasters and add each to the raster stack
  datalist = list()
  indx= 0
  for ( ColDate in dailyColNames)
  {
    indx = indx + 1
    nrdTotals <- WMD.rain[, c("X", "Y",  ColDate)]
    colnames(nrdTotals) <- c("x", "y",  ColDate)
    
    spg <- nrdTotals
    coordinates(spg) <- ~ x + y
    proj4string(spg) = HARNSP17ft
    
    spg <- spTransform(spg, HARNSP17ft)
    
    lwcRAS_1 <- rasterize(spg, lwcRAS,  ColDate, fun = mean)
    datestr <-gsub("X","",ColDate)
    #gageByDay <-subset(LWC_Rain_G_data_2003to2018, DAILY_DATE == as.Date(datestr) & is.na(LWC_Rain_G_data_2003to2018$CODE))
    gageByDay <-subset(LWC_Rain_G_data_2003to2018, DAILY_DATE == as.Date(datestr) )
    names(gageByDay)
    gageByDay$CODE <- NULL
    gageByDay$AGENCY <- NULL
    avgGage <-dcast(gageByDay,STATION+XCOORD+YCOORD ~ DAILY_DATE,fun=mean,value.var = "VALUE")
    coordinates(avgGage) <- c("XCOORD", "YCOORD")
    projection(avgGage) <- HARNSP17ft
    
    RainGage.df <-data.frame(DailyDate=datestr,avgGage["STATION"],avgGage[datestr],
                             extract(lwcRAS_1, avgGage, 
                                     method='simple'))
    names(RainGage.df)
    names(RainGage.df) <- c("DailyDate","RainGage","XCOORD","YCOORD","OPTIONAL","Rainfall","XCRD","YCRD","TRUEFALSE","NRD")
    RainGage.df[c("OPTIONAL","XCRD","YCRD","TRUEFALSE")]<- NULL
    datalist[[indx]] <- RainGage.df
 #   rasStack <- stack(rasStack, lwcRAS_1)
    # plot (lwcRAS_1,main= ColDate)
  }
  
  oneYear = do.call(rbind,datalist)
  one_file <- paste("//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/LWC_NRDvsGageLY",
                year,".csv",sep = "")
  write_csv( oneYear, one_file)
  noNA <-na.omit(oneYear[c("Rainfall", "NRD")])
#  noNA <-subset(noNA,RainGage < 1)
#  noNA <-subset(noNA,NRD < 1)
#  noNA <-subset(noNA,RainGage >=0.01)
#  noNA <-subset(noNA,NRD >=0.01)
#  noNA$roundNRD <- round(noNA$NRD,digits=2) 
#  noNA$roundRG <- round(noNA$RainGage,digits=2) 
#  noNA = subset(noNA,noNA$roundNRD > .25)
#  noNA = subset(noNA,noNA$roundRG > 0.25 )
#  noNA = subset(noNA,(round(RainGage-NRD,3)< -.005 )| (round(RainGage-NRD,3)> 0.005 ))
#  filename = paste("//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD_Files/LEC_NRD",year,"%03d.png", sep = "" )
#  png(  file = filename, width = 3000,height = 3000,units = "px",  res=300)
  plot(x=noNA$logRain,y=noNA$logNRD)
#  plot(x=noNA$roundRG,y=noNA$roundNRD)
#  model <-lm(noNA$roundRG ~ noNA$roundNRD)
  noNA <-subset(noNA,Rainfall >=0.001)
  noNA <-subset(noNA,NRD >=0.001)
  model <-lm(noNA$Rainfall ~ noNA$NRD)
  noNA$logNRD <- log(noNA$NRD)
  noNA$logRain <- log(noNA$NRD)
  model.exp <- lm(noNA$logRain~noNA$logNRD)
  fit_lm = lm(log(noNA$Rainfall)~log(noNA$NRD))
  plot (fit_lm,which=1)
#  abline(model)
  summary(model)
  summary(model.exp)
  model.exp$coefficients
  ggplot(noNA, aes(x=Rainfall,y=NRD)) + 
    geom_point() +
    geom_smooth(method="auto", se=TRUE)
  lines(noNA$Rainfall,noNA$Rainfall^model.exp$coefficients[2])
  #  legend("topright", bty="n", legend=paste("Year = ",year, "R2 =",format(summary(model)$adj.r.squared, digits=4)))
  #  dev.off()
  #  names(rasStack) <- dailyColNames
  #  options(scipen = 10000)
  
  #  print(    levelplot(      rasStack,      contour = FALSE,      par.settings = myTheme,
  #      at = geomSeries(base = 2, max = 520),      layout = c(5, 1),      scales = list(x = list(rot = 90))
  #    ) + layer(sp.polygons(WMDbnd))
  #  )
  dev.off()
}