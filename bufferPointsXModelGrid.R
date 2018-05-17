
myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))

basePAth <-  "//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/"

LWC_NRD1999 <- read.csv(paste0(basePAth, "LWC_NRD1999.csv"))
NRD1999 <- na.omit(melt(LWC_NRD1999, id = c("Pixel_id", "X", "Y")))
PixelCoords <- NRD2016[c("Pixel_id", "X", "Y")]
PixelPoints<-SpatialPointsDataFrame(coords = PixelCoords[, c("X", "Y")],
                       data = PixelCoords,proj4string = HARNSP17ft)
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
# Create a buffer polygon around pixels
#-------------------------------------------------

PixelPolys <-gBuffer(PixelPoints, width=12000, id=PixelPoints$Pixel_id, byid=TRUE )
