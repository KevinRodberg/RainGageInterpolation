library(RODBC)
library(RODBCext)
library(reshape2)
library(data.table)
library(raster)
library(dplyr)
library(rasterVis)

# Connect to DBHdryo's Database Instance (wrep) as PUB/PUB
channel <- odbcConnect("wrep", uid="pub", pwd="pub", believeNRows=FALSE)

# Check that connection is working (Optional)
odbcGetInfo(channel)
rm(nrdAnnualSums)
fileTots <- paste("Y:/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/LWC_NRD_Totals.csv")

# fileTots <- paste("Z:/DistrictAreaProj/LEC/Data/LEC_NRD_ETTotals.csv")

for (year in 2000:2017){
  print(paste("The year is", year))

  # ET data can be queried changing NRD_TIME_SERIES to NRD_ET_TIME_SERIES_VW 
  # and TSTYPEID = 3 (Rainfall) to TSTYPEID = 5 (Reference ET)
  # SQL sums data grouped daily so 15 min Rainfall and Daily ET summaries are consistent
  # Switch all references to Rain and RET
  
  #  file <- paste("Z:/DistrictAreaProj/LEC/Data/LEC_NRD_ET",year,".csv",sep="")
  file <- paste("Y:/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/LWC_NRD",year,".csv",sep="")
  
  source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/NexRad/queryTotalNRDbyParams.R")
  dataframe<-queryTotalNRD(year,NRDtype,MapArea)
    
 # dataframe <- sqlQuery(channel, sqlQ)
  names(dataframe)
  #--- switch Rain and RET  
  colnames(dataframe) <- c("Pixel_id", "X", "Y","Date", "Rain")
  #  "Mean" is an arbitrary function yet needed to pivot the data
  nrd <- dcast(dataframe, Pixel_id+X+Y ~ Date, mean )
  # Calculate Pixel sums for each year
  #--- switch Rain and RET  
  nrdSum <- aggregate(dataframe[c("Rain")],FUN=sum, by= list(dataframe$Pixel_id))
  nrdMax <- aggregate(dataframe[c("Rain")],FUN=max, by= list(dataframe$Date))
  cName <- paste("Total",toString(year),sep="")
  #--- switch Rain and RET  
  nrd[,cName]<- nrdSum$Rain
  # create new dataframe during 1st year of processing
  # append Annual Rainfall totals by pixel
  if (!exists('nrdAnnualSums')){
    nrdAnnualSums <- subset(nrd, select = c(1,2,3,ncol(nrd)))
  } else {
    nrdAnnualSums[,cName] <-nrdSum$Rain
  }
  # Export 
  write.table(nrd, file,  sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="0.00")
}
# Clear raster Stack
rasStack <-stack()

# Calculate raster extents
xmin = floor(min(nrd[c('X')]))
xmax = ceiling(max(nrd[c('X')]))
ymin = floor(min(nrd[c('Y')]))
ymax = ceiling(max(nrd[c('Y')]))

# calculate number of rows and columns
# 6561.679 = feet in 2 Kilometer pixel spacing
rasRows <- (ymax-ymin)/6561.679
rasCols <- (xmax-xmin)/6561.679

# define raster template
lecRAS <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=NA)

#fileTots <- paste("Z:/DistrictAreaProj/LEC/Data/LEC_NRD_Totals.csv")
#rm(nrdAnnualSums)
#nrdAnnualSums <-read.table(fileTots, header=TRUE, sep=",",  na.strings="0.00")

annualColNames = names(nrdAnnualSums[,-(1:3)])
# create points from each column and convert 
# them to rasters and add each to the raster stack
for (TotCol in annualColNames)
{
  testCname <-TotCol
  nrdTotals <- nrdAnnualSums[,c("X","Y",testCname)]
  colnames(nrdTotals) <- c("x", "y",TotCol)
  
  spg <- nrdTotals
  coordinates(spg) <- ~ x + y
  lecRAS_1 <-rasterize(spg,lecRAS,TotCol,fun=mean)
  rasStack <- stack(rasStack, lecRAS_1)
}
names(rasStack)<-annualColNames
levelplot(rasStack,contour=FALSE,par.settings = RdBuTheme, at= seq(10,100,length=10))

# Number of years equals number of columns in nrdAnnualSum - 3 columns for Pixel_Id, X, Y
nrdAnnualSums <-within(nrdAnnualSums, rm(Mean))
numberOfYrs = ncol(nrdAnnualSums)-3
annualColNames = names(nrdAnnualSums[-c(1:3)])

# Create a final column summary as a mean of all annual totals or POR Annual Rainfall or Refernce ET 
if (numberOfYrs > 1) {
  nrdAnnualSums$Mean <- rowSums(nrdAnnualSums[,-c(1:3)])/numberOfYrs
} else {
  nrdAnnualSums$Mean <-nrdAnnualSums[,-c(1:3)]
}
write.table(nrdAnnualSums, fileTots,  sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="0.00")
