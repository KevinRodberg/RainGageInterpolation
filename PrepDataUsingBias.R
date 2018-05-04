library(reshape2)
library(readr)
library(dplyr)
library(data.table)
basePAth <-"//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/"
LWC_Rain_03to18 <- read_csv(paste0(basePAth,"LWC_Rain_G_data_2003to2018.csv"))
LWC_NRD2007 <- read_csv(paste0(basePAth,"LWC_NRD2007.csv"))
LWC_NRD2008 <- read_csv(paste0(basePAth,"LWC_NRD2008.csv"))
LWC_NRD2009 <- read_csv(paste0(basePAth,"LWC_NRD2009.csv"))
NRD2007<-melt(LWC_NRD2007,id=c("Pixel_id", "X","Y"))
NRD2008<-melt(LWC_NRD2008,id=c("Pixel_id", "X","Y"))
NRD2009<-melt(LWC_NRD2009,id=c("Pixel_id", "X","Y"))
LWC_NRDvsGage08 <- read_csv(paste0(basePAth,"LWC_NRDvsGageLY2008.csv"))
LWC_NRDvsGage2000_16 <- read_excel(paste0(basePAth,"Compared_/LWC_NRDvsGageLY2000-2016.xlsx"))

RvsG <-as.data.frame(LWC_NRDvsGage2000_16[LWC_NRDvsGage2000_16$Rainfall <20 &
                         as.Date(LWC_NRDvsGage2000_16$DailyDate) >= as.Date("01/01/2007","%m/%d/%Y") &
                        as.Date(LWC_NRDvsGage2000_16$DailyDate) <= as.Date("12/31/2009","%m/%d/%Y"),])
names(RvsG)
RvsG[c("newNRD","biasNRD")]<-NA
RainStats<- data.frame(gage=character(),coefficient=numeric(),nrows=numeric())
stationList = unique(RvsG$RainGage)
for (gage in stationList){
  if ((nrow(na.omit(RvsG[RvsG$RainGage==gage & RvsG$NRD >.25,-c(14,15)])) )>0 ){
    nobsSigRain <-nrow(na.omit(RvsG[RvsG$RainGage==gage & RvsG$NRD >.25,-c(14,15)]))
    nobsNRD <-nrow(RvsG[RvsG$RainGage==gage & !is.na(RvsG$NRD),])
    nobsRain <-nrow(RvsG[RvsG$RainGage==gage & !is.na(RvsG$Rainfall),])
    meanNRD <-mean(as.numeric(RvsG[RvsG$RainGage==gage & !is.na(RvsG$Rainfall),11]))
    meanRain <-mean(RvsG[RvsG$RainGage==gage & !is.na(RvsG$Rainfall),10])
    bias = meanRain/meanNRD
    fit0<-lm(Rainfall~0+NRD,data=RvsG[RvsG$RainGage==gage & RvsG$NRD > .25,])
    summary(fit0)
    cat(paste(gage,",",fit0$coefficients))
    cat(paste(",",nrow(na.omit(RvsG[RvsG$RainGage==gage ,-c(14,15)])) ),"\n")
    
    RvsG[RvsG$RainGage==gage &
           RvsG$NRD > .25 & !is.na(RvsG$NRD),14]<- (RvsG[RvsG$RainGage==gage & RvsG$NRD > .25 &
                                                                               !is.na(RvsG$NRD),11])*fit0$coefficients[1]
    RvsG[RvsG$RainGage==gage &
           RvsG$NRD <= .25 & !is.na(RvsG$NRD),14]<- (RvsG[RvsG$RainGage==gage & RvsG$NRD <= .25 &
                                                           !is.na(RvsG$NRD),11])*1.0
    
    RvsG[RvsG$RainGage==gage & !is.na(RvsG$NRD),15]<- na.omit((RvsG[RvsG$RainGage==gage  & !is.na(RvsG$NRD),11])*bias)
    
    ValsList = c(gage,as.numeric(unlist(fit0$coefficients[1])), as.numeric(nobsSigRain), nobsNRD, nobsRain )
    record = as.data.frame(t(ValsList))
    RainStats<-rbind.data.frame(RainStats,record)
  } 
}
names(RainStats)<-c("RainGage","coefficient","SigRainObs","NRDobs","RainObs")
names(RvsG)
newNRDSum<-aggregate(RvsG$newNRD,list(RvsG$YEAR,RvsG$RainGage,RvsG$XCOORD,RvsG$YCOORD),FUN=sum)
names(newNRDSum) = c("YEAR","RainGage","XCOORD","YCOORD","newNRD")
NRDSum<-aggregate(RvsG$NRD,list(RvsG$YEAR,RvsG$RainGage,RvsG$XCOORD,RvsG$YCOORD),FUN=sum)
names(NRDSum) = c("YEAR","RainGage","XCOORD","YCOORD","NRD")
GageSum<-aggregate(RvsG$Rainfall,list(RvsG$YEAR,RvsG$RainGage,RvsG$XCOORD,RvsG$YCOORD),FUN=sum)
names(GageSum) = c("YEAR","RainGage","XCOORD","YCOORD","Rainfall")
BiasSum<-aggregate(RvsG$biasNRD,list(RvsG$YEAR,RvsG$RainGage,RvsG$XCOORD,RvsG$YCOORD),FUN=sum)
names(BiasSum) = c("YEAR","RainGage","XCOORD","YCOORD","biasNRD")

compareSums<-merge(merge(merge(newNRDSum,NRDSum,all=TRUE),GageSum,all=TRUE),BiasSum,all=TRUE)

combineData <- plyr::join(compareSums, RainStats, by='RainGage',type ='left', match='all')

min(as.numeric(combineData$coefficient),na.rm=TRUE)
mean(as.numeric(combineData$coefficient),na.rm=TRUE)
max(as.numeric(combineData$coefficient),na.rm=TRUE)
csvFile <-paste0(basePAth,"newNRDstats.csv")
fwrite(na.omit(combineData),csvFile)
gage <- "S79_R"

testColl <-na.omit(RvsG[RvsG$RainGage==gage ,]) 
ggplot(testColl,aes(Rainfall,newNRD)) + stat_smooth(method=lm) + geom_point(aes(Rainfall,NRD),color="red") + geom_point(aes(Rainfall,newNRD)) 
nrow(na.omit(RvsG[RvsG$RainGage==gage & RvsG$NRD> .25,]))


ggplot(testColl) + geom_point(aes(Rainfall,NRD))
