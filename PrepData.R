library(reshape2)
library(readr)
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

RvsG <-LWC_NRDvsGage2000_16[LWC_NRDvsGage2000_16$Rainfall <20 &
                         as.Date(LWC_NRDvsGage2000_16$DailyDate) >= as.Date("01/01/2007","%m/%d/%Y") &
                        as.Date(LWC_NRDvsGage2000_16$DailyDate) <= as.Date("12/31/2009","%m/%d/%Y"),]
names(RvsG)
gage <- "L005"
RvsG[c("newNRD")]<-NA
stationList = unique(RvsG$RainGage)
for (gage in stationList){
  if ((nrow(na.omit(RvsG[RvsG$RainGage==gage & RvsG$NRD >.25,-14])) )>0 ){
    fit0<-lm(Rainfall~0+NRD,data=RvsG[RvsG$RainGage==gage & RvsG$NRD > .25,])
    summary(fit0)
    cat(paste(gage,",",fit0$coefficients))
    cat(paste(",",nrow(na.omit(RvsG[RvsG$RainGage==gage ,])) ),"\n")
    
    RvsG[RvsG$RainGage==gage &
           RvsG$NRD > .25 & !is.na(RvsG$NRD),14]<- (RvsG[RvsG$RainGage==gage & RvsG$NRD > .25 &
                                                                               !is.na(RvsG$NRD),11])*fit0$coefficients[1]
    RvsG[RvsG$RainGage==gage &
           RvsG$NRD <= .25 & !is.na(RvsG$NRD),14]<- (RvsG[RvsG$RainGage==gage & RvsG$NRD <= .25 &
                                                           !is.na(RvsG$NRD),11])*1.0
    
  } 
}
testColl <-na.omit(RvsG[RvsG$RainGage==gage & RvsG$NRD >.25,]) 
testColl <-na.omit(RvsG[RvsG$RainGage==gage ,]) 
#ggplot(testColl,aes(Rainfall,NRD)) + stat_smooth(method=lm,formula = fit0$call$formula) + geom_point(aes(Rainfall,NRD))
ggplot(testColl,aes(Rainfall,newNRD)) + stat_smooth(method=lm) + geom_point(aes(Rainfall,NRD),color="red") + geom_point(aes(Rainfall,newNRD)) 
nrow(na.omit(RvsG[RvsG$RainGage==gage & RvsG$NRD> .25,]))


ggplot(testColl) + geom_point(aes(Rainfall,NRD))
