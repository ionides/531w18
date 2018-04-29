#setwd("C:/Users/yjin9/Documents/531project/final")
require(dplyr)

dengue<-read.csv("dengue_raw_data.csv")
dengue$wkid<-(dengue$Year-2012)*52+dengue$Week
dengue[dengue$Year==2015,]$wkid<-dengue[dengue$Year==2015,]$wkid+1
dengue[dengue$District=="Petalng",]$District<-"Petaling"
dengue[dengue$District=="Petalilng",]$District<-"Petaling"

by_cat<-group_by(dengue,wkid, District)
weekly_case<-summarise(by_cat,count=sum(Total_Cases),Week=mean(Week),Year=mean(Year))
petaling<-weekly_case[weekly_case$District=="Petaling",]
petaling$month<-floor(petaling$Week/4.345)+1
petaling[petaling$month==13,]$month<-12
petaling<-petaling[3:95,]

# Monthly in mm
rainfall<-read.csv("precipitation_raw_data.csv")
rainfall<-select(rainfall,pr,X.Year,Month)
rainfall<-rainfall[rainfall$X.Year>2011,]

data<-left_join(petaling,rainfall,by=c("Year"="X.Year","month"="Month"))
data$wkid<-data$wkid-77
write.csv(data,file="data_combined.csv")