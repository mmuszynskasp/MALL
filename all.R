rm(list=ls())

library("MortalitySmooth")
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(readr)


out.dir <- "C://Users//Magdalena//demography//philosopyineq//data"
deaths.dir <- "C://Users//Magdalena//demography//philosopyineq//data//Deaths_1x1"
exposures.dir <- "C://Users//Magdalena//demography//philosopyineq//data//Exposures_1x1"
flt.dir <- "C://Users//Magdalena//demography//philosopyineq//data//fltper_1x1"
mlt.dir <- "C://Users//Magdalena//demography//philosopyineq//data//mltper_1x1"


######prepare data from HMD: exposures and deaths for single age groups
countrylist <- c("AUS", "AUT","BEL","BGR","CAN","HRV","CHE","CZE","DEUTNP","DNK","ESP","EST","FIN","FRATNP","GRC","HUN","IRL",
                 "ISL","ISR","ITA","JPN", "LTU", "LUX","LVA","NLD", "NOR","NZL_NP","POL","PRT","RUS","SVK","SVN","SWE",    
                 "GBR_NP","USA")

setwd(deaths.dir)
i=1
deaths <- read_table(paste0(countrylist[i],".Deaths_1x1.txt"),col_names=TRUE,skip=2)
Country <- countrylist[i]
deaths <- cbind(deaths,Country)
setwd(out.dir)
write.table(deaths, file="Deaths.txt",sep=",",row.names=FALSE,col.names=TRUE)

setwd(exposures.dir)
exposures <- read_table(paste0(countrylist[i],".Exposures_1x1.txt"),col_names=TRUE,skip=2)
exposures <- cbind(exposures,Country)
setwd(out.dir)
write.table(exposures, file="Exposures.txt",sep=",",row.names=FALSE,col.names=TRUE)

setwd(flt.dir)
flt <- read_table(paste0(countrylist[i],".fltper_1x1.txt"),col_names=TRUE,skip=2)
Sex <- "Female"
flt <- cbind(flt,Country,Sex)
setwd(out.dir)
write.table(flt, file="lt.txt",sep=",",row.names=FALSE,col.names=TRUE)

setwd(mlt.dir)
mlt <- read_table(paste0(countrylist[i],".mltper_1x1.txt"),col_names=TRUE,skip=2)
Sex <- "Male"
mlt <- cbind(mlt,Country,Sex)
setwd(out.dir)
write.table(mlt, file="lt.txt",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

for (i in 2:length(countrylist)){
  setwd(deaths.dir)
  deaths <- read_table(paste0(countrylist[i],".Deaths_1x1.txt"),col_names=TRUE,skip=2)
  Country <- countrylist[i]
  deaths <- cbind(deaths,Country)
  
  setwd(exposures.dir)
  exposures <- read_table(paste0(countrylist[i],".Exposures_1x1.txt"),col_names=TRUE,skip=2)
  exposures <- cbind(exposures,Country)
  
  setwd(flt.dir)
  flt <- read_table(paste0(countrylist[i],".fltper_1x1.txt"),col_names=TRUE,skip=2)
  Sex <- "Female"
  flt <- cbind(flt,Country,Sex)
  
  setwd(mlt.dir)
  mlt <- read_table(paste0(countrylist[i],".mltper_1x1.txt"),col_names=TRUE,skip=2)
  Sex <- "Male"
  mlt <- cbind(mlt,Country,Sex)
  
  setwd(out.dir)
  write.table(deaths, file="Deaths.txt",sep=",",row.names=FALSE,col.names=FALSE, append=TRUE)
  write.table(exposures, file="Exposures.txt",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
  write.table(flt, file="lt.txt",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
  write.table(mlt, file="lt.txt",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
}


#####################################################################################################################
#######smoothing with P-splines ages 10-109 
rm(list=ls())

library("MortalitySmooth")
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(readr)


out.dir <- "C://Users//Magdalena//demography//philosopyineq//data"
setwd(out.dir)
deaths <- read.table("Deaths.txt", sep=",", header=TRUE) 
exposures <- read.table("Exposures.txt", sep=",", header=TRUE)

ages <-10:105 #we start at age 10 as recommended by Horiuchi et al. (2013)
years <- 1900:2020
sex <- c("Female","Male")

countrylist <- sort(unique(deaths$Country))

################   #prepare data function  by country,sex,year
yearsc <- function(exposures,countryj){
  expcountry <- exposures %>%
    filter(Country==countrylist[countryj], Year %in%years) 
  yearscountry <- sort(unique(expcountry$Year))
  return(yearscountry)
}


dataprep <- function(deaths,exposures,countryj,yearscountry,yeari,sexk){
  year <- yearscountry[yeari]
  
  mydata <- deaths %>%
    filter(Country==countrylist[countryj], Year==year, Age%in%ages) %>%
    pivot_longer(cols = Female:Total,
                 names_to = "Sex", values_to = "Deaths") %>%
    left_join(exposures %>%
                filter(Country==countrylist[countryj], Year==year, Age%in%ages)%>%
                pivot_longer(cols = Female:Total,
                             names_to = "Sex", values_to = "Exposures"))  %>%
    mutate(Age=as.numeric(Age)) 
  return(mydata)
}


#smoothing   
modesmooth <- function(mydata,sexk){
  mydata <- mydata %>%
    filter(Sex==sex[sexk])
  
  countrym <- unique(mydata$Country)
  sexm <- unique(mydata$Sex)
  yearm <- unique(mydata$Year)
  
  if (is.na(mydata$Deaths[1])){
    outdata <- cbind(countrym,sexm,yearm,NA)
  } else{
    mydata <- mydata %>%
      mutate(Exposures=ifelse(Exposures==0,0.1,Exposures))
    
    fit1Df <- Mort1Dsmooth(x=mydata$Age, y=mydata$Deaths, offset=log(mydata$Exposures))
    
    mx <- exp(fit1Df$logmortality)
    qx <- mx/(1+0.5*mx) #easy for now
    lx <- matrix(100000,nrow=nrow(mx),ncol=ncol(mx))
    dx <- matrix(100000,nrow=nrow(mx),ncol=ncol(mx))
    
    for (i in 2:length(lx)){
      lx[i] <- lx[i-1]*(1-qx[i-1])
      dx[i-1] <- lx[i-1]-lx[i]
    }
    dx[length(lx)] <- lx[length(lx)]
    
    age <- mydata$Age
    dx <- as.vector(dx)
    disdx <-as_tibble(cbind(age,dx)) 
    mode1 <- disdx$age[dx==max(dx)]
    if (mode1==100){
      mode <- NA} else{
        nm1 <- disdx$dx[disdx$age==(as.numeric(mode1)-1)]
        np1 <- disdx$dx[disdx$age==(as.numeric(mode1)+1)]
        
        mode <- mode1+(max(dx)-nm1)/((max(dx)-nm1)+(max(dx)-np1))
      }
    outdata <- cbind(countrym,sexm,yearm,mode)
    colnames(outdata) <- c("country","sex","year","mode")
  }
  return(outdata)
}


#first file
countryj <- 1
countryjy <- yearsc(exposures=exposures,countryj=countryj)

#first country, first year
yeari <- 1
sexk <- 1
datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yearscountry=countryjy, yeari=yeari)
modesm <- modesmooth(mydata=datap,sexk=sexk)

write.table(modesm, file="mode.txt", sep=",", row.names=FALSE)

#males,1st country
sexk <- 2
datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yearscountry=countryjy,yeari=yeari)
modesm <- modesmooth(mydata=datap,sexk=sexk)

write.table(modesm, file="mode.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)

##rest years for 1st country
for (yeari in 2:length(countryjy)){
  for (sexk in 1:2){
    datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
    modesm <- modesmooth(mydata=datap,sexk=sexk)
    write.table(modesm, file="mode.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  }
}

#all other countries
for (countryj in c(2:11,14:length(countrylist))){  
  countryjy <- yearsc(exposures=exposures,countryj=countryj)
  for (yeari in 1:length(countryjy)){
    for (sexk in 1:2){
      datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
      modesm <- modesmooth(mydata=datap,sexk=sexk)
      write.table(modesm, file="mode.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    }
  }
}

countryj <- 12   #Finland, problem for males in 1941
countryjy <- yearsc(exposures=exposures,countryj=countryj)
sexk <- 2  
for (yeari in c(1:41,46:length(countryjy))){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  modesm <- modesmooth(mydata=datap,sexk=sexk)
  write.table(modesm, file="mode.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}

sexk <- 1
for (yeari in c(1:length(countryjy))){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  modesm <- modesmooth(mydata=datap,sexk=sexk)
  write.table(modesm, file="mode.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}



countryj <- 13   #France, problem for males in 1915
countryjy <- yearsc(exposures=exposures,countryj=countryj)
sexk <- 2  
for (yeari in c(1:15,18:length(countryjy))){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  modesm <- modesmooth(mydata=datap,sexk=sexk)
  write.table(modesm, file="mode.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}

sexk <- 1
for (yeari in c(1:length(countryjy))){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  modesm <- modesmooth(mydata=datap,sexk=sexk)
  write.table(modesm, file="mode.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}




#####################################################################################################################
#######smoothing with P-splines ages 5-15 
rm(list=ls())

library("MortalitySmooth")
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(readr)


out.dir <- "C://Users//Magdalena//demography//philosopyineq//data"
setwd(out.dir)
deaths <- read.table("Deaths.txt", sep=",", header=TRUE) 
exposures <- read.table("Exposures.txt", sep=",", header=TRUE)

ages <-5:15 #as in Ebeling
years <- 1900:2020
sex <- c("Female","Male")

countrylist <- sort(unique(deaths$Country))

################   #prepare data function  by country,sex,year
yearsc <- function(exposures,countryj){
  expcountry <- exposures %>%
    filter(Country==countrylist[countryj], Year %in%years) 
  yearscountry <- sort(unique(expcountry$Year))
  return(yearscountry)
}


dataprep <- function(deaths,exposures,countryj,yearscountry,yeari,sexk){
  year <- yearscountry[yeari]
  
  mydata <- deaths %>%
    filter(Country==countrylist[countryj], Year==year, Age%in%ages) %>%
    pivot_longer(cols = Female:Total,
                 names_to = "Sex", values_to = "Deaths") %>%
    left_join(exposures %>%
                filter(Country==countrylist[countryj], Year==year, Age%in%ages)%>%
                pivot_longer(cols = Female:Total,
                             names_to = "Sex", values_to = "Exposures"))  %>%
    mutate(Age=as.numeric(Age)) 
  return(mydata)
}


#smoothing   
minsmooth <- function(mydata,sexk){
  mydata <- mydata %>%
    filter(Sex==sex[sexk])
  
  countrym <- unique(mydata$Country)
  sexm <- unique(mydata$Sex)
  yearm <- unique(mydata$Year)
  
  if (is.na(mydata$Deaths[1])){
    outdata <- cbind(countrym,sexm,yearm,NA)
  } else{
    mydata <- mydata %>%
      mutate(Exposures=ifelse(Exposures==0,0.1,Exposures)) %>%
      mutate(Deaths=ifelse(Deaths==0,0.1,Deaths))
    
    fit1Df <- Mort1Dsmooth(x=mydata$Age, y=mydata$Deaths, offset=log(mydata$Exposures))
    
    mx <- exp(fit1Df$logmortality)
    qx <- mx/(1+0.5*mx) #easy for now
    lx <- matrix(100000,nrow=nrow(mx),ncol=ncol(mx))
    dx <- matrix(100000,nrow=nrow(mx),ncol=ncol(mx))
    
    for (i in 2:length(lx)){
      lx[i] <- lx[i-1]*(1-qx[i-1])
      dx[i-1] <- lx[i-1]-lx[i]
    }
    
    age <- mydata$Age
    dx <- as.vector(dx)
    disdx <-as_tibble(cbind(age,dx)) 
    amin1 <- disdx$age[dx==min(dx)]
    if (amin1==min(age)){
      amin <- NA} else{
        nm1 <- disdx$dx[disdx$age==(as.numeric(amin1)-1)]
        np1 <- disdx$dx[disdx$age==(as.numeric(amin1)+1)]
        amin <- amin1+(min(dx)-nm1)/((min(dx)-nm1)+(min(dx)-np1))
        
      }
    outdata <- cbind(countrym,sexm,yearm,amin)
    colnames(outdata) <- c("country","sex","year","agemin")
  }
  return(outdata)
}


#first file
countryj <- 1
countryjy <- yearsc(exposures=exposures,countryj=countryj)

#first country, first year
yeari <- 1
sexk <- 1
datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yearscountry=countryjy, yeari=yeari)
minsm <- minsmooth(mydata=datap,sexk=sexk)

write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE)

#males,1st country
sexk <- 2
datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yearscountry=countryjy,yeari=yeari)
minsm <- minsmooth(mydata=datap,sexk=sexk)

write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)

##rest years for 1st country
for (yeari in 2:length(countryjy)){
  for (sexk in 1:2){
    datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
    minsm <- minsmooth(mydata=datap,sexk=sexk)
    write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  }
}

#all other countries
for (countryj in c(2,4:11,14:33, 35:length(countrylist))){  
  countryjy <- yearsc(exposures=exposures,countryj=countryj)
  for (yeari in 1:length(countryjy)){
    for (sexk in 1:2){
      datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
      minsm <- minsmooth(mydata=datap,sexk=sexk)
      write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    }
  }
}

countryj <- 12   #Finland, problem for males in 1941
countryjy <- yearsc(exposures=exposures,countryj=countryj)
sexk <- 2  
for (yeari in c(1:41,46:length(countryjy))){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  minsm <- minsmooth(mydata=datap,sexk=sexk)
  write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}

sexk <- 1
for (yeari in c(1:length(countryjy))){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  minsm <- minsmooth(mydata=datap,sexk=sexk)
  write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}



countryj <- 13   #France, problem for males in 1915
countryjy <- yearsc(exposures=exposures,countryj=countryj)
sexk <- 2  
for (yeari in c(1:15,18:length(countryjy))){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  minsm <- minsmooth(mydata=datap,sexk=sexk)
  write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}

sexk <- 1
for (yeari in c(1:length(countryjy))){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  minsm <- minsmooth(mydata=datap,sexk=sexk)
  write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}


countryj <- 34 #sweden
countryjy <- yearsc(exposures=exposures,countryj=countryj)
sexk <- 2
for (yeari in 1:length(countryjy)){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  minsm <- minsmooth(mydata=datap,sexk=sexk)
  write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}

ages <-3:15 
sexk <- 1
for (yeari in 1:length(countryjy)){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  minsm <- minsmooth(mydata=datap,sexk=sexk)
  write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}


countryj <- 3 #belgium
countryjy <- yearsc(exposures=exposures,countryj=countryj)
ages <-5:15 
sexk <- 2
for (yeari in 1:length(countryjy)){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  minsm <- minsmooth(mydata=datap,sexk=sexk)
  write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}

ages <-3:15 
sexk <- 1
for (yeari in 1:length(countryjy)){
  datap <- dataprep(deaths=deaths,exposures=exposures,countryj=countryj,yeari=yeari,yearscountry=countryjy)
  minsm <- minsmooth(mydata=datap,sexk=sexk)
  write.table(minsm, file="agemin.txt", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
}



#########################
rm(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(readr)

out.dir <- "C://Users//Magdalena//demography//philosopyineq//data"
setwd(out.dir)
mode <- read.table("mode.txt", sep=",", header=TRUE) 
lt <- read.table("lt.txt", sep=",", header=TRUE) %>%
  mutate(Age=ifelse(Age=="110+",110,Age))


countrylist <- sort(unique(mode$country))
ages <-20:100
years <- 1900:2019
sex <- c("Female","Male")

MALL <- lt %>%
  mutate(sex=Sex,country=Country,year=Year) %>%
  left_join(mode) %>%
  filter(Age==floor(mode)) %>%
  mutate(MALL=mode-(1.25*ex-mode+floor(mode))) %>%
  filter(!is.na(MALL)) %>%
  select(Year,Country,Sex,mode,MALL)


Iprem <- lt %>%
  left_join(lt %>%
              dplyr::group_by(Year,Country,Sex) %>%
              dplyr::summarise(dx0=sum(dx,na.rm=TRUE)), by = c("Year", "Country", "Sex")) %>%
  left_join(MALL, by = c("Year", "Country", "Sex")) %>%
  filter(Age==floor(MALL)) %>%
  mutate(Iprem=(100000-lx)/100000+(mode-floor(mode))*dx/100000) %>% # how many deaths prior to MALL, takes into account number of deaths from lower age interval to the actuall MALL
  select(Year,Country,Sex,mode,MALL,Iprem)  %>%
  left_join(lt %>%  #Iprem 0-10
              filter(as.numeric(Age)<10, Country %in% countrylist, Year %in% years ) %>%
              dplyr::group_by(Year,Country,Sex) %>%
              dplyr::summarise(dx10=sum(dx,na.rm=TRUE)) %>% 
              mutate(Iprem010=dx10/100000)) %>%
  mutate(Iprem10=Iprem-Iprem010) %>%
  select(-c(dx10,Iprem010)) 

Depprembase <- lt %>%
  left_join(lt %>%
              left_join(MALL, by = c("Year", "Country", "Sex")) %>%
              filter(Age==floor(MALL)) %>%
              mutate(dx.MALL=dx) %>%
              select(Year,Country,Sex,dx.MALL) %>%
              left_join(lt %>%
                          left_join(MALL, by = c("Year", "Country", "Sex")) %>% 
                          filter(Age==floor(MALL)) %>%
                          mutate(Tx.MALL=Tx) %>%
                          select(Year,Country,Sex,Tx.MALL,MALL), by = c("Year", "Country", "Sex"))%>%
              mutate(Tx.MALL=Tx.MALL-0.25*dx.MALL),  by = c("Year", "Country", "Sex")) %>% #correction in Tx for the time between exact MALL and actual MALL
  mutate(ex.lost=(Tx-Tx.MALL+ax*dx)/(lx-ax*dx),   #mean number of years lost between age at death and MALL, in the middle of the age interval or at x+ax
         dx.ex.lost= ex.lost*dx) %>%
  filter(as.numeric(Age)<MALL) 

Depprem <- Depprembase %>%
  dplyr::group_by(Country,Sex,Year) %>%
  dplyr::summarise(edagger=sum(dx.ex.lost)/100000) %>%
  left_join(Depprembase %>%
              filter(as.numeric(Age)>=10) %>%
              left_join(Depprembase %>%
                          filter(as.numeric(Age)==10) %>%
                          mutate(l10=lx) %>% 
                          select(Year,Country,Sex,l10), by = c("Year", "Country", "Sex"))%>%
              dplyr::group_by(Country,Sex,Year) %>%
              dplyr::summarise(edagger10=sum(dx.ex.lost)/l10)) %>%
  distinct() %>%
  left_join(Iprem) %>%
  mutate(SM=mode-MALL) %>%
  pivot_longer(edagger:SM,names_to="Stat_name",values_to = "Stat_value")

#################################################################################################
##################plots
selcountries <- c("BEL","SWE","ITA","POL","USA","JPN")

bothMALL <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="MALL") %>%
  ggplot(aes(x=Year, y=Stat_value, Group = Country, color=Country)) + 
  geom_smooth(size = 0.7,se=F)+
  #  scale_color_manual(values=c("#0F8B8D","#EC9A29","#A8201A", "black","red","blue"))+
  facet_grid(vars(Sex))+
  theme_minimal() +
  theme(legend.position="bottom")+
  scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
  scale_y_continuous(name="Years", limits=c(60,90))+  
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))

ggsave(bothMALL, file="C:/Users/Magdalena/demography/philosopyineq/figures/bothMALL.pdf",width = 16, height = 20, units = "cm")

####
bothmode <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="mode") %>%
  ggplot(aes(x=Year, y=Stat_value, Group = Country, color=Country)) + 
  geom_smooth(size = 0.7,se=F)+
  #  scale_color_manual(values=c("#0F8B8D","#EC9A29","#A8201A", "black","red","blue"))+
  facet_grid(vars(Sex))+
  theme_minimal() +
  theme(legend.position="bottom")+
  scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
  scale_y_continuous(name="Years", limits=c(60,100))+  
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))

bothSM <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="SM") %>%
  ggplot(aes(x=Year, y=Stat_value, Group = Country, color=Country)) + 
  geom_smooth(size = 0.7,se=F)+
  #  scale_color_manual(values=c("#0F8B8D","#EC9A29","#A8201A", "black","red","blue"))+
  facet_grid(vars(Sex))+
  theme_minimal() +
  theme(legend.position="bottom")+
  scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
  scale_y_continuous(name="Years", limits=c(5,10))+  
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))


both2 <- ggarrange(bothmode,bothSM,nrow=1,ncol=2,common.legend = TRUE, legend="bottom")

ggsave(both2, file="C:/Users/Magdalena/demography/philosopyineq/figures/both.pdf",width = 25, height = 25, units = "cm")


#########cor mode, s(m+)
cormodesm <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="mode"|Stat_name=="SM") %>%
  pivot_wider(names_from = Stat_name, values_from = Stat_value) %>%
  ggplot(aes(x=mode, y=SM,Group = Country, color=Country)) + 
  facet_grid(vars(Sex))+
  geom_point(size=0.3)+
  labs(y= "S(M+)", x = "M")+
  #  geom_smooth(se=F,method='lm',size=0.5)+
  theme_minimal() +
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))

ggsave(cormodesm, file="C:/Users/Magdalena/demography/philosopyineq/figures/cormodesm.pdf",width = 12, height = 15, units = "cm")



cormallI <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="MALL"|Stat_name=="Iprem") %>%
  pivot_wider(names_from = Stat_name, values_from = Stat_value) %>%
  ggplot(aes(x=MALL, y=Iprem, Group = Country, color=Country)) + 
  facet_grid(vars(Sex))+
  geom_point(size=0.7)+
  #  geom_smooth(se=F,method='lm',size=0.5)+
  theme_minimal() +
  labs(y= "I", x = "MALL")+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))+
  ggtitle("MALL and  I")

cormallI10 <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="MALL"|Stat_name=="Iprem10") %>%
  pivot_wider(names_from = Stat_name, values_from = Stat_value) %>%
  ggplot(aes(x=MALL, y=Iprem10, Group = Country, color=Country)) + 
  facet_grid(vars(Sex))+
  geom_point(size=0.7)+
  #  geom_smooth(se=F,method='lm',size=0.5)+
  theme_minimal() +
  labs(y= "I*", x = "MALL")+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))+
  ggtitle("MALL and I*")

bothI <- ggarrange(cormallI,cormallI10,nrow=1,ncol=2,common.legend = TRUE, legend="bottom")

ggsave(bothI, file="C:/Users/Magdalena/demography/philosopyineq/figures/cormallI.pdf",width = 25, height = 25, units = "cm")



cortest <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="MALL"|Stat_name=="Iprem") %>%
  pivot_wider(names_from = Stat_name, values_from = Stat_value) %>%
  dplyr::group_by(Sex,Country) %>%
  dplyr::summarise(corc=cor(MALL,Iprem))



cortest10 <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="MALL"|Stat_name=="Iprem10") %>%
  pivot_wider(names_from = Stat_name, values_from = Stat_value) %>%
  filter(!is.na(Iprem10)) %>%
  dplyr::group_by(Sex,Country) %>%
  dplyr::summarise(corc=cor(MALL,Iprem10))


####
mystats <- c("Iprem","Iprem10")
Iprem <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name %in% mystats) %>%
  mutate(Statistic=Stat_name,
         Statistic=recode(Statistic,"Iprem"="I","Iprem10"="I*")) %>%
  ggplot(aes(x=Year, y=Stat_value, Group = Country, color=Country)) + 
  geom_smooth(size = 0.7,se=F)+
  labs(y= "Incidence", x = "Years")+
  #  scale_color_manual(values=c("#0F8B8D","#EC9A29","#A8201A", "black","red","blue"))+
  facet_grid(vars(Sex),vars(Statistic))+
  theme_minimal() +
  scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))

ggsave(Iprem, file="C:/Users/Magdalena/demography/philosopyineq/figures/Iprem.pdf",width = 25, height = 20, units = "cm")


mystats <- c("edagger","edagger10")

edagger <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name %in% mystats) %>%
  mutate(Statistic=Stat_name,
         Statistic=recode(Statistic,"edagger"="YPLL","edagger10"="YPLL*")) %>%
  ggplot(aes(x=Year, y=Stat_value, Group = Country, color=Country)) + 
  geom_smooth(size = 0.7,se=F)+
  labs(y= "Years", x = "Years")+
  #  scale_color_manual(values=c("#0F8B8D","#EC9A29","#A8201A", "black","red","blue"))+
  facet_grid(vars(Sex),vars(Statistic))+
  theme_minimal() +
  scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))

ggsave(edagger, file="C:/Users/Magdalena/demography/philosopyineq/figures/edagger.pdf",width = 25, height = 20, units = "cm")


corIdep10 <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="Iprem10"|Stat_name=="edagger10") %>%
  pivot_wider(names_from = Stat_name, values_from = Stat_value) %>%
  ggplot(aes(x=Iprem10, y=edagger10, Group = Country, color=Country)) + 
  labs(x= "I*", y = "YPLL*")+
  facet_grid(vars(Sex))+
  geom_point(size=0.7)+
  #  geom_smooth(se=F,method='lm',size=0.5)+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))

corIdep <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="Iprem"|Stat_name=="edagger") %>%
  pivot_wider(names_from = Stat_name, values_from = Stat_value) %>%
  ggplot(aes(x=Iprem, y=edagger, Group = Country, color=Country)) + 
  facet_grid(vars(Sex))+
  geom_point(size=0.7)+
  #    geom_smooth(se=F,method='lm',size=0.5)+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))+
  labs(y= "YPLL", x = "I")


bothIDep <- ggarrange(corIdep,corIdep10,nrow=1,ncol=2,common.legend = TRUE, legend="bottom")


ggsave(bothIDep, file="C:/Users/Magdalena/demography/philosopyineq/figures/corIdep.pdf",width = 25, height = 25, units = "cm")


corIdep10 <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="Iprem10"|Stat_name=="edagger10") %>%
  pivot_wider(names_from = Stat_name, values_from = Stat_value) %>%
  filter(!is.na(Iprem10)) %>%
  dplyr::group_by(Country,Sex)%>%
  dplyr::summarise(corc=cor(Iprem10,edagger10))



corIdep <- Depprem %>% 
  filter(Country=="USA", Stat_name=="Iprem"|Stat_name=="edagger", Sex=="Female") %>%
  pivot_wider(names_from = Stat_name, values_from = Stat_value)

###########separating age premature adult plot
out.dir <- "C://Users//Magdalena//demography//philosopyineq//data"
setwd(out.dir)
amin <- read.table(file="agemin.txt", sep=",", header=TRUE) %>%
  mutate(agemin=ifelse(agemin>13,NA,agemin),
         agemin=ifelse(agemin<8,NA,agemin))


selcountries <- c("BEL","SWE","ITA","POL","USA","JPN")

bothagemin <- amin %>% 
  filter(country %in% selcountries) %>%
  ggplot(aes(x=year, y=agemin, Group = country, color=country)) + 
  geom_smooth(size = 0.7,se=F)+
  #  scale_color_manual(values=c("#0F8B8D","#EC9A29","#A8201A", "black","red","blue"))+
  facet_grid(vars(sex))+
  theme_minimal() +
  theme(legend.position="bottom")+
  scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
  # scale_y_continuous(name="Years", limits=c(60,90))+  
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))

ggsave(bothagemin, file="C:/Users/Magdalena/demography/philosopyineq/figures/bothagemin.pdf",width = 16, height = 20, units = "cm")



#################################################################################################
##################plots
selcountries <- c("BEL","SWE","ITA","POL","USA","JPN")

bothMALL <- Depprem %>% 
  filter(Country %in% selcountries, Stat_name=="MALL") %>%
  ggplot(aes(x=Year, y=Stat_value, Group = Country, color=Country)) + 
  geom_smooth(size = 0.7,se=F)+
  #  scale_color_manual(values=c("#0F8B8D","#EC9A29","#A8201A", "black","red","blue"))+
  facet_grid(vars(Sex))+
  theme_minimal() +
  theme(legend.position="bottom")+
  scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
  scale_y_continuous(name="Years", limits=c(60,90))+  
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.y = element_text(angle=0))

ggsave(bothMALL, file="C:/Users/Magdalena/demography/philosopyineq/figures/bothMALL.pdf",width = 16, height = 20, units = "cm")



######################################################
#### statistics for the paper
selcountries2 <- c("BEL","SWE","ITA","USA")

fordiff <- Depprem %>%
  filter(Stat_name=="MALL", Country %in% selcountries2) %>%
  dplyr::group_by(Sex,Year) %>%
  dplyr::summarise(min=min(as.numeric(Stat_value)),max=max(as.numeric(Stat_value))) %>%
  mutate(maxmin=max-min) %>%
  mutate()



males <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% selcountries), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="none")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("I"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% selcountries), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="none")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(2,20))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("YPLL"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% selcountries), aes(x=Year, y=I10p, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="none")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0.3,0.5))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("I(10+)"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% selcountries), aes(x=Year, y=edagger10, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="none")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(2,10))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("YPLL(10+)"), nrow=2,ncol=2,common.legend = TRUE, legend="bottom")

ggsave(males, file="C:/Users/Magdalena/demography/philosopyineq/figures/males.pdf",width = 25, height = 25, units = "cm")



females <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% selcountries), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="none")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("I"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% selcountries), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="none")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(2,20))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("YPLL"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% selcountries), aes(x=Year, y=I10p, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="none")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0.3,0.5))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("I(10+)"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% selcountries), aes(x=Year, y=edagger10, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="none")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(2,10))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("YPLL(10+)"), nrow=2,ncol=2,common.legend = TRUE, legend="bottom")

ggsave(females, file="C:/Users/Magdalena/demography/philosopyineq/figures/females.pdf",width = 25, height = 25, units = "cm")





##############################all countries - to be corrected labels and ranges
Country <- sort(unique(Depprem$Country))
Country[Country=="FRATNP"] <- "FRA"
countryloc <- c("O","W","W","E","O","W","E","O","N","S","E","N","W","O","S","E","E","O","N","O","S","O","E","O","E","W","N","O","E","S","E","E","E","N","O")
countries <- as_tibble(cbind(Country,countryloc))

w <- countries %>% filter(countryloc=="W")
e <- countries %>% filter(countryloc=="E")
n <- countries %>% filter(countryloc=="N")
s <- countries %>% filter(countryloc=="S")
o <- countries %>% filter(countryloc=="O")

malesMALL <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% w$Country), aes(x=Year, y=MALL, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(60,90))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% e$Country), aes(x=Year, y=MALL, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(55,90))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% n$Country), aes(x=Year, y=MALL, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(60,90))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% s$Country), aes(x=Year, y=MALL, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(60,90))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(malesMALL, file="C:/Users/Magdalena/demography/philosopyineq/figures/MALLmales.pdf",width = 25, height = 25, units = "cm")



femalesMALL<- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% w$Country), aes(x=Year, y=MALL, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(60,90))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% e$Country), aes(x=Year, y=MALL, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(55,90))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% n$Country), aes(x=Year, y=MALL, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(60,90))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% s$Country), aes(x=Year, y=MALL, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(60,90))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(femalesMALL, file="C:/Users/Magdalena/demography/philosopyineq/figures/MALLfemales.pdf",width = 25, height = 25, units = "cm")




malesI <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% w$Country), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% e$Country), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% n$Country), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% s$Country), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(malesI, file="C:/Users/Magdalena/demography/philosopyineq/figures/Ipremmales.pdf",width = 25, height = 25, units = "cm")




femalesI <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% w$Country), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% e$Country), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% n$Country), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% s$Country), aes(x=Year, y=Iprem, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0.3,0.7))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(femalesI, file="C:/Users/Magdalena/demography/philosopyineq/figures/Ipremfemales.pdf",width = 25, height = 25, units = "cm")



femalesedagger <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% w$Country), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(3,22))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% e$Country), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(3,22))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% n$Country), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(3,22))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% s$Country), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(3,22))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(femalesedagger, file="C:/Users/Magdalena/demography/philosopyineq/figures/femalesdagger.pdf",width = 25, height = 25, units = "cm")





malesedagger <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% w$Country), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(3,22))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% e$Country), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(3,22))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% n$Country), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(3,22))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% s$Country), aes(x=Year, y=edagger, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(3,22))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(malesedagger, file="C:/Users/Magdalena/demography/philosopyineq/figures/malesdagger.pdf",width = 25, height = 25, units = "cm")


femalesIshare <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% w$Country), aes(x=Year, y=Icontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0,0.5))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% e$Country), aes(x=Year, y=Icontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.5))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% n$Country), aes(x=Year, y=Icontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.5))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% s$Country), aes(x=Year, y=Icontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.5))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(femalesIshare, file="C:/Users/Magdalena/demography/philosopyineq/figures/Icontr010females.pdf",width = 25, height = 25, units = "cm")


malesIshare <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% w$Country), aes(x=Year, y=Icontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0,0.5))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% e$Country), aes(x=Year, y=Icontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.5))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% n$Country), aes(x=Year, y=Icontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.5))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% s$Country), aes(x=Year, y=Icontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.5))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(malesIshare, file="C:/Users/Magdalena/demography/philosopyineq/figures/Icontr010males.pdf",width = 25, height = 25, units = "cm")


femalesedshare <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% w$Country), aes(x=Year, y=edcontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% e$Country), aes(x=Year, y=edcontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% n$Country), aes(x=Year, y=edcontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.7))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Female",Country %in% s$Country), aes(x=Year, y=edcontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.7))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(femalesedshare, file="C:/Users/Magdalena/demography/philosopyineq/figures/edcontr010females.pdf",width = 25, height = 25, units = "cm")


malesedshare <- ggarrange(
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% w$Country), aes(x=Year, y=edcontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    scale_y_continuous(name="Years", limits=c(0,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("West"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% e$Country), aes(x=Year, y=edcontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    #   theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.7))+  
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("East"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% n$Country), aes(x=Year, y=edcontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.7))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("North"),
  ggplot(data=Depprem %>% filter(Sex=="Male",Country %in% s$Country), aes(x=Year, y=edcontr010, Group = Country, color=Country)) + 
    geom_smooth(size = 0.7,se=F)+
    theme_minimal() +
    theme(legend.position="bottom")+
    scale_y_continuous(name="Years", limits=c(0,0.7))+ 
    scale_x_continuous(name="Year",limits=c(1920,2020), breaks=seq(from=1920,to=2020,by=20))+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.y = element_text(angle=0))+
    ggtitle("South"),  
  nrow=2,ncol=2)

ggsave(malesedshare, file="C:/Users/Magdalena/demography/philosopyineq/figures/edcontr010males.pdf",width = 25, height = 25, units = "cm")


