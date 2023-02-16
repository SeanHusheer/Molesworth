######This code provides data summaries, species lists and plot locations for Molesworth plots measured between 2006, 2007 and 2008, and re-measured between Jan and March 2016.
#####The 2006 survey consisted of the HFI plots established in 1989, the 2007 survey 80 randomly located plots, and the 2008 survey paired fenced plots.

###Written by Sean Husheer October 2015  #######
rm(list=ls()) #Clear the workspace of all objects so R has a clean start
setwd("/home/sean/Documents/Molesworth") # code for LINUX - sets working directory where data is stored
RPackages = "/home/sean/.RPackages" ; .libPaths(RPackages)
##remove.packages(c("mapdata", "bartMachineJARs"))
## install.packages(c("bartMachine"),  RPackages, repos = "https://cran.stat.auckland.ac.nz",  dependencies = TRUE, type = "source")
##install.packages(c("bartMachineJARs", "bartMachine"), lib=RPackages, repos=c("http://glmmadmb.r-forge.r-project.org/repos",   getOption("repos")),  type="source") ###Need "R2admb" library(glmmADMB); 
## update.packages(lib.loc = RPackages, repos = "https://cran.stat.auckland.ac.nz",  ask = FALSE, clean = TRUE, dependencies = TRUE)
options(width = 120)
library("glmmTMB"); library(buildmer); library(DHARMa); library(Hmisc); library(rdddr); library(plotrix)
library(rgrass7); library(rgeos); library(maptools); library(spatstat); library(rgdal); library(maps); library(raster) ##library(geoR) ###Spatial libraries
library(xtable); library(gridExtra); library(plotrix); library(jpeg); library(gplots); library(stargazer); library(tables); library(plyr); library(car); library(latex2exp) ###library(RColorBrewer) ### 
library(pgirmess); library(nlme); library (lme4); library(arm); library(reshape); library(contrast); library(effects); library(corrgram)  ##Mixed Model Packages.  ARM from Gelman and Hill - Applied Regression Models.  
library(vegan); library(R2jags) ; library(mgcv); library(fso); library(diverse); library(codyn); library(dplyr)  ###Libraries for multivariate vegetation
source("/home/sean/Documents/Molesworth/Data/glmmLatex.r") ###Import functon to use xtable for glmm summary tables. function = glmmTMB.latex
s.e. <- function(x) sqrt(var(na.omit(x))/length(na.omit(x))); l.s.d <-  function(x) (sqrt(var(na.omit(x))/length(na.omit(x))))*(qt(0.95,(length(na.omit(x))))) ## SEM  and the more conservative least significant difference
###load("/home/sean/Documents/Molesworth/Data/Objects/HFIcommonMeans.Table.Rdata") ## Object = hfi.common.means (includes lmer coefficents by spp), lmer.commonspp.coe
options(width = 200)
#########################################################################################################################################################################
  #####################################################################################################################################################################                                #################################################        Recce Data    ########################################################
         ################################################################################################################################################################################################################################   Cover Data from RECCE  plots     ######################################################################
             #######################################                  2016 data                           #############################################

## system("ls /home/sean/Documents/Molesworth2016/Data/Molesworth2016random >> /home/sean/Documents/Molesworth2016/plotlist.txt")
## random.plots.entered <- read.table("/home/sean/Documents/Molesworth2016/plotlist.txt", header=FALSE, skip = 0, colClasses="character", stringsAsFactors=FALSE)
random.plots.entered <- dir("/home/sean/Documents/Molesworth/Data/Molesworth2016random") # better than above two lines
omit.plots <- c("B-2A.csv", "P-2A.csv") ;  random.plots.entered <- random.plots.entered[-which(random.plots.entered %in% omit.plots)] ##i <- "A-1.csv"
for(i in random.plots.entered) {
temp.dat <- read.csv(paste("/home/sean/Documents/Molesworth/Data/Molesworth2016random/", i, sep=""), sep = ",", header=TRUE, skip = 0, stringsAsFactors=FALSE)
temp.dat <- temp.dat[duplicated(temp.dat$Code)==FALSE, ] ###removes the second duplicate if there is one
temp.dat$Plot <- gsub(".csv", "",i); temp.dat <- temp.dat[, c(1, 20, 2:19)]
write.table(temp.dat, file = paste("/home/sean/Documents/Molesworth/Data/WithName/", i, sep=""), sep = ",", na = "0", dec = ".", row.names = FALSE, col.names = FALSE)
}     

paired.plots.entered <- dir("/home/sean/Documents/Molesworth/Data/Molesworth2016paired") # better than above two lines
for(i in paired.plots.entered) {
temp.dat <- read.csv(paste("/home/sean/Documents/Molesworth/Data/Molesworth2016paired/", i, sep=""), sep = ",", header=TRUE, skip = 0, stringsAsFactors=FALSE)
temp.dat <- temp.dat[duplicated(temp.dat$Code)==FALSE, ] ###removes the second duplicate if there is one
temp.dat$Plot <- gsub(".csv", "",i); temp.dat <- temp.dat[, c(1, 20, 2:19)]
write.table(temp.dat, file = paste("/home/sean/Documents/Molesworth/Data/WithName/", i, sep=""), sep = ",", na = "0", dec = ".", row.names = FALSE, col.names = FALSE)
}     

recce.2016 <- do.call(rbind,lapply(paste("/home/sean/Documents/Molesworth/Data/WithName/", c(random.plots.entered, paired.plots.entered), sep=""), read.csv, header=FALSE))
names(recce.2016) <- c("code",  "plot", "1<.1", "1<.3", "1<1",  "2<.1", "2<.3", "2<1",  "3<.1", "3<.3", "3<1",  "4<.1", "4<.3", "4<1",  "5<.1", "5<.3", "5<1",  "6<.1", "6<.3", "6<1")

spp.list <- read.csv("~/Documents/Molesworth/Data/Corrections/MolesworthSpeciesListNZFS.csv", sep = ",", header=TRUE, skip = 5, colClasses="character")
names(spp.list) <- c("code", "family", "species", "native", "form","annual", "CommonName", "Book", "Page","Status", "Notes")
spp.list[spp.list$BioStatus== "Indigenous Endemic" | spp.list$BioStatus== "Indigenous Non-Endemic", "BioStatus"] <- "Native"

options(useFancyQuotes = FALSE)
recce.2016.corrections <- read.csv("~/Documents/Molesworth/Data/Corrections/SpeciesCodeCorrections2016.csv", sep = ",", header=TRUE, skip = 5, colClasses="character", stringsAsFactors=FALSE)
recce.2016.corrections$OldCode <- toupper(recce.2016.corrections$OldCode); recce.2016.corrections$NewCode  <- toupper(recce.2016.corrections$NewCode) 

text.code <- paste("recce.2016[recce.2016$plot==", dQuote(recce.2016.corrections$Plot), " & recce.2016$code==", dQuote(recce.2016.corrections$OldCode),
                    ",  ", dQuote("code")," ]  <-   ", dQuote(recce.2016.corrections$NewCode), sep="")
write.table(text.code, file = "/home/sean/Documents/Molesworth/batch.R", quote = FALSE, sep = ",", na = "", dec = ".", row.names = FALSE, col.names = FALSE)#
source("/home/sean/Documents/Molesworth/batch.R")
recce.2016$plot <- as.character(recce.2016$plot); recce.2016$code <- as.character(recce.2016$code) ##to prepare for merge 
recce.2016 <- merge(recce.2016, spp.list[, c("code", "species", "family", "native", "form", "annual")], by= "code", all.x=TRUE, all.y=FALSE)
spp.to.check <- recce.2016[is.na(recce.2016$species)==TRUE,] ##find which codes don't match with species list
#write.table(spp.to.check[ order(spp.to.check$plot),c(1:2)], file = "~/Documents/Molesworth2016/SpeciesToCheck.txt", row.names = FALSE, col.names = FALSE)
recce.2016[recce.2016$native=="Indigenous Endemic"  | recce.2016$native=="Indigenous Non-Endemic", "native" ] <- "Native"

non.spp.cover <- c("Vegetation","Moss","Litter","Bare Ground","Rock", "VEGETATION", "MOSS","LITTER","BARE GROUND","ROCK", "LICHEN")

## ###Code to include tiers higher than 1 m
## recce.high <- read.csv("~/Documents/Molesworth2016/Data/RandomPlotCoverScores2016.csv", sep = ",", header=TRUE, skip = 5) ##, stringsAsFactors=FALSE) ### Import the higher tiers.

## names(recce.high) <- c("plot", "code", "1<.1", "1<.3", "1<1", "1<2", "1<5", "1>5", "2<.1", "2<.3", "2<1", "2<2", "2<5", "2>5", "3<.1", "3<.3", "3<1", "3<2", "3<5", "3>5", "4<.1", "4<.3", "4<1", "4<2", "4<5", "4>5", "5<.1", "5<.3", "5<1", "5<2", "5<5", "5>5",  "6<.1", "6<.3", "6<1", "6<2", "6<5", "6>5")
## recce.high$code <- toupper(recce.high$code) ##prepare for merge
## recce.high[c("1<2", "1<5", "1>5", "2<2", "2<5", "2>5", "3<2", "3<5", "3>5", "4<2", "4<5", "4>5", "5<2", "5<5", "5>5", "6<2", "6<5", "6>5")][is.na(recce.high[c("1<2", "1<5", "1>5", "2<2", "2<5", "2>5", "3<2", "3<5", "3>5", "4<2", "4<5", "4>5", "5<2", "5<5", "5>5", "6<2", "6<5", "6>5")])] <- 0 ##change na to zero
##recce.high <- recce.high[-which(recce.high$Code %in% non.spp.cover),]  ##remove unusable data
## unfactor <- c("1<2", "1<5", "1>5", "2<2", "2<5", "2>5", "3<2", "3<5", "3>5", "4<2", "4<5", "4>5", "5<2", "5<5", "5>5", "6<2", "6<5", "6>5")
## recce.high[,unfactor] <- lapply(unfactor, function(x) as.numeric(as.character(recce.high[,x])))
## recce.2016 <- (merge(recce.2016, recce.high[, c("plot", "code", "1<2", "1<5", "1>5", "2<2", "2<5", "2>5", "3<2", "3<5", "3>5", "4<2", "4<5", "4>5", "5<2", "5<5", "5>5", "6<2", "6<5", "6>5")], by=c("plot", "code"), all.x=TRUE, all.y=TRUE))

## for (i in 1:6){
## assign(paste("cover", i, sep="."), rowSums(recce.2016[ ,c(paste(i, "<.1", sep=""), paste(i, "<.3", sep=""), paste(i,"<1", sep=""), paste(i,"<2", sep=""), paste(i,"<5", sep=""), paste(i,">5", sep=""))], na.rm=TRUE))
## }
## ###Produce a data frame that has sums of all height tiers
## year <- rep("2016", times=length(recce.2016$plot))
## recce.16 <- as.data.frame(cbind(year, recce.2016$plot, recce.2016$code, cover.1, cover.2,  cover.3,  cover.4,  cover.5, cover.6),  stringsAsFactors = FALSE)
## unfactor <- c("year", "cover.1", "cover.2",  "cover.3",  "cover.4",  "cover.5", "cover.6"); names(recce.16) <- c("year", "plot", "code", "sub.1","sub.2","sub.3","sub.4","sub.5","cover") 

year <- rep("2016", times=length(recce.2016$plot))
recce.16 <- cbind(year, recce.2016[, c("plot", "species","1<.1", "2<.1", "3<.1", "4<.1", "5<.1", "6<.1", "family",  "native",  "form", "annual")])     
sub.plot <- c("sub.1","sub.2","sub.3","sub.4","sub.5","sub.6")
names(recce.16) <- c("year", "plot", "species", sub.plot, "family",  "native",  "form", "annual")     

###remove duplicates  i <-  "A-1"
for (i in unique(recce.16$plot)){
temp.dat <- recce.16[recce.16$plot==i, ]
write.table(temp.dat[duplicated(temp.dat$species)==FALSE, ], file = paste("/home/sean/Documents/Molesworth/Data/WithName/", i, ".csv", sep=""), sep = ",", na = "0", dec = ".", row.names = FALSE, col.names = FALSE)
}     
recce.016 <- do.call(rbind,lapply(paste("/home/sean/Documents/Molesworth/Data/WithName/", unique(recce.16$plot), ".csv", sep=""), read.csv, header=FALSE))
names(recce.016) <- names(recce.16)


##################2007 random plot cover score data
recce.2007 <- read.csv("~/Documents/Molesworth/Data/Molesworth2007Cover.csv", sep = ",", header=TRUE, skip = 5, stringsAsFactors=FALSE)
recce.2007$plot <- paste(recce.2007$Transect, recce.2007$Plot, sep="-")
names(recce.2007) <- c("Transect", "Plot", "code", "1<.1", "1<.3", "1<1", "1<2", "1<5", "1>5", "2<.1", "2<.3", "2<1", "2<2", "2<5", "2>5", "3<.1", "3<.3", "3<1", "3<2", "3<5", "3>5", "4<.1", "4<.3", "4<1", "4<2", "4<5", "4>5", "5<.1", "5<.3", "5<1", "5<2", "5<5", "5>5",  "6<.1", "6<.3", "6<1", "6<2", "6<5", "6>5", "browser",  "BrowserIntensity", "plot")  
## for (i in 1:6){
## assign(paste("cover", i, sep="."), rowSums(recce.2007[ ,c(paste(i, "<.1", sep=""), paste(i, "<.3", sep=""), paste(i,"<1", sep=""), paste(i,"<2", sep=""), paste(i,"<5", sep=""), paste(i,">5", sep=""))], na.rm=TRUE))
## }


year <- rep("2007", times=length(recce.2007$plot))
recce.07 <- cbind(year, recce.2007[, c("plot", "code","1<.1", "2<.1", "3<.1", "4<.1", "5<.1", "6<.1")])     
##  as.data.frame(cbind(year, recce.2007$plot, recce.2007$code, cover.1, cover.2,  cover.3,  cover.4,  cover.5, cover.6),  stringsAsFactors = FALSE)
## unfactor <- c("year", "cover.1", "cover.2",  "cover.3",  "cover.4",  "cover.5", "cover.6")
## recce.07[,unfactor] <- lapply(unfactor, function(x) as.integer(as.character(recce.07[,x])))
## names(recce.07) <- c("year", "plot", "code", "sub.1", "sub.2", "sub.3", "sub.4", "sub.5", "cover")## Six plot sizes from o.25m2 to full 20 x20.
recce.07 <- recce.07[-which(recce.07$code %in% non.spp.cover),]  ##remove unusable data

recce.2007.corrections <- read.csv("~/Documents/Molesworth/Data/Corrections/RandomPlotSpeciesCodeCorrections2007.csv", sep = ",", header=TRUE, skip = 5, colClasses="character", stringsAsFactors=FALSE)
recce.2007.corrections$OldCode <- toupper(recce.2007.corrections$OldCode); recce.2007.corrections$NewCode  <- toupper(recce.2007.corrections$NewCode) 

text.code <- paste("recce.07[recce.07$plot==", dQuote(recce.2007.corrections$Plot), " & recce.07$code==", dQuote(recce.2007.corrections$OldCode),
                   ",  ", dQuote("code")," ]  <-   ", dQuote(recce.2007.corrections$NewCode), sep="")
write.table(text.code, file = "/home/sean/Documents/Molesworth/batch.R", quote = FALSE, sep = ",", na = "", dec = ".", row.names = FALSE, col.names = FALSE)#
source("/home/sean/Documents/Molesworth/batch.R")

recce.2007.correctionsTWO <- read.csv("~/Documents/Molesworth/Data/Corrections/SpeciesCodeCorrections2007edited.csv", sep = ",", header=TRUE, skip = 5, colClasses="character", stringsAsFactors=FALSE)
recce.2007.correctionsTWO$OldCode <- toupper(recce.2007.correctionsTWO$OldCode); recce.2007.correctionsTWO$NewCode  <- toupper(recce.2007.correctionsTWO$NewCode) 
text.code <- paste("recce.07[recce.07$plot==", dQuote(recce.2007.correctionsTWO$Plot), " & recce.07$code==", dQuote(recce.2007.correctionsTWO$OldCode),
                   ",  ", dQuote("code")," ]  <-   ", dQuote(recce.2007.correctionsTWO$NewCode), sep="")
write.table(text.code, file = "/home/sean/Documents/Molesworth/batch.R", quote = FALSE, sep = ",", na = "", dec = ".", row.names = FALSE, col.names = FALSE)#
source("/home/sean/Documents/Molesworth/batch.R")

recce.07 <- merge(recce.07, spp.list[, c("code", "species", "family", "native", "form", "annual")], by= "code", all.x=TRUE, all.y=FALSE)
spp.to.check.07 <- recce.07[is.na(recce.07$species)==TRUE,] ##find which codes don't match with species list
recce.07 <- recce.07[is.na(recce.07$species)==FALSE,] ##remove a couple of blank species
names(recce.07) <-  c("code", "year", "plot",  sub.plot, "species", "family",  "native",  "form",    "annual") 

###remove duplicates  i <-  "A-1"
for (i in unique(recce.07$plot)){
temp.dat <- recce.07[recce.07$plot==i, ]
write.table(temp.dat[duplicated(temp.dat$species)==FALSE, ], file = paste("/home/sean/Documents/Molesworth/Data/WithName/", i, ".csv", sep=""), sep = ",", na = "0", dec = ".", row.names = FALSE, col.names = FALSE)
}     
recce.007 <- do.call(rbind,lapply(paste("/home/sean/Documents/Molesworth/Data/WithName/", unique(recce.07$plot), ".csv", sep=""), read.csv, header=FALSE))
names(recce.007) <- names(recce.07)


##################2008 paired plot cover score data
recce.2008 <- read.csv("~/Documents/Molesworth/Data/Molesworth2008Cover.csv", sep = ",", header=TRUE, skip = 5, stringsAsFactors=FALSE)
recce.2008$Plot <- paste(recce.2008$Transect, recce.2008$Plot, sep="-")
recce.2008$Plot <- sub(" unfenced", "-Un", recce.2008$Plot); recce.2008$Plot <- sub(" fenced", "-Fe", recce.2008$Plot)

names(recce.2008) <- c("Transect", "Plot", "code", "1<.1", "1<.3", "1<1", "1<2", "1<5", "1>5", "2<.1", "2<.3", "2<1", "2<2", "2<5", "2>5", "3<.1", "3<.3", "3<1", "3<2", "3<5", "3>5", "4<.1", "4<.3", "4<1", "4<2", "4<5", "4>5", "5<.1", "5<.3", "5<1", "5<2", "5<5", "5>5",  "6<.1", "6<.3", "6<1", "6<2", "6<5", "6>5", "browser",  "BrowserIntensity", "Browser2Species", "Browser2Intensity", "plot")  
## for (i in 1:6){
## assign(paste("cover", i, sep="."), rowSums(recce.2008[ ,c(paste(i, "<.1", sep=""), paste(i, "<.3", sep=""), paste(i,"<1", sep=""), paste(i,"<2", sep=""), paste(i,"<5", sep=""), paste(i,">5", sep=""))], na.rm=TRUE))
## }

year <- rep("2008", times=length(recce.2008$plot))
recce.08 <- cbind(year, recce.2008[, c("plot", "code","1<.1", "2<.1", "3<.1", "4<.1", "5<.1", "6<.1")])     

recce.08 <- recce.08[-which(recce.08$code %in% non.spp.cover),]  ##remove unusable data

recce.2008.corrections <- read.csv("~/Documents/Molesworth/Data/Corrections/PairedPlotSpeciesCodeCorrections2008.csv", sep = ",", header=TRUE, skip = 5, colClasses="character", stringsAsFactors=FALSE)
recce.2008.corrections$OldCode <- toupper(recce.2008.corrections$OldCode); recce.2008.corrections$NewCode  <- toupper(recce.2008.corrections$NewCode) 
recce.2008.corrections$Plot <- sub("Un", "-Un", recce.2008.corrections$Plot); recce.2008.corrections$Plot <- sub("Fe", "-Fe", recce.2008.corrections$Plot)

text.code <- paste("recce.08[recce.08$plot==", dQuote(recce.2008.corrections$Plot), " & recce.08$code==", dQuote(recce.2008.corrections$OldCode),
                   ",  ", dQuote("code")," ]  <-   ", dQuote(recce.2008.corrections$NewCode), sep="")
write.table(text.code, file = "/home/sean/Documents/Molesworth/batch.R", quote = FALSE, sep = ",", na = "", dec = ".", row.names = FALSE, col.names = FALSE)#
source("/home/sean/Documents/Molesworth/batch.R")

recce.08 <- merge(recce.08, spp.list[, c("code", "species", "family", "native", "form", "annual")], by= "code", all.x=TRUE, all.y=FALSE)
recce.08 <- recce.08[is.na(recce.08$species)==FALSE,] ##remove a blank species
spp.to.check.08 <- recce.08[is.na(recce.08$species)==TRUE,] ##find which codes don't match with species list
names(recce.08) <-  c("code", "year", "plot",  sub.plot, "species", "family",  "native",  "form",    "annual") 

for (i in unique(recce.08$plot)){
temp.dat <- recce.08[recce.08$plot==i, ]
write.table(temp.dat[duplicated(temp.dat$species)==FALSE, ], file = paste("/home/sean/Documents/Molesworth/Data/WithName/", i, ".csv", sep=""), sep = ",", na = "0", dec = ".", row.names = FALSE, col.names = FALSE)
}     
recce.008 <- do.call(rbind,lapply(paste("/home/sean/Documents/Molesworth/Data/WithName/", unique(recce.08$plot), ".csv", sep=""), read.csv, header=FALSE))
names(recce.008) <- names(recce.08)

##############################################################   1989 Cover Data from HFI plots     ######################################################################

hfi.1989 <- read.csv("~/Documents/Molesworth/Data/Molesworth1989-2006-HFI.csv", sep = ",", header=TRUE, skip = 5)
names(hfi.1989)[1:4] <- c("Year",  "Plot", "Spp", "Code"); hfi.1989$Code <- toupper(hfi.1989$Code) 
hfi.1989.corrections <- read.csv("~/Documents/Molesworth/Data/Corrections/MolesworthSpeciesCodeCorrectionsHFI1989.csv", sep = ",", header=TRUE, skip = 5)
text.code <- paste("hfi.1989[hfi.1989$Plot==", dQuote(hfi.1989.corrections$Plot), " &  hfi.1989$Year==", dQuote(hfi.1989.corrections$Year), " & hfi.1989$Code==", dQuote(hfi.1989.corrections$OldCode),
                   ",  ", dQuote("Code")," ]  <-   ", dQuote(hfi.1989.corrections$NewCode), sep="")
write.table(text.code, file = "/home/sean/Documents/Molesworth/batch.R", quote = FALSE, sep = ",", na = "", dec = ".", row.names = FALSE, col.names = FALSE)#
source("/home/sean/Documents/Molesworth/batch.R")
names(hfi.1989)[names(hfi.1989) == "Code"] <- "code"  ##; names(hfi.1989) <-  c("year", "plot", "code", "Tier1")
hfi.1989 <- (merge(hfi.1989[, c("Year", "Plot", "code", "Tier1")], spp.list[, c("code", "species")], by="code", all.x=TRUE, all.y=FALSE))
hfi.1989$cover <- 1; hfi.1989[which(hfi.1989$Tier1 > 1 & hfi.1989$Tier1 <= 5),"cover"] <- 2; hfi.1989[which(hfi.1989$Tier1 > 5 & hfi.1989$Tier1 <= 25),"cover"] <- 3
hfi.1989[which(hfi.1989$Tier1 > 25 & hfi.1989$Tier1 <= 50),"cover"] <- 4; hfi.1989[which(hfi.1989$Tier1 > 50 & hfi.1989$Tier1 <= 75),"cover"] <- 5;
hfi.1989[hfi.1989$Tier1 >=75,"cover"] <- 6; hfi.1989 <- hfi.1989[, c("Year","Plot","species","cover")]; names(hfi.1989) <- c("year","plot","species","cover")
hfi.1989$plot <- sub("S", "Saxton", hfi.1989$plot); hfi.1989$plot <- sub("M", "Molesworth", hfi.1989$plot)

locs.1989 <- read.csv("~/Documents/Molesworth/Data/SiteData2016.csv", sep = ",", header=TRUE, skip = 5)
locs.1989$plot <- paste(locs.1989$Line, locs.1989$Plot, sep="" )
locs.1989 <- locs.1989[locs.1989$plot %in% unique(hfi.1989$plot), c("plot", "East", "North")]
locs.1989 <- locs.1989[is.na(locs.1989$East)==FALSE, ]; names(locs.1989) <- c("plot", "x", "y")


##############################################################   1987 PNA  Cover Data from RECCE plots     ######################################################################
locs.1987 <- as.data.frame(read.csv("~/Documents/Molesworth/Data/Molesworth1987ReccesNVS/Coordinates.csv", sep = ",", header=TRUE, skip = 0, colClasses="character")[, c("Plot","EastingMG","NorthingMG")], stringsAsFactors=FALSE)
names(locs.1987) <- c("plot", "x", "y"); locs.1987$x <- as.numeric(locs.1987$x); locs.1987$y <- as.numeric(locs.1987$y); locs.1987 <- locs.1987[is.na(locs.1987$x)==FALSE, ]
coordinates(locs.1987) <- c("x", "y"); proj4string(locs.1987) <- CRS("+init=epsg:27200") #georeference to NZMG ##CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #georeference to for GE
locs.1987 <- spTransform(locs.1987, CRS("+init=epsg:2193")); locs.1987.spatial.points <- locs.1987

recce.locs.Tier1 <- as.data.frame(read.csv("~/Documents/Molesworth/Data/MolesworthRecceTier1/Coordinates.csv", sep = ",", header=TRUE, skip = 0, colClasses="character")[, c("Plot", "AbsoluteCoordXEastLong", "AbsoluteCoordYNorthLat", "MapProjection")], stringsAsFactors=FALSE); names(recce.locs.Tier1) <-  c("plot", "x", "y", "MapProjection")
recce.locs.Tier1 <- recce.locs.Tier1[recce.locs.Tier1$MapProjection=="NZTM2000",]
recce.locs.Tier1 <- recce.locs.Tier1[is.na(recce.locs.Tier1$x)==FALSE & recce.locs.Tier1$x!="", c("plot", "x", "y")]
recce.locs.Tier1$x <- as.numeric(recce.locs.Tier1$x);  recce.locs.Tier1$y <- as.numeric(recce.locs.Tier1$y)  
coordinates(recce.locs.Tier1) <- c("x", "y") #save cordinates in prep for spatial conversion 
proj4string(recce.locs.Tier1) <- CRS("+init=epsg:2193") #georeference to NZTM

############Recce data from 1987    Use only cover score from the ground tier = tier 6
dat.1987 <- as.data.frame(read.csv("/home/sean/Documents/Molesworth/Data/Molesworth1987ReccesNVS/TaxonCategoryValue.csv",
                                   sep = ",", header=TRUE, skip = 0)[, c("Plot","PreferredSpeciesName","Tier","Category")], stringsAsFactors=FALSE)
names(dat.1987) <- c("plot", "species", "tier", "cover")  ###Tier 6 = 0-0.1 m, Tier 5 = 0.1 - 0.3 m, Tier 4 = 0.3 - 2 m, Tier 3 = 2 - 5 m, Tier 2 = 5-12 m, Tier 1 = >12 m
dat.1987[grep("P", dat.1987$cover), "cover"] <- 1; dat.1987$cover <- as.integer(dat.1987$cover)  ##length(unique(dat.1987$plot))
##dat.1987.tiers.summed <- as.data.frame.table(tapply(dat.1987$cover, list(dat.1987$plot, dat.1987$species), sum, na.rm=TRUE), stringsAsFactors=FALSE); names(dat.1987) <- c("plot", "species", "sum.cover")
dat.1987$year <-  1987
dat.1987 <- dat.1987[is.na(dat.1987$cover)==FALSE & dat.1987$tier==6, c("year", "plot", "species", "cover")]; dat.1987 <- dat.1987[which(dat.1987$plot %in% locs.1987$plot),]  ###select plots found within Molesworth boundary
dat.1987 <- dat.1987[is.na(dat.1987$species)==FALSE,]
dat.1987$species <- as.character(dat.1987$species); dat.1987$plot <- as.character(dat.1987$plot) ##to prepare for corrections 
recce.1987.corrections <- read.csv("~/Documents/Molesworth/Data/Corrections/SpeciesCodeCorrections1987.csv", sep = ",", header=TRUE, skip = 5, colClasses="character", stringsAsFactors=FALSE)
text.code <- paste("dat.1987[dat.1987$species==", dQuote(recce.1987.corrections$OldCode), ",  ", dQuote("species")," ]  <-   ", dQuote(recce.1987.corrections$NewCode), sep="")
write.table(text.code, file = "/home/sean/Documents/Molesworth/batch.R", quote = FALSE, sep = ",", na = "", dec = ".", row.names = FALSE, col.names = FALSE)#
source("/home/sean/Documents/Molesworth/batch.R")
dat.1987 <- dat.1987[dat.1987$species!="NonVascular", ]


##############Tier 1 cover data Use only cover score from the ground tier = tier 6
dat.tier1.a <- as.data.frame(read.csv("/home/sean/Documents/Molesworth/Data/MolesworthRecceTier1/LUCAS Vascular Recce (27268) 20.csv")[, c("Plot","PreferredSpeciesName", "Cover.class.Tier.1", "Cover.class.Tier.2", "Cover.class.Tier.3", "Cover.class.Tier.4", "Cover.class.Tier.5", "Cover.class.Tier.6A")], sep = ",", header=TRUE, skip = 0, stringsAsFactors=FALSE)
dat.tier1.b <- as.data.frame(read.csv("/home/sean/Documents/Molesworth/Data/MolesworthRecceTier1/LUCAS Vascular Recce (24227) 20.csv")[, c("Plot","PreferredSpeciesName", "Cover.class.Tier.1", "Cover.class.Tier.2", "Cover.class.Tier.3", "Cover.class.Tier.4", "Cover.class.Tier.5", "Cover.class.Tier.6A")], sep = ",", header=TRUE, skip = 0, stringsAsFactors=FALSE)
dat.tier1 <- rbind(dat.tier1.a, dat.tier1.b); dat.tier1$Cover.class.Tier.6A <- as.character(dat.tier1$Cover.class.Tier.6A); dat.tier1[,c(3:8)]  <- sapply(dat.tier1[,c(3:8)], as.numeric)
##dat.tier1$sum.cover <- rowSums(dat.tier1[,c(3:8)], na.rm=TRUE); dat.tier1 <- dat.tier1[, c("Plot", "PreferredSpeciesName", "sum.cover")]; names(dat.tier1) <- c("plot", "species", "sum.cover")

##Which threatened species in 1987 recces dat.1987[dat.1987$species %in% threatened.species, ]


dat.tier1$year <- 2012  ##arbirtray measurement year
dat.tier1 <- dat.tier1[,c("year", "Plot","PreferredSpeciesName","Cover.class.Tier.6A")] ; names(dat.tier1) <- c("year", "plot", "species", "cover") 
dat.tier1$species <- as.character(dat.tier1$species); dat.tier1$plot <- as.character(dat.tier1$plot) ##to prepare for corrections 
recce.DOCtier1.corrections <- read.csv("~/Documents/Molesworth/Data/Corrections/DOCtier1SpeciesCodeCorrections2007.csv", sep = ",", header=TRUE, skip = 5, colClasses="character", stringsAsFactors=FALSE)
text.code <- paste("dat.tier1[dat.tier1$species==", dQuote(recce.DOCtier1.corrections$OldCode), ",  ", dQuote("species")," ]  <-   ", dQuote(recce.DOCtier1.corrections$NewCode), sep="")
write.table(text.code, file = "/home/sean/Documents/Molesworth/batch.R", quote = FALSE, sep = ",", na = "", dec = ".", row.names = FALSE, col.names = FALSE)#
source("/home/sean/Documents/Molesworth/batch.R")

## for (i in unique(dat.tier1$plot)){  ###Check for duplicates
## print(dat.tier1[duplicated(dat.tier1[dat.tier1$plot== i, "species"])==TRUE, ]) }

recces <- rbind(dat.tier1, dat.1987); recces$species <- gsub("^([^ ]* [^ ]*) .*$", "\\1", recces$species)
recces <-  recces[-(grep("species", recces$species)),]  ##; recces <-recces[which(recces$plot %in% locs.1987$plot), ];
recces[recces$species== "Craspedia 'woolly", "species"] <- "Craspedia lanata" # could also be C.unifora or C. incana
recces[is.na(recces$cover)==TRUE, "cover"] <- 1    ;  recces$plot <- as.factor(recces$plot)[drop=TRUE]


###################################################### Wraight cover data WAIRAU GRASSLAND 1959-1960 and  1972-1973

wairau.plots <- c("F1 1 (F1 1_TRA)", "F1 2 (F1 2_TRA)", "F1 3 (F1 3_TRA)", "F2 1 (F2 1_TRA)", "F2 2 (F2 2_TRA)", "F3 1 (F3 1_TRA)", "F3 2 (F3 2_TRA)", "F4 1 (F4 1_TRA)", "F4 2 (F4 2_TRA)", "F5 1 (F5 1_TRA)")
wairau.locs <- as.data.frame(read.csv("~/Documents/Molesworth/Data/WraightWairauSurveys/Coordinates.csv", sep = ",", header=TRUE, skip = 0)) ##, stringsAsFactors=FALSE) ### Import the higher tiers.
wairau.locs <- wairau.locs[wairau.locs$Plot %in% wairau.plots, c("Plot","Project", "AbsoluteCoordXEastLong", "AbsoluteCoordYNorthLat")]
names(wairau.locs) <- c("plot", "project", "x", "y"); wairau.locs <- wairau.locs[wairau.locs$project == "WAIRAU GRASSLAND 1972-1973", c("plot", "x", "y")] 
wairau.dat <- as.data.frame(read.csv("~/Documents/Molesworth/Data/WraightWairauSurveys/TaxonSimpleValue.csv", sep = ",", header=TRUE, skip = 0)) ##, stringsAsFactors=FALSE) ### Import the higher tiers.
wairau.dat <- wairau.dat[wairau.dat$Plot %in% wairau.plots, c("Plot","Project", "PreferredSpeciesName", "Value")]
names(wairau.dat) <- c("plot", "project", "species", "value")            
wairau.dat <-  wairau.dat[-which(wairau.dat$species %in% c("(Unknown)", "Lichen species", "Moss species")), ]                           
wairau.dat$species <- as.character(wairau.dat$species); wairau.dat$plot <- as.character(wairau.dat$plot) ##to prepare for corrections 
wairau.corrections <- read.csv("~/Documents/Molesworth/Data/Corrections/SpeciesCodeCorrectionsWraight60.csv", sep = ",", header=TRUE, skip = 5, colClasses="character", stringsAsFactors=FALSE)
text.code <- paste("wairau.dat[wairau.dat$species==", dQuote(wairau.corrections$OldCode), ",  ", dQuote("species")," ]  <-   ", dQuote(wairau.corrections$NewCode), sep="")
write.table(text.code, file = "/home/sean/Documents/Molesworth/batch.R", quote = FALSE, sep = ",", na = "", dec = ".", row.names = FALSE, col.names = FALSE)#
source("/home/sean/Documents/Molesworth/batch.R")
wairau.dat$year <- 1960; wairau.dat[wairau.dat$project=="WAIRAU GRASSLAND 1972-1973", "year"] <- 1972 ##December 1959 to Feb 1960, November 1972 to Feb 1973.
wairau.dat$cover <- 1; wairau.dat[which(wairau.dat$value > 1 & wairau.dat$value <= 5),"cover"] <- 2; wairau.dat[which(wairau.dat$value > 5 & wairau.dat$value <= 25),"cover"] <- 3
wairau.dat[which(wairau.dat$value > 25 & wairau.dat$value <= 50),"cover"] <- 4; wairau.dat[which(wairau.dat$value > 50 & wairau.dat$value <= 75),"cover"] <- 5;
wairau.dat[wairau.dat$value >=75,"cover"] <- 6; wairau.dat <- wairau.dat[, c("year","plot","species","cover")]  

##################  Moore releve data
recce.1952 <- as.data.frame(read.csv("~/Documents/Molesworth/Data/Moore1952releve.csv", sep = ",", header=TRUE, skip = 5)) ##, stringsAsFactors=FALSE) ### Import the higher tiers.
names(recce.1952) <- c("species", paste("Moore-", seq(1:24), sep=""), "comment", "table.order")
recce.1952 <- merge(recce.1952, spp.list[, c("species", "family", "native", "form")], by = "species", all.x=TRUE, all.y=FALSE)     
recce.1952 <- melt(recce.1952[,c(1:25)], id=c("species"), na.rm=TRUE); recce.1952$year <- 1952; names(recce.1952) <- c("species", "plot", "cover", "year")

##############################################################   Change in species composition  ######################################################################
######################### All species ground tier cover scores from 1952 to 2016 - 64 years of change ###################################################### 

recce.dat <- rbind(recce.016[, c("year", "plot", "species",  "sub.6")], recce.008[, c("year", "plot", "species",  "sub.6")], recce.007[, c("year", "plot", "species",  "sub.6")])
names(recce.dat) <- names(recces) ## prepare for bind
recce.dat <- rbind(recce.dat, hfi.1989, recces, recce.1952[,names(recce.dat)], wairau.dat)  ##bind all data together
recce.dat[recce.dat$cover==0, "cover"] <- 1  ###zeros recorded but species probably presnet in the plot
recce.dat <- merge(recce.dat, spp.list[,c("species", "family", "native","form")], by ="species", all.x=TRUE, all.y=FALSE)

####

no.list <- sort(as.character(unique(recce.dat[is.na(recce.dat$family)==TRUE, "species"]))) ##Which species are not in the NZFS spplist???
nvs.spp <- read.csv("/home/sean/Documents/Molesworth/Data/CurrentNVSNames.csv", sep = ",", header=TRUE, skip = 0, stringsAsFactors=FALSE)
nvs.spp <- nvs.spp[nvs.spp$SpeciesName %in% no.list, c("SpeciesName", "Family", "BioStatus", "GrowthForm")] ; names(nvs.spp) <- c("species", "family", "native", "form")
recce.dat.list <- recce.dat[-which(recce.dat$species %in% no.list),]
recce.dat.no.list <- merge(recce.dat[recce.dat$species %in% no.list, c("species", "year", "plot", "cover")], nvs.spp, by="species", all.x=TRUE, all.y=TRUE)
recce.dat.no.list[recce.dat.no.list$species=="Pimelea nitens" |  recce.dat.no.list$species=="Pimelea mesoa", "family"] <- "Thymelaeaceae"
recce.dat.no.list[recce.dat.no.list$species=="Pimelea nitens" |  recce.dat.no.list$species=="Pimelea mesoa", "native"] <- "Indigenous Endemic"
recce.dat.no.list[recce.dat.no.list$species=="Pimelea nitens" |  recce.dat.no.list$species=="Pimelea mesoa", "form"] <- "Forb"
recce.dat <- rbind(recce.dat.list, recce.dat.no.list)##; recce.dat[recce.dat$native=="Unknown","native"] <- "Indigenous Endemic"
##m1 <-lme(cover  ~  year + species, random=~year|species, data=recce.dat, na.action=na.omit)  ##




################################################################################################################################################################################
 ################################# Import Polygons to make a grid of envirnomental data  #############################################################################
##system("ogrinfo /home/sean/Documents/Molesworth2016/WinterSummerGrazing/Grazing2016Megerd.shp")
#dem <- readGDAL("/home/sean/Documents/Molesworth2016/Grass/MolesworthDEM.tif")##import using the gdal package
dem <- readGDAL("/home/sean/Documents/Molesworth/Data/Spatial/CroppedDEM.tif")
terrace <- readGDAL("/home/sean/Documents/Molesworth/Data/Spatial/Terrace400.tif") ### image(terrace[terrace$band1==1,]) selects only terrace 
temperature <- readGDAL("/home/sean/Documents/Molesworth/Data/Spatial/MolesworthTemperature.tif")  
lcdb.polygons <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/MolesworthLCDB2012Polygon.shp", "MolesworthLCDB2012Polygon"); proj4string(lcdb.polygons) <- CRS("+init=epsg:2193") 
lenz.soil.polygons <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/SoilLENZ.shp", "SoilLENZ"); proj4string(lenz.soil.polygons) <- CRS("+init=epsg:2193") 
soil.polygons <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/SoilPolygons.kml", "SoilPolygons")
soil.polygons <- spTransform(soil.polygons, CRS("+init=epsg:2193")) #georeference to NZTM
rain.polygons <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/RainfallPolygons.kml", "RainfallPolygons")
rain.polygons <- spTransform(rain.polygons, CRS("+init=epsg:2193")) #georeference to NZTM
grazing2016.polygons <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/Grazing2016Megerd.shp",  "Grazing2016Megerd")
grazing2016.polygons <- spTransform(grazing2016.polygons, CRS("+init=epsg:2193")) #georeference to NZTM
fenced.polygons <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/FencedPolygons.kmz", "FencedPolygons")
fenced.polygons <- spTransform(fenced.polygons, CRS("+init=epsg:2193")) #georeference to NZTM
grazing1987.polygons <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/Grazing_areas_summer_winter_split.shp",  "Grazing_areas_summer_winter_split")
grazing1987.polygons <- spTransform(grazing1987.polygons, CRS("+init=epsg:2193")) #georeference to NZTM
oversow.polygons <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/OversowingPolygon.kml", "OversowingPolygon")
oversow.polygons <- spTransform(oversow.polygons, CRS("+init=epsg:2193")) #georeference to NZTM
oversow.polygons.improved <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/ImprovedLandPolygons.kml", "ImprovedLandPolygons")
oversow.polygons.improved <- spTransform(oversow.polygons.improved, CRS("+init=epsg:2193")) #georeference to NZTM
oversow.polygons <- rbind(oversow.polygons, oversow.polygons.improved) #### Join two spatial data frames 




cell.space <- 100 ##grid spacing
grid <- as.data.frame(cbind(rep(seq(from=xmin(dem), to=xmax(dem), by=cell.space), times=length(seq(from=ymin(dem), to=ymax(dem), by=cell.space))), rep(seq(from=ymin(dem), to=ymax(dem), by=cell.space), each=length(seq(from=xmin(dem), to=xmax(dem), by=cell.space))))) 
names(grid) <- c("x", "y"); coordinates(grid) <- c("x", "y"); proj4string(grid) <- CRS("+init=epsg:2193"); 
grid.alt <- as.data.frame(cbind(round(extract(raster(dem), grid), digits=1), coordinates(grid)), stringsAsFactors=FALSE); names(grid.alt) <- c("dem.altitude", "x", "y")
grid.terrace <- as.data.frame(cbind(extract(raster(terrace), grid), coordinates(grid))); names(grid.terrace) <- c("terrace", "x", "y")
grid.terrace[is.na(grid.terrace)==TRUE, "terrace"] <-  0; grid.terrace$terrace <-  as.character(grid.terrace$terrace);
grid.terrace[grid.terrace$terrace=="1", "terrace"] <- "Terrace"; grid.terrace[grid.terrace$terrace=="0", "terrace"] <- "Slope"
grid.temp <- as.data.frame(cbind(extract(raster(temperature), grid), coordinates(grid)), stringsAsFactors=FALSE); names(grid.temp) <- c("temp", "x", "y"); 
grid.lcdb <- as.data.frame(cbind(over(grid, lcdb.polygons[,1]), coordinates(grid)[,1:2])); names(grid.lcdb) <- c("lcdb.class", "x", "y")
grid.lenz.soil <- as.data.frame(cbind(over(grid, lenz.soil.polygons[,1]), coordinates(grid)[,1:2])); names(grid.lenz.soil) <- c("lenz.soil.class", "x", "y")
grid.soil <- as.data.frame(cbind(over(grid, soil.polygons[,1]), coordinates(grid)[,1:2])); names(grid.soil) <- c("soil.class", "x", "y")  ###Soil classes from MOlesworth book
grid.rain <- as.data.frame(cbind(over(grid, rain.polygons[,1]), coordinates(grid)[,1:2])); names(grid.rain) <- c("rain.class", "x", "y")
grid.grazing1987 <- as.data.frame(cbind(over(grid, grazing1987.polygons[3]), coordinates(grid)[,1:2]), stringsAsFactors=FALSE); names(grid.grazing1987) <- c("grazing1987", "x", "y")##Extract Summer, Winter dataframe from the shape file
grid.grazing1987$grazing1987 <- as.character(grid.grazing1987$grazing1987); grid.grazing1987[is.na(grid.grazing1987)==TRUE, "grazing1987"] <- "Ungrazed"
grid.grazing2016 <- as.data.frame(cbind(over(grid, grazing2016.polygons), coordinates(grid)[,1:2]), stringsAsFactors=FALSE); names(grid.grazing2016) <- c("grazing2016", "x", "y")##Extract Summer, Winter dataframe from the shape file
grid.grazing2016[is.na(grid.grazing2016)==TRUE, "grazing2016"]<-  0; grid.grazing2016$grazing2016 <-  as.numeric(as.character(grid.grazing2016$grazing2016))
grid.grazing2016[grid.grazing2016$grazing2016>=1, "grazing20162016"] <- "Grazed"; grid.grazing2016[grid.grazing2016$grazing2016==0, "grazing2016"] <- "Ungrazed"
grid.fenced <- as.data.frame(cbind(over(grid, fenced.polygons[1]), coordinates(grid)[,1:2]), stringsAsFactors=FALSE); names(grid.fenced) <- c("area.fenced", "x", "y")##Extract Summer, Winter dataframe from the shape file
grid.fenced$area.fenced <- as.character(grid.fenced$area.fenced); grid.fenced[is.na(grid.fenced)==TRUE, "area.fenced"] <- "Unfenced"
grid.oversow <- as.data.frame(cbind(over(grid, oversow.polygons[,1]), coordinates(grid)[,1:2])); names(grid.oversow) <- c("oversow", "x", "y")
grid.oversow$oversow <-  as.character(grid.oversow$oversow); grid.oversow[is.na(grid.oversow$oversow)==FALSE, "oversow"] <- "yes"; grid.oversow[is.na(grid.oversow$oversow)==TRUE, "oversow"] <- "no"


## ##Correct NZMG to NZTM - Some plots had gps set to NZMG in Jan 2016
## NZMGlocs <- read.csv("/home/sean/Documents/Molesworth2016/NZMGlocs.csv", header=TRUE, skip = 0, stringsAsFactors=FALSE)
## coordinates(NZMGlocs) <- c("east", "north"); proj4string(NZMGlocs) <- CRS("+init=epsg:27200") #georeference to NZMG ##CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #georeference to for GE
## NZMGlocs <- spTransform(NZMGlocs, CRS("+init=epsg:2193"))
## write.table(NZMGlocs, file = "/home/sean/Documents/Molesworth2016/NZTMlocs.csv", sep = ",", na = "0", dec = ".", row.names = FALSE, col.names = TRUE)

locs.2016 <- as.data.frame(read.csv("~/Documents/Molesworth/Data/Spatial/SiteData2016.csv", sep = ",", header=TRUE, skip = 5, colClasses="character")[, c("Line", "Plot", "East", "North")], stringsAsFactors=FALSE)
locs.2016 <- locs.2016[is.na(locs.2016$East)==FALSE & locs.2016$East!="",]
locs.2016$plot <- paste(locs.2016$Line, locs.2016$Plot, sep="-"); 
locs.2016$plot <- sub("Fe", "-Fe", locs.2016$plot); locs.2016$plot <- sub("Un", "-Un", locs.2016$plot)
no.recce.plots <- locs.2016[c(grep("Saxton",locs.2016$plot), grep("Molesworth",locs.2016$plot)), "plot"]
locs.2016 <- locs.2016[-which(locs.2016$plot %in% no.recce.plots), c("plot", "East", "North")]; names(locs.2016) <- c("plot", "x", "y")
locs.1987 <- as.data.frame(read.csv("~/Documents/Molesworth/Data/Molesworth1987ReccesNVS/Coordinates.csv", sep = ",", header=TRUE, skip = 0, colClasses="character")[, c("Plot","EastingMG","NorthingMG")], stringsAsFactors=FALSE)
names(locs.1987) <- c("plot", "x", "y"); locs.1987$x <- as.numeric(locs.1987$x); locs.1987$y <- as.numeric(locs.1987$y); locs.1987 <- locs.1987[is.na(locs.1987$x)==FALSE, ]
coordinates(locs.1987) <- c("x", "y"); proj4string(locs.1987) <- CRS("+init=epsg:27200") #georeference to NZMG ##CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #georeference to for GE
locs.1987 <- spTransform(locs.1987, CRS("+init=epsg:2193"))  #Convert to NZTM; locs.1987.spatial.points <- locs.1987
locs.1987 <-  as.data.frame(locs.1987)
moore.locs <- data.frame(read.csv("/home/sean/Documents/Molesworth/Data/Spatial/MolesworthReccce1952.csv", sep = ",", header=TRUE, skip = 5), stringsAsFactors=FALSE)
moore.locs <- moore.locs[is.na(moore.locs$x)==FALSE,]; moore.locs$plot <- paste("Moore-", moore.locs$plot, sep="")

locs <- rbind(locs.2016, as.data.frame(recce.locs.Tier1), locs.1989, locs.1987, wairau.locs, moore.locs[,c("plot","x","y")]) ### binds random, paired, doc tier1, 1987 PNA and moore plots. 

locs[,c("x","y")] <- lapply(c("x","y"), function(x) as.integer(locs[,x])); locs <- locs[duplicated(locs$plot)==FALSE,  ] ##remove duplicates from tier 1 data
coordinates(locs) <- c("x", "y") #save cordinates in prep for spatial conversion - need to convert NZMG for 2009 into NZTM
proj4string(locs) <- CRS("+init=epsg:2193") 
boundary <- readOGR("/home/sean/Documents/Molesworth/Data/Spatial/MolesworthBoundary/Molesworth.shp", "Molesworth"); boundary <- spTransform(boundary, CRS("+init=epsg:2193")) #georeference to NZTM from NZMG
locs <- locs[boundary,]  ##subset locs only found within Molesworth boundary
##molesworth.1987.locs <- locs.1987[boundary,]

locs.terrace <- as.data.frame(cbind(locs$plot,extract(raster(terrace), coordinates(locs)))); names(locs.terrace) <- c("plot", "terrace")
locs.terrace[is.na(locs.terrace$terrace)==TRUE, "terrace"] <- 0 ##Calls a couple of mssing values 0 for slope. Terrace is 1.
fenced.area.ha <- (sum(sapply(slot(fenced.polygons, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area"))))/10000 #How much area all polygons in ha
oversow.area.ha <- (sum(sapply(slot(oversow.polygons, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area"))))/10000 #
boundary.area.ha <- (sum(sapply(slot(boundary, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area"))))/10000 #How much area all polygons in ha

locs <-  (as.data.frame(cbind(locs$plot, coordinates(locs), round(extract(raster(dem), locs), digits=1), (extract(raster(temperature), locs))/10, over(locs, lcdb.polygons[,1]),  over(locs, lenz.soil.polygons[,1]), over(locs, soil.polygons[,1]), over(locs, rain.polygons[,1]), over(locs, oversow.polygons[,1]), over(locs, grazing1987.polygons[3]), over(locs, grazing2016.polygons), over(locs, fenced.polygons[1])), stringsAsFactors=FALSE)) # 
names(locs) <- c("plot", "x", "y", "dem.altitude", "temp", "lcdb", "lenz.soil", "soil.class","rain.class", "oversow", "grazing1987", "grazing2016", "fenced")
locs$fenced <- as.character(locs$fenced); locs[grep("-Fe", locs$plot), "fenced"] <- "fenced"; locs[is.na(locs$fenced)==FALSE, "fenced"] <- "Fenced"; locs[is.na(locs$fenced)==TRUE, "fenced"] <- "Unfenced"
locs <- merge(locs, locs.terrace, by="plot", all.x=TRUE, all.y=TRUE)

## locs.1987[, c("x", "y", "dem.altitude", "temp")] <- sapply(locs.1987[, c("x", "y", "dem.altitude", "temp")], as.numeric)
## locs.1987[, c("lcdb", "lenz.soil", "soil.class","rain.class", "oversow", "grazing1987", "grazing2016", "fenced")] <- sapply(locs.1987[, c("lcdb", "lenz.soil", "soil.class","rain.class", "oversow", "grazing1987", "grazing2016", "fenced")], as.character)
## locs.1987[is.na(locs.1987$oversow)==FALSE, "oversow" ] <- "Oversown"; locs.1987[is.na(locs.1987$oversow)==TRUE, "oversow" ] <- "No"
## locs.1987[is.na(locs.1987$grazing2016)==FALSE, "grazing2016" ] <- "Grazed"; locs.1987[is.na(locs.1987$grazing2016)==TRUE, "grazing2016" ] <- "Ungrazed"
## locs.1987[is.na(locs.1987$grazing1987)==TRUE, "grazing1987" ] <- "Ungrazed"
## locs.1987[is.na(locs.1987$fenced)==FALSE, "fenced" ] <- "Fenced"; locs.1987[is.na(locs.1987$fenced)==TRUE, "fenced" ] <- "Unfenced"
## dca.site.2 <- merge(dca.site, locs.1987, by="plot"); dca.dat <- merge(dca.site.2, env.1987[,c("plot","terrace", "water.deficit")], by="plot", all.X=TRUE, all.Y=FALSE)
## dca.dat[dca.dat$terrace=="", "terrace"] <- "Slope";  dca.dat$water.deficit <- as.numeric(dca.dat$water.deficit)
## ##dca.dat$x <- dca.dat$x + runif(1, min = 1, max = 50); dca.dat$y <- dca.dat$y + runif(1, min = 1, max = 50) ##Add imprecision to grid reference data
## ##Extract envirnomental data for recce plots.
## env.1987 <- as.data.frame(read.csv("~/Documents/Molesworth2016/RecceTerraceWater.csv", sep = ",", header=TRUE, skip = 0, colClasses="character"), stringsAsFactors=FALSE)


## #############################################################                Ordination of all species                  #############################################

recce.dat$plot.yr <- paste(recce.dat$plot, recce.dat$year, sep=".")
recce.dat$species <- as.character(recce.dat$species)
recce.dat <- recce.dat[recce.dat$plot %in% locs$plot, ] ## only plots within the boundary of Molesworth
#recce.dat[recce.dat$cover>0,] <- 1 ## Make presence absence

nplot <- length(table(recce.dat$plot.yr)) # what number of plots
nspp  <- length(table(recce.dat$species))  # what number of species

comp <- matrix(0, nrow=nplot, ncol=nspp) # produce a matix of length plots and width species
rownames(comp) <- unique(recce.dat$plot.yr) # give row names plots 
colnames(comp) <- unique(recce.dat$species) # give column names species

plot <-recce.dat$plot.yr  # matrix py takes the form of py from data frame dat
species <- recce.dat$species
cover <- recce.dat$cover
        
for(i in 1:nrow(recce.dat)) {
  r <- which(rownames(comp)==plot[i])
  c <- which(colnames(comp)==species[i])
  comp[r, c] <- cover[i]
}

##.94 + 1.1 + 3.7 + 1.8

 #########################################################  Ordination for Most common species      ######################################################################

comp.common <-as.data.frame(comp); comp.common$year <-  as.numeric(substr(rownames(comp.common), regexpr("\\.",rownames(comp.common))+1, regexpr("\\.", rownames(comp.common))+5))
comp.common <- comp.common[comp.common$year==2016, c(1:length(names(comp.common))-1)]
comp.common$plot <-  as.character(substr(rownames(comp.common), 1, regexpr("\\.",rownames(comp.common))-1))
random.plots <- unique(recce.2007$plot); random.plots <- random.plots[-which(random.plots %in% c("B-2A", "P-2A"))]
comp.common <- comp.common[comp.common$plot %in% random.plots, c(1:length(names(comp.common))-1)]
for(i in names(comp.common)) {comp.common[comp.common[names(comp.common)==i] != 0, names(comp.common)==i] <- 1} ## produce a matrix with 1 for occurance
col.sums.comp.common <- as.data.frame(colSums(comp.common)); col.sums.comp.common$spp <- rownames(col.sums.comp.common)

common.spp <- sort(col.sums.comp.common[(col.sums.comp.common[1] > 32)==TRUE, "spp"])  ##Produce a list of common species. 20 for species occuring in more than 10 plots in 2016 random plot survey

recce.common <- recce.dat[recce.dat$species %in% common.spp, ]

nplot <- length(table(recce.common$plot.yr)) # what number of plots
nspp  <- length(table(recce.common$species))  # what number of species

comp <- matrix(0, nrow=nplot, ncol=nspp) # produce a matix of length plots and width species
rownames(comp) <- unique(recce.common$plot.yr) # give row names plots 
colnames(comp) <- unique(recce.common$species) # give column names species

plot <-recce.common$plot.yr  # matrix py takes the form of py from data frame dat
species <- recce.common$species
cover <- recce.common$cover
        
for(i in 1:nrow(recce.common)) {
  r <- which(rownames(comp)==plot[i])
  c <- which(colnames(comp)==species[i])
  comp[r, c] <- cover[i]
}

comp.common.spp <- comp



#####################################################################################################################################################################################
     ##########################################################        nMDS  Analaysis                                   ##############################################################


## mds.comp  <- metaMDS(comp, k=2,  distance = "euclidean", trymax = 100,  autotransform = FALSE) ###Reduces dimensions to two with a Eculidean distance matrix. Less stress than Bray
## ##Non default number of iterations (trymax) required for repeated success. 

## nmds.comp.run <- as.data.frame(mds.comp$points[,c(1,2)]);  nmds.comp.run$plot  <- substr(rownames(nmds.comp.run), 1, nchar(rownames(nmds.comp.run))-5)
## ##Data frame with nMDS of common species
## nmds.comm  <- merge(nmds.comp.run, locs[ , c("plot", "x", "y", "dem.altitude", "oversow", "fenced", "terrace")])

## fit.nmds <- envfit(mds.comp ~  nmds.comm$x * nmds.comm$y * nmds.comm$dem.altitude * nmds.comm$oversow * nmds.comm$fenced * nmds.comm$terrace, na.rm = TRUE) ###Two altutudes missing

## nmds.tab <-  round(as.data.frame(cbind(fit.nmds$factors$r, fit.nmds$factors$pvals)), digits=3)
## nmds.tab$env.var <- c("Temperature", "O2 COD ratio", "Total solids", "Time")
## names(nmds.voc.tab)  <-  c("R2", "P", "Factor")

## nmds.voc.arrows <-  as.data.frame(fit.nmds.voc$factors$centroids)
## nmds.voc.arrows$var <- substr(rownames(nmds.voc.arrows), regexpr("dat", rownames(nmds.voc.arrows))+4, regexpr(")", rownames(nmds.voc.arrows))-1)
## nmds.voc.arrows$var.level <- substr(rownames(nmds.voc.arrows), regexpr(")", rownames(nmds.voc.arrows))+1, nchar(rownames(nmds.voc.arrows)))

## nmds.voc <- as.data.frame(mds.voc$species[,c(1,2)]) ; nmds.voc$voc <- rownames(nmds.voc)
## names(nmds.voc) <- c("nmds1", "nmds2", "name"); nmds.voc <- nmds.voc[ ,c("name", "nmds1", "nmds2")]
## nmds.voc  <- merge(nmds.voc, name.list[, c("name", "Jan.Grouping", "Colour", "GertyGrouping", "NewColour")], by="name", all.x=TRUE, all.y=FALSE)

## nmds.voc$index.num  <- order(nmds.voc$Jan.Grouping, nmds.voc$name)
## index.string  <- paste(", ", nmds.voc$index.num, ": ", nmds.voc$name, sep="")
## index.string  <- paste("Figure 1. Non-metric multidimensional scaling diagram of mg TS of 69 VOCs. A Bray distance matrix was used to arrange data prior to ordination. Stress of ", round(mds.voc$stress, digits=3), "indicates a reliable represenation of all VOCs in two dimensions (Clarke, 1993). The envfit function in the R library vegan was used to fit environmental factors into the ordination. Temperature (red) and O2:COD ratio (blue) factors were related to the composition of VOCs. Time and Total Solids were not (grey). Each compound is coloured by categories of Hydrocarbon #000000,  Sulphur compounds #009E73, Aldehydes, ketones, alcohols, nitrile #56B4E9, Furans and pyrroles #D55E00 and Aromatics #E69F00. Each compound is numbered ", paste(index.string[order(as.numeric(substr(index.string, 3, regexpr(":", index.string)-1)))], collapse=""))
## write.xlsx(index.string, "/home/sean/Documents/VOC/TabulatedResults.xlsx", sheetName = "Codes",  col.names = FALSE, row.names = FALSE, append = TRUE)


#####################################################################################################################################################################################
##########################################################        DCA Analaysis                                        ##############################################################

ord <- decorana(comp)

dca.site <- as.data.frame(ord$rproj[, c(1:2)])
dca.site$year <-  as.numeric(substr(rownames(dca.site), regexpr("\\.",rownames(dca.site))+1, regexpr("\\.", rownames(dca.site))+5))
dca.site$plot <- substr(rownames(dca.site), 1, regexpr("\\.", rownames(dca.site))-1)
dca.spp <-  as.data.frame(ord$cproj[, c(1:2)])
dca.spp$species <- rownames(dca.spp)
dat <- merge(dca.site, locs, by="plot", all.x=TRUE, all.y=TRUE)
dat <- dat[is.na(dat$DCA1)==FALSE, ]
dat <- dat[is.na(dat$DCA1)==FALSE & is.na(dat$dem.altitude)==FALSE, c("plot", "DCA1", "DCA2", "year", "fenced", "x", "y", "dem.altitude", "terrace", "oversow", "grazing1987")]  ##

dat$oversow <-  as.character(dat$oversow) ; dat[is.na(dat$oversow)==FALSE, "oversow"] <- "oversown"; dat[is.na(dat$oversow)==TRUE, "oversow"] <- "natural"
dat$terrace <-  as.character(dat$terrace) ; dat[dat$terrace=="1", "terrace"] <- "Terrace"; dat[dat$terrace=="0", "terrace"] <- "Slope"
dat[dat$terrace == "Slope", "terrace"]  <- 0; dat[dat$terrace == "Terrace", "terrace"]  <- 1; dat$terrace  <- as.integer(as.character(dat$terrace))
dat[dat$fenced == "Unfenced", "fenced"]  <- 0; dat[dat$fenced == "Fenced", "fenced"]  <- 1; dat$fenced  <- as.integer(as.character(dat$fenced))
dat[dat$oversow == "natural", "oversow"]  <- 0; dat[dat$oversow == "oversown", "oversow"]  <- 1; dat$oversow  <- as.integer(as.character(dat$oversow))

dat[dat$year==2007 |  dat$year==2008, "x"] <- dat[dat$year==2007 | dat$year==2008, "x"] + 1.5
dat[dat$year==2007 | dat$year==2008, "y"] <- dat[dat$year==2007 | dat$year==2008, "y"] + 1.5
dat[duplicated(dat$x)==TRUE, "x"] <-  dat[duplicated(dat$x)==TRUE, "x"] + 1.5  ##Make changes for things duplicated once
dat[duplicated(dat$y)==TRUE, "y"] <-  dat[duplicated(dat$y)==TRUE, "y"] + 1.5  
dat[duplicated(dat$x)==TRUE, "x"] <-  dat[duplicated(dat$x)==TRUE, "x"] - 5  ##Make changes for things duplicated twice
dat[duplicated(dat$y)==TRUE, "y"] <-  dat[duplicated(dat$y)==TRUE, "y"] - 5  
dat[duplicated(dat$x)==TRUE, "x"] <-  dat[duplicated(dat$x)==TRUE, "x"] - 2.5  ##Make changes for things duplicated twice
dat[duplicated(dat$y)==TRUE, "y"] <-  dat[duplicated(dat$y)==TRUE, "y"] - 2.5  
dat.locs <- dat[, c("plot", "year", "x", "y")]; coordinates(dat.locs) <- c("x", "y"); proj4string(dat.locs) <- CRS("+init=epsg:2193") #georeference to NZTM 

dca.site.tab  <-  dca.site
dca.site.tab$plot.yr  <- rownames(dca.site.tab)

comp.site  <- as.data.frame(comp)
comp.site$plot.yr  <- rownames(comp.site)
dca.site.tab <- merge(dca.site.tab, comp.site, by="plot.yr", all.x=TRUE, all.y=TRUE)

common.table.mean  <- data.frame(species = common.spp, mean.1952= NA, mean.1960= NA, mean.1972= NA, mean.1987= NA, mean.1989= NA, mean.1995= NA, mean.2001= NA, mean.2006=NA,  mean.2007= NA, mean.2008= NA, mean.2012= NA, mean.2016= NA) 
common.table.sem  <- data.frame(species = common.spp, sem.1952= NA, sem.1960= NA, sem.1972= NA, sem.1987= NA, sem.1989= NA, sem.1995= NA, sem.2001= NA, sem.2006=NA,  sem.2007= NA, sem.2008= NA, sem.2012= NA, sem.2016= NA) 


for (i in common.spp) {
    common.table.mean[common.table.mean$species == i , 2:13]  <-  as.character(round(tapply(dca.site.tab[ , i], dca.site.tab[ , "year"], mean),digits=2))
    common.table.sem[common.table.sem$species == i , 2:13]  <-  as.character(round(tapply(dca.site.tab[ , i], dca.site.tab[ , "year"], s.e.),digits=2)) 
    }

common.table  <- cbind(common.table.mean, common.table.sem[, -1])
common.table  <-  common.table[ , c("species", "mean.1952", "sem.1952",  "mean.1960", "sem.1960", "mean.1972", "sem.1972", "mean.1987", "sem.1987",
                                "mean.1989", "sem.1989", "mean.1995", "sem.1995", "mean.2001", "sem.2001", "mean.2006", "sem.2006",
                               "mean.2007", "sem.2007", "mean.2008", "sem.2008", "mean.2012", "sem.2012", "mean.2016", "sem.2016") ]

common.table <-  merge(common.table, spp.list[, c("species", "family", "native", "form")], by="species", all.x=TRUE, all.y=FALSE)
common.table <- common.table[duplicated(common.table$species)==FALSE, ]
common.table[common.table$native!="Exotic","native"] <- "Native"; common.table[common.table$form %in%  c("Shrub", "Tree", "Vine","SubShrub", "Treefern"), "form"] <- "Woody"
###Add native and form classes to the table

common.table.recce <- rbind( c("Native herbaceous plants", rep("", length(names(common.table))-3)),                      
                      common.table[common.table$native=="Native" & common.table$form=="Forb", names(common.table)[1:26]],
                     c("Exotic herbaceous plants", rep("", length(names(common.table))-3)),                      
                      common.table[common.table$native=="Exotic" & common.table$form=="Forb", names(common.table)[1:26]],
                     c("Native grasses, tussocks, sedges and rushes", rep("", length(names(common.table))-3)),                      
                      common.table[common.table$native=="Native" & common.table$form=="Graminoid", names(common.table)[1:26]],
                     c("Exotic grasses, sedges and rushes", rep("", length(names(common.table))-3)),                      
                      common.table[common.table$native=="Exotic" & common.table$form=="Graminoid", names(common.table)[1:26]],
                     c("Native woody plants", rep("", length(names(common.table))-3)),                      
                      common.table[common.table$native=="Native" & common.table$form=="Woody", names(common.table)[1:26]],
                     c("Exotic woody plants", rep("", length(names(common.table))-3)),                      
                      common.table[common.table$native=="Exotic" & common.table$form=="Woody", names(common.table)[1:26]])

   ##############     ##############     ##############     ##############     ##############     ##############     ##############  
##############For statements and single results table in main manuscript
#####Produce a data frame with statistics quoted in text in the results section. Included in a list, and saved later. 
 dca.site.tab$plot.yr  <-  rownames(dca.site.tab)

### Four introduced plants associated with grazing in valley bottoms had low DCA axis 1 scores...
    four.exotics  <- dca.spp[dca.spp$species %in% c("Agrostis capillaris", "Anthoxanthum odoratum", "Rosa rubiginosa", "Trifolium repens"), ]
    four.exotics.range  <- c(paste0("range = ", min(four.exotics$DCA1), "--", round(max(four.exotics$DCA1), digits=3)))

#### Hypochaeris radicata and Hieraciinae exotics had high DCA axis 1 an 2 scores, related to their increase in cover between 2007 and 2016.......
    five.daisies  <- dca.spp[dca.spp$species %in% c("Hypochaeris radicata", "Hieracium pollichiae", "Pilosella caespitosa", "Pilosella officinarum", "Pilosella praealta"), ]
    five.daisies.range  <- c(paste0("DCA 1 range = ", round(min(five.daisies$DCA1), digits=2), "--", round(max(five.daisies$DCA1), digits=3)), c(paste0("DCA 2 range = ", round(min(five.daisies$DCA2), digits=2), "--", round(max(five.daisies$DCA2), digits=3))))


####Some grass and tussock species increased in cover and biomass, and had accordingly high DCA scores....
    five.grasses  <- dca.spp[dca.spp$species %in% c("Anthosachne solandri", "Chionochloa flavescens", "Festuca novae-zelandiae", "Koeleria novozelandica", "Poa colensoi"), ]
    five.grasses.DCArange  <- c(paste0("DCA 1 range = ", round(min(five.grasses$DCA1), digits=2), "--", round(max(five.grasses$DCA1), digits=3)), c(paste0("DCA 2 range = ", round(min(five.grasses$DCA2), digits=2), "--", round(max(five.grasses$DCA2), digits=3))))

### Other widespread grass species showed less  indication of change, but still had high DCA axis 2 scores....
    three.grasses  <- dca.spp[dca.spp$species %in% c("Anthoxanthum odoratum", "Deyeuxia avenoides", "Rytidosperma setifolium"), ]
    three.grasses.DCArange  <- c(paste0("DCA 1 range = ", round(min(three.grasses$DCA1), digits=2), "--", round(max(three.grasses$DCA1), digits=3)), c(paste0("DCA 2 range = ", round(min(three.grasses$DCA2), digits=2), "--", round(max(three.grasses$DCA2), digits=3))))

### Some  woody species increased in cover between 2007 and 2016 by up to 40%,  had high DCA axis 2 scores, and included low-statured native species....
    four.shrubs  <- dca.spp[dca.spp$species %in% c("Acrothamnus colensoi", "Ozothamnus vauvilliersii", "Pimelea oreophila", "Rosa rubiginosa"), ]
    four.shrubs.DCArange  <- c(paste0("DCA 1 range = ", round(min(four.shrubs$DCA1), digits=2), "--", round(max(four.shrubs$DCA1), digits=3)), c(paste0("DCA 2 range = ", round(min(four.shrubs$DCA2), digits=2), "--", round(max(four.shrubs$DCA2), digits=3))))

noted.spp  <-  c("Acrothamnus colensoi", "Anthoxanthum odoratum", "Anthosachne solandri", "Chionochloa flavescens", "Deyeuxia avenoides", "Festuca novae-zelandiae", "Hypochaeris radicata", "Hieracium pollichiae", "Koeleria novozelandica", "Ozothamnus vauvilliersii",
                 "Pilosella caespitosa", "Pilosella officinarum", "Pilosella praealta", "Pimelea oreophila", "Poa colensoi", "Rosa rubiginosa", "Rytidosperma setifolium")

noted.spp.dat  <-  data.frame(species=noted.spp, mean.cover.1952=NA, mean.cover.1987=NA, mean.cover.07.08 = NA, mean.cover.2016=NA, change.cover.52.16 = NA, change.cover.87.16 = NA, change.cover.0708.16 = NA)

## for (i in noted.spp) { ##
## i  <-  "Hypochaeris radicata"
## temp.dat  <- data.frame(plot.yr=rownames(comp.site), spp.i = comp.site[, i]); names(temp.dat)   <-  c("plot.yr", "spp.i") 
## temp.dat  <- merge(dca.site, temp.dat, by = "plot.yr")
## mean.1952  <- round(mean(temp.dat[temp.dat$year==1952, "spp.i"], na.rm=TRUE), digits=1); mean.1987  <- round(mean(temp.dat[temp.dat$year==1987, "spp.i"], na.rm=TRUE), digits=1)
## mean.07.08  <- round(mean(temp.dat[temp.dat$year %in% c(2007,2008), "spp.i"], na.rm=TRUE), digits=1); mean.2016  <- round(mean(temp.dat[temp.dat$year==2016, "spp.i"], na.rm=TRUE), digits=1)
## mean.1987  <- round(mean(temp.dat[temp.dat$year==1987, "spp.i"], na.rm=TRUE), digits=1); mean.2016  <- round(mean(temp.dat[temp.dat$year==2016, "spp.i"], na.rm=TRUE), digits=1)
## temp.mean.change.52  <- round(((mean(temp.dat[temp.dat$year==2016, "spp.i"], na.rm=TRUE) - mean(temp.dat[temp.dat$year==1952, "spp.i"], na.rm=TRUE)) / mean(temp.dat[temp.dat$year==1952, "spp.i"], na.rm=TRUE)) * 100, digits=1)
## temp.mean.change.87  <- round(((mean(temp.dat[temp.dat$year==2016, "spp.i"], na.rm=TRUE) - mean(temp.dat[temp.dat$year==1987, "spp.i"], na.rm=TRUE)) / mean(temp.dat[temp.dat$year==1987, "spp.i"], na.rm=TRUE)) * 100, digits=1)
## temp.mean.change.0708  <- round(((mean(temp.dat[temp.dat$year==2016, "spp.i"], na.rm=TRUE) - mean(temp.dat[temp.dat$year %in% c(2007,2008), "spp.i"], na.rm=TRUE)) / mean(temp.dat[temp.dat$year==%in% c(2007,2008), "spp.i"], na.rm=TRUE)) * 100, digits=1)
## temp.paired  <- merge(temp.dat[temp.dat$year %in% c(2007,2008),  ], temp.dat[temp.dat$year==2016, c("plot", "spp.i")], by="plot")
## names(temp.paired)  <- c("plot", "plot.yr", "DCA1", "DCA2", "establish.year", "cover.0708", "cover.16")
## temp.paired$change  <- (temp.paired$cover.16 - temp.paired$cover.0708)/ temp.paired$cover.0708

## ##; temp.paired$perc.change  <- temp.paired$change  
## ##t.test(temp.dat[temp.dat$year==2016, "spp.i"],  temp.dat[temp.dat$year==1987, "spp.i"],  alternative = "two.sided", paired = FALSE, var.equal = FALSE)
## noted.spp.dat[noted.spp.dat$species==i,  c("mean.cover.1952", "mean.cover.1987", "mean.cover.07.08", "mean.cover.2016", "change.cover.87.16")]  <- c(mean.1952, mean.1987, mean.07.08, mean.2016, temp.mean.change)
## }


####A list of the 32 common species in 80 random plots in 2016 (occured in m more than 32)
common.spp.list  <-  c("Aciphylla aurea","Acrothamnus colensoi", "Agrostis capillaris","Anisotome aromatica","Anthosachne solandri","Anthoxanthum odoratum","Celmisia gracilenta","Celmisia spectabilis","Chionochloa flavescens",
                        "Deyeuxia avenoides","Dracophyllum rosmarinifolium","Festuca novae-zelandiae","Gaultheria depressa","Hieracium pollichiae","Hypochaeris radicata","Koeleria novozelandica","Leucopogon fraseri","Luzula rufa",
                        "Lycopodium fastigiatum","Olearia cymbifolia","Ozothamnus vauvilliersii","Pilosella caespitosa","Pilosella officinarum","Pilosella praealta","Pimelea oreophila","Poa colensoi","Raoulia subsericea",
                        "Rosa rubiginosa","Rumex acetosella","Rytidosperma setifolium","Trifolium repens","Wahlenbergia albomarginata")
### Import HFI data for site variables and tabulation
hfi.common <- read.csv("/home/sean/Documents/Molesworth/Data/HFISpeciesComposition.csv", sep = ",", header=TRUE, skip = 0)
############### Produce a table with change in %cover and %HFI
change.table  <- data.frame(species = common.spp.list, Cover=NA, c.se=NA, ranef.cov=NA, t.146.cover=NA, P.cover=NA,  HFI = NA, h.se=NA, ranef.hfi=NA, t.146.hfi=NA, P.hfi=NA)  
change.table <-  merge(change.table, spp.list[duplicated(spp.list$species)==FALSE, c("species", "family", "native", "form")], by="species", all.x=TRUE, all.y=FALSE)
##change.table$HFI  <- as.integer(change.table$HFI)


####First cover... 
hfi.site.dat  <- hfi.common[, c("plot.yr", "plot", "year")] ##Bind with hfi common because of ease of access to environmental variables ,:
hfi.site.dat$plot.yr  <- paste0(substr(hfi.site.dat$plot.yr, 1, nchar(hfi.site.dat$plot.yr)-5), ".",
                  substr(hfi.site.dat$plot.yr, nchar(hfi.site.dat$plot.yr)-3, nchar(hfi.site.dat$plot.yr)))
## "x", "y", "soil.class", "rain.class",  "oversow.year", "terrace", "water.deficit", "temperature", "lcdb", "lenz.soil"
## "dem.altitude", "dem.slope", "grazing.areas", "fenced.areas", "summer.winter.grazing.areas""fenced.year"

###Convert cover spp x plot matrix to long format
cover.dat  <- merge(melt(data = dca.site.tab[, 4:37], id.vars = c("plot", "year"), variable.name = "species",  value.name = "cover"),
                    dat[, c("plot", "year", "x", "y", "dem.altitude", "fenced", "terrace", "oversow")], by = c("plot", "year"), all.x=TRUE, all.y=FALSE) 
## ###PLots measured repeatedly, so random effect. Specie drawn
## ##m.1.gauu <- glmmTMB(cover ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 + (1|plot) + (1|species), data=cover.dat[cover.dat$year %in% c(2007,2008, 2016), ], family="gaussian")
## ##m.1 <- lmer(cover ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 + (1|plot) + (year|species), data=cover.dat[cover.dat$year %in% c(2007,2008, 2016), ]) ###Does not converge

####Save instead of run to save time
 ## m.1 <- glmmTMB(cover ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 + (1|plot) + (year|species), data=cover.dat[cover.dat$year %in% c(2007,2008, 2016), ], family=genpois) 
 ## m.null <- glmmTMB(cover ~ 1 +  (1|plot) + (year|species), data=cover.dat[cover.dat$year %in% c(2007,2008, 2016), ], family="genpois")
 ## save(m.1, m.null, file="/home/sean/Documents/Molesworth/Data/Objects/glmmTMBCover.Rdata")
load("/home/sean/Documents/Molesworth/Data/Objects/glmmTMBCover.Rdata")
###POisson model fits better, and converges with the inclusion of specie varying by year.
###all of the effects of X1 (year - i.e. all of the predictors that we would get from a linear model using y ~ X1) vary across the groups/levels defined by X2 (species). Note x1 numeric, x2 categorical

ranef.spp  <- as.data.frame(ranef(m.1)$cond$species)*100*(2016-2008)
recce.spp.R2  <- MuMIn::r.squaredLR(m.1, m.null)[1] ##recce.spp.cover  <- MuMIn::r.squaredGLMM(m.1, m.null). #### r.squaredGLMM
glmm.model.spp.recce  <-  glmmTMB.latex(m.1)

###### Zero inflation term did not improve model
###m.1.zero <- glmmTMB(cover ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 + (1|plot) + (year|species), data=cover.dat[cover.dat$year %in% c(2007,2008, 2016), ], ziformula=~1, family="genpois") 

for (i in common.spp.list) { ##    i  <-  "Chionochloa flavescens"
## temp.dat  <- data.frame(plot.yr=rownames(comp.site), spp.i = comp.site[, i]); names(temp.dat)   <-  c("plot.yr", "spp.i") 
## temp.dat  <- merge(hfi.site.dat[, c("plot.yr", "plot", "year")], temp.dat, by = "plot.yr")
## temp.dat  <- merge(temp.dat[temp.dat$year %in%  c(2007, 2008), c("plot", "spp.i")], temp.dat[temp.dat$year == 2016, c("plot", "spp.i")], by = "plot", all.x=TRUE, all.y=FALSE)
temp.dat  <- merge(dca.site.tab[dca.site.tab$year %in%  c(2007, 2008), c("plot", i)], dca.site.tab[dca.site.tab$year == 2016, c("plot", i)], by = "plot", all.x=TRUE, all.y=FALSE)    
names(temp.dat)  <- c("plot", "i.0708", "i.2016");  temp.dat$change  <-  ((temp.dat$i.2016  -  temp.dat$i.0708)/ (temp.dat$i.0708)) *100  
      temp.dat[temp.dat$i.0708 == 0 & temp.dat$i.2016 == 1, "change"] <-  100
      temp.dat[temp.dat$i.0708 == 0 & temp.dat$i.2016 >= 2, "change"] <-  200
      temp.dat[temp.dat$i.0708 == 1 & temp.dat$i.2016 == 0, "change"] <-  -100
      temp.dat[temp.dat$i.0708 >= 2 & temp.dat$i.2016 == 0, "change"] <-  -200
      temp.dat[temp.dat$i.0708 == 0 & temp.dat$i.2016 == 0, "change"] <-  0    
##  temp.dat[is.finite(temp.dat$change) == FALSE, "change"]  <-  0
change.table[change.table$species == i, "Cover"]  <-   round(mean(temp.dat$change, na.rm=TRUE), digits=0)
change.table[change.table$species == i, "c.se"]  <-   round(s.e.(temp.dat$change), digits=0)
temp.ttest  <- t.test(temp.dat$i.0708, temp.dat$i.2016,  alternative = "two.sided", paired = FALSE, var.equal = FALSE)
    change.table[change.table$species == i, "ranef.cov"]  <- round(ranef.spp[rownames(ranef.spp) == i, "year"], digits=1)
    change.table[change.table$species == i, "t.146.cover"]  <- round(abs(temp.ttest$statistic), digits=3)
    change.table[change.table$species == i, "P.cover"]  <- round(abs(temp.ttest$p.value), digits=3)
    }

###Second HFI . Absolute mean changed in number of intercepts for each species in each plot established in 2007 and 2008. i.e random and paired fenced. 

###Fix up some name changes that occured between surveys.
hfi.common[hfi.common$ELYSOL != 0, "ANTSOL"]  <- hfi.common[hfi.common$ELYSOL != 0, "ELYSOL"] ##Was Elymus solandri now Anthosachne solandri
hfi.common[hfi.common$ELYAFR != 0, "ANTSOL"]  <- hfi.common[hfi.common$ELYAFR != 0, "ELYAFR"] ##Was Elymus species (aff. rectisetus) not on NZPCN
hfi.common[hfi.common$BROHOR != 0, "KOENOV"]  <- hfi.common[hfi.common$BROHOR != 0, "BROHOR"] ##Bromus hordeaceus mistaken for Koeleria novozelandica? Deyeuxia avenoides?
hfi.common[hfi.common$CYNCRI != 0, "KOENOV"]  <- hfi.common[hfi.common$CYNCRI != 0, "CYNCRI"] ##Bromus hordeaceus mistaken for Koeleria novozelandica? Deyeuxia avenoides?
hfi.common[hfi.common$CHIFSB != 0, "CHIFLA"]  <- hfi.common[hfi.common$CHIFSB != 0, "CHIFSB"] ##Change code used for Chionochloa flavescens
hfi.common[hfi.common$HIECAE != 0, "PILCAE"]  <- hfi.common[hfi.common$HIECAE != 0, "HIECAE"] ##Hieracium to Pilosella. H. caespitosa.
hfi.common[hfi.common$HIEPIL != 0, "PILOFF"]  <- hfi.common[hfi.common$HIEPIL != 0, "HIEPIL"] ##Hieracium pilosella to Pilosella officinarum.
hfi.common[hfi.common$HIEPRA != 0, "PILPRA"]  <- hfi.common[hfi.common$HIEPRA != 0, "HIEPRA"] ##Hieracium pra to Pilosella praealta.

nvs.spp <- read.csv("/home/sean/Documents/Molesworth/Data/CurrentNVSNames.csv", sep = ",", header=TRUE, skip = 0, stringsAsFactors=FALSE)
nvs.grasses  <- nvs.spp[nvs.spp$NVSCode %in% names(hfi.common)[28:338] & nvs.spp$GrowthForm == "Graminoid", c("SpeciesName", "NVSCode")] ###For checking grass spp
nvs.common  <- nvs.spp[nvs.spp$SpeciesName %in% common.spp, c("SpeciesName", "NVSCode")]
nvs.common  <- rbind(nvs.common, c("Pilosella praealta", "PILPRA"))
nvs.common  <-  nvs.common[order(nvs.common$SpeciesName),]
hfi.common.spp  <- hfi.common[ , c("ACIAUR","ACRCOL","AGRCAP","ANIARO","ANTSOL","ANTODO","CELGRA","CELSPE","CHIFLA","DEYAVE","DRAROS","FESNOV","GAUDEP",
                   "HIEPOL","HYPRAD","KOENOV","LEUFRA","LUZRUF","LYCFAS","OLECYM","OZOVAU",
                   "PILCAE","PILOFF","PILPRA","PIMORE","POACOL","RAOSUB","ROSRUB","RUMACE","RYTSET","TRIREP","WAHALB")]
names(hfi.common.spp)  <- common.spp.list
hfi.common  <- cbind(hfi.common[, c("plot.yr", "plot", "year", "native.n", "n", "exotic.n", "fenced.year")], hfi.common.spp)                 

###Convert cover spp x plot matrix to long format
hfi.dat  <- merge(melt(data = hfi.common[, c(2:3, 8:39)], id.vars = c("plot", "year"), variable.name = "species",  value.name = "hfi"),
                    dat[, c("plot", "year", "x", "y", "dem.altitude", "fenced", "terrace", "oversow")], by = c("plot", "year"), all.x=TRUE, all.y=FALSE) 
####Save instead of run to save time
## m.1 <- glmmTMB(hfi ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 + (1|plot) + (year|species), data=hfi.dat[hfi.dat$year %in% c(2007,2008, 2016), ], family="gaussian")
## m.fixed.slope  <- glmmTMB(hfi ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 + (1|plot) + (1|species), data=hfi.dat[hfi.dat$year %in% c(2007,2008, 2016), ], family="gaussian")
## m.null <- glmmTMB(hfi ~ 1 +  (1|plot) + (1|species), data=hfi.dat[hfi.dat$year %in% c(2007,2008, 2016), ], family="gaussian")  ###Note. WOuld not fit with year|species. so a model for fixed slope used
## save(m.1, m.fixed.slope, m.null, file="/home/sean/Documents/Molesworth/Data/Objects/glmmTMBHFI.Rdata")
load("/home/sean/Documents/Molesworth/Data/Objects/glmmTMBHFI.Rdata")


###m.1.pois <- glmmTMB(hfi ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 + (1|plot) + (year|species), data=hfi.dat[hfi.dat$year %in% c(2007,2008, 2016), ], family="genpois") 
####POisson model fits better, but not significantly better .   anova(m.1.gaus, m.1)
ranef.spp  <- as.data.frame(ranef(m.1)$cond$species)*100*(2016-2008)
hfi.spp.R2  <- MuMIn::r.squaredLR(m.fixed.slope, m.null)[1] ##recce.spp.cover  <- MuMIn::r.squaredGLMM(m.1, m.null). #### r.squaredGLMM
glmm.model.hfi  <-  glmmTMB.latex(m.1)

####
for (i in common.spp.list) { ##i  <-  "Chionochloa flavescens"
temp.dat  <- merge(hfi.common[hfi.common$year %in%  c(2007, 2008), c("plot", i)], hfi.common[hfi.common$year == 2016, c("plot", i)], by = "plot", all.x=TRUE, all.y=FALSE)
names(temp.dat)  <- c("plot", "i.0708", "i.2016");  temp.dat$change  <-  ((temp.dat$i.2016  -  temp.dat$i.0708)/ (temp.dat$i.0708)) * 100 
###temp.dat[is.finite(temp.dat$change) == FALSE, "change"]  <-  0
     temp.dat[temp.dat$i.0708 == 0 & temp.dat$i.2016 == 1, "change"] <-  100
      temp.dat[temp.dat$i.0708 == 0 & temp.dat$i.2016 >= 2, "change"] <-  200
      temp.dat[temp.dat$i.0708 == 1 & temp.dat$i.2016 == 0, "change"] <-  -100
      temp.dat[temp.dat$i.0708 >= 2 & temp.dat$i.2016 == 0, "change"] <-  -200
      temp.dat[temp.dat$i.0708 == 0 & temp.dat$i.2016 == 0, "change"] <-  0    
change.table[change.table$species == i, "HFI"]  <-   round(mean(temp.dat$change, na.rm=TRUE), digits=0)
change.table[change.table$species == i, "h.se"]  <-   round(s.e.(temp.dat$change), digits=0)
temp.ttest  <- t.test(temp.dat$i.0708, temp.dat$i.2016,  alternative = "two.sided", paired = FALSE, var.equal = FALSE)
change.table[change.table$species == i, "ranef.hfi"]  <- round(ranef.spp[rownames(ranef.spp) == i, "year"], digits=1)
change.table[change.table$species == i, "t.146.hfi"]  <- round(abs(temp.ttest$statistic), digits=3)
change.table[change.table$species == i, "P.hfi"]  <- round(abs(temp.ttest$p.value), digits=3)
}

## change.table[change.table$P.hfi <= 0, "P.hfi"]  <- "$<$0.001" 
## change.table[change.table$HFI <= 0.1, "HFI"]  <- "$<$1"
## change.table[change.table$h.se <= 0.1, "h.se"]  <- "$<$1" 

change.table[change.table$native!="Exotic","native"] <- "Native"; change.table[change.table$form %in%  c("Shrub", "Tree", "Vine","SubShrub", "Treefern"), "form"] <- "Woody"

change.table <- rbind(c("Native herbaceous plants", rep("", length(names(change.table))-3)),                      
                       change.table[change.table$native=="Native" & change.table$form=="Forb", names(change.table)[1:12]],
          c("Exotic herbaceous plants", rep("", length(names(change.table))-3)),                      
                     change.table[change.table$native=="Exotic" & change.table$form=="Forb", names(change.table)[1:12]],
          c("Native grasses, tussocks, sedges and rushes", rep("", length(names(change.table))-3)),                      
                      change.table[change.table$native=="Native" & change.table$form=="Graminoid", names(change.table)[1:12]],
          c("Exotic grasses, sedges and rushes", rep("", length(names(change.table))-3)),                      
                      change.table[change.table$native=="Exotic" & change.table$form=="Graminoid", names(change.table)[1:12]],
          c("Native woody plants", rep("", length(names(change.table))-3)),                      
                      change.table[change.table$native=="Native" & change.table$form=="Woody", names(change.table)[1:12]],
          c("Exotic woody plants", rep("", length(names(change.table))-3)),                      
                      change.table[change.table$native=="Exotic" & change.table$form=="Woody", names(change.table)[1:12]])


text.stats  <- list(four.exotics.range, five.daisies.range, five.grasses.DCArange)

save(change.table, glmm.model.spp.recce, recce.spp.R2, glmm.model.hfi, hfi.spp.R2 ,file="/home/sean/Documents/Molesworth/SppTablesStats.Rdata")




     ############################################    ############################################    ############################################
             ############################################  Statistical Modelling   ############################################
## Following reviewers comments additional models were trialled. 

##dat <-  dat[dat$year!=2012, ]

## m1 <- gls(DCA1  ~ (year + x + y)^2 + dem.altitude + terrace + oversow,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## m1.f <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow +fenced,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)

## m1.tmb <- gl(DCA1  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2 + terrace + oversow +  fenced + (1|plot), correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)

## m1.tmb <- glmmTMB(cbind(native, exotic)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2 + terrace + oversow +  fenced + (1|plot),
##                   data=native.dat, family=binomial)

## m1.cor <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow,  data=dat)
## m2 <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow + grazing1987,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## m3 <- gls(DCA1  ~ year + fenced + x + y + dem.altitude + terrace + oversow + grazing1987,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## anova(m1, m1.cor, m2, m3)  ## M1 the winner


## library("sjstats")


## m1 <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)

## m1.interactions <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow + year * oversow,  data=dat)
## m1.interactions.cor <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow + year * oversow,  correlation = corExp(form = ~ x + y, nugget = TRUE),  data=dat)

## m1.interactions.xy <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow + year * oversow + year * x + year * y,  data=dat)

## m1.interactions.cor <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow + year * oversow, correlation = corExp(form = ~ x + y, nugget = TRUE),  data=dat)
## m2.interactions.cor <- gls(DCA2  ~ year + dem.altitude, correlation = corExp(form = ~ x + y, nugget = TRUE),  data=dat)

## fit1ml <- update(m1.interactions.cor, . ~ ., method = "ML")
## fit2ml <- update(m2.interactions.cor, . ~ ., method = "ML")

## fit3ml <- update(m1.interactions.xy.cor, . ~ ., method = "ML")
## anova(fit1ml, fit2ml)#, fit3ml)

## m2 <- gls(DCA2  ~ scale(year) + scale(x) + scale(y) + scale(dem.altitude) + terrace + oversow,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)

## ## m2.lme <- lme(DCA1  ~  scale(x) + scale(y) + terrace  + scale(dem.altitude) + oversow + fenced,  random=~1|year, data=dat)
## m1.lm <- lm(DCA1  ~  (scale(year) + scale(x) + scale(y)  + scale(dem.altitude) + terrace + oversow + fenced)^2, data=dat)
## m1.gls <- gls(DCA1  ~  (scale(year) + scale(x) + scale(y)  + scale(dem.altitude) + terrace + oversow + fenced)^2, correlation = corExp(form = ~ x + y), data=dat)
## m1.lme  <- lme(DCA1  ~  (scale(year) + scale(x) + scale(y)  + scale(dem.altitude) + terrace + oversow + fenced)^2,  random=~1|plot, data=dat)
## m1.lme.s  <- lme(DCA1  ~  (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced))^2,  random=~1|plot, correlation = corExp(form = ~ x + y), data=dat) 

## dat$pos <- numFactor(dat$x, dat$y) ###Check not duplicated dat$pos[duplicated(dat$pos)==TRUE]
## dat$group <- factor(rep(1, nrow(dat))) ###Dummy grouping variable


## ###Remove scale from categorical varibales
## ##m.1 <- glmmTMB(DCA1 ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced))^2 +(1|plot), data=dat) #exp(pos + 0 | group)
## ####Save instead of run to save time
 ## m.1 <- glmmTMB(DCA1 ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 +(1|plot), data=dat) #exp(pos + 0 | group)
 ## m.null <- glmmTMB(DCA1 ~ 1 +  (1|plot), data=dat, family = gaussian)
 ## save(m.1, m.null, file="/home/sean/Documents/Molesworth/Data/Objects/glmmTMBDCA1.Rdata")
load("/home/sean/Documents/Molesworth/Data/Objects/glmmTMBDCA1.Rdata")

##recce.table.dca1  <- glmmTMB.latex(m.1)
##R2.1  <- (cor(dat$hfi, dat$model.fit))^2
recce.dca1.R2  <- MuMIn::r.squaredGLMM(m.1, m.null)
glmm.model.DCA1  <-  glmmTMB.latex(m.1)

##m.1.step  <- buildglmmTMB(DCA1 ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced))^2 +(1|plot), data=dat)

##1 - sum((dat$hfi - dat$model.fit)^2)/sum((dat$hfi-mean(dat$hfi))^2)

###

## ##m.hfi.null <- lme(distance ~ 1, ~ 1 | Subject, data = Orthodont)
## m.hfi.null <- glmmTMB(DCA1 ~ 1 +  (1|plot), data=dat, family = gaussian)

##MuMIn::r.squaredGLMM(m.1, m.hfi.null)

## ####Spatial model of 2016 data so that Moran's I can be calculated    
## ##m.spatial  <- glmmTMB(hfi  ~ scale(x) + scale(y) + scale(dem.altitude) +  (1|plot) , data=dat.spatial, family = nbinom2)
## ## To fit the model, a numFactor and a dummy grouping variable must be added to the dataset:
## dat$pos <- numFactor(dat$x, dat$y)
## dat$group <- factor(rep(1, nrow(dat)))
## m.spatial.f <- glmmTMB(DCA1 ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced))^2  +  exp(pos + 0 | group), data=dat)
## simulationOutput <- simulateResiduals(fittedModel = m.spatial.f, plot = F)
## testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat.spatial$x, y = dat.spatial$y, plot=FALSE)

## ####Save instead of run to save time
## m.2 <- glmmTMB(DCA2 ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 +(1|plot), data=dat) #exp(pos + 0 | group)
## m.null <- glmmTMB(DCA2 ~ 1 +  (1|plot), data=dat, family = gaussian)
## save(m.2, m.null, file="/home/sean/Documents/Molesworth/Data/Objects/glmmTMBDCA2.Rdata")
load("/home/sean/Documents/Molesworth/Data/Objects/glmmTMBDCA2.Rdata")

##m.2.step  <- buildglmmTMB(DCA2 ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced))^2 +(1|plot), data=dat) #

## simulationOutput <- simulateResiduals(fittedModel = m.1, plot = F)
## plot(simulationOutput)
## testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat$x, y = dat$y, plot=FALSE)
glmm.model.DCA2  <-  glmmTMB.latex(m.2) ##recce.table.dca2  <- glmmTMB.latex(m.2)
recce.dca2.R2  <- MuMIn::r.squaredGLMM(m.2, m.null)


##qqnorm(m.1)

## vario2 <- Variogram(m2.lme.cor, form = ~x + y, resType = "pearson")
## plot(vario2, smooth = TRUE, ylim = c(0, 1.5))


## library(gstat)
## vario <- variogram(DCA1~1, data=dat, locations= ~x+y, cutoff=1000)
## plot(vario$dist, vario$gamma)

## library(ape)
## dca.dists <- as.matrix(dist(cbind(dat$x, dat$y)))

## dca.dists.inv <- 1/dca.dists
## diag(dca.dists.inv) <- 0
 
## ##dca.dists.inv[1:5, 1:5]
## Moran.I(dat$DCA1, dca.dists.inv)
## Moran.I(dat$DCA2, dca.dists.inv)



## acf(residuals(m2.lme.cor, retype="normalized"))
## pacf(residuals(m2.lme, retype="normalized"))
## plot(residuals(m2.lme, retype="normalized")~dat$y)




## simulationOutput <- simulateResiduals(fittedModel = m2.lme, plot = F)
## plot(simulationOutput)
## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)
## testSpatialAutocorrelation(simulationOutput = simulationOutput, x = native.dat$x, y= native.dat$y)


## m1.cor <- gls(DCA2  ~ year + x + y + dem.altitude + terrace + oversow,  data=dat)
## m2 <- gls(DCA2  ~ year + x + y + dem.altitude + terrace + oversow + grazing1987,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## m3 <- gls(DCA2  ~ year + fenced + x + y + dem.altitude + terrace + oversow + grazing1987,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## anova(m1, m1.cor, m2, m3)  ## M1.2 the winner



## ##m.dca1 <- as.data.frame(round(summary(fit1ml)$tTable, digits=3))
## m.dca1 <- as.data.frame(round(summary(m1)$tTable, digits=3))
## m.dca1[,names(m.dca1)] <- lapply(names(m.dca1), function(x) as.character(m.dca1[,x]))
## m.dca1[ m.dca1=="0"] <-  "$<$0.001"
## m.dca2 <- as.data.frame(round(summary(m2)$tTable, digits=3))
## m.dca2[,names(m.dca2)] <- lapply(names(m.dca2), function(x) as.character(m.dca2[,x]))
## m.dca2[m.dca2=="0"] <-  "$<$0.001"
## common.recce.table <- rbind(rep("", length(names(m.dca1))), m.dca1, rep("", length(names(m.dca1))), m.dca2)
## save(common.recce.table, common.spp, file="/home/sean/Documents/Molesworth2020/Analysis/Objects/CommonRecceTable.Rdata")



###############RGB spatial prediction for common species composition
s.e. <- function(x) sqrt(var(na.omit(x))/length(na.omit(x))); l.s.d <-  function(x) (sqrt(var(na.omit(x))/length(na.omit(x))))*(qt(0.95,(length(na.omit(x))))) ## SEM  and the more conservative least significant difference

dat.locs[dat.locs$year %in% c(1952, 1960, 1972), "year"] <- 1952; dat.locs[dat.locs$year %in% c(1989, 1995, 2001, 2006, 2007, 2008, 2012, 2016), "year"] <- 2016 ## Pool into three survey groups


text.colours  <- expand.grid(year = c(1952, 1987, 2016), ValleyEle = c("high", "low"), predict.m1=NA, predict.m2=NA,  red=NA, green=NA, blue=NA) 



####### Colours for figs 1 and 2 - combined functional types
rgb2hex <- function(r,g,b) sprintf('#%s',paste(as.hexmode(c(r,g,b)),collapse = ''))
dca1.native.woody  <- mean(dca.spp[dca.spp$species %in% common.table[common.table$form == "Woody" & common.table$native=="Native", "species"], "DCA1"])
dca2.native.woody  <- mean(dca.spp[dca.spp$species %in% common.table[common.table$form == "Woody" & common.table$native=="Native", "species"], "DCA2"])
dca1.native.herb  <- mean(dca.spp[dca.spp$species %in% common.table[common.table$form %in% c("Graminoid", "Forb") & common.table$native=="Native", "species"], "DCA1"])
dca2.native.herb  <- mean(dca.spp[dca.spp$species %in% common.table[common.table$form %in% c("Graminoid", "Forb") & common.table$native=="Native", "species"], "DCA2"])
dca1.exotic.woody  <- mean(dca.spp[dca.spp$species %in% common.table[common.table$form == "Woody" & common.table$native=="Exotic", "species"], "DCA1"])
dca2.exotic.woody  <- mean(dca.spp[dca.spp$species %in% common.table[common.table$form == "Woody" & common.table$native=="Exotic", "species"], "DCA2"])
dca1.exotic.herb  <- mean(dca.spp[dca.spp$species %in% common.table[common.table$form %in% c("Graminoid", "Forb") & common.table$native=="Exotic", "species"], "DCA1"])
dca2.exotic.herb  <- mean(dca.spp[dca.spp$species %in% common.table[common.table$form %in% c("Graminoid", "Forb") & common.table$native=="Exotic", "species"], "DCA2"])

native.woody.col.rgb  <- c(1- (dca1.native.woody / max(dca.spp$DCA1)), 1 -(((max(dca.spp$DCA1) - dca1.native.woody)) / (max(dca.spp$DCA1))), 1- (dca2.native.woody / max(dca.spp$DCA2)))
native.herb.col.rgb  <- c(1- (dca1.native.herb / max(dca.spp$DCA1)), 1 -(((max(dca.spp$DCA1) - dca1.native.herb)) / (max(dca.spp$DCA1))), 1- (dca2.native.herb / max(dca.spp$DCA2)))
exotic.woody.col.rgb  <- c(1- (dca1.exotic.woody / max(dca.spp$DCA1)), 1 -(((max(dca.spp$DCA1) - dca1.exotic.woody)) / (max(dca.spp$DCA1))), 1- (dca2.exotic.woody / max(dca.spp$DCA2)))
exotic.herb.col.rgb  <- c(1- (dca1.exotic.herb / max(dca.spp$DCA1)), 1 -(((max(dca.spp$DCA1) - dca1.exotic.herb)) / (max(dca.spp$DCA1))), 1- (dca2.exotic.herb / max(dca.spp$DCA2)))


nat.herb.col  <- rgb2hex(as.integer(native.herb.col.rgb[1]*255), as.integer(native.herb.col.rgb[2]*255), as.integer(native.herb.col.rgb[3]*255))
exot.herb.col  <- rgb2hex(as.integer(exotic.herb.col.rgb[1]*255), as.integer(exotic.herb.col.rgb[2]*255), as.integer(exotic.herb.col.rgb[3]*255))
nat.woody.col  <- rgb2hex(as.integer(native.woody.col.rgb[1]*255), as.integer(native.woody.col.rgb[2]*255), as.integer(native.woody.col.rgb[3]*255))
exot.woody.col  <- "#F40A02"  ## rgb2hex(as.integer(exotic.woody.col.rgb[1]*255), as.integer(exotic.woody.col.rgb[2]*255), as.integer(exotic.woody.col.rgb[3]*255)) ###Didnt work

## nat.herb.col  <- "#3aeb13"
## nat.wood.col  <- "#d1d487"
## exot.herb.col  <- "#cf0e0e"
## exot.wood.col  <- "#8c3424"
save(nat.woody.col, exot.woody.col, nat.herb.col, exot.herb.col, native.woody.col.rgb, exotic.woody.col.rgb, native.herb.col.rgb, exotic.herb.col.rgb, file="Data/Objects/SpeciesColourCodesGraphics.rdata")  ###For HFI graphs



pdf("/home/sean/Documents/Molesworth/Graphs/CommonSpeciesMap.pdf")

par(mfrow = c(3, 1), cex = 0.6, mar = c(0, 0, 0, 0), oma = c(0.5, 2.5, 2.5, 0.5))
for (i in c(1952, 1987, 2016)){  ## ##    i <-  1952  ##sort(unique(dat$year)
grid.predict <-  as.data.frame(cbind(grid.alt, grid.terrace[1], grid.oversow[1], grid.fenced[1], rep(i, length(grid.alt[,1])))); names(grid.predict) <-  c("dem.altitude", "x", "y", "terrace", "oversow", "fenced", "year")
grid.predict <- grid.predict[is.na(grid.predict$dem.altitude)==FALSE,]
grid.predict[grid.predict$fenced!="Unfenced", "fenced"]  <- "Fenced"
   
grid.predict[grid.predict$terrace == "Slope", "terrace"]  <- 0; grid.predict[grid.predict$terrace == "Terrace", "terrace"]  <- 1;
grid.predict$terrace  <- as.integer(as.character(grid.predict$terrace))
grid.predict[grid.predict$fenced == "Unfenced", "fenced"]  <- 0; grid.predict[grid.predict$fenced == "Fenced", "fenced"]  <- 1;
grid.predict$fenced  <- as.integer(as.character(grid.predict$fenced))
grid.predict[grid.predict$oversow == "no", "oversow"]  <- 0; grid.predict[grid.predict$oversow == "yes", "oversow"]  <- 1;
grid.predict$oversow  <- as.integer(as.character(grid.predict$oversow))
grid.predict$plot  <- NA

predict.m1 <- predict(m.1, grid.predict, type="response", re.form=NA, se.fit = FALSE); predict.m2 <- predict(m.2, grid.predict, type="response", re.form=NA, se.fit = FALSE)
predict.dca <- cbind(grid.predict[, c("x", "y")], predict.m1, predict.m2)
predict.dca$red <- 1- abs(predict.dca$predict.m1 / max(predict.dca$predict.m1, na.rm=TRUE)) ##red band. Increasing red with DCA 1
predict.dca$green <- 1 - abs((max(predict.dca$predict.m1, na.rm=TRUE)) - predict.dca$predict.m1) / max(predict.dca$predict.m1, na.rm=TRUE)  ##green band. Decreasing green with DCA 1
predict.dca$blue <-1 -  abs(predict.dca$predict.m2 / max(predict.dca$predict.m2, na.rm=TRUE)) ##blue band. INcreasing with DCA 2
dca.r <- predict.dca[is.na(predict.dca$predict.m1)==FALSE, c("x", "y", "red", "green", "blue")]
##coordinates(dca.r) <- c("x", "y"); proj4string(dca.r) <-  CRS("+init=epsg:2193")##NZTM
dca.r <- rasterFromXYZ(dca.r, crs = "+init=epsg:2193")
##grid.check <- grid.predict[grid.predict$year==2016,]; coordinates(grid.check) <- c("x", "y"); proj4string(grid.check) <- CRS("+init=epsg:2193") #georeference to NZTM
##writeRaster(dca.r, filename=paste("/home/sean/Documents/Molesworth2020/Analysis/Graphs/CommonSpeciesMap", i,".tif", sep=""), bandorder="BIL", overwrite=TRUE) }

##pdf(paste("/home/sean/Documents/Molesworth2016/Graphs/CommonSpeciesMap",i,".pdf", sep=""))
plotRGB(dca.r, r="red", g="green", b="blue", scale=1.1, interpolate=TRUE, asp=1, box=FALSE, axes=FALSE, labels= FALSE, frame.plot=FALSE) #xlim=c(min(1570000), max(1640000)) 
plot(dat.locs[dat.locs$year==i,], add=TRUE, pch=19, cex=0.5)
text(1590000, 5350000, i)##, col.main="black")#,  sub="Predictive model from recce data 1952-2016", col.sub="grey", cex.lab=0.5)
#text(1590000, 5353000, "Common species DCA ordination", col="black", cex=0.75) 
##dev.off()
assign(paste("predict.dca1.", i, sep=""), c(mean(predict.m1), s.e.(predict.m1), mean(predict.m2), s.e.(predict.m2), i))
##mtext("Change in species composition", outer = TRUE)


## Lower Clarence across River from Quail Flat 460 m 1638182 5330094 
text.colours[text.colours$year == i & text.colours$ValleyEle == "low",
               c("predict.m1", "predict.m2", "red", "green", "blue") ]  <- predict.dca[predict.dca$x==1638182 & predict.dca$y == 5330094,
                                                                                       c("predict.m1", "predict.m2", "red", "green", "blue") ]

###Nearby St Bernard at 2233 m. About 6 km away 1635382 5335094
text.colours[text.colours$year == i & text.colours$ValleyEle == "high",
               c("predict.m1", "predict.m2", "red", "green", "blue") ]  <- predict.dca[predict.dca$x==1635382 & predict.dca$y == 5335094,
                                                                                       c("predict.m1", "predict.m2", "red", "green", "blue") ]

## #####R1 at 688 m
## text.colours[text.colours$year == i & text.colours$ValleyEle == "low",
##              c("predict.m1", "predict.m2", "red", "green", "blue") ] <-  predict.dca[predict.dca$x==1593482 & predict.dca$y == 5325794,
##                                                                                      c("predict.m1", "predict.m2", "red", "green", "blue") ]
## #### G-1 at 1464 m
## text.colours[text.colours$year == i & text.colours$ValleyEle == "high",
##                c("predict.m1", "predict.m2", "red", "green", "blue") ]  <- predict.dca[predict.dca$x==1593496  & predict.dca$y == 5325826,
###                                                                                        c("predict.m1", "predict.m2", "red", "green", "blue") ]
}

dev.off()

##### Low plots are 460 m. High plots 2230 grid.predict[grid.predict$dem.altitude <= min(grid.predict$dem.altitude) +10, ]

###R-1 688 m  1601670 5307040, G-1 1464 1593496 5325826  grid.predict[grid.predict$dem.altitude >= max(grid.predict$dem.altitude) - 10, ]


ord.fit.dca <- envfit(dat[, c("DCA1", "DCA2")]  ~ dat$year  + dat$x + dat$y + dat$dem.altitude + dat$terrace + dat$oversow +dat$fenced)
dca.arrows <-  as.data.frame(cbind((ord.fit.dca$vectors)$arrows, (ord.fit.dca$vectors)$r, (ord.fit.dca$vectors)$pvals))### Extract ordfit vectors
dca.arrows$env.var <- rownames(dca.arrows); 
dca.arrows[dca.arrows$env.var == "dat$year", "env.var"] <- "Year"
dca.arrows[dca.arrows$env.var == "dat$x", "env.var"] <- "East"
dca.arrows[dca.arrows$env.var == "dat$y", "env.var"] <- "North"
dca.arrows[dca.arrows$env.var == "dat$dem.altitude", "env.var"] <- "Altitude"
dca.arrows[dca.arrows$env.var == "dat$terrace", "env.var"] <- "Terrace"
dca.arrows[dca.arrows$env.var == "dat$oversow", "env.var"] <- "Oversow"
dca.arrows[dca.arrows$env.var == "dat$fenced", "env.var"] <- "Fenced"
names(dca.arrows) <- c("DCA1", "DCA2", "r2", "P", "env.var")
## ord.fit.factors <- as.data.frame(ord.fit.dca$vectors[1])
## ##ord.fit.factors <-  cbind( rbind((ord.fit.factors["dat.terraceTerrace",] - ord.fit.factors["dat.terraceSlope",]) , (ord.fit.factors["dat.oversowyes",] - ord.fit.factors["dat.oversowno",])),   as.data.frame(ord.fit.dca$factors[2]), c(.001,.001), c("Terrace", "Oversown")) ;
## names(ord.fit.factors) <- names(dca.arrows)
## dca.arrows <-  rbind(dca.arrows, ord.fit.factors)
dca.arrows$position <- c(4,2,4,4,1,3,2) ##Specify direction of offset


##dca.arrows$position <- c(3,3,3,4,4) ##Specify direction of offset

evals <- ord$evals ##Extract eignevalues

dca.site.colours <- (cbind(dca.site, 1- (dca.site$DCA1 / max(dca.site$DCA1)), 1 -(((max(dca.site$DCA1) - dca.site$DCA1)) / (max(dca.site$DCA1))), 1- (dca.site$DCA2 / max(dca.site$DCA2))))
names(dca.site.colours) <- c("DCA1", "DCA2", "year", "plot", "red", "green", "blue")


lambda.dca2  <- as.character(round((evals["DCA2"]/sum(evals)*100), digits=1))



pdf("/home/sean/Documents/Molesworth/Graphs/CommonSppDCASite.pdf")
plot(dca.site$DCA2 ~ dca.site$DCA1, ylim=c(min(-1.5), max(ceiling(max(dca.site$DCA2)))), xlim=c(min(-1.5), max(ceiling(max(dca.site$DCA1)))), bty="l", type="n", cex.axis = 1.25, cex.lab = 0.9,
     ylab = bquote("DCA 2 scores (" ~ lambda ~ "  = " ~  .(round((evals["DCA2"]/sum(evals)*100), digits=1))~ "%)"), ##sep="")), #expression(paste(plain(sin) * phi, "  and  ", plain(cos) * phi)),
         ##TeX(paste("DCA 2 scores ($\\lambda  = ",  round((evals["DCA2"]/sum(evals)*100), digits=1) ,"%)", sep="")),
     xlab = bquote("DCA 1 scores (" ~ lambda ~ "  = " ~  .(round((evals["DCA1"]/sum(evals)*100), digits=1))~ "%)"))
         ##TeX(paste("DCA 1 scores ($\\lambda  = ",  round((evals["DCA1"]/sum(evals)*100), digits=1) ,"%)", sep="")))
points(dca.site$DCA2 ~ dca.site$DCA1, pch=19, cex = .8, type="p", col= rgb(dca.site.colours[,c(5:7)]))
## points(dca.spp$DCA2 ~ dca.spp$DCA1, pch=19, cex = 0.5, type="p")
## pointLabel(dca.spp$DCA2 ~ dca.spp$DCA1, labels  = dca.spp$species, cex = 0.75, col = "black", offset =0.1, vfont = c("serif", "italic"),  method = c("SANN"), allowSmallOverlap = FALSE, doPlot = TRUE)
text(dca.arrows$DCA2 ~ dca.arrows$DCA1, labels = dca.arrows$env.var, cex = 1, col = "black", pos = dca.arrows$position, offset = .5)
arrows(0, 0, dca.arrows$DCA1, dca.arrows$DCA2, length = 0.2, angle = 15, code = 2)

for (i in 1:length(dca.arrows$env.var)){    ##i <- 1   
text(dca.arrows[i,"DCA2"]-0.15 ~ dca.arrows[i, "DCA1"], labels = bquote(italic(R)^2 ~ "=" ~ .(round(dca.arrows[i, "r2"], digits=2)) ~ "P =" ~ .(round(dca.arrows[i, "P"], digits=3))), cex = 0.65, col = "black", pos = dca.arrows[i, "position"])
}
##text(dca.arrows$DCA2 ~ dca.arrows$DCA1, labels = paste(expression(italic(R)^2, " = ", round(dca.arrows[i, "r2"], digits=2))), cex = 1, col = "black", pos = c(2,1,1,4,4), offset = .5)   

dev.off()

##################################################   Graph change by year  ##############################################################################

dca.year <- as.data.frame(cbind(tapply(dat$DCA1, dat$year, mean), tapply(dat$DCA1, dat$year, s.e.), tapply(dat$DCA2, dat$year, mean), tapply(dat$DCA2, dat$year, s.e.))); dca.year$year <- rownames(dca.year)
names(dca.year) <- c("dca1", "dca1.se", "dca2", "dca2.se", "year") 
dca.year$red <- 1- abs(dca.year$dca1 / max(predict.dca$predict.m1, na.rm=TRUE)) ##red band. Increasing red with DCA 1
dca.year$green <- 1 - abs((max(predict.dca$predict.m1, na.rm=TRUE)) - dca.year$dca1) / max(predict.dca$predict.m1, na.rm=TRUE)  ##green band. Decreasing green with DCA 1
dca.year$blue <-1 -  abs(dca.year$dca2 / max(predict.dca$predict.m2, na.rm=TRUE)) 

fit.raw.1 <- as.data.frame(cbind(dat$year, predict(m.1, dat, type="response", re.form=NA, se.fit = FALSE)));
names(fit.raw.1) <- c("year", "fitted");
fit.raw.2 <- as.data.frame(cbind(dat$year, predict(m.2, dat, type="response", re.form=NA, se.fit = FALSE))); names(fit.raw.2) <- c("year", "fitted")
fitted <- as.data.frame(cbind(tapply(fit.raw.1$fitted, fit.raw.1$year, mean), tapply(fit.raw.1$fitted, fit.raw.1$year, s.e.),
                              tapply(fit.raw.2$fitted, fit.raw.2$year, mean), tapply(fit.raw.2$fitted, fit.raw.2$year, s.e.)))
fitted$year <- rownames(fitted); names(fitted) <- c("fitted1", "fitted1.se", "fitted2", "fitted2.se", "year") 
fitted$red <- 1 - abs(fitted$fitted1 / max(predict.dca$predict.m1, na.rm=TRUE)) ##red band. Increasing red with DCA 1
fitted$green <- 1 - abs((max(predict.dca$predict.m1, na.rm=TRUE) - fitted$fitted1) / max(predict.dca$predict.m1, na.rm=TRUE))  ##green band. Decreasing green with DCA 1
fitted$blue <-1 -  abs(fitted$fitted2 / max(predict.dca$predict.m1, na.rm=TRUE)) 

for (i in sort(unique(dat$year))){  
grid.predict <-  as.data.frame(cbind(grid.alt, grid.terrace[1], grid.oversow[1], grid.fenced[1], rep(i, length(grid.alt[,1]))));
names(grid.predict) <-  c("dem.altitude", "x", "y", "terrace", "oversow", "fenced", "year")

grid.predict[grid.predict$terrace == "Slope", "terrace"]  <- 0; grid.predict[grid.predict$terrace == "Terrace", "terrace"]  <- 1; grid.predict$terrace  <- as.integer(as.character(grid.predict$terrace))
grid.predict[grid.predict$fenced == "Unfenced", "fenced"]  <- 0; grid.predict[grid.predict$fenced == "Fenced", "fenced"]  <- 1; grid.predict$fenced  <- as.integer(as.character(grid.predict$fenced))
grid.predict[grid.predict$oversow == "no", "oversow"]  <- 0; grid.predict[grid.predict$oversow == "yes", "oversow"]  <- 1; grid.predict$oversow  <- as.integer(as.character(grid.predict$oversow))
grid.predict$plot  <- NA
grid.predict <- grid.predict[is.na(grid.predict$dem.altitude)==FALSE,]

predict.m1 <- predict(m.1, grid.predict, type="response", re.form=NA, se.fit = FALSE);
predict.m2 <- predict(m.2, grid.predict, type="response", re.form=NA, se.fit = FALSE);
predict.dca <- cbind(grid.predict[, c("x", "y")], predict.m1, predict.m2)
assign(paste("predict.dca1.", i, sep=""), c(mean(predict.m1, na.rm=TRUE), s.e.(predict.m1), quantile(predict.m1, na.rm=TRUE)["25%"], quantile(predict.m1, na.rm=TRUE)["75%"],
                                            mean(predict.m2, na.rm=TRUE), s.e.(predict.m2), quantile(predict.m2, na.rm=TRUE)["25%"], quantile(predict.m2, na.rm=TRUE)["75%"], i))
}

predicted <- as.data.frame(rbind(predict.dca1.1952, predict.dca1.1960,predict.dca1.1972, predict.dca1.1987, predict.dca1.2007, predict.dca1.2016))
names(predicted) <- c("predicted1", "predicted1.se", "predicted1.25", "predicted1.75", "predicted2", "predicted2.se", "predicted2.25", "predicted2.75", "year")      
predicted$red <- 1 - abs(predicted$predicted1 / max(predict.dca$predict.m1, na.rm=TRUE)) ##red band. Increasing red with DCA 1
predicted$green <- 1 - abs((max(predict.dca$predict.m1, na.rm=TRUE) - predicted$predicted1) / max(predict.dca$predict.m1, na.rm=TRUE))  ##green band. Decreasing green with DCA 1
predicted$blue <-1 -  abs(predicted$predicted2 / max(predict.dca$predict.m2, na.rm=TRUE)) 

predicted$red.25 <- 1 - abs(predicted$predicted1.25 / max(predict.dca$predict.m1, na.rm=TRUE)) ##red band. 
predicted$green.25 <- 1 - abs((max(predict.dca$predict.m1, na.rm=TRUE) - predicted$predicted1.25) / max(predict.dca$predict.m1, na.rm=TRUE))  ##green band. 
predicted$blue.25 <-1 -  abs(predicted$predicted2.25 / max(predict.dca$predict.m2, na.rm=TRUE)) 

predicted$red.75 <- 1 - abs(predicted$predicted1.75 / max(predict.dca$predict.m1, na.rm=TRUE)) ##red band. 
predicted$green.75 <- 1 - abs((max(predict.dca$predict.m1, na.rm=TRUE) - predicted$predicted1.75) / max(predict.dca$predict.m1, na.rm=TRUE))  ##green band. 
predicted$blue.75 <-1 -  abs(predicted$predicted2.75 / max(predict.dca$predict.m2, na.rm=TRUE)) 


pdf("/home/sean/Documents/Molesworth/Graphs/CommonSppDCASpecies.pdf")
plot(dca.site$DCA2 ~ dca.site$DCA1, ylim=c(min(-1), max(ceiling(max(dca.site$DCA2)))), xlim=c(min(-1), max(ceiling(max(dca.site$DCA1)))), bty="l", type="n", cex.axis = 1.25, cex.lab = 0.9,
     ylab = bquote("DCA 2 scores (" ~ lambda ~ "  = " ~  .(round((evals["DCA2"]/sum(evals)*100), digits=1))~ "%)"),
         ##TeX(paste("DCA 2 scores ($\\lambda  = ",  round((evals["DCA2"]/sum(evals)*100), digits=1) ,"%)", sep="")),
     xlab = bquote("DCA 1 scores (" ~ lambda ~ "  = " ~  .(round((evals["DCA1"]/sum(evals)*100), digits=1))~ "%)"))
         ##TeX(paste("DCA 1 scores ($\\lambda  = ",  round((evals["DCA1"]/sum(evals)*100), digits=1) ,"%)", sep="")))
points(dca.spp$DCA2 ~ dca.spp$DCA1, pch=19, cex = 0.5, type="p")
pointLabel(dca.spp$DCA2 ~ dca.spp$DCA1, labels  = dca.spp$species, cex = 0.75, col = "black", offset =0.1, vfont = c("serif", "italic"),  method = c("SANN"), allowSmallOverlap = FALSE, doPlot = TRUE)
points(predicted$predicted2 ~ predicted$predicted1, pch=19, cex = 3, type="p" , col= rgb(predicted[,c(10:12)]))
points(predicted$predicted2 ~ predicted$predicted1, pch=19, cex = 1, type="l" , col= "black", lwd=1)
text(predicted$predicted2[c(1,6)] ~ predicted$predicted1[c(1,6)], labels = predicted$year[c(1,6)], cex = 1, col =  "black", pos = c(4,2), offset =c(.5, 1.5))
arrows(predicted$predicted1, predicted$predicted2, predicted$predicted1 + predicted$predicted1.se, length=0.025, angle=90, code=3, lwd = 0.5) ## error bar right
arrows(predicted$predicted1, predicted$predicted2, predicted$predicted1 - predicted$predicted1.se, length=0.025, angle=90, code=3, lwd = 0.5) ## error bar left
arrows(predicted$predicted1, predicted$predicted2,predicted$predicted1, predicted$predicted2 + predicted$predicted2.se, length=0.025, angle=90, code=3, lwd = 0.5) ## error bar up
arrows(predicted$predicted1, predicted$predicted2, predicted$predicted1,predicted$predicted2 - predicted$predicted2.se, length=0.025, angle=90, code=3, lwd = 0.5) ## error bar down
dev.off()


DCA1.change <- round(((predicted$predicted1[2]- predicted$predicted1[1])/(range(dca.site$DCA1)[2] - range(dca.site$DCA1)[1]))*100, digits=1)
DCA2.change <- round(((predicted$predicted2[2]- predicted$predicted2[1])/(range(dca.site$DCA2)[2] - range(dca.site$DCA2)[1]))*100, digits=1)
DCA1.range <- round((range(dca.site$DCA1)[2] - range(dca.site$DCA1)[1]), digits=1); DCA2.range <- round((range(dca.site$DCA2)[2] - range(dca.site$DCA2)[1]), digits=1)

save(DCA1.change, DCA2.change, DCA1.range, DCA2.range, common.table.recce, predicted, file="/home/sean/Documents/Molesworth/Data/Objects/RecceStatistics.Rdata")





























###############   Delete from here. 





pdf("/home/sean/Documents/Molesworth/Graphs/TableAsPlotDCA1.pdf")
plot(dca.spp$DCA1 ~ c(1:length(dca.spp$species)) , bty="l", cex.axis = 1.25, cex.lab = 0.9, ylab = "DCA axis 2 scores", xaxt = "n", xlab="", pch=19, cex=1.1)
axis(1, at=c(1:length(dca.spp$species)),  labels = FALSE) #
text(1:length(dca.spp$species), par("usr")[3] - 0.15, srt = 45, adj = 0.9,  labels = dca.spp$species, xpd = TRUE, cex=0.75) ## 
dev.off()

pdf("/home/sean/Documents/Molesworth/Graphs/TableAsPlotDCA2.pdf")
plot(dca.spp$DCA2 ~ c(1:length(dca.spp$species)) , bty="l", cex.axis = 1.25, cex.lab = 0.9, ylab = "DCA axis 2 scores", xaxt = "n", xlab="", pch=19, cex=1.1)
axis(1, at=c(1:length(dca.spp$species)),  labels = FALSE) #
text(1:length(dca.spp$species), par("usr")[3] - 0.15, srt = 45, adj = 0.9,  labels = dca.spp$species, xpd = TRUE, cex=0.75) ## 
dev.off()


load("/home/sean/Documents/Molesworth/Data/Objects/RandomRecceTable.Rdata") ##Object = random.recce.table.dcamodel, random.recce.table.dcameans, random.common.tabl
load("/home/sean/Documents/Molesworth/Data/Objects/HFIcommonMeans.Table.Rdata") ## Object = hfi.common.means (includes lmer coefficents by spp), lmer.commonspp.coe

pdf("/home/sean/Documents/Molesworth/Graphs/TableAsPlotHFI.pdf")
plot(hfi.common.means$mean ~ c(1:length(hfi.common.means$species)) , bty="l", cex.axis = 1.25, cex.lab = 0.9,
     ylab = "Mean HFI", xaxt = "n", xlab="", pch=19, cex=1.1)
axis(1, at=c(1:length(hfi.common.means$species)),  labels = FALSE) #
text(1:length(hfi.common.means$species), par("usr")[3] - 0.15, srt = 45, adj = 0.9,  labels = hfi.common.means$species, xpd = TRUE, cex=0.75) ## 
dev.off()





## ####################################################################################################################################################################################################
## ######################################################                         Paired Plot DCA                                       ###############################################################

## paired.plots <- as.character(unique(recce.dat[recce.dat$year==2008, "plot"]))
## recce.paired <- recce.dat[recce.dat$plot %in% paired.plots,]
## ##recce.dat[recce.dat$species %in% common.spp, ]

## nplot <- length(table(recce.paired$plot.yr)) # what number of plots
## nspp  <- length(table(recce.paired$species))  # what number of species

## comp <- matrix(0, nrow=nplot, ncol=nspp) # produce a matix of length plots and width species
## rownames(comp) <- unique(recce.paired$plot.yr) # give row names plots 
## colnames(comp) <- unique(recce.paired$species) # give column names species

## plot <-recce.paired$plot.yr  # matrix py takes the form of py from data frame dat
## species <- recce.paired$species
## cover <- recce.paired$cover
        
## for(i in 1:nrow(recce.paired)) {
##   r <- which(rownames(comp)==plot[i])
##   c <- which(colnames(comp)==species[i])
##   comp[r, c] <- cover[i]
## }

## ord <- decorana(comp)
## dca.site <- as.data.frame(ord$rproj[, c(1:2)])
## dca.site$year <-  as.numeric(substr(rownames(dca.site), regexpr("\\.",rownames(dca.site))+1, regexpr("\\.", rownames(dca.site))+5))
## dca.site$site <- substr(rownames(dca.site), 1, 2)
## dca.site$plot <- substr(rownames(dca.site), 1, 4)
## dca.site$fenced <- "Fenced" ; dca.site[grep("Un", rownames(dca.site)), "fenced"] <- "Unfenced"


## fenced.year.dca1 <- round(as.data.frame(cbind(tapply(dca.site$DCA1, list(dca.site$year, dca.site$fenced),  mean), tapply(dca.site$DCA1, list(dca.site$year, dca.site$fenced),  s.e.))), digits=2)
## fenced.year.dca2 <- round(as.data.frame(cbind(tapply(dca.site$DCA2, list(dca.site$year, dca.site$fenced),  mean), tapply(dca.site$DCA2, list(dca.site$year, dca.site$fenced),  s.e.))), digits=2)

## fenced.recce.table <- rbind(c("DCA axis 1", rep("",length(names(fenced.year.dca1))-1)), fenced.year.dca1, c("DCA axis 2", rep("",length(names(fenced.year.dca2))-1)), fenced.year.dca2)[,c(1,3,2,4)] 
## names(fenced.recce.table) <- c("Fenced", "Fenced.sem", "Unfenced", "Unfenced.sem")

## ##m1 <- aov(DCA1 ~ plot + fenced*year, data=dca.site)  ##randomised block. Plot fixed effect 
## ##m1 <- aov(DCA1 ~ plot + fenced + year + Error(plot/fenced), data=dca.site)  ##fenced fixed effect, plot random effect. Need balanced model to use aov
## m2 <- lme(DCA1 ~ fenced + year, random=~1|plot, data=dca.site)  ##fenced fixed effect, plot random effect. 
## #m2a <- lme(DCA1 ~ fenced + year, random=~year|plot, data=dca.site)  ##fenced fixed effect, plot random effect. 
## m2.2 <- lme(DCA2 ~ fenced + year, random=~1|plot, data=dca.site)  ##fenced fixed effect, plot random effect. 

## ## m2.aov <- aov(DCA1 ~ fenced + year + Error(plot), data=dca.site) ##plot is a **random** error component of variance
## ## m2.lm <- lm(DCA1 ~ plot + fenced + year , data=dca.site) ##plot is fixed effect 

## ## coef(m2)[3]  ## fencedUnfenced coefficient - higher than average
## ## coef(m2)[4]  ## Year coefficient - lower over time than average

## paired.dca1 <- as.data.frame(round(summary(m2)$tTable, digits=3)); paired.dca2 <- round(summary(m2.2)$tTable, digits=3)
## paired.recce.table <- rbind(c("DCA axis 1", rep("",length(names(paired.dca1))-1)), paired.dca1, c("DCA axis 2", rep("",length(names(paired.dca1))-1)), paired.dca2)
## save(paired.recce.table, fenced.recce.table, file="/home/sean/Documents/Molesworth2020/Analysis/Objects/PairedRecceTable.Rdata")

   ####################################################################################################################################################################
####################################################       Random  Plot DCA - Plots established on lines in 2007     ####################################################

random.plots <- as.character(unique(recce.dat[recce.dat$year==2007, "plot"]))
recce.random <- recce.dat[recce.dat$plot %in% random.plots,]
##recce.dat[recce.dat$species %in% common.spp, ]

nplot <- length(table(recce.random$plot.yr)) # what number of plots
nspp  <- length(table(recce.random$species))  # what number of species

comp <- matrix(0, nrow=nplot, ncol=nspp) # produce a matix of length plots and width species
rownames(comp) <- unique(recce.random$plot.yr) # give row names plots 
colnames(comp) <- unique(recce.random$species) # give column names species

plot <-recce.random$plot.yr  # matrix py takes the form of py from data frame dat
species <- recce.random$species
cover <- recce.random$cover
        
for(i in 1:nrow(recce.random)) {
  r <- which(rownames(comp)==plot[i])
  c <- which(colnames(comp)==species[i])
  comp[r, c] <- cover[i]
}

ord <- decorana(comp)
dca.site <- as.data.frame(ord$rproj[, c(1:2)])
dca.site$year <-  as.numeric(substr(rownames(dca.site), regexpr("\\.",rownames(dca.site))+1, regexpr("\\.",rownames(dca.site))+ 5))
dca.site$plot <- as.character(substr(rownames(dca.site), 1, regexpr("\\.",rownames(dca.site))-1))
dca.spp <-  as.data.frame(ord$cproj[, c(1:2)]) ; dca.spp$species <- rownames(dca.spp)

random.year.dca1 <- round(as.data.frame(cbind(tapply(dca.site$DCA1, dca.site$year,  mean), tapply(dca.site$DCA1, dca.site$year,  s.e.))), digits=2)
random.year.dca2 <- round(as.data.frame(cbind(tapply(dca.site$DCA2, dca.site$year,  mean), tapply(dca.site$DCA2, dca.site$year,  s.e.))), digits=2)
random.recce.table.dcameans <- rbind(c("DCA axis 1", rep("",length(names(random.year.dca1))-1)), random.year.dca1, c("DCA axis 2", rep("",length(names(random.year.dca2))-1)), random.year.dca2)
names(random.recce.table.dcameans) <- c("dca", "dca.sem")

## m1 <- lme(DCA1 ~ year, random=~1|plot, data=dca.site)  ##fenced fixed effect, plot random effect. 
## ##m1a <- lme(DCA1 ~ year, random=~year|plot, data=dca.site)  ##fenced fixed effect, plot random effect. 
## m2 <- lme(DCA2 ~ year, random=~1|plot, data=dca.site)  ##fenced fixed effect, plot random effect. 
## ##m2a <- lme(DCA2 ~ year, random=~year|plot, data=dca.site)  ##fenced fixed effect, plot random effect. 

## random.dca1 <- as.data.frame(round(summary(m1)$tTable, digits=3)); random.dca2 <- round(summary(m2)$tTable, digits=3)
## random.recce.table.dcamodel <- rbind(c("DCA axis 1", rep("",length(names(random.dca1))-1)), random.dca1, c("DCA axis 2", rep("",length(names(random.dca1))-1)), random.dca2)
## random.recce.table.dcamodel[,names(random.recce.table.dcamodel)] <- lapply(names(random.recce.table.dcamodel), function(x) as.character(random.recce.table.dcamodel[,x]))
## random.recce.table.dcamodel[random.recce.table.dcamodel=="0"] <-  "$<$0.001"

######################Tabulate random recce raw means for common species

##recce.random.2007 <-  recce.random[recce.random$year== 2007 & recce.random$species %in% common.spp, ] ### only has the common species occuring in more than 10 plots in 2016

comp.2007 <- as.data.frame(comp[rownames(comp) %in% rownames(dca.site[dca.site$year==2007,]), ])
comp.2007.mean <- as.data.frame(apply(comp.2007, 2, mean)); comp.2007.mean$species <- rownames(comp.2007.mean); names(comp.2007.mean) <- c("mean.2007", "species")
comp.2007.se <- as.data.frame(apply(comp.2007, 2, s.e.)); comp.2007.se$species <- rownames(comp.2007.se); names(comp.2007.se) <- c("se.2007", "species")
comp.2007 <- merge(comp.2007.mean, comp.2007.se, by="species")
comp.2016 <- as.data.frame(comp[rownames(comp) %in% rownames(dca.site[dca.site$year==2016,]), ])
comp.2016.mean <- as.data.frame(apply(comp.2016, 2, mean)); comp.2016.mean$species <- rownames(comp.2016.mean); names(comp.2016.mean) <- c("mean.2016", "species")
comp.2016.se <- as.data.frame(apply(comp.2016, 2, s.e.)); comp.2016.se$species <- rownames(comp.2016.se); names(comp.2016.se) <- c("se.2016", "species")
comp.2016 <- merge(comp.2016.mean, comp.2016.se, by="species")
random.common.table <- merge(comp.2007, comp.2016, by="species"); random.common.table <- random.common.table[random.common.table$species %in% common.spp, ]
random.common.table <-  merge(random.common.table, dca.spp, by="species", all.x=TRUE, all.y=FALSE)
random.common.table[2:7] <- round(random.common.table[2:7], digits=2); 
random.common.table <-  merge(random.common.table, spp.list[, c("species", "family", "native", "form")], by="species", all.x=TRUE, all.y=FALSE)
random.common.table <- random.common.table[duplicated(random.common.table$species)==FALSE, ]
random.common.table[random.common.table$native!="Exotic","native"] <- "Native"; random.common.table[random.common.table$form %in%  c("Shrub", "Tree", "Vine","SubShrub", "Treefern"), "form"] <- "Woody"
###Add native and form classes to the table

random.common.table <- rbind( c("Native herbaceous plants", rep("", length(names(random.common.table))-1)),
random.common.table[random.common.table$native=="Native" & random.common.table$form=="Forb", c("species", "mean.2007", "se.2007", "mean.2016", "se.2016", "DCA1", "DCA2", "family")],
c("Exotic herbaceous plants", rep("", length(names(random.common.table))-1)),
random.common.table[random.common.table$native=="Exotic" & random.common.table$form=="Forb", c("species", "mean.2007", "se.2007", "mean.2016", "se.2016", "DCA1", "DCA2", "family")],
c("Native grasses, tussocks, sedges and rushes", rep("", length(names(random.common.table))-1)),
random.common.table[random.common.table$native=="Native" & random.common.table$form=="Graminoid", c("species", "mean.2007", "se.2007", "mean.2016", "se.2016", "DCA1", "DCA2", "family")],
c("Exotic grasses, sedges and rushes", rep("", length(names(random.common.table))-1)),
random.common.table[random.common.table$native=="Exotic" & random.common.table$form=="Graminoid", c("species", "mean.2007", "se.2007", "mean.2016", "se.2016", "DCA1", "DCA2", "family")],
c("Native woody plants", rep("", length(names(random.common.table))-1)),
random.common.table[random.common.table$native=="Native" & random.common.table$form=="Woody", c("species", "mean.2007", "se.2007", "mean.2016", "se.2016", "DCA1", "DCA2", "family")],
c("Exotic woody plants", rep("", length(names(random.common.table))-1)),
random.common.table[random.common.table$native=="Exotic" & random.common.table$form=="Woody", c("species", "mean.2007", "se.2007", "mean.2016", "se.2016", "DCA1", "DCA2", "family")]
)

save(random.recce.table.dcameans, random.common.table, text.colours, file="/home/sean/Documents/Molesworth/Data/Objects/RandomRecceTable.Rdata") ##random.recce.table.dcamodel,

####################################################################################################################################################################################################
################################################################         Native and form change 1952-2016  #########################################################################################

recce.dat[recce.dat$native %in% c("Indigenous Endemic", "Indigenous Non-Endemic"),"native"] <- "Native" ; recce.dat[recce.dat$form %in%  c("Shrub", "Tree", "Vine","SubShrub", "Treefern"), "form"] <- "Woody"
native.dat <- as.data.frame(table(paste(recce.dat$plot, recce.dat$year, sep="."), recce.dat$native)) ; names(native.dat) <- c("plot.yr", "native", "freq")
native.dat$year <-  as.numeric(substr(native.dat$plot.yr, regexpr("\\.", native.dat$plot.yr)+1, regexpr("\\.", native.dat$plot.yr)+ 5))
native.dat$plot <- as.character(substr(native.dat$plot.yr, 1, regexpr("\\.", native.dat$plot.yr)-1))
native.dat <- merge(native.dat[native.dat$native=="Native", c("plot.yr","plot", "year", "freq") ], native.dat[native.dat$native=="Exotic",c("plot.yr", "freq") ], by="plot.yr", all.x=TRUE, all.y=TRUE)
names(native.dat) <- c("plot.yr", "plot", "year", "native",  "exotic" )


## m1 <-    glm(cbind(native.dat$native, native.dat$exotic) ~ year, data=native.dat, family = binomial)  ##fenced fixed effect, site random effect. Need balanced model to use aov
## dev.off()
## plot(m1)
## sum(residuals(m1, type = "deviance")^2)/m1$df.residual 


## testDispersion(m1)
## simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
## plot(simulationOutput)
## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)

####GLMER fit is better

## m1  <- NULL
## m1 <- glmer(cbind(native.dat$native, native.dat$exotic) ~ year +(1|plot), data=native.dat, family = binomial)  ##

##  library(DHARMa)
## simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
## plot(simulationOutput)
## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)



## Over-dispersion in a GLM arises when a model restricts the variance of the response variable, and the data exhibits greater variance than the model restriction allows. This is highly relevant with count data where a Poisson error familiy is included in a glm model. As reviewer 2 points out there are diagnostic tests that can be used. We do not use count data. We use close to normally distributed data for gls modelling of DCA scores. We use a binomial model for comparing the proportion of native species, where overdispersion may occur.  Our models allow for the variance parameter to be fit to the data. In contrast to count data where a Poisson distibution, where variance is proportional to the mean, because as variance increases with an increasing mean. Where possible we have checked for overdispersion as suggested, and find it not relevant to the analysis we present. Of greater concern to us was model fit, independence of residuals and heterogeneity of variance (rather than overdispersion of variance). In other words, our model fit was not perfect. We consider our models adequately fit data. To reassure highly informed readers (which the reviewers clearly are), we have included a section on model diagnostics in supplimentary material. 


## m2 <- glm(native.dat$native ~ year, data=native.dat)  ##fenced fixed effect, site random effect. Need balanced model to use aov
## m3 <- glm(native.dat$exotic ~ year, data=native.dat)  ##fenced fixed effect, site random effect. Need balanced model to use aov

count.dat <- recce.dat; count.dat$cover <- 1  ##there are no nas for cover in recce.dat, so converting to 1 for cover allows a spp count per plot per class.  
count.dat$plot <-  as.character(count.dat$plot)
count.dat.1987 <- count.dat[count.dat$year==1987,] ##Have to summarise 1987 data seperately because only in this survey do some plots have exotics missing
native.dat.1987 <- as.data.frame.table(tapply(count.dat.1987$cover, list(count.dat.1987$year, count.dat.1987$plot, count.dat.1987$native), sum, na.rm=TRUE))
names(native.dat.1987) <- c("year", "plot", "native", "freq"); native.dat.1987[is.na(native.dat.1987$freq)==TRUE, "freq" ] <- 0

count.dat <- count.dat[count.dat$year!=1987,]
native.dat <- as.data.frame.table(tapply(count.dat$cover, list(count.dat$year, count.dat$plot, count.dat$native), sum, na.rm=TRUE))
names(native.dat) <- c("year", "plot", "native", "freq"); 
native.dat <- native.dat[is.na(native.dat$freq)==FALSE, ]  ##Remove plots with na which are not surveyed in some years.
native.dat <- rbind(native.dat, native.dat.1987); native.dat$year <- as.integer(as.character(native.dat$year))
native.dat <-  merge(native.dat[native.dat$native=="Native", c("year", "plot", "freq")], native.dat[native.dat$native=="Exotic", c("year", "plot", "freq")], by=c("year", "plot")) ##For binomial model
names(native.dat) <- c("year", "plot", "native", "exotic")
#m1 <- gls(freq ~ year + native, data=native.dat,)  ##


################################################################################################################################################################################################
        #######################################################  Spatial Prediction of Form and Natives  ################################################################
#####Native to Exotic Ratio
native.dat <- merge(native.dat, locs, by="plot", all.x=TRUE, all.y=TRUE)
native.dat <- native.dat[is.na(native.dat$native)==FALSE & is.na(native.dat$dem.altitude)==FALSE, c("plot", "native", "exotic", "year", "fenced", "x", "y", "dem.altitude", "terrace", "oversow", "grazing1987")]  ##
native.dat$oversow <-  as.character(native.dat$oversow) ; native.dat[is.na(native.dat$oversow)==FALSE, "oversow"] <- "yes"; native.dat[is.na(native.dat$oversow)==TRUE, "oversow"] <- "no"
native.dat$terrace <-  as.character(native.dat$terrace) ; native.dat[native.dat$terrace=="1", "terrace"] <- "Terrace"; native.dat[native.dat$terrace=="0", "terrace"] <- "Slope"
native.dat[native.dat$terrace=="1", "terrace"] <- "Terrace"; native.dat[native.dat$terrace=="0", "terrace"] <- "Slope"
native.dat[native.dat$year==2007 |  native.dat$year==2008, "x"] <- native.dat[native.dat$year==2007 | native.dat$year==2008, "x"] + 1.5
native.dat[native.dat$year==2007 | native.dat$year==2008, "y"] <- native.dat[native.dat$year==2007 | native.dat$year==2008, "y"] + 1.5
native.dat[duplicated(native.dat$x)==TRUE, "x"] <-  native.dat[duplicated(native.dat$x)==TRUE, "x"] + 1.5  ##Make changes for things duplicated once
native.dat[duplicated(native.dat$y)==TRUE, "y"] <-  native.dat[duplicated(native.dat$y)==TRUE, "y"] + 1.5  
native.dat[duplicated(native.dat$x)==TRUE, "x"] <-  native.dat[duplicated(native.dat$x)==TRUE, "x"] - 5  ##Make changes for things duplicated twice
native.dat[duplicated(native.dat$y)==TRUE, "y"] <-  native.dat[duplicated(native.dat$y)==TRUE, "y"] - 5  
native.dat[duplicated(native.dat$x)==TRUE, "x"] <-  native.dat[duplicated(native.dat$x)==TRUE, "x"] - 2.5  ##Make changes for things duplicated thrice
native.dat[duplicated(native.dat$y)==TRUE, "y"] <-  native.dat[duplicated(native.dat$y)==TRUE, "y"] - 2.5  

native.dat$terrace <-  as.character(native.dat$terrace) ;
native.dat[native.dat$terrace == "Slope", "terrace"]  <- 0; native.dat[native.dat$terrace == "Terrace", "terrace"]  <- 1;
native.dat$terrace  <- as.integer(as.character(native.dat$terrace))
native.dat$fenced <-  as.character(native.dat$fenced) ;
native.dat[native.dat$fenced == "Unfenced", "fenced"]  <- 0; native.dat[native.dat$fenced == "Fenced", "fenced"]  <- 1;
native.dat$fenced  <- as.integer(as.character(native.dat$fenced))
native.dat$oversow <-  as.character(native.dat$oversow) ;
native.dat[native.dat$oversow == "no", "oversow"]  <- 0; native.dat[native.dat$oversow == "yes", "oversow"]  <- 1;
native.dat$oversow  <- as.integer(as.character(native.dat$oversow))

## ####Save instead of run to save time
## m.1 <- glmmTMB(cbind(native, exotic) ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 +  (1|plot),
##                 data=native.dat, family=binomial) #exp(pos + 0 | group)
## m.null <- glmmTMB(cbind(native, exotic) ~ 1 +  (1|plot), data=native.dat, family = binomial)
## save(m.1, m.null, file="/home/sean/Documents/Molesworth/Data/Objects/glmmTMBNative.Rdata")
load("/home/sean/Documents/Molesworth/Data/glmmTMBNative.Rdata")





## Error in .subset2(x, i, exact = exact) : no such index at level 1
## m.1.step <- buildglmmTMB(cbind(native, exotic) ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced))^2 +(1|plot),
##                 data=native.dat, family=binomial) #exp(pos + 0 | group)


recce.table.native  <- glmmTMB.latex(m.1)

recce.native.R2  <- MuMIn::r.squaredGLMM(m.1, m.null)



## #############Statement in results:
## #### There has also been a commensurate increase in the proportion of woody species present in plots of $\approx$7\% from in 1952 to

## new.hfi.yrs  <- expand.grid(year = c(1952, 1989, 2016), x = mean(hfi$x), y = mean(hfi$y), dem.altitude = mean(hfi$dem.altitude), fenced=c(0,1),
##                             terrace=c(0,1), oversow=c(0,1),  plot=NA)

## new.hfi.yrs$rich  <-  predict(m.1, new.hfi.yrs, type="response", re.form=NA, se.fit = FALSE)
## native.rich.yr  <- as.data.frame.table(tapply(new.hfi.yrs$rich, list(new.hfi.yrs$year, new.hfi.yrs$terrace), mean))
## names(native.rich.yr)  <- c("year", "terrace", "proportion.native") 



## ## simulationOutput <- simulateResiduals(fittedModel = m.1, plot = F)
## ## plot(simulationOutput)
## ## testDispersion(simulationOutput)
## ## testZeroInflation(simulationOutput)
## ## testSpatialAutocorrelation(simulationOutput = simulationOutput, x = native.dat$x, y= native.dat$y)






## ##  m1 <- glm(cbind(native, exotic)  ~ year + x + y + dem.altitude + terrace + oversow + year + fenced , data=native.dat, family=binomial)
## ## ##m2 <- glm(freq  ~ native  + year + x + y + dem.altitude + terrace , data=native.dat, family=poisson)
## ## stepAIC(m1) ##keep em all

## ##Use library DHARMa
## m1 <- glmer(cbind(native, exotic)  ~ scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  (1|plot), data=native.dat, family=binomial)
## m1.tmb <- glmmTMB(cbind(native, exotic)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2 + terrace + oversow +  fenced + (1|plot),
##                   data=native.dat, family=binomial)








## m1.native <- lmer(native  ~ scale(year) + scale(x) + scale(y) + scale(dem.altitude) + terrace + oversow +  fenced +  (1|plot), data=native.dat)
## m1.native <- lm(native  ~ scale(year) + scale(x) + scale(y) + scale(dem.altitude) + terrace + oversow +  fenced , data=native.dat)

## simulationOutput <- simulateResiduals(fittedModel = m1.native, plot = F)
## plot(simulationOutput)
## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)
## testSpatialAutocorrelation(simulationOutput = simulationOutput, x = native.dat$x, y= native.dat$y)

## m1.exotic <- glmer(exotic  ~ scale(year) + scale(x) + scale(y) + scale(dem.altitude) + terrace + oversow +  fenced +  (1|plot), data=native.dat, family=poisson)
## simulationOutput <- simulateResiduals(fittedModel = m1.exotic, plot = F)
## plot(simulationOutput)
## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)
## testSpatialAutocorrelation(simulationOutput = simulationOutput, x = native.dat$x, y= native.dat$y)

## ## m3.tmb.native <- glmmTMB(native  ~ scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  (1|plot), data=native.dat, family=poisson) #
## m1.tmb.native <- glmer(native  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2 +  (1|plot), data=native.dat, family=poisson) ## 
## ## m4.tmb.native <- lmer(native  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2 +  (1|plot), data=native.dat) ## 
## ## m2.tmb.native <- glmmTMB(native  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2 + terrace + oversow +  fenced +  (1|plot), data=native.dat, family=poisson)

## ##anova(m1.tmb.native, m2.tmb.native, m3.tmb.native, m4.tmb.native)


## ##Use library(AER); dispersiontest(m1.tmb.native)

## simulationOutput <- simulateResiduals(fittedModel = m1.tmb.native, plot = F)
## testOverdispersion(simulationOutput)
## plot(simulationOutput)

## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)
## testSpatialAutocorrelation(simulationOutput = simulationOutput, x = native.dat$x, y= native.dat$y)


## m1.tmb.exotic <- glmmTMB(exotic  ~ scale(year) * scale(x) * scale(y) * scale(dem.altitude) +  (1|plot), data=native.dat, family=poisson)
## m2.tmb.native <- glmmTMB(exotic  ~ scale(year) + scale(x) + scale(y) + scale(dem.altitude) + terrace + oversow +  fenced +  (1|plot), data=native.dat, family=poisson)

## simulationOutput <- simulateResiduals(fittedModel = m1.tmb.exotic, plot = F)
## plot(simulationOutput)
## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)
## testSpatialAutocorrelation(simulationOutput = simulationOutput, x = native.dat$x, y= native.dat$y)




## ## m1.cor <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow,  data=dat)
## ## m2 <- gls(DCA1  ~ year + x + y + dem.altitude + terrace + oversow + grazing1987,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## ## m3 <- gls(DCA1  ~ year + fenced + x + y + dem.altitude + terrace + oversow + grazing1987,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## ## anova(m1, m1.cor, m2, m3)  ## M1 the winner

## ##m2 <- gls(DCA2  ~ year + x + y + dem.altitude + terrace + oversow,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## ## m1.cor <- gls(DCA2  ~ year + x + y + dem.altitude + terrace + oversow,  data=dat)
## ## m2 <- gls(DCA2  ~ year + x + y + dem.altitude + terrace + oversow + grazing1987,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## ## m3 <- gls(DCA2  ~ year + fenced + x + y + dem.altitude + terrace + oversow + grazing1987,  correlation = corExp(form = ~ x + y, nugget = TRUE), data=dat)
## ## anova(m1, m1.cor, m2, m3)  ## M1.2 the winner

## native.table <- as.data.frame(round(coef(summary(m1))[,c(1:2,4)], digits=3)); 
## native.table[,names(native.table)] <- lapply(names(native.table), function(x) as.character(native.table[,x]))
## native.table[native.table=="0"] <-  "$<$0.001"


## native.dat$ratio <- (native.dat$native/ (native.dat$native + native.dat$exotic)); native.dat[native.dat$ratio=="Inf", "ratio"] <- native.dat[native.dat$ratio=="Inf", "native"] ###make ratio value for plots with no exotics the number of natives
## native.table.means <- as.data.frame(rbind( tapply(native.dat$ratio, native.dat$year, mean), tapply(native.dat$ratio, native.dat$year, s.e.)))

## ##native.predicted.means <- rbind(predict.native.1952, predict.native.2016); names(native.predicted.means) <- c("ratio.mean", "ratio.se", "year")

## native.means <- as.data.frame(cbind(names(native.table.means), round(t(native.table.means), digits=3))); names(native.means) <- c("year", "raw.mean", "raw.se")
## ##native.means[,names(native.means)] <- lapply(names(native.means), function(x) as.character(native.means[,x]))
## ##native.means[native.means=="0"] <-  "$<$0.001"

## save(native.means, native.table, file="/home/sean/Documents/Molesworth2020/Analysis/Objects/ProportionNativeTable.Rdata") ###


############### Spatial prediction for modelled ratio of native to exotic

grid.fenced[grid.fenced$area.fenced!= "Unfenced", "area.fenced"] <- "Fenced"

origin.box <- c(1635000, 1635000, 1640000, 1640000, 5302000, 5297000, 5297000, 5302000)
increase.box <- c(0,0,0,0,5000,5000,5000,5000)
boxes <- as.data.frame(rbind(origin.box, origin.box + increase.box, origin.box + increase.box*2, origin.box + increase.box*3, origin.box + increase.box*4))
names(boxes) <-  c("x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4")
rownames(boxes) <- seq(1:5) ; boxes$colours <- heat.colors(5); boxes$greys <- c("grey90", "grey80", "grey70", "grey60", "grey50")

pdf("/home/sean/Documents/Molesworth/Graphs/NativeSpeciesMap.pdf", width=7, height = 3.5)
par(mfrow = c(1, 2), cex = 0.8, mar = c(0, 0, 0, 0), oma = c(0.5, 2.5, 2.5, 0.5))
for (i in c(1952, 2016)){  ## ##   i <-  1952 sort(unique(native.dat$year))

grid.predict <-  as.data.frame(cbind(grid.alt, grid.terrace[1], grid.oversow[1], grid.fenced[1], rep(i, length(grid.alt[,1]))));
names(grid.predict) <-  c("dem.altitude", "x", "y", "terrace", "oversow", "fenced", "year")

grid.predict[grid.predict$terrace == "Slope", "terrace"]  <- 0; grid.predict[grid.predict$terrace == "Terrace", "terrace"]  <- 1; grid.predict$terrace  <- as.integer(as.character(grid.predict$terrace))
grid.predict[grid.predict$fenced == "Unfenced", "fenced"]  <- 0; grid.predict[grid.predict$fenced == "Fenced", "fenced"]  <- 1; grid.predict$fenced  <- as.integer(as.character(grid.predict$fenced))
grid.predict[grid.predict$oversow == "no", "oversow"]  <- 0; grid.predict[grid.predict$oversow == "yes", "oversow"]  <- 1; grid.predict$oversow  <- as.integer(as.character(grid.predict$oversow))
grid.predict$plot  <- NA
grid.predict <- grid.predict[is.na(grid.predict$dem.altitude)==FALSE,]

predict.m1 <- predict(m.1, grid.predict, type="response", re.form=NA, se.fit = FALSE);
predict.native <- cbind(grid.predict[, c("x", "y")], predict.m1)
native.r <-  rasterFromXYZ(predict.native, crs = "+init=epsg:2193")
image(native.r, asp=1, xlim=c(min(1570000), max(1640000)), box=FALSE, axes=FALSE, labels= FALSE, frame.plot=FALSE) 
#plot(dat.locs[dat.locs$year==i,], add=TRUE, pch=19, cex=0.5)

text(1585000, 5350000, i, cex=1.5)##, col.main="black")#,  sub="Predictive model from recce data 1952-2016", col.sub="grey", cex.lab=0.5)
mean.1 <- (round(mean(predict.m1), digits=2))*100; se.1 <- round(s.e.(predict.m1)*100, digits=2)
text(1587000, 5345000, eval(substitute(paste("Native = ", mean.1, "\u00B1", se.1 , "%", sep=""))))
## for (j in boxes$colours){ polygon(boxes[boxes$colours== j, c(1:4)], boxes[boxes$colours== j, c(5:8)], col=j) } 
## text(c(rep(1630000,3)), c(5300000, 5310000, 5320000), c("100% exotic", "50%-50%", "100% native"))
##text(1590000, 5353000, "Common species DCA ordination", col="black", cex=0.75) 
##assign(paste("predict.native.", i, sep=""), c(mean(predict.m1), s.e.(predict.m1), i))
##mtext("Proportion of native species", outer = TRUE)
#legend.gradient(c(1620000, 1620000, 1625000, 1625000), c(5315000, 5295000, 5295000, 5315000), cols = heat.colors, limits = c(0, 1),  title = "Legend")
}
plot(oversow.polygons, add=TRUE) ##plot(fenced.polygons, add=TRUE)
for (j in boxes$colours){ polygon(boxes[boxes$colours== j, c(1:4)], boxes[boxes$colours== j, c(5:8)], col=j) } 
 text(c(rep(1634000,3)), c(5300000, 5310000, 5320000), c("100% \n exotic", "50%-50%", "100% \n native"), cex=.95, adj=1)

dev.off()



#####Woody to herbaceous ratio
recce.dat$woody <-  "Herbaceous";
recce.dat[recce.dat$form== "Woody" & is.na(recce.dat$form)==FALSE, "woody"] <- "Woody"
form.dat <- as.data.frame(table(paste(recce.dat$plot, recce.dat$year, sep="."), recce.dat$woody)); names(form.dat) <- c("plot.yr", "woody", "freq")
form.dat$year <-  as.numeric(substr(form.dat$plot.yr, regexpr("\\.", form.dat$plot.yr)+1, regexpr("\\.", form.dat$plot.yr)+ 5))
form.dat$plot <- as.character(substr(form.dat$plot.yr, 1, regexpr("\\.", form.dat$plot.yr)-1))
form.dat <- merge(form.dat[form.dat$woody=="Herbaceous", c("plot.yr", "year",  "plot", "freq")], form.dat[form.dat$woody=="Woody", c("plot.yr", "freq") ], by="plot.yr", all.x=TRUE, all.y=TRUE)
names(form.dat) <- c("plot.yr", "year", "plot", "herb", "woody")

##m1 <- glm(cbind(form.dat$woody, form.dat$herb) ~ year, data=form.dat, family = binomial)  ##fenced fixed effect, site random effect. Need balanced model to use aov
form.dat <- merge(form.dat, locs, by="plot", all.x=TRUE, all.y=TRUE)
form.dat <- form.dat[is.na(form.dat$herb)==FALSE & is.na(form.dat$dem.altitude)==FALSE, c("plot", "herb", "woody", "year", "fenced", "x", "y", "dem.altitude", "terrace", "oversow", "grazing1987")]  ##
form.dat$oversow <-  as.character(form.dat$oversow) ; form.dat[is.na(form.dat$oversow)==FALSE, "oversow"] <- "yes"; form.dat[is.na(form.dat$oversow)==TRUE, "oversow"] <- "no"
form.dat$terrace <-  as.character(form.dat$terrace) ; form.dat[form.dat$terrace=="1", "terrace"] <- "Terrace"; form.dat[form.dat$terrace=="0", "terrace"] <- "Slope"
form.dat[form.dat$terrace=="1", "terrace"] <- "Terrace"; form.dat[form.dat$terrace=="0", "terrace"] <- "Slope"

form.dat$terrace <-  as.character(form.dat$terrace);
form.dat[form.dat$terrace == "Slope", "terrace"]  <- 0; form.dat[form.dat$terrace == "Terrace", "terrace"]  <- 1;
form.dat$terrace  <- as.integer(as.character(form.dat$terrace))
form.dat$fenced <-  as.character(form.dat$fenced) ;
form.dat[form.dat$fenced == "Unfenced", "fenced"]  <- 0; form.dat[form.dat$fenced == "Fenced", "fenced"]  <- 1;
form.dat$fenced  <- as.integer(as.character(form.dat$fenced))
form.dat$oversow <-  as.character(form.dat$oversow) ;
form.dat[form.dat$oversow == "no", "oversow"]  <- 0; form.dat[form.dat$oversow == "yes", "oversow"]  <- 1;
form.dat$oversow  <- as.integer(as.character(form.dat$oversow))

## ####Save instead of run to save time
## m.1 <- glmmTMB(cbind(woody, herb) ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 +(1|plot),
##                 data=form.dat, family=binomial) #exp(pos + 0 | group)
## m.null <- glmmTMB(cbind(woody, herb) ~ 1 +  (1|plot), data=form.dat, family = binomial)
## save(m.1, m.null, file="/home/sean/Documents/Molesworth2020/Analysis/Objects/glmmTMBWoody.Rdata")
load("/home/sean/Documents/Molesworth/Data/Objects/glmmTMBWoody.Rdata")
recce.table.woody  <- glmmTMB.latex(m.1)
recce.woody.R2  <- MuMIn::r.squaredGLMM(m.1, m.null)


## Error in .subset2(x, i, exact = exact) : no such index at level 1

##  m.1.step <- buildglmmTMB(cbind(woody, herb) ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced))^2 +(1|plot),
##                 data=form.dat, family=binomial) #




#############Statement in results:
#### There has also been a commensurate increase in the proportion of woody species present in plots of $\approx$4\% from in 1952 to $\approx$15\% in 15% in 2016 

new.form.yrs  <- expand.grid(year = c(1952, 1989, 2016), x = mean(form.dat$x), y = mean(form.dat$y), dem.altitude = mean(form.dat$dem.altitude), fenced=c(0,1),
                            terrace=c(0,1), oversow=c(0,1),  plot=NA)

new.form.yrs$rich  <-  predict(m.1, new.form.yrs, type="response", re.form=NA, se.fit = FALSE)

fenced.slope.oversown  <- new.form.yrs[new.form.yrs$fenced==1 & new.form.yrs$terrace==0 & new.form.yrs$oversow==0, ]

not.fenced.slope.oversown  <- new.form.yrs[new.form.yrs$fenced==0 & new.form.yrs$terrace==1 & new.form.yrs$oversow==1, ]

save(text.stats, text.colours, fenced.slope.oversown, not.fenced.slope.oversown,  file="/home/sean/Documents/Molesworth/Data/Objects/TextStats.Rdata")



## new

## woody.rich.yr  <- as.data.frame.table(tapply(new.hfi.yrs$rich, new.hfi.yrs$year, mean))
## names(woody.rich.yr)  <- c("year", "proportion.woody") 




## simulationOutput <- simulateResiduals(fittedModel = m.1, plot = F)
## plot(simulationOutput)
## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)
## testSpatialAutocorrelation(simulationOutput = simulationOutput, x = native.dat$x, y= native.dat$y)


## m1 <- glm(cbind(woody, herb)  ~ year + x + y + dem.altitude + terrace + oversow  , data=form.dat, family=binomial)
## ##m2 <- glm(freq  ~ native  + year + x + y + dem.altitude + terrace , data=native.dat, family=poisson)
## ##stepAIC(m1) ## remove fenced

## form.table <- as.data.frame(round(coef(summary(m1))[,c(1:2,4)], digits=3)); 
## form.table[,names(form.table)] <- lapply(names(form.table), function(x) as.character(form.table[,x]))
## form.table[form.table=="0"] <-  "$<$0.001"

## form.dat$ratio <- (form.dat$woody/ (form.dat$woody + form.dat$herb))
## form.table.means <- as.data.frame(rbind( tapply(form.dat$ratio, form.dat$year, mean), tapply(form.dat$ratio, form.dat$year, s.e.)))

###############RGB spatial prediction for modeled ratio of native to exotic

greys <- c("grey100", "grey90","grey80", "grey70", "grey60", "grey50", "grey40", "grey30")  ##greys <- grey.colors(5,.55,.95)

pdf("/home/sean/Documents/Molesworth/Data/Graphs/ShrubMap.pdf")
par(mfrow = c(2, 1), cex = 0.6, mar = c(0, 0, 0, 0), oma = c(0.5, 2.5, 2.5, 0.5))
for (i in c(1952, 2016)){  ## ##    sort(unique(form.dat$year)) i <-  1952

grid.predict <-  as.data.frame(cbind(grid.alt, grid.terrace[1], grid.oversow[1], grid.fenced[1], rep(i, length(grid.alt[,1]))));
names(grid.predict) <-  c("dem.altitude", "x", "y", "terrace", "oversow", "fenced", "year")

grid.predict[grid.predict$terrace == "Slope", "terrace"]  <- 0; grid.predict[grid.predict$terrace == "Terrace", "terrace"]  <- 1; grid.predict$terrace  <- as.integer(as.character(grid.predict$terrace))
grid.predict[grid.predict$fenced == "Unfenced", "fenced"]  <- 0; grid.predict[grid.predict$fenced == "Fenced", "fenced"]  <- 1; grid.predict$fenced  <- as.integer(as.character(grid.predict$fenced))
grid.predict[grid.predict$oversow == "no", "oversow"]  <- 0; grid.predict[grid.predict$oversow == "yes", "oversow"]  <- 1; grid.predict$oversow  <- as.integer(as.character(grid.predict$oversow))
grid.predict$plot  <- NA
grid.predict <- grid.predict[is.na(grid.predict$dem.altitude)==FALSE,]


## predict.native <- cbind(grid.predict[, c("x", "y")], predict.m1)
## grid.predict <-  as.data.frame(cbind(grid.alt, grid.terrace[1], grid.oversow[1], rep(i, length(grid.alt[,1])))); names(grid.predict) <-  c("dem.altitude", "x", "y", "terrace", "oversow", "year")
## grid.predict <- grid.predict[is.na(grid.predict$dem.altitude)==FALSE,]
predict.m1 <- predict(m.1, grid.predict, type="response", re.form=NA, se.fit = FALSE);
predict.form <- cbind(grid.predict[, c("x", "y")], predict.m1)
form.r <-  rasterFromXYZ(predict.form, crs = "+init=epsg:2193")
image(form.r, col=greys,  asp=1, axes=FALSE,  frame.plot=FALSE)
##image(form.r, col=greys,  breaks=seq(.1,.5,.05), asp=1, box=FALSE, axes=FALSE, labels= FALSE, frame.plot=FALSE)##breaks=c(.1,.2,.3,.4,.5,.75))##, col =  terrain.colors[6]) ##rainbow, heat.colors,, terrain.colors )
#plot(dat.locs[dat.locs$year==i,], add=TRUE, pch=19, cex=0.5)
plot(oversow.polygons, add=TRUE, lwd=0.25)
plot(boundary, add=TRUE, lwd=0.25)
text(1592000, 5355000, i)##, col.main="black")#,  sub="Predictive model from recce data 1952-2016", col.sub="grey", cex.lab=0.5)
text(1590000, 5345000, paste("Woody =", (round(mean(predict.m1), digits=2))*100, "% +", round(s.e.(predict.m1)*100, digits=2), "%"))
for (j in boxes$colours){ polygon(boxes[boxes$colours== j, c(1:4)], boxes[boxes$colours== j, c(5:8)], col=boxes[boxes$colours== j, "greys"])} 
text(c(rep(1630000,3)), c(5300000, 5310000, 5320000), c(" <20% woody", " <30% woody", " >50% woody"))

##text(1590000, 5353000, "Common species DCA ordination", col="black", cex=0.75) 
##assign(paste("predict.woody.", i, sep=""), c(mean(predict.m1), s.e.(predict.m1), i))
}
mtext("Proportion of woody species - relev\'es", outer = TRUE, cex=0.6, vfont=c("serif", "plain"))
dev.off()



## #woody.predicted.means <- rbind(predict.woody.1952, predict.woody.1987, predict.woody.2007, predict.woody.2008, predict.woody.2012, predict.woody.2016); names(woody.predicted.means) <- c("ratio.mean", "ratio.se", "year")
## woody.means <- as.data.frame(cbind(names(form.table.means), round(t(form.table.means), digits=2))); names(woody.means) <- c("year", "raw.mean", "raw.se")
## ## woody.means[,names(woody.means)] <- lapply(names(woody.means), function(x) as.character(woody.means[,x]))
## ## woody.means[woody.means=="0"] <-  "$<$0.001"
## save(woody.means, form.table, file="/home/sean/Documents/Molesworth2020/Analysis/Objects/ProportionWoodyTable.Rdata") ###

### writeRaster(terrace, "/home/sean/Documents/Molesworth2020/Analysis/Graphs/TerraceRaster.tif", overwrite=TRUE)
writeGDAL(terrace, "/home/sean/Documents/Molesworth/Graphs/TerraceRaster.tif", drivername = "GTiff")


pdf("/home/sean/Documents/Molesworth/Graphs/TerraceMap.pdf")
par(cex = 0.6, mar = c(0, 0, 0, 0), oma = c(0.5, 2.5, 2.5, 0.5))
image(terrace, col=c("grey90", "grey40"))
plot(dat.locs, add=TRUE, pch=19, cex=0.25)
mtext(expression("Terraces (<400 m from watercourse and <20" * degree *" slope)"), outer = TRUE)
dev.off()





  #######################################  #######################################  #######################################  #######################################
  #### Compare form

dat  <- recce.dat[recce.dat$native %in% c("Native",  "Exotic"), ]
dat  <- as.data.frame.table(tapply(dat$cover, list(dat$plot, dat$year, dat$native, dat$woody), length))
names(dat)  <-  c("plot", "year", "native", "woody", "n.spp")
dat  <- dat[is.na(dat$n.spp) == FALSE, ]

dat$form  <-  NA
dat[dat$native ==  "Native" & dat$woody=="Herbaceous",  "form"] <-  "Native herbaceous"
dat[dat$native == "Native" & dat$woody=="Woody",  "form"]     <-  "Native woody"
dat[dat$native == "Exotic" & dat$woody=="Herbaceous",  "form"]     <-  "Exotic herbaceous"
dat[dat$native == "Exotic" & dat$woody=="Woody", "form"]     <-  "Exotic woody"
dat$hex.col  <-  NA

## nat.herb.col  <- "#3aeb13"
## nat.wood.col  <- "#d1d487"
## exot.herb.col  <- "#cf0e0e"
## exot.wood.col  <- "#8c3424"

dat[dat$form == "Native herbaceous",  "hex.col"]     <-  nat.herb.col 
dat[dat$form == "Native woody",  "hex.col"]     <-  nat.woody.col 
dat[dat$form == "Exotic herbaceous",  "hex.col"]     <-  exot.herb.col
dat[dat$form == "Exotic woody",  "hex.col"]     <- exot.woody.col   #### "#0d0000"


dat  <-  merge(dat, native.dat[ ,c("plot", "year", "fenced", "x", "y", "dem.altitude", "terrace", "oversow")],
                     by=c("plot", "year"), all.x=FALSE, all.y=FALSE)

dat$native  <-  as.character(dat$native); dat$woody  <-  as.character(dat$woody)
dat[dat$native == "Exotic", "native"]  <-  0; dat[dat$native == "Native", "native"]  <-  1 ##Recode to 0 and 1 for exotic and native
dat[dat$woody == "Herbaceous", "woody"]  <-  0; dat[dat$woody == "Woody", "woody"]  <-  1 ##Recode to 0 and 1 for herb and native
dat$native  <-  as.integer(dat$native) ; dat$woody  <-  as.integer(dat$woody)
dat$year  <-  as.integer(as.character(dat$year)); dat$plot  <-  as.character(dat$plot)

## m.spp.1 <- glmmTMB(n.spp  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  as.factor(terrace) + as.factor(oversow) + as.factor(fenced))^2 +
##                        (1|plot), data=dat, family=poisson)
###m.spp.2 <- glmmTMB(n.spp  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  as.factor(terrace) + as.factor(oversow) + as.factor(fenced) + as.factor(woody) + as.factor(native))^2 + (1|plot),  data=dat, family=gaussian)  #####Didnt converge

m.spp <- glmmTMB(n.spp  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  as.factor(terrace) + as.factor(oversow) + as.factor(fenced) + as.factor(woody) + as.factor(native))^2 +   (1|plot), data=dat, family=nbinom1)  #### Best model

## Convergence failure. Reducing terms and retrying... The failure was:
## m.spp.step <- buildglmmTMB(n.spp  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced) + scale(woody) + scale(native))^2 + (1|plot),  ziformula = ~0, data=dat, family=truncated_nbinom1)


###anova(m.spp.1,  m.spp.2, m.spp)


## simulationOutput <- simulateResiduals(fittedModel = m.spp, plot = F)
## plot(simulationOutput)
## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)

m.null <- glmmTMB(n.spp ~ 1 +  (1|plot), data=dat, family = nbinom1)

recce.n.spp.R2  <- MuMIn::r.squaredGLMM(m.spp, m.null) ###Cant used truncated negative binomial for R2 calc
glmm.model.recce.n.spp  <- glmmTMB.latex(m.spp)

new.dat  <- expand.grid(year = sort(unique(dat$year)), x = mean(dat$x), y = mean(dat$y), oversow=0, fenced=0, terrace=0,  dem.altitude = mean(dat$dem.altitude),
                            woody = c(0,1),  native = c(0,1), plot=NA)
new.dat$form  <-  NA
new.dat[new.dat$native == 1 & new.dat$woody==0,  "form"]     <-  "Native herbaceous"
new.dat[new.dat$native == 1 & new.dat$woody==1,  "form"]     <-  "Native woody"
new.dat[new.dat$native == 0 & new.dat$woody==0,  "form"]     <-  "Exotic herbaceous"
new.dat[new.dat$native == 0 & new.dat$woody==1,  "form"]     <-  "Exotic woody"

new.dat$hex.col  <-  NA
new.dat[new.dat$form == "Native herbaceous",  "hex.col"]     <-  nat.herb.col ##Light green 
new.dat[new.dat$form == "Native woody",  "hex.col"]     <-  nat.woody.col ##"#d1d487"  ##"#165e06" ### browny yellow
new.dat[new.dat$form == "Exotic herbaceous",  "hex.col"]     <- exot.herb.col##  "#cf0e0e"  ##red
new.dat[new.dat$form == "Exotic woody",  "hex.col"]     <- exot.woody.col  ###"#8c3424" ## "#0d0000"  


new.dat$rich  <-  predict(m.spp, new.dat, type="response", re.form=NA, se.fit = FALSE)
new.dat$rich.se  <-  predict(m.spp, new.dat, type="response", re.form=NA, se.fit = TRUE)$se.fit

## total.sum  <- as.data.frame.table(tapply(dat.total$hfi, dat.total$year, mean))
## names(total.sum)  <-  c("year", "hfi"); total.sum$year  <- as.integer(as.character(total.sum$year)); total.sum$hfi  <- as.integer(as.character(total.sum$hfi))



##################Include HFI in estimates

hfi.spp <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/HFISpecies.csv", sep = ",", header=TRUE, skip = 0, stringsAsFactors=FALSE)  
hfi.spp.add <- as.data.frame(cbind("B2", "RUMACE", 2016, (hfi.spp[hfi.spp$plot=="B2" & hfi.spp$species=="Rumex acetosella", -c(1:3)]))); names(hfi.spp.add) <- names(hfi.spp) ##add a rum ace so that all plots have at least one spp
hfi.spp <-  rbind(hfi.spp, hfi.spp.add)
hfi.spp$woody  <- "Herbaceous"; hfi.spp[hfi.spp$form %in% c("Shrub", "SubShrub", "Vine", "Tree"), "woody"]  <- "Woody"
hfi.spp$oversow  <- 0; hfi.spp[hfi.spp$oversow.year %in% c("1958A", "1957A", "1963A", "1964B", "1963B", "1953A"), "oversow"]  <- 1
hfi.spp$fenced  <- 1; hfi.spp[is.na(hfi.spp$fenced.year)==TRUE, "fenced"]  <- 0 
hfi.dat  <- hfi.spp[, c("plot", "species", "year", "hfi", "native", "woody", "x", "y", "dem.altitude", "oversow", "fenced", "terrace", "fenced.year")]
hfi.dat  <-  hfi.dat[hfi.dat$plot %nin% c(unique(hfi.dat$plot)[grep("Molesworth", unique(hfi.dat$plot))], unique(hfi.dat$plot)[grep("Saxton", unique(hfi.dat$plot))]), ]

hfi  <- as.data.frame.table(tapply(hfi.dat$hfi, list(hfi.dat$plot, hfi.dat$year, hfi.dat$native, hfi.dat$woody), length))
names(hfi)  <-  c("plot", "year", "native", "woody", "n.spp")
hfi  <- hfi[is.na(hfi$n.spp) == FALSE, ]



hfi$form  <-  NA
hfi[hfi$native ==  "Native" & hfi$woody=="Herbaceous",  "form"] <-  "Native herbaceous"
hfi[hfi$native == "Native" & hfi$woody=="Woody",  "form"]     <-  "Native woody"
hfi[hfi$native == "Exotic" & hfi$woody=="Herbaceous",  "form"]     <-  "Exotic herbaceous"
hfi[hfi$native == "Exotic" & hfi$woody=="Woody", "form"]     <-  "Exotic woody"
hfi$hex.col  <-  NA
hfi[hfi$form == "Native herbaceous",  "hex.col"]     <- nat.herb.col ##  "#3aeb13"
hfi[hfi$form == "Native woody",  "hex.col"]     <-  nat.woody.col ##"d1d487"
hfi[hfi$form == "Exotic herbaceous",  "hex.col"]     <-  exot.herb.col ###  "#cf0e0e"
hfi[hfi$form == "Exotic woody",  "hex.col"]     <- exot.woody.col ##  "#8c3424"



native.dat[native.dat$terrace == "Slope", "terrace"]  <- 0; native.dat[native.dat$terrace == "Terrace", "terrace"]  <- 1;
native.dat$terrace  <- as.integer(as.character(native.dat$terrace))




hfi  <-  merge(hfi, unique(hfi.dat[ ,c("plot", "year", "fenced", "x", "y", "dem.altitude", "terrace", "oversow")]),
                     by=c("plot", "year"), all.x=TRUE, all.y=FALSE)

hfi$native  <-  as.character(hfi$native); hfi$woody  <-  as.character(hfi$woody)
hfi[hfi$native == "Exotic", "native"]  <-  0; hfi[hfi$native == "Native", "native"]  <-  1 ##Recode to 0 and 1 for exotic and native
hfi[hfi$woody == "Herbaceous", "woody"]  <-  0; hfi[hfi$woody == "Woody", "woody"]  <-  1 ##Recode to 0 and 1 for herb and native
hfi$native  <-  as.integer(hfi$native) ; hfi$woody  <-  as.integer(hfi$woody)
hfi$year  <-  as.integer(as.character(hfi$year)); hfi$plot  <-  as.character(hfi$plot)

hfi[hfi$terrace == "Slope", "terrace"]  <- 0; hfi[hfi$terrace == "Terrace", "terrace"]  <- 1;  hfi$terrace  <- as.integer(as.character(hfi$terrace))


m.hfi <- glmmTMB(n.spp  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  as.factor(terrace) + as.factor(oversow) + as.factor(fenced) + as.factor(woody) + as.factor(native))^2 +
                     (1|plot),  data=hfi[hfi$year > 2005, ], family=truncated_nbinom1) #family=poisson, data=hfi) +  terrace + oversow + fenced +
m.hfi.1 <- glmmTMB(n.spp  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) +  as.factor(terrace) + as.factor(oversow) + as.factor(fenced) + as.factor(woody) + as.factor(native))^2 +
                       (1|plot),  ziformula = ~0, data=hfi[hfi$year > 2005, ], family=poisson)
##m.hfi.2 <- glmmTMB(n.spp  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced) + scale(woody) + scale(native))^2 +  (1|plot),  data=hfi[hfi$year > 2005, ], family=truncated_poisson) ##ziformula = ~1,


## m.hfi.step <- buildglmmTMB(glmmTMB(n.spp  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced) + scale(woody) + scale(native))^2 +  (1|plot),  data=hfi[hfi$year > 2005, ], family=truncated_nbinom1)) 



anova(m.hfi, m.hfi.1)##, m.hfi.2)

simulationOutput <- simulateResiduals(fittedModel = m.hfi, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

m.null <- glmmTMB(n.spp ~ 1 +  (1|plot), data=hfi[hfi$year > 2005, ], family = poisson)

hfi.n.spp.R2  <- MuMIn::r.squaredGLMM(m.hfi.1, m.null)

glmm.model.hfi.n.spp  <-  glmmTMB.latex(m.hfi)


save(glmm.model.DCA1, recce.dca1.R2, glmm.model.DCA2, recce.dca2.R2,
     glmm.model.recce.n.spp, recce.n.spp.R2, glmm.model.hfi.n.spp, hfi.n.spp.R2,
     file = "/home/sean/Documents/Molesworth2020/Analysis/RecceGlmmTMB.rdata")


new.hfi  <- expand.grid(year = c(1989, 2016), x = mean(hfi$x), y = mean(hfi$y), dem.altitude = mean(hfi$dem.altitude),
                            oversow=0, fenced=0, terrace=0, woody = c(0,1),  native = c(0,1), plot=NA)
new.hfi$form  <-  NA
new.hfi[new.hfi$native == 1 & new.hfi$woody==0,  "form"]     <-  "Native herbaceous"
new.hfi[new.hfi$native == 1 & new.hfi$woody==1,  "form"]     <-  "Native woody"
new.hfi[new.hfi$native == 0 & new.hfi$woody==0,  "form"]     <-  "Exotic herbaceous"
new.hfi[new.hfi$native == 0 & new.hfi$woody==1,  "form"]     <-  "Exotic woody"

new.hfi$hex.col  <-  NA
new.hfi[new.hfi$form == "Native herbaceous",  "hex.col"]     <-  nat.herb.col ##"#3aeb13" ##Light green 
new.hfi[new.hfi$form == "Native woody",  "hex.col"]     <-  nat.woody.col ## "#d1d487"  ##"#165e06" ### browny yellow
new.hfi[new.hfi$form == "Exotic herbaceous",  "hex.col"]     <- exot.herb.col ## "#cf0e0e"  ##red
new.hfi[new.hfi$form == "Exotic woody",  "hex.col"]     <-  exot.woody.col ## "#8c3424" ## "#0d0000"  


new.hfi$rich  <-  predict(m.hfi, new.hfi, type="response", re.form=NA, se.fit = FALSE)
new.hfi$rich.se  <-  predict(m.hfi, new.hfi, type="response", re.form=NA, se.fit = TRUE)$se.fit



  #################  #################  #################  #################  #################  #################
   #################  #################            Plot results               #################  #################
pdf("/home/sean/Documents/Molesworth/Graphs/PredictedFormRecce.pdf")
plot(n.spp ~ year, data=dat, type="n", bty="l", ylab=expression("Number of species (plot"^~2~")"), xaxt = "n", xlab="", ylim=(c(0,25)))#, xlim=(c(0,600)), axes=FALSE) ##main="Observed",
axis(1, at=sort(unique(dat$year)), labels = sort(unique(dat$year)))#c(substr(month.abb[11:12], 1,1), rep(substr(month.abb, 1,1), 2)), cex=.5)
##    lines(hfi ~ year, data=total.sum, lwd=1, col = "black",  lty=1)

for (i in unique(dat$form)) {  #    i   <-  "Native herbaceous"
    x.yr  <- new.dat[new.dat$form==i, "year"]
    Con.Low  <-  new.dat[new.dat$form==i, "rich"] - new.dat[new.dat$form==i, "rich.se"]
    Con.High  <-  new.dat[new.dat$form==i, "rich"] + new.dat[new.dat$form==i, "rich.se"]
    hex.colour  <- new.dat[new.dat$form==i,"hex.col"]
    ####Plot recce data from 1952 - 2016
    lines(rich ~ year, data=new.dat[new.dat$form==i,], lwd=2, col = hex.colour ,  lty=3)
    polygon(c(x.yr,rev(x.yr)),c(Con.Low,rev(Con.High)), col= hex_add_alpha(hex.colour, alpha = 0.1), border=NA) ## col=GISTools::add.alpha(hex.colour, alpha = 0.1)
    points(jitter(n.spp, amount=.7) ~ jitter(year, amount=.7), data=dat[dat$form==i,], pch=19, cex=0.25, col= dat[dat$form==i,"hex.col"])
    text(1966, (mean(new.dat[new.dat$form==i & new.dat$year==1987, "rich"]))*0.85, i, cex=.75)    

    ###HFI data from 2006 to 2016
    lines(rich ~ year, data=new.hfi[new.hfi$form==i,], lwd=2, col = hex.colour ,  lty=1)
    x.yr  <- new.hfi[new.hfi$form==i, "year"]
    Con.Low  <-  new.hfi[new.hfi$form==i, "rich"] - new.hfi[new.hfi$form==i, "rich.se"]
    Con.High  <-  new.hfi[new.hfi$form==i, "rich"] + new.hfi[new.hfi$form==i, "rich.se"]
    hex.colour  <- new.hfi[new.hfi$form==i,"hex.col"]
    polygon(c(x.yr,rev(x.yr)),c(Con.Low,rev(Con.High)),col=hex_add_alpha(hex.colour, alpha = 0.1),border=NA)

    ##text(2013, mean(new.dat[new.dat$form==i & new.dat$year==2016, "all.hfi"])+8, paste("All classes", i), cex=.75)
    
}
dev.off()






############################################################## ############################################################## ##############################################################
##################################################################     Pairwise Correlations for Most common species      ######################################################################

pairwise.common <-as.data.frame(comp); pairwise.common$year <-  as.numeric(substr(rownames(pairwise.common), regexpr("\\.20",rownames(pairwise.common))+1, regexpr("\\.20", rownames(pairwise.common))+5))
pairwise.common.2016 <- pairwise.common[pairwise.common$year==2016, c(1:length(names(pairwise.common))-1)]
rownames(pairwise.common.2016) <- substr(rownames(pairwise.common.2016), 1, regexpr("\\.20", rownames(pairwise.common.2016))-1)
pairwise.common.2016 <- pairwise.common.2016[order(rownames(pairwise.common.2016)),]

##occurance.matrix <- as.data.frame(matrix(0, nrow(pairwise.common.2016), ncol(pairwise.common.2016))); rownames(occurance.matrix) <- rownames(pairwise.common.2016); names(occurance.matrix) <- names(pairwise.common.2016)
occurance.matrix <- pairwise.common.2016
for(i in names(occurance.matrix)) {occurance.matrix[occurance.matrix[names(occurance.matrix)==i] != 0, names(occurance.matrix)==i] <- 1} ## produce a matrix with 1 for occurance
col.sums.occurance.matrix <- as.data.frame(colSums(occurance.matrix));  col.sums.occurance.matrix$spp <- rownames(col.sums.occurance.matrix)
common.spp.hfi <- sort(col.sums.occurance.matrix[col.sums.occurance.matrix[,1]>=20,2]) ###The species occuring in 20 or greater hfi plots

pairwise.common.2007 <- pairwise.common[pairwise.common$year==2007, c(1:length(names(pairwise.common))-1)]
rownames(pairwise.common.2007) <- substr(rownames(pairwise.common.2007), 1, regexpr("\\.20", rownames(pairwise.common.2007))-1)
pairwise.common.2007 <- pairwise.common.2007[order(rownames(pairwise.common.2007)),]

pairwise.change <- ((pairwise.common.2016[, common.spp.hfi] - pairwise.common.2007[, common.spp.hfi])) 
##pairwise.change[is.na(pairwise.change)] <- 0

pairwise.common.correlation <- cor(pairwise.change)

for (i in 1:ncol(pairwise.common.correlation)) {pairwise.common.correlation[pairwise.common.correlation[i,]==1, i] <- NA}

hfi.spp.correlations <- as.data.frame(apply(pairwise.common.correlation, 1, FUN=mean, na.rm=TRUE)); names(hfi.spp.correlations) <- "mean"; hfi.spp.correlations$species <- rownames(hfi.spp.correlations)

corr.stats <- cbind(c("min", "max", "mean"), round(c(min(hfi.spp.correlations$mean), max(hfi.spp.correlations$mean), mean(hfi.spp.correlations$mean)), digits=2)) 

save(corr.stats, common.spp, file="/home/sean/Documents/Molesworth2020/Analysis/Objects/HfiSppCorrelations.Rdata")

## pdf(paste("/home/sean/Documents/Molesworth2020/Analysis/Graphs/PairedHFIpeciesComparisons.pdf", sep="" )) ##i <- "Yes"
## corrgram(pairwise.common.correlation, use="pairwise.complete.obs", lower.panel=panel.shade, upper.panel=NULL,
##          main="Pairwise correlations of species biomass" , label.srt=60)
## dev.off()

library(ellipse)
library(RColorBrewer)
## my_colors <- brewer.pal(5, "Spectral")
## my_colors <- colorRampPalette(my_colors)(length(common.spp))
 
## # Order the correlation matrix
## ord <- order(pairwise.common.correlation[1, ])
## data_ord <- data[ord, ord]
## plotcorr(pairwise.common.correlation , col=my_colors[pairwise.common.correlation*length(my_colors)+length(my_colors)] , mar=c(1,1,1,1)  )





## plotcorr(pairwise.common.correlation, type = "lower",  label.pos = c(0.5, 0.5), cex.labels = 0.5, order=TRUE, col.regions = brewer.pal(n = 11, name = "RdYlGn"),  outer.labels = TRUE)

spp.labels  <- paste(substr(names(pairwise.change), 1,1), ". ", substr(names(pairwise.change), regexpr(" ", names(pairwise.change))+1, nchar(names(pairwise.change))), sep="")#[1:20]
spp.labels.twolines  <- paste(substr(names(pairwise.change), 1, regexpr(" ", names(pairwise.change))-1), "\n",
                              substr(names(pairwise.change), regexpr(" ", names(pairwise.change))+1, nchar(names(pairwise.change))), sep="")#[1:20]


pdf(paste("/home/sean/Documents/Molesworth2020/Analysis/Graphs/PairedHFIpeciesComparisons.pdf", sep="" )) ##i <- "Yes"

corrgram(cor(pairwise.change), labels= spp.labels, label.pos = c(0.5, .6), type = "cor", cex.labels = 0.35, order=TRUE, upper.panel = NULL,
         main="Pairwise correlations of species biomass" ,   col.regions=colorRampPalette(c("red", "orange", "darkkhaki", "darkgreen")),
oma=c(7, 7, 2, 2), outer.labels=list(bottom=list(labels=spp.labels,cex=0.5,srt=60), left=list(labels=spp.labels.twolines,cex=0.2,srt=0,adj=c(1.5,0))))


dev.off()




##################### Compare HFI with recce scores for reviewer 1


hfi.2007 <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/Molesworth2007HFI.csv", sep = ",", header=TRUE, skip = 5, stringsAsFactors=FALSE)
hfi.2007  <- hfi.2007[, c("Transect", "Plot", "Intercept", "Species", 	"Tier1")]
hfi.2007$plot  <-  paste(hfi.2007$Transect, hfi.2007$Plot, sep="") 
hfi.2007  <- as.data.frame.table(tapply(hfi.2007$Tier1, list(hfi.2007$plot, hfi.2007$Species), sum))
names(hfi.2007)  <- c("plot", "spp", "hfi")
hfi.2007  <-  hfi.2007[is.na(hfi.2007$hfi)==FALSE, ]
##hfi.2007[hfi.2007$hfi >=2, "hfi"]  <-  2 * hfi.2007[hfi.2007$hfi >=2, "hfi"]  ###Convert to 100 intercepts from 50

## hfi.2007$cover.hfi  <- 1; hfi.2007[which(hfi.2007$hfi > 1 & hfi.2007$hfi <= 5),"cover.hfi"] <- 2; hfi.2007[which(hfi.2007$hfi > 5 & hfi.2007$hfi <= 25),"cover.hfi"] <- 3
## hfi.2007[which(hfi.2007$hfi > 25 & hfi.2007$hfi <= 50),"cover.hfi"] <- 4; hfi.2007[which(hfi.2007$hfi > 50 & hfi.2007$hfi <= 75),"cover.hfi"] <- 5; hfi.2007[hfi.2007$hfi >=75,"cover.hfi"] <- 6;

hfi.2007$cover.hfi  <- 1; hfi.2007[which(hfi.2007$hfi > 3 & hfi.2007$hfi <= 10),"cover.hfi"] <- 2; hfi.2007[which(hfi.2007$hfi > 10 & hfi.2007$hfi <= 15),"cover.hfi"] <- 3
hfi.2007[which(hfi.2007$hfi > 15 & hfi.2007$hfi <= 25),"cover.hfi"] <- 4; hfi.2007[which(hfi.2007$hfi > 25 & hfi.2007$hfi <= 40),"cover.hfi"] <- 5; hfi.2007[hfi.2007$hfi >=40,"cover.hfi"] <- 6;



cover.2007 <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/Molesworth2007Cover.csv", sep = ",", header=TRUE, skip = 5, stringsAsFactors=FALSE)
cover.2007$cover  <- rowSums(cover.2007[,  c("Subplot4.Tier1", "Subplot4.Tier2", "Subplot4.Tier3", "Subplot4.Tier4", "Subplot4.Tier5",  "Subplot4.Tier6")], na.rm = TRUE)
cover.2007$plot  <-  paste(cover.2007$Transect, cover.2007$Plot, sep="") 
cover.2007  <- cover.2007[ , c("plot", "Species", "cover")]
names(cover.2007)  <-  c("plot", "spp", "cover")

dat.2007  <- merge(cover.2007, hfi.2007, by=c("plot", "spp"))
dat.2007  <-  dat.2007[dat.2007$cover >=1 & dat.2007$cover <=6, ]


cor(dat.2007$cover, dat.2007$cover.hfi)


pdf("/home/sean/Documents/Molesworth2020/Analysis/Graphs/HFIvsCover.pdf")
plot(cover.hfi ~ cover, data = dat.2007, ylim=c(min(0), max(6)), xlim=c(min(0), max(6)),
     bty="l", type="n", cex.axis = 1.25, cex.lab = 0.9,
     ylab = "Cover from HFI",     xlab = "Cover from recce scores")
points(jitter(dat.2007$cover) ~ jitter(dat.2007$cover.hfi), pch=19, cex = .15, type="p", col= "black")
dev.off()

pdf("/home/sean/Documents/Molesworth2020/Analysis/Graphs/PairCorrelation.pdf")
pairs(dat[, c("fenced",  "x", "y", "dem.altitude",  "terrace", "oversow")])
dev.off()

xtable(cor(dat[, c("fenced",  "x", "y", "dem.altitude",  "terrace", "oversow")]))


###### Numbers of spp found in HFI compared to Recce. Sub plot 3 is 2 x 2 . SUbplot 4 is 5 x 5.

hfi.2007 <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/Molesworth2007HFI.csv", sep = ",", header=TRUE, skip = 5, stringsAsFactors=FALSE)
hfi.2007  <- hfi.2007[, c("Transect", "Plot", "Intercept", "Species", 	"Tier1")]
hfi.2007$plot  <-  paste(hfi.2007$Transect, hfi.2007$Plot, sep="")

hfi.spp  <- unique(hfi.2007[, c("plot", "Species")])
hfi.spp  <- hfi.spp[hfi.spp$Species %nin% c("ACASP", "AGRSP", "CARSP",  "EPISP", "GRASP", "LICHEN", "MOSS", "RAOSP", "RYTSP", "POASP", "UNCSP"), ]
hfi.spp$hfi  <- 1



cover.2007 <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/Molesworth2007Cover.csv", sep = ",", header=TRUE, skip = 5, stringsAsFactors=FALSE)
cover.2007$plot  <-  paste(cover.2007$Transect, cover.2007$Plot, sep="") 

cover.2007  <- cover.2007[cover.2007$Species %nin% c("BARE GROUND", "GRASP", "ROCK", "UNCSP", "VEGETATION", "EPISP", "LICHEN", "LITTER", "RYTSP", "AGRSP", "CARSP", "ISOSP", "POASP", "MOSS"), ]

cover.spp  <- unique(cover.2007[is.na(cover.2007$Subplot4.Tier1)==FALSE, c("plot", "Species")])
cover.spp$recce  <- 1

spp.dat  <- merge(hfi.spp, cover.spp, by=c("plot", "Species"), all.x=TRUE, all.y=TRUE)

### HOw many recces missed spp found in HFI?
nrow(spp.dat[is.na(spp.dat$recce)==TRUE,])

### HOw many hfi missed spp found in recce?
nrow(spp.dat[is.na(spp.dat$hfi)==TRUE,])

##Summarise number of spp found for recce and hfi in each plot
spp.n  <- cbind(as.data.frame.table(tapply(spp.dat[is.na(spp.dat$hfi)==FALSE, "hfi"], spp.dat[is.na(spp.dat$hfi)==FALSE, "plot"], sum)),
     as.data.frame.table(tapply(spp.dat[is.na(spp.dat$recce)==FALSE, "recce"], spp.dat[is.na(spp.dat$recce)==FALSE, "plot"], sum))[, 2])
##spp.n[spp.n$plot == "N2", ]
names(spp.n)  <- c("plot", "hfi.n", "recce.n")
mean(spp.n$recce.n / spp.n$hfi.n)
sd(spp.n$recce.n / spp.n$hfi.n)
mean(spp.n$recce.n)
mean(spp.n$hfi.n)


Box.A  <- 9 * 2 * 2.4 
Box.B  <- 9 * 2.5 * 1.8 
Box.C  <-  5.8 * 1.5 * 1.6
Box.D  <- 3 * 1.8 * 1.8
Box.E  <- 4.5 * 2.3 * 3.7

(Box.E / (Box.A + Box.B + Box.C + Box.D + Box.E))  * 1350  ##(Box.A + Box.B + Box.C + Box.D + Box.E) 
