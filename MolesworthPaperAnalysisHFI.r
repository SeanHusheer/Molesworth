######This code provides data summaries, species lists and plot locations for Molesworth plots measured between 2006, 2007 and 2008, and re-measured between Jan and March 2016.
#####The 2006 survey consisted of the HFI plots established in 1989, the 2007 survey 80 randomly located plots, and the 2008 survey paired fenced plots.

###Written by Sean Husheer May 2021###############
#########################################################################################################################################################################
  #####################################################################################################################################################################
             ############   ###   ############   ###                  HFI Data      ############   ###   ############   ###   
#### SummaryHFIwithTiers.csv  is from SummariseHFI.r. Corrected for 50 intercepts
rm(list=ls()); setwd("/home/sean/Documents/Molesworth"); RPackages = "/home/sean/.RPackages"; .libPaths(RPackages)
RPackages = "/home/sean/.RPackages" # save typing later
.libPaths(RPackages) 
 options(width = 60) ### for desplay of summary
##options(width = 600)  ### For display of data frame
source("/home/sean/Documents/Molesworth2020/Analysis/glmmLatex.r") ###Import functon to use xtable for glmm summary tables. function = glmmTMB.latex

##install.packages("GISTools", RPackages, repos = "https://cran.stat.auckland.ac.nz",  dependencies = TRUE)
##update.packages(lib.loc = RPackages, repos = "https://cran.stat.auckland.ac.nz",  ask = FALSE, clean = TRUE, dependencies = TRUE)

dat <- read.csv("/home/sean/Documents/Molesworth/Data/SummaryHFIwithTiers.csv", , sep = ",", header=TRUE, skip = 0)#
packages.hfi  <- c("tools", "GISTools", "Hmisc", "sampSurf",  "sp", "sf", "glmmTMB", "DHARMa", "rgrass7", "rgeos", "maptools", "spatstat", "rgdal", "maps", "raster", "geoR", "xtable", "gridExtra", "plotrix", "jpeg", "gplots", "stargazer", "tables", "plyr", "car", "latex2exp", "pgirmess", "nlme", "lme4", "arm", "reshape", "contrast", "effects", "corrgram", "vegan", "R2jags", "mgcv", "fso", "diverse", "codyn", "dplyr", "rdddr"); 
invisible(lapply(packages.hfi, function(x) require(x, character.only = T, quietly = T)))

s.e. <- function(x) sqrt(var(na.omit(x))/length(na.omit(x))); l.s.d <-  function(x) (sqrt(var(na.omit(x))/length(na.omit(x))))*(qt(0.95,(length(na.omit(x)))))


load("Data/SpeciesColourCodesGraphics.rdata")


dat$form  <-  NA
dat[dat$native ==  "Native" & dat$woody=="Herbaceous",  "form"]     <-  "Native herbaceous"
dat[dat$native == "Native" & dat$woody=="Woody",  "form"]     <-  "Native woody"
dat[dat$native == "Exotic" & dat$woody=="Herbaceous",  "form"]     <-  "Exotic herbaceous"
dat[dat$native == "Exotic" & dat$woody=="Woody", "form"]     <-  "Exotic woody"
dat$hex.col  <-  NA
dat[dat$form == "Native herbaceous",  "hex.col"]     <-  nat.herb.col ## "#3aeb13"
dat[dat$form == "Native woody",  "hex.col"]     <- nat.woody.col ###  "#165e06"
dat[dat$form == "Exotic herbaceous",  "hex.col"]     <-   exot.herb.col ##  ##"#cf0e0e"
dat[dat$form == "Exotic woody",  "hex.col"]     <- exot.woody.col ## "#0d0000"


#######################################################################################################################################################################
      ##############################################################  Richness - Count                      ########################################################

m1.tmb <- glmmTMB(hfi  ~ scale(year) * scale(x) * scale(y) * scale(dem.altitude) +  (1|plot), data=dat, family=poisson)

dat[dat$native == "Exotic", "native"]  <-  0; dat[dat$native == "Native", "native"]  <-  1 ##Recode to 0 and 1 for exotic and native
dat[dat$woody == "Herbaceous", "woody"]  <-  0; dat[dat$woody == "Woody", "woody"]  <-  1 ##Recode to 0 and 1 for herb and native
dat$native  <-  as.integer(dat$native) ; dat$woody  <-  as.integer(dat$woody)

dat.total  <- as.data.frame.table(tapply(dat$hfi, list(dat$plot, dat$year), sum)); names(dat.total)  <-  c("plot", "year", "hfi")
dat.spatial  <- merge(dat.total[dat.total$year==2016, ], unique(dat[, c("plot", "x","y", "dem.altitude", "terrace", "fenced")]), all.x=TRUE, all.y=TRUE)

dat.total$un50  <- as.data.frame.table(tapply(dat$un50, list(dat$plot, dat$year), sum))[, 3]
dat.total$ov50  <- as.data.frame.table(tapply(dat$ov50, list(dat$plot, dat$year), sum))[, 3]

dat.total  <- merge(dat.total, unique(dat[, c("plot", "x","y", "dem.altitude", "terrace", "fenced")]), all.x=TRUE, all.y=TRUE)
dat.total  <-  dat.total[is.na(dat.total$hfi)==FALSE, ]
dat.total$year  <-  as.integer(as.character(dat.total$year))

####Spatial model of 2016 data so that Moran's I can be calculated    
m.spatial  <- glmmTMB(hfi  ~ scale(x) + scale(y) + scale(dem.altitude) +  (1|plot) , data=dat.spatial, family = nbinom2)
## To fit the model, a numFactor and a dummy grouping variable must be added to the dataset:
dat.spatial$pos <- numFactor(dat.spatial$x, dat.spatial$y)
dat.spatial$group <- factor(rep(1, nrow(dat.spatial)))
m.spatial.f <- glmmTMB(hfi ~ scale(x) + scale(y) + scale(dem.altitude) +  exp(pos + 0 | group), data=dat.spatial)
simulationOutput <- simulateResiduals(fittedModel = m.spatial.f, plot = F)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat.spatial$x, y = dat.spatial$y, plot=FALSE)

m.tot <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2  +  (1|plot),
                  data=dat.total, ziformula = ~0, family = nbinom2) ###No zero inflation
## m.tot.1 <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + terrace +fenced)^2  +  (1|plot),
##                   data=dat.total, ziformula = ~0, family = nbinom2) ###No zero inflation
## m.tot.3 <- glmmTMB(hfi  ~ scale(year) + scale(x) + scale(y) + scale(dem.altitude)  +  (1|plot),
##                   data=dat.total, ziformula = ~0, family = nbinom2) ###No zero inflation
## m.tot.4 <- glmmTMB(hfi  ~ scale(year) + scale(x) + scale(y) + scale(dem.altitude) ,    data=dat.total, ziformula = ~0, family = nbinom2) ###No zero inflation
## anova(m.tot.1, m.tot.2, m.tot.3, m.tot.4)

m.un50 <- glmmTMB(un50  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2  +  (1|plot),
                  data=dat.total, ziformula = ~0, family = nbinom2) ###No zero inflation
m.ov50 <- glmmTMB(un50  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2  +  (1|plot),
                  data=dat.total, ziformula = ~0, family = nbinom2) ###No zero inflation


## m.tot.p <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2  +  (1|plot),
##                   data=dat.total, ziformula = ~0, family = poisson) ###No zero inflation
## m.tot.g <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2  +  (1|plot),
##                   data=dat.total, ziformula = ~0, family = gaussian) ###No zero inflation
## anova(m.tot, m.tot.p, m.tot.g)

## simulationOutput <- simulateResiduals(fittedModel = m.tot, plot = F)
## plot(simulationOutput)


dat$terrace <-  as.character(dat$terrace) ;
dat[dat$terrace == "Slope", "terrace"]  <- 0; dat[dat$terrace == "Terrace", "terrace"]  <- 1;
dat$terrace  <- as.integer(as.character(dat$terrace))
dat$fenced <-  as.character(dat$fenced) ;
dat[dat$fenced == "unfenced", "fenced"]  <- 0; dat[dat$fenced == "fenced", "fenced"]  <- 1;
dat$fenced  <- as.integer(as.character(dat$fenced))
dat$oversow <-  as.character(dat$oversow) ;
dat[dat$oversow == "natural", "oversow"]  <- 0; dat[dat$oversow == "oversown", "oversow"]  <- 1;
dat$oversow  <- as.integer(as.character(dat$oversow))


##m1.tmb <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(woody) + scale(native))^2, data=dat, family = nbinom2)
m.hfi <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(oversow) + as.factor(fenced) + as.factor(native) + as.factor(woody))^2
                 +  (1|plot),  data=dat, family = nbinom2) ###No zero inflation ziformula = ~1,
## m.hfi.1 <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced)+scale(native) +scale(woody))^2  +
##                        (1|plot) ,    data=dat,  family = nbinom1) ###No zero inflation ziformula = ~1,

## m.hfi.2 <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced) +scale(native) +scale(woody))^2  +
##                        (1|plot) ,  data=dat,  family = poisson) ###No zero inflation ziformula = ~1,

## m.hfi.3 <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced) +scale(native) +scale(woody))^2  +
##                        (1|plot) ,   data=dat, ziformula = ~1,  family=truncated_nbinom1) ###No zero inflation ziformula = ~1,

## m.hfi.4 <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(oversow) + scale(fenced) +scale(native) +scale(woody))^2  +
##                        (1|plot) ,  data=dat, ziformula = ~1, family=truncated_poisson) ###No zero inflation ziformula = ~1,

## anova(m.hfi, m.hfi.1, m.hfi.2, m.hfi.3, m.hfi.4)


## dat$model.fit <- predict(m.hfi, data=dat)
## R2.1  <- (cor(dat$hfi, dat$model.fit))^2
##R2.2  <-  1 - sum((dat$hfi - dat$model.fit)^2)/sum((dat$hfi-mean(dat$hfi))^2)

##MuMIn::r.squaredLR(m.hfi)  <- lme(distance ~ 1, ~ 1 | Subject, data = Orthodont)
m.hfi.null <- glmmTMB(hfi ~ 1 +  (1|plot), data=dat, family = nbinom2)

simulationOutput <- simulateResiduals(fittedModel = m.hfi, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

hfi.biomass.R2  <- MuMIn::r.squaredGLMM(m.hfi, m.hfi.null)
glmm.model.hfi.biomass  <- glmmTMB.latex(m.hfi)

## glmmTMB.latex <- function(glmm.table) {
 
## glmm.table  <- as.data.frame(round(summary(m.hfi)$coefficients$cond, digits=3))
## glmm.table  <- glmm.table[order(abs(glmm.table$Estimate), decreasing = TRUE), ]
## rownames(glmm.table)  <-  gsub("scale", "", rownames(glmm.table))
## rownames(glmm.table)  <-  gsub("\\(", "", rownames(glmm.table))
## rownames(glmm.table)  <-  gsub("\\)", "", rownames(glmm.table))
## rownames(glmm.table)  <-  gsub("dem.", "", rownames(glmm.table))
## rownames(glmm.table)  <-  toTitleCase(rownames(glmm.table))
## rownames(glmm.table)  <-gsub("^(x)", "Easting (m)", rownames(glmm.table)) ###first character in string only
## rownames(glmm.table)  <-gsub("^(y)", "Northing (m)", rownames(glmm.table)) 
## rownames(glmm.table)  <-  gsub(":", " x ", rownames(glmm.table))
## names(glmm.table)  <- c("Estimate", "Estimate Std Error", "Z value", "P value")  

## return(xtable(glmm.table))
## }


m0.tmb <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(woody) + as.factor(native))^2 +  (1|plot),
                  data=dat, ziformula = ~ 1, family = nbinom2)

m.un50 <- glmmTMB(un50  ~ (scale(year) + scale(dem.altitude) + as.factor(woody) + as.factor(native))^3 + scale(x) + scale(y)  +  (1|plot),
                  data=dat, ziformula = ~0, family = nbinom2) ###

m.ov50 <- glmmTMB(ov50  ~ (scale(year) + scale(dem.altitude) + as.factor(woody) + as.factor(native))^3 + scale(x) + scale(y)  +  (1|plot),
                  data=dat, ziformula = ~0, family = nbinom2) ###

new.hfi.dat  <- expand.grid(year = sort(unique(dat$year)), x = mean(dat$x), y = mean(dat$y),  dem.altitude = mean(dat$dem.altitude),
                            oversow=0, terrace=0, fenced=0, woody = c(0,1),  native = c(0,1), plot=NA)
new.hfi.dat$form  <-  NA
new.hfi.dat[new.hfi.dat$native == 1 & new.hfi.dat$woody==0,  "form"]     <-  "Native herbaceous"
new.hfi.dat[new.hfi.dat$native == 1 & new.hfi.dat$woody==1,  "form"]     <-  "Native woody"
new.hfi.dat[new.hfi.dat$native == 0 & new.hfi.dat$woody==0,  "form"]     <-  "Exotic herbaceous"
new.hfi.dat[new.hfi.dat$native == 0 & new.hfi.dat$woody==1,  "form"]     <-  "Exotic woody"

new.hfi.dat$hex.col  <-  NA
new.hfi.dat[new.hfi.dat$form == "Native herbaceous",  "hex.col"]     <-  nat.herb.col ## "#3aeb13" ##Light green 
new.hfi.dat[new.hfi.dat$form == "Native woody",  "hex.col"]     <-  nat.woody.col ## "#d1d487"  ##"#165e06" ### browny yellow
new.hfi.dat[new.hfi.dat$form == "Exotic herbaceous",  "hex.col"]     <- exot.herb.col ##  "#cf0e0e"  ##red
new.hfi.dat[new.hfi.dat$form == "Exotic woody",  "hex.col"]     <- exot.woody.col ## "#8c3424" ## "#0d0000"  


new.hfi.dat$all.hfi  <-  predict(m.hfi, new.hfi.dat, type="response", re.form=NA, se.fit = FALSE)
new.hfi.dat$all.hfi.se  <-  predict(m.hfi, new.hfi.dat, type="response", re.form=NA, se.fit = TRUE)$se.fit
new.hfi.dat$un50.hfi  <-  predict(m.un50, new.hfi.dat, type="response", re.form=NA)
new.hfi.dat$un50.hfi.se  <-  predict(m.un50, new.hfi.dat, type="response", re.form=NA, se.fit = TRUE)$se.fit
new.hfi.dat$ov50.hfi  <-  predict(m.ov50, new.hfi.dat, type="response", re.form=NA)
new.hfi.dat$ov50.hfi.se  <-  predict(m.ov50, new.hfi.dat, type="response", re.form=NA, se.fit = TRUE)$se.fit

total.sum  <- as.data.frame.table(tapply(dat.total$hfi, dat.total$year, mean))
names(total.sum)  <-  c("year", "hfi"); total.sum$year  <- as.integer(as.character(total.sum$year)); total.sum$hfi  <- as.integer(as.character(total.sum$hfi))

pdf("/home/sean/Documents/Molesworth/Graphs/PredictedHFIform.pdf")
plot(hfi ~ year, data=dat, type="n", bty="l", ylab=expression("HFI (plot"^~2~")"), xaxt = "n", xlab="", ylim=(c(0,200)))#, xlim=(c(0,600)), axes=FALSE) ##main="Observed",
axis(1, at=sort(unique(dat$year)), labels = sort(unique(dat$year)))#c(substr(month.abb[11:12], 1,1), rep(substr(month.abb, 1,1), 2)), cex=.5)
##    lines(hfi ~ year, data=total.sum, lwd=1, col = "black",  lty=1)

for (i in unique(dat$form)) {  #    i   <-  "Native herbaceous"
    x.yr  <- new.hfi.dat[new.hfi.dat$form==i, "year"]
    Con.Low  <-  new.hfi.dat[new.hfi.dat$form==i, "all.hfi"] - new.hfi.dat[new.hfi.dat$form==i, "all.hfi.se"]
    Con.High  <-  new.hfi.dat[new.hfi.dat$form==i, "all.hfi"] + new.hfi.dat[new.hfi.dat$form==i, "all.hfi.se"]
    hex.colour  <- new.hfi.dat[new.hfi.dat$form==i,"hex.col"]
    polygon(c(x.yr,rev(x.yr)),c(Con.Low,rev(Con.High)), col= hex_add_alpha(hex.colour, alpha = 0.1), border=NA)

    lines(all.hfi ~ year, data=new.hfi.dat[new.hfi.dat$form==i,], lwd=2, col = hex.colour ,  lty=1)
    points(hfi ~ jitter(year, amount=.7), data=dat[dat$form==i,], pch=19, cex=0.25, col= hex.colour)
    text(2013, (mean(new.hfi.dat[new.hfi.dat$form==i & new.hfi.dat$year==2016, "all.hfi"])+8)*0.95, i, cex=.75)
    
}
dev.off()

## pdf("/home/sean/Documents/Molesworth2020/Analysis/Graphs/PredictedHFIformUnder50.pdf")
## plot(hfi ~ year, data=dat, type="n", bty="l", ylab=expression("HFI (plot"^~2~")"), xaxt = "n", xlab="", ylim=(c(0,200)))#, xlim=(c(0,600)), axes=FALSE) ##main="Observed",
## axis(1, at=sort(unique(dat$year)), labels = sort(unique(dat$year)))#c(substr(month.abb[11:12], 1,1), rep(substr(month.abb, 1,1), 2)), cex=.5)
## ##    lines(hfi ~ year, data=total.sum, lwd=1, col = "black",  lty=1)

## for (i in unique(dat$form)) {  #    i   <-  "Native herbaceous"
##     x.yr  <- new.hfi.dat[new.hfi.dat$form==i, "year"]
##     Con.Low  <-  new.hfi.dat[new.hfi.dat$form==i, "all.hfi"] - new.hfi.dat[new.hfi.dat$form==i, "all.hfi.se"]
##     Con.High  <-  new.hfi.dat[new.hfi.dat$form==i, "all.hfi"] + new.hfi.dat[new.hfi.dat$form==i, "all.hfi.se"]
##     hex.colour  <- new.hfi.dat[new.hfi.dat$form==i,"hex.col"]
##     polygon(c(x.yr,rev(x.yr)),c(Con.Low,rev(Con.High)),col=GISTools::add.alpha(hex.colour, alpha = 0.1),border=NA)

##     ## lines(all.hfi ~ year, data=new.hfi.dat[new.hfi.dat$form==i,], lwd=2, col = hex.colour ,  lty=1)
##     ## lines(ov50.hfi ~ year, data=new.hfi.dat[new.hfi.dat$form==i,], lwd=1, col = hex.colour,  lty=2) 
##     lines(un50.hfi ~ year, data=new.hfi.dat[new.hfi.dat$form==i,], lwd=1, col = hex.colour,  lty=3) 

##     points(hfi ~ jitter(year, amount=.7), data=dat[dat$form==i,], pch=19, cex=0.25, col= dat[dat$form==i,"hex.col"])
## ##    text(2013, mean(new.hfi.dat[new.hfi.dat$form==i & new.hfi.dat$year==2016, "all.hfi"])+8, paste("All classes", i), cex=.75)
## }
## dev.off()













## dat$pos <- numFactor(dat$x, dat$y)
## dat$group <- factor(rep(1, nrow(dat)))

## m.spat.tmb <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(woody) + scale(native))^2  +  (1|plot),
##                   data=dat, ziformula = ~0, family = nbinom2) ###Include spatial covariance matrix


## simulationOutput <- simulateResiduals(fittedModel = m.hfi, plot = F)
## plot(simulationOutput)
## testDispersion(simulationOutput)
## testZeroInflation(simulationOutput)

## m1 <-lme(log(hfi+1)  ~ year * woody * native, random = ~ 1|plot, data=dat) ##, family=poisson)  ## method="ML"
## m2 <-lme(log(un50+1)  ~ year * woody * native, random = ~ 1|plot, data=dat) ##, family=poisson)  ## method="ML"
## m3 <-lme(log(ov50+1)  ~ year * woody * native, random = ~ 1|plot, data=dat) ##, family=poisson)  ## method="ML"














## for(i in survey.list) {
##     i  <-  "haphazard.plot.list"
## title <- paste(toupper(substr(i,1,1)), substr(i, 2, regexpr(".plot", i)[1]-1), " Plots", sep="")
## temp.plot.list <- plot.list[plot.list[,1]==i,2]
## hfi.temp <- dat[which(dat$plot %in% temp.plot.list),]
## n.site.mean <- as.data.frame.table(tapply(hfi.temp$n, hfi.temp$year, mean))
## n.site.sem <- as.data.frame.table(tapply(hfi.temp$n, list(hfi.temp$fenced.year, hfi.temp$year), s.e.))
## ##    n.site.sem <- as.data.frame.table(tapply(hfi.temp$n, list(hfi.temp$fenced.year, hfi.temp$year), s.e.))
## n <- (cbind(n.site.mean, n.site.sem))[,c(1:3,6)]
## names(n) <- c("fenced", "year", "n.mean", "n.sem"); n$year <- as.numeric(as.character(n$year))
## pdf(paste("/home/sean/Documents/Molesworth2020/Analysis/Graphs/SpeciesCount", sub(" ","", title),  ".pdf", sep="" ))
## plot(n$n.mean ~ n$year, ylim= c(min(0), max(round_any(max(n$n.mean), 1, ceiling))), xlim=c(min(min(n$year)),  max(max(n$year))), type="n", bty="l", ylab = "Mean Species Richness", xlab = "Year", main=paste(title, "\nCounts of plot species richness", sep=""))
## points(n[n$fenced=="Fenced", "n.mean"] ~ n[n$fenced=="Fenced", "year"], pch=19, cex=0.5, type="b")
## points(n[n$fenced=="Unfenced", "n.mean"] ~ n[n$fenced=="Unfenced", "year"], pch="O", cex=0.5, type="b", lty=2)
## arrows(n$year, (n$n.mean + n$n.sem), n$year, (n$n.mean - n$n.sem), length=0.025, angle=90, code=3, lwd = 0.25)
## dev.off()
## }
#######################################################################################################################################################################
##############################################################      Proportion of native species               ########################################################

dat <- read.csv("/home/sean/Documents/Molesworth/Data/ProportionSummaryHFI.csv", sep = ",", header=TRUE, skip = 0)#

dat[dat$terrace == "Slope", "terrace"]  <- 0; dat[dat$terrace == "Terrace", "terrace"]  <- 1; dat$terrace  <- as.integer(as.character(dat$terrace))
dat[dat$fenced == "unfenced", "fenced"]  <- 0; dat[dat$fenced == "fenced", "fenced"]  <- 1; dat$fenced  <- as.integer(as.character(dat$fenced))
dat[dat$oversown == "natural", "oversown"]  <- 0; dat[dat$oversown == "oversown", "oversown"]  <- 1; dat$oversown  <- as.integer(as.character(dat$oversown))

m.hfi <- glmmTMB(cbind(hfi.native, hfi.exotic) ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(fenced) + as.factor(oversown))^2 +
                     (1|plot),   data=dat, family = binomial)
## m.hfi <- glmmTMB(cbind(hfi.native, hfi.exotic)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace))^2 +  (1|plot),
##                 data=dat, family = binomial)

## m3.tmb <- glmmTMB(cbind(hfi.native, hfi.exotic) ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace))^2 +  (1|plot),
##                 data=dat, family = binomial, ziformula = ~ 1)

## anova(m1.tmb, m2.tmb, m3.tmb)  

## simulationOutput <- simulateResiduals(fittedModel = m1.tmb, plot = F)
## plot(simulationOutput)

m.hfi.null <- glmmTMB(cbind(hfi.native, hfi.exotic) ~ 1 +  (1|plot), data=dat, family = binomial)

hfi.prop.native.R2  <- MuMIn::r.squaredGLMM(m.hfi, m.hfi.null)
glmm.model.hfi.prop.native  <- glmmTMB.latex(m.hfi)


temp.dat  <- data.frame(temp.j = rep(c("terrace", "fenced", "oversown"), each=2), level.j = c("Slope", "Terrace", "Unfenced", "Fenced", "Natural", "Oversown"),
           num.j = rep(c(0,1), 3), col.j  = rep(c("#f2ea0a", "#0762f5"), 3))

pdf("/home/sean/Documents/Molesworth/Graphs/PredictedProportionNative.pdf")
par(mfrow=c(3,1), cex=1, mar=c(0,5,2,2)+0.1) #pty="s",  The default mar is ‘c(5, 4, 4, 2) +0.1’. ‘c(bottom, left, top, right)’

for (j in c("terrace", "fenced", "oversown")) { ##j  <-     "terrace"

new.hfi.dat  <- expand.grid(year = sort(unique(dat$year)), x = mean(dat$x), y = mean(dat$y),  dem.altitude = mean(dat$dem.altitude),
                            temp = c(0,1) , plot=NA)

names(new.hfi.dat)[names(new.hfi.dat) == "temp"]  <-  j
      
## new.hfi.dat[ , toTitleCase(j)]  <-  NA
## new.hfi.dat[new.hfi.dat$terrace == 0,  "Terrace"]     <-  "Slope"
## new.hfi.dat[new.hfi.dat$terrace == 1,  "Terrace"]     <-  "Terrace"
## new.hfi.dat$hex.col  <-  NA
## new.hfi.dat[new.hfi.dat$terrace == 0,  "hex.col"]     <-  "#f2ea0a" ##An deep yellow
## new.hfi.dat[new.hfi.dat$terrace == 1,  "hex.col"]     <-  "#0762f5" ###An intense blue

formula.i <- as.formula(paste("cbind(hfi.native, hfi.exotic)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(", j, "))^2 +  (1|plot)", sep=""))

m.hfi <- glmmTMB(formula.i, data=dat, family = binomial)                
                
new.hfi.dat$hfi  <-  predict(m.hfi, new.hfi.dat, type="response", re.form=NA, se.fit = FALSE)
new.hfi.dat$hfi.se  <-  predict(m.hfi, new.hfi.dat, type="response", re.form=NA, se.fit = TRUE)$se.fit

plot(hfi ~ year, data=new.hfi.dat, type="n", bty="l", ylab="Proportion native", xaxt = "n", xlab="", ylim=(c(0.2,1.1)))#, xlim=(c(0,600)), axes=FALSE) ##main="Observed",
axis(1, at=sort(unique(dat$year)), labels = sort(unique(dat$year)))#c(substr(month.abb[11:12], 1,1), rep(substr(month.abb, 1,1), 2)), cex=.5)

for (i in c(0,1)) {  #
   
    x.yr  <- new.hfi.dat[new.hfi.dat[, j] == i, "year"]
    Con.Low  <-  new.hfi.dat[new.hfi.dat[, j]==i, "hfi"] - new.hfi.dat[new.hfi.dat[, j]==i, "hfi.se"]
    Con.High  <-  new.hfi.dat[new.hfi.dat[, j]==i, "hfi"] + new.hfi.dat[new.hfi.dat[, j]==i, "hfi.se"]
    hex.colour  <- temp.dat[temp.dat$temp.j == j & temp.dat$num.j == i, "col.j"]
    polygon(c(x.yr,rev(x.yr)),c(Con.Low,rev(Con.High)), col=hex_add_alpha(hex.colour, alpha = 0.1), border=NA) ##col=GISTools::add.alpha(hex.colour, alpha = 0.1)

if (i == 0) { 
y.pos  <- new.hfi.dat[new.hfi.dat$year == 1989 & new.hfi.dat[, j]==i, "hfi"] - 0.15
} else {
y.pos  <- new.hfi.dat[new.hfi.dat$year == 1989 & new.hfi.dat[, j]==i, "hfi"] + 0.15
}

        
    lines(hfi ~ year, data=new.hfi.dat[new.hfi.dat[, j]==i,], lwd=2, col = hex.colour ,  lty=1)
    text(1992, y.pos, temp.dat[temp.dat$temp.j == j & temp.dat$num.j==i, "level.j"])
   ## points(hfi ~ jitter(year, amount=.7), data=dat[dat$terrace==i,], pch=19, cex=0.25, col= dat[dat$terrace==i,"hex.col"])
}
    }
dev.off()

#######################################################################################################################################################################
##############################################################      Proportion of woody species               ########################################################

m.hfi <- glmmTMB(cbind(hfi.woody, hfi.herbaceous) ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + as.factor(terrace) + as.factor(fenced) + as.factor(oversown))^2 +
                     (1|plot),   data=dat, family = binomial)

simulationOutput <- simulateResiduals(fittedModel = m.hfi, plot = F)
plot(simulationOutput)

m.hfi.null <- glmmTMB(cbind(hfi.woody, hfi.herbaceous) ~ 1 +  (1|plot), data=dat, family = binomial)

hfi.prop.woody.R2  <- MuMIn::r.squaredGLMM(m.hfi, m.hfi.null)
glmm.model.hfi.prop.woody  <- glmmTMB.latex(m.hfi)



save(glmm.model.hfi.prop.woody, hfi.prop.woody.R2, glmm.model.hfi.prop.native, hfi.prop.native.R2,  glmm.model.hfi.biomass, hfi.biomass.R2,
     file = "/home/sean/Documents/Molesworth/Data/Objects/HFIglmmTMB.rdata")



pdf("/home/sean/Documents/Molesworth/Graphs/PredictedProportionWoody.pdf")

par(mfrow=c(2,1), cex=1, mar=c(0,5,2,2)+0.1) #pty="s",  The default mar is ‘c(5, 4, 4, 2) +0.1’. ‘c(bottom, left, top, right)’

for (j in c("terrace", "fenced")) { ##j  <-     "terrace"

new.hfi.dat  <- expand.grid(year = sort(unique(dat$year)), x = mean(dat$x), y = mean(dat$y),  dem.altitude = mean(dat$dem.altitude),
                            temp = c(0,1) , plot=NA)

names(new.hfi.dat)[names(new.hfi.dat) == "temp"]  <-  j
      
formula.i <- as.formula(paste("cbind(hfi.woody, hfi.herbaceous)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(", j, "))^2 +  (1|plot)", sep=""))

m.hfi <- glmmTMB(formula.i, data=dat, family = binomial)                
                
new.hfi.dat$hfi  <-  predict(m.hfi, new.hfi.dat, type="response", re.form=NA, se.fit = FALSE)
new.hfi.dat$hfi.se  <-  predict(m.hfi, new.hfi.dat, type="response", re.form=NA, se.fit = TRUE)$se.fit

plot(hfi ~ year, data=new.hfi.dat, type="n", bty="l", ylab="Proportion woody", xaxt = "n", xlab="", ylim=(c(0,0.4)))#, xlim=(c(0,600)), axes=FALSE) ##main="Observed",
axis(1, at=sort(unique(dat$year)), labels = sort(unique(dat$year)))#c(substr(month.abb[11:12], 1,1), rep(substr(month.abb, 1,1), 2)), cex=.5)

for (i in c(0,1)) {  #
   
    x.yr  <- new.hfi.dat[new.hfi.dat[, j] == i, "year"]
    Con.Low  <-  new.hfi.dat[new.hfi.dat[, j]==i, "hfi"] - new.hfi.dat[new.hfi.dat[, j]==i, "hfi.se"]
    Con.High  <-  new.hfi.dat[new.hfi.dat[, j]==i, "hfi"] + new.hfi.dat[new.hfi.dat[, j]==i, "hfi.se"]
    hex.colour  <- temp.dat[temp.dat$temp.j == j & temp.dat$num.j == i, "col.j"]
    polygon(c(x.yr,rev(x.yr)),c(Con.Low,rev(Con.High)),col=hex_add_alpha(hex.colour, alpha = 0.1), border=NA) ##GISTools::add.alpha(hex.colour, alpha = 0.1)

if (i == 0) { 
y.pos  <- new.hfi.dat[new.hfi.dat$year == 1989 & new.hfi.dat[, j]==i, "hfi"] - 0.015
} else {
y.pos  <- new.hfi.dat[new.hfi.dat$year == 1989 & new.hfi.dat[, j]==i, "hfi"] + 0.015
}

        
    lines(hfi ~ year, data=new.hfi.dat[new.hfi.dat[, j]==i,], lwd=2, col = hex.colour ,  lty=1)
    text(1992, y.pos, temp.dat[temp.dat$temp.j == j & temp.dat$num.j==i, "level.j"])
   ## points(hfi ~ jitter(year, amount=.7), data=dat[dat$terrace==i,], pch=19, cex=0.25, col= dat[dat$terrace==i,"hex.col"])
}
    }
dev.off()







#### Cut from here
























## m1.tmb <- glmmTMB(cbind(hfi.woody, hfi.herbaceous)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(fenced))^2 +  (1|plot),
##                 data=dat, family = binomial)
## m.hfi <- glmmTMB(cbind(hfi.woody, hfi.herbaceous)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(fenced))^2 +  (1|plot),
##                 data=dat, family = binomial)
## m3.tmb <- glmmTMB(cbind(hfi.woody, hfi.herbaceous) ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(fenced))^2 +  (1|plot),
##                 data=dat, family = binomial, ziformula = ~ 1)

## anova(m1.tmb, m.hfi, m3.tmb)  

## simulationOutput <- simulateResiduals(fittedModel = m.hfi, plot = F)
## plot(simulationOutput)

## new.hfi.dat  <- expand.grid(year = sort(unique(dat$year)), x = mean(dat$x), y = mean(dat$y),  dem.altitude = mean(dat$dem.altitude),
##                             fenced = c(0,1), plot=NA)
## new.hfi.dat$Fenced  <-  NA
## new.hfi.dat[new.hfi.dat$fenced == 0,  "Fenced"]     <-  "Cattle present"
## new.hfi.dat[new.hfi.dat$fenced == 1,  "Fenced"]     <-  "Cattle excluded"

## new.hfi.dat$hex.col  <-  NA
## new.hfi.dat[new.hfi.dat$fenced == 0,  "hex.col"]     <-  "#f2ea0a" ##An deep yellow
## new.hfi.dat[new.hfi.dat$fenced == 1,  "hex.col"]     <-  "#0762f5" ###An intense blue

## new.hfi.dat$hfi  <-  predict(m.hfi, new.hfi.dat, type="response", re.form=NA, se.fit = FALSE)
## new.hfi.dat$hfi.se  <-  predict(m.hfi, new.hfi.dat, type="response", re.form=NA, se.fit = TRUE)$se.fit

## plot(hfi ~ year, data=new.hfi.dat, type="n", bty="l", ylab="Proportion woody", xaxt = "n", xlab="", ylim=(c(0,0.1)))#, xlim=(c(0,600)), axes=FALSE) ##main="Observed",
## axis(1, at=sort(unique(dat$year)), labels = sort(unique(dat$year)))#c(substr(month.abb[11:12], 1,1), rep(substr(month.abb, 1,1), 2)), cex=.5)

## for (i in unique(dat$fenced)) {  #
   
##     x.yr  <- new.hfi.dat[new.hfi.dat$fenced== i, "year"]
##     Con.Low  <-  new.hfi.dat[new.hfi.dat$fenced==i, "hfi"] - new.hfi.dat[new.hfi.dat$fenced==i, "hfi.se"]
##     Con.High  <-  new.hfi.dat[new.hfi.dat$fenced==i, "hfi"] + new.hfi.dat[new.hfi.dat$fenced==i, "hfi.se"]
##     hex.colour  <- new.hfi.dat[new.hfi.dat$fenced==i,"hex.col"]
##     polygon(c(x.yr,rev(x.yr)),c(Con.Low,rev(Con.High)),col=GISTools::add.alpha(hex.colour, alpha = 0.1),border=NA)
##     lines(hfi ~ year, data=new.hfi.dat[new.hfi.dat$fenced==i,], lwd=2, col = hex.colour ,  lty=1)
    

## }






















m1.tmb <- glm(cbind(hfi.woody, hfi.herbaceous)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2 , data=dat, family = binomial)




m1.tmb <- glmmTMB(cbind(hfi.woody, hfi.herbaceous)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(terrace) + scale(fenced))^2 +  (1|plot),
                  data=dat, family = binomial)
simulationOutput <- simulateResiduals(fittedModel = m1.tmb, plot = F)
plot(simulationOutput)





m1 <- glm(c(hfi.native, hfi.exotic)  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude))^2, data=dat, family = binomial)



m.hfi <- glmmTMB(hfi  ~ (scale(year) + scale(dem.altitude) + scale(woody) + scale(native))^3 + scale(x) + scale(y)  +  (1|plot),
                  data=dat, ziformula = ~0, family = nbinom2) ###No zero inflation
m0.tmb <- glmmTMB(hfi  ~ (scale(year) + scale(x) + scale(y) + scale(dem.altitude) + scale(woody) + scale(native))^2 +  (1|plot),
                  data=dat, ziformula = ~ 1, family = nbinom2)












survey.list <- c("haphazard.plot.list", "random.plot.list", "paired.plot.list")
haphazard.plot.list <- cbind("haphazard.plot.list", c(unique(as.character(dat[dat$year==2006, "plot"])), "Saxton4"))
random.plot.list <- cbind("random.plot.list", unique(as.character(dat[dat$year==2007, "plot"])))
paired.plot.list <- cbind("paired.plot.list", unique(as.character(dat[dat$year==2008, "plot"])))
plot.list <- rbind(haphazard.plot.list, random.plot.list, paired.plot.list)




for(i in survey.list) {
title <- paste(toupper(substr(i,1,1)), substr(i, 2, regexpr(".plot", i)[1]-1), " Plots", sep="")
temp.plot.list <- plot.list[plot.list[,1]==i,2]
hfi.temp <- hfi.dat[which(hfi.dat$plot %in% temp.plot.list),]

m.native <- glm(cbind(native.n,  exotic.n) ~ year * fenced.year, data=hfi.temp, family=binomial)
table <- as.data.frame(round(coef(summary(m.native))[,c(1:2,4)], digits=3))
table[,names(table)] <- lapply(names(table), function(x) as.character(table[,x]))
table[table=="0"] <-  "$<$0.001"
assign(paste("native.table.", substr(i, 1, nchar(i)-10), sep=""), table)
   
native.site.mean <- as.data.frame.table(tapply(hfi.temp$prop.native, list(hfi.temp$fenced.year, hfi.temp$year), mean))
native.site.sem <- as.data.frame.table(tapply(hfi.temp$prop.native, list(hfi.temp$fenced.year, hfi.temp$year), s.e.))
native <- (cbind(native.site.mean, native.site.sem))[,c(1:3,6)]
names(native) <- c("fenced", "year", "native.mean", "native.sem"); native$year <- as.numeric(as.character(native$year))

pdf(paste("/home/sean/Documents/Molesworth/Graphs/NativeSpecies", sub(" ","", title),  ".pdf", sep="" ))
  
plot(native$native.mean ~ native$year, ylim= c(min(0), max(round_any(max(native$native.mean), 1, ceiling))), xlim=c(min(min(native$year)),  max(max(native$year))), type="n", bty="l", ylab = "Proportion of native species", xlab = "Year", main=paste(title, "\nProportion of native species", sep=""))
points(native[native$fenced=="Fenced", "native.mean"] ~ native[native$fenced=="Fenced", "year"], pch=19, cex=0.5, type="b")
points(native[native$fenced=="Unfenced", "native.mean"] ~ native[native$fenced=="Unfenced", "year"], pch="O", cex=0.5, type="b", lty=2)
arrows(native$year, (native$native.mean + native$native.sem), native$year, (native$native.mean - native$native.sem), length=0.025, angle=90, code=3, lwd = 0.25)
dev.off()
}
save(native.table.haphazard, native.table.random, native.table.paired, file = "/home/sean/Documents/Molesworth/Data/Objects/HFINative.Table.Rdata") ##


hfi.2016 <- hfi.dat[hfi.dat$year==2016, c(1:24)]
m1 <-  glm(cbind(native.n, exotic.n) ~ fenced.year + dem.altitude + x + y + water.deficit, data = hfi.2016, family=binomial)
m1 <-  glm(prop.native ~ dem.altitude + x + y, data = hfi.2016)
m1.alt <-  glm(cbind(native.n, exotic.n) ~ fenced.year + dem.altitude, data = hfi.2016, family=binomial)
m1.x <-  glm(cbind(native.n, exotic.n) ~ fenced.year + x, data = hfi.2016, family=binomial)
m1.y <-  glm(cbind(native.n, exotic.n) ~ fenced.year + y, data = hfi.2016, family=binomial)
m1.def <-  glm(cbind(native.n, exotic.n) ~ fenced.year * water.deficit, data = hfi.2016, family=binomial)


x.alt= seq(0,2000,by=10); x.def=seq(0,200); x.alt.high= seq(1000,2000,by=5); x.x = seq(min(hfi.2016$x), max(hfi.2016$x), by=100); x.y = seq(min(hfi.2016$y), max(hfi.2016$y), by=100)

plot(hfi.2016$prop.native ~ hfi.2016$dem.altitude, ylim= c(min(0), max(1)), xlim=c(min(800),  max(1400)), type="n", bty="l", ylab = expression("Proportion native-exotic species"), xlab = "Altitude (m a.s.l.)", main="")
points(hfi.2016[hfi.2016$fenced.year=="Fenced", "prop.native"] ~ hfi.2016[hfi.2016$fenced.year=="Fenced", "dem.altitude"], pch=19, cex=0.75, col="green") 
lines(x.alt, predict(m1.alt, list(dem.altitude = x.alt, fenced.year=rep("Fenced", length(x.alt))), type="response", se.fit = FALSE), lwd=2, col="green",  lty=1)
points(hfi.2016[hfi.2016$fenced.year=="Unfenced", "prop.native"] ~ hfi.2016[hfi.2016$fenced.year=="Unfenced", "dem.altitude"], pch="O", cex=0.75, col="brown") 
lines(x.alt, predict(m1.alt, list(dem.altitude = x.alt, fenced.year=rep("Unfenced", length(x.alt))), type="response", se.fit = FALSE), lwd=2, col="brown",  lty=2)

plot(hfi.2016$prop.native ~ hfi.2016$x, ylim= c(min(0), max(1)), xlim=c(min(min(hfi.2016$x)),  max(max(hfi.2016$x))), type="n", bty="l", ylab = expression("Proportion native-exotic species"), xlab = "Easting (NZTM)", main="")
points(hfi.2016[hfi.2016$fenced.year=="Fenced", "prop.native"] ~ hfi.2016[hfi.2016$fenced.year=="Fenced", "x"], pch=19, cex=0.75, col="green") 
lines(x.x, predict(m1.x, list(x = x.x, fenced.year=rep("Fenced", length(x.x))), type="response", se.fit = FALSE), lwd=2, col="green",  lty=1)
points(hfi.2016[hfi.2016$fenced.year=="Unfenced", "prop.native"] ~ hfi.2016[hfi.2016$fenced.year=="Unfenced", "x"], pch="O", cex=0.75, col="brown") 
lines(x.x, predict(m1.x, list(x = x.x, fenced.year=rep("Unfenced", length(x.x))), type="response", se.fit = FALSE), lwd=2, col="brown",  lty=2)

plot(hfi.2016$prop.native ~ hfi.2016$y, ylim= c(min(0), max(1)), xlim=c(min(min(hfi.2016$y)),  max(max(hfi.2016$y))), type="n", bty="l", ylab = expression("Proportion native-exotic species"), xlab = "Northing (NZTM)", main="")
points(hfi.2016[hfi.2016$fenced.year=="Fenced", "prop.native"] ~ hfi.2016[hfi.2016$fenced.year=="Fenced", "y"], pch=19, cex=0.75, col="green") 
lines(x.y, predict(m1.y, list(y = x.y, fenced.year=rep("Fenced", length(x.y))), type="response", se.fit = FALSE), lwd=2, col="green",  lty=1)
points(hfi.2016[hfi.2016$fenced.year=="Unfenced", "prop.native"] ~ hfi.2016[hfi.2016$fenced.year=="Unfenced", "y"], pch="O", cex=0.75, col="brown") 
lines(x.y, predict(m1.y, list(y = x.y, fenced.year=rep("Unfenced", length(x.y))), type="response", se.fit = FALSE), lwd=2, col="brown",  lty=2)

plot(hfi.2016$prop.native ~ hfi.2016$water.deficit, ylim= c(min(0), max(1)), xlim=c(min(0),  max(max(hfi.2016$water.deficit))), type="n", bty="l", ylab = expression("Proportion native-exotic species"), xlab = "Water deficit (mm)")
points(hfi.2016[hfi.2016$fenced.year=="Fenced", "prop.native"] ~ hfi.2016[hfi.2016$fenced.year=="Fenced", "water.deficit"], pch=19, cex=0.75, col="green") 
lines(x.def, predict(m1.def, list(water.deficit = x.def, fenced.year=rep("Fenced", length(x.def))), type="response", se.fit = FALSE), lwd=2, col="green",  lty=1)
points(hfi.2016[hfi.2016$fenced.year=="Unfenced", "prop.native"] ~ hfi.2016[hfi.2016$fenced.year=="Unfenced", "water.deficit"], pch="O", cex=0.75, col="brown") 
lines(x.def, predict(m1.def, list(water.deficit = x.def, fenced.year=rep("Unfenced", length(x.def))), type="response", se.fit = FALSE), lwd=2, col="brown",  lty=2)


#######################################################################################################################################################################
 ##############################################               Biomass in three tiers                          ########################################################




#######################################################################################################################################################################
##############################################################               Biomass                           ########################################################

hfi.dat <- read.csv("/home/sean/Documents/Molesworth/Data/HFISpeciesComposition.csv", sep = ",", header=TRUE, skip = 0)
rownames(hfi.dat) <- hfi.dat$plot.yr; comp <- hfi.dat[,c(25:334)]  ### produce a frame for compositionn data
survey.list <- c("haphazard.plot.list", "random.plot.list", "paired.plot.list")
haphazard.plot.list <- cbind("haphazard.plot.list", c(as.character(hfi.dat[hfi.dat$year==2006, "plot"]), "Saxton4"))
random.plot.list <- cbind("random.plot.list", as.character(hfi.dat[hfi.dat$year==2007, "plot"]))
paired.plot.list <- cbind("paired.plot.list", as.character(hfi.dat[hfi.dat$year==2008, "plot"]))
plot.list <- rbind(haphazard.plot.list, random.plot.list, paired.plot.list)



pdf("/home/sean/Documents/Molesworth/Graphs/Biomass2016.pdf")
par(mfrow=c(1,2), cex=1, mar=c(4,5,4,2)+0.1) #pty="s",  The default mar is ‘c(5, 4, 4, 2) +0.1’. ‘c(bottom, left, top, right)’
biomass.random <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/RandomPlotBiomass.csv", sep = ",", header=TRUE, skip = 0)
biomass.random <-  merge(biomass.random, hfi.dat[hfi.dat$year==2016, c("plot", "fenced.year", "x", "y", "soil.class","rain.class", "oversow.year","terrace", "water.deficit", "temperature", "lcdb", "lenz.soil", "dem.altitude", "dem.slope" )], by="plot", all.X=TRUE, all.Y=FALSE)
boxplot2(biomass.random$biomass.2016 ~ biomass.random$fenced, top=TRUE, ylim= c(min(0), max(round_any(max(biomass.random$biomass.2016), 1, ceiling))), frame.plot = FALSE,  outline=FALSE, ylab =  expression("Biomass 2016 (g m"^~2 ~")")) 
mtext("a) Random plots", side = 3, line=2, adj = 0, cex=1)

biomass.paired <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/PairedPlotBiomass2016.csv", sep = ",", header=TRUE, skip = 0)
biomass.paired <-  merge(biomass.paired, hfi.dat[hfi.dat$year==2016, c("plot", "fenced.year", "x", "y", "soil.class","rain.class", "oversow.year","terrace", "water.deficit", "temperature", "lcdb", "lenz.soil", "dem.altitude", "dem.slope" )], by="plot", all.X=TRUE, all.Y=FALSE)
boxplot2(biomass.paired$biomass.2016 ~ biomass.paired$fenced, top=TRUE, ylim= c(min(0), max(round_any(max(biomass.paired$biomass.2016), 1, ceiling))), frame.plot = FALSE,  outline=FALSE, ylab =  expression("Biomass 2016 (g m"^~2 ~")")) 
mtext("b) Paired plots", side = 3, line=2, adj = 0, cex=1)
dev.off()

pdf("/home/sean/Documents/Molesworth/Graphs/BiomassBaseline.pdf")
par(mfrow=c(1,2), cex=1, mar=c(4,5,4,2)+0.1) #pty="s",  The default mar is ‘c(5, 4, 4, 2) +0.1’. ‘c(bottom, left, top, right)’
biomass.random <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/RandomPlotBiomassBaseline.csv", sep = ",", header=TRUE, skip = 0)
biomass.random <-  merge(biomass.random, hfi.dat[hfi.dat$year==2016, c("plot", "fenced.year", "x", "y", "soil.class","rain.class", "oversow.year","terrace", "water.deficit", "temperature", "lcdb", "lenz.soil", "dem.altitude", "dem.slope" )], by="plot", all.X=TRUE, all.Y=FALSE)
boxplot2(biomass.random$biomass.2007 ~ biomass.random$fenced, top=TRUE, ylim= c(min(0), max(round_any(max(biomass.random$biomass.2007), 1, ceiling))), frame.plot = FALSE,  outline=FALSE, ylab =  expression("Biomass 2007 (g m"^~2 ~")")) 
mtext("a) Random plots 2007", side = 3, line=2, adj = 0, cex=1)

biomass.paired <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/PairedPlotBiomassBaseline.csv", sep = ",", header=TRUE, skip = 0)
biomass.paired <-  merge(biomass.paired, hfi.dat[hfi.dat$year==2016, c("plot", "fenced.year", "x", "y", "soil.class","rain.class", "oversow.year","terrace", "water.deficit", "temperature", "lcdb", "lenz.soil", "dem.altitude", "dem.slope" )], by="plot", all.X=TRUE, all.Y=FALSE)
boxplot2(biomass.paired$biomass.2008 ~ biomass.paired$fenced, top=TRUE, ylim= c(min(0), max(round_any(max(biomass.paired$biomass.2008), 1, ceiling))), frame.plot = FALSE,  outline=FALSE, ylab =  expression("Biomass 2008 (g m"^~2 ~")")) 
mtext("b) Paired plots 2008", side = 3, line=2, adj = 0, cex=1)
dev.off()


############# Mixed models on HFI
hfi.sum <- read.csv( "/home/sean/Documents/Molesworth2020/Analysis/Data/RandomPlotHFI.csv", sep = ",", header=TRUE, skip = 0)
hfi.sum.2007 <- (hfi.sum[, c("plot","native","form","annual","hfi.2007")]); hfi.sum.2007$year <- 2007;names(hfi.sum.2007) <- c("plot","native","form","annual","hfi", "year")
hfi.sum.2016 <- (hfi.sum[, c("plot","native","form","annual","hfi.2016")]); hfi.sum.2016$year <- 2016; names(hfi.sum.2016) <- c("plot","native","form","annual","hfi", "year")
hfi.sum <- rbind(hfi.sum.2007, hfi.sum.2016); hfi.sum[hfi.sum$form=="Tree", "form"] <- "Shrub"
hfi.sum <- merge(hfi.sum, (hfi.dat[hfi.dat$year==2016, c("plot", "fenced.year", "x","y","oversow.year","terrace","water.deficit","temperature","lcdb","lenz.soil","dem.altitude","dem.slope")]), by="plot", all.X=TRUE, all.Y=FALSE)


pdf(paste("/home/sean/Documents/Molesworth/Graphs/RandomPlotHFIforms.pdf", sep="" ))
plot(hfi.sum$hfi ~ hfi.sum$year, ylim= c(min(0), max(100)), xlim=c(min(2006),  max(2017)), type="n", bty="l", ylab = "HFI", xlab = "Year")
for (i in unique(hfi.sum$form)){
points(c("2007", "2016"), tapply(hfi.sum[hfi.sum$form == i & hfi.sum$fenced.year=="Fenced", "hfi"], hfi.sum[hfi.sum$form == i & hfi.sum$fenced.year=="Fenced", "year"], mean), pch=19, cex=0.5, type="b", lty=1)
points(c("2007", "2016"), tapply(hfi.sum[hfi.sum$form == i & hfi.sum$fenced.year=="Unfenced", "hfi"], hfi.sum[hfi.sum$form == i & hfi.sum$fenced.year=="Unfenced", "year"], mean), pch=19, cex=0.5, type="b", lty=2)
text(2016.6, mean(hfi.sum[hfi.sum$form == i & hfi.sum$fenced.year=="Unfenced" & hfi.sum$year==2016, "hfi"]), paste("Unfenced\n", i, sep=""), cex = 0.8, col = "black",  pos = 1, offset = -0.5)
}; dev.off()



pdf(paste("/home/sean/Documents/Molesworth/Graphs/RandomPlotHFIannuals.pdf", sep="" )) ##i <- "Yes"
plot(hfi.sum$hfi ~ hfi.sum$year, ylim= c(min(0), max(80)), xlim=c(min(2006),  max(2017)), type="n", bty="l", ylab = "HFI", xlab = "Year", main="Paired plot biomass (HFI summed)")
upper.se <- mean(hfi.sum[hfi.sum$annual == "Yes" & hfi.sum$fenced.year=="Fenced" & hfi.sum$year==2016, "hfi"]) + s.e.(hfi.sum[hfi.sum$annual == "Yes" & hfi.sum$fenced.year=="Fenced" & hfi.sum$year==2016, "hfi"])
lower.se <- mean(hfi.sum[hfi.sum$annual == "Yes" & hfi.sum$fenced.year=="Fenced" & hfi.sum$year==2016, "hfi"]) - s.e.(hfi.sum[hfi.sum$annual == "Yes" & hfi.sum$fenced.year=="Fenced" & hfi.sum$year==2016, "hfi"])
arrows(2016, upper.se, 2016, lower.se, length=0.025, angle=90, code=3, lwd = 0.25)
for (i in unique(hfi.sum$annual)){
points(c("2007", "2016"), tapply(hfi.sum[hfi.sum$annual == i & hfi.sum$fenced.year=="Fenced", "hfi"], hfi.sum[hfi.sum$annual == i & hfi.sum$fenced.year=="Fenced", "year"], mean), pch=19, cex=0.5, type="b", lty=1)
points(c("2007", "2016"), tapply(hfi.sum[hfi.sum$annual == i & hfi.sum$fenced.year=="Unfenced", "hfi"], hfi.sum[hfi.sum$annual == i & hfi.sum$fenced.year=="Unfenced", "year"], mean), pch=19, cex=0.5, type="b", lty=2)
text(2016.6, mean(hfi.sum[hfi.sum$annual == "Yes" & hfi.sum$fenced.year=="Fenced" & hfi.sum$year==2016, "hfi"]), "Fenced\nannuals", cex = 0.8, col = "black",  pos = 3, offset = 0.5)
text(2016.6, mean(hfi.sum[hfi.sum$annual == "" & hfi.sum$fenced.year=="Unfenced" & hfi.sum$year==2016, "hfi"]), "Unfenced\nperenials", cex = 0.8, col = "black",  pos = 4, offset = -1)
text(2016.6, mean(hfi.sum[hfi.sum$annual == "Yes" & hfi.sum$fenced.year=="Unfenced" & hfi.sum$year==2016, "hfi"]), "Unfenced\nannuals", cex = 0.8, col = "black",  pos = 4, offset = -1)
text(2007, mean(hfi.sum[hfi.sum$annual == "" & hfi.sum$fenced.year=="Fenced" & hfi.sum$year==2007, "hfi"]), "Fenced perrenials", cex = 0.8, col = "black",  pos = 1, offset = .5)

}
dev.off()

##m2 <-glmer(hfi  ~   year*form*native + (1|plot), data=hfi.sum, family=poisson)  ##

m1 <-lme(log(hfi)  ~   year*form*native -1,  random = ~ 1|plot, method="ML", data=hfi.sum)  ##
#m2 <-lme(log(hfi)  ~   year*form*native + fenced.year + dem.altitude + x + y, random = ~ 1|plot, method="ML", data=hfi.sum)  ## Srepwise procedure simplified model. No significan intercept - so removed
## m2 <-lme(log(hfi)  ~   year*form + native + fenced.year, random = ~ 1|plot, data=hfi.sum)  ##
## m3 <-lme(log(hfi)  ~   year + form + native + fenced.year, random = ~ 1|plot, data=hfi.sum)  ##
## m4 <-lme(log(hfi)  ~   year*native , random = ~ 1|plot, data=hfi.sum)  ##

hfi.form.table.random <- as.data.frame(round(coef(summary(m1))[,c(1:2,4:5)], digits=3)); 
hfi.form.table.random[,names(hfi.form.table.random)] <- lapply(names(hfi.form.table.random), function(x) as.character(hfi.form.table.random[,x]))
hfi.form.table.random[hfi.form.table.random=="0"] <-  "$<$0.001"
hfi.form.table.random$factors <- c("Year", "Forb", "Grasses", "Shrubs", "Native", "Grasses x year", "Shrubs x year", "Native x year", "Native grasses", "Native shrubs", "Native grasses x year", "Native shrubs x year") 
##    c("Intercept. Fenced, exotic forbs", "Year", "Grasses", "Shrubs", "Native", "Unfenced x year", "Grasses x year", "Shrubs x year", "Native x year", "Native grasses", "Native shrubs", "Native grasses x year", "Native shrubs x year") 
hfi.form.table.random <- hfi.form.table.random[,c(5,1:4)]

pdf(paste("/home/sean/Documents/Molesworth/Graphs/RandomPlotHFIformNative.pdf", sep="" ))
plot(hfi.sum$hfi ~ hfi.sum$year, ylim= c(min(0), max(80)), xlim=c(min(2006),  max(2017)), type="n", bty="l", ylab = "HFI", xlab = "Year")
for (i in unique(hfi.sum$form)){
points(c("2007", "2016"), tapply(hfi.sum[hfi.sum$form == i & hfi.sum$native=="Native", "hfi"], hfi.sum[hfi.sum$form == i & hfi.sum$native=="Native", "year"], mean), pch=19, cex=0.5, type="b", lty=1)
points(c("2007", "2016"), tapply(hfi.sum[hfi.sum$form == i & hfi.sum$native=="Exotic", "hfi"], hfi.sum[hfi.sum$form == i & hfi.sum$native=="Exotic", "year"], mean),  pch=19, cex=0.5, type="b", lty=2)
text(2016.6, mean(hfi.sum[hfi.sum$form == i & hfi.sum$native=="Exotic" & hfi.sum$year==2016, "hfi"]), paste("Exotic\n", i, sep=""), cex = 0.8, col = "black",  pos = 1, offset = -0.5)
text(2006.4, mean(hfi.sum[hfi.sum$form == i & hfi.sum$native=="Native" & hfi.sum$year==2007, "hfi"]),  i, cex = 0.8, col = "black",  pos = 1, offset = -0.5)
}; dev.off()

##Paired plots
hfi.paired.sum <- read.csv( "/home/sean/Documents/Molesworth2020/Analysis/Data/PairedPlotHFI.csv", sep = ",", header=TRUE, skip = 0)
hfi.paired.sum.2008 <- (hfi.paired.sum[, c("plot","native","form","annual","hfi.2008")]); hfi.paired.sum.2008$year <- 2008; names(hfi.paired.sum.2008) <- c("plot","native","form","annual","hfi", "year")
hfi.paired.sum.2016 <- (hfi.paired.sum[, c("plot","native","form","annual","hfi.2016")]); hfi.paired.sum.2016$year <- 2016; names(hfi.paired.sum.2016) <- c("plot","native","form","annual","hfi", "year")
hfi.paired.sum <- rbind(hfi.paired.sum.2008, hfi.paired.sum.2016); hfi.paired.sum[hfi.paired.sum$form=="Tree", "form"] <- "Shrub"; hfi.paired.sum$form <- as.character(hfi.paired.sum$form)
hfi.paired.sum <- merge(hfi.paired.sum, (hfi.dat[hfi.dat$year==2016, c("plot", "fenced.year", "x","y","oversow.year","terrace","water.deficit","temperature","lcdb","lenz.soil","dem.altitude","dem.slope")]), by="plot", all.X=TRUE, all.Y=FALSE)

hfi.paired.sum$plot <-  as.character(hfi.paired.sum$plot); hfi.paired.sum$fenced <- substr(hfi.paired.sum$plot, (nchar(hfi.paired.sum$plot))-1, nchar(hfi.paired.sum$plot))
hfi.paired.sum$plot <- substr(hfi.paired.sum$plot, 1, (nchar(hfi.paired.sum$plot))-3)




## rownames(npp.table.random) <- NULL; names(npp.table.random) <- c("Random plot NPP", "sem", "P"); 
## npp.table.random$Effect <- c("Intercept - Fenced", "Unfenced", "Altitude (m.a.s.l)", "Year$\\times$ Unfenced")    
## table <- xtable(npp.table.random[,c(4,1:3)], label= "RandomNPP", caption = "Results from a GLM on NPP from randomly located plots established in 2007 and remeasured in 2016.") 
## align(table) <-   c("l", "l","r @{$\\pm$}", "l","l")
## print(table, floating = TRUE, include.rownames=FALSE, sanitize.text.function = function(x){x}, scalebox = 0.85)





pdf(paste("/home/sean/Documents/Molesworth/Graphs/PairedPlotHFIforms.pdf", sep="" ))
plot(hfi.paired.sum$hfi ~ hfi.paired.sum$year, ylim= c(min(0), max(100)), xlim=c(min(2008),  max(2017)), type="n", bty="l", ylab = "HFI", xlab = "Year")
upper.se <- mean(hfi.paired.sum[hfi.paired.sum$form == "Graminoid" & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2016, "hfi"]) + s.e.(hfi.paired.sum[hfi.paired.sum$form == "Graminoid" & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2016, "hfi"])
lower.se <- mean(hfi.paired.sum[hfi.paired.sum$form == "Graminoid" & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2016, "hfi"]) - s.e.(hfi.paired.sum[hfi.paired.sum$form == "Graminoid" & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2016, "hfi"])
arrows(2016, upper.se, 2016, lower.se, length=0.025, angle=90, code=3, lwd = 0.25)
for (i in c("Forb", "Graminoid", "Shrub")){
points(c("2008", "2016"), tapply(hfi.paired.sum[hfi.paired.sum$form == i & hfi.paired.sum$fenced.year=="Fenced", "hfi"], hfi.paired.sum[hfi.paired.sum$form == i & hfi.paired.sum$fenced.year=="Fenced", "year"], mean), pch=19, cex=0.5, type="b", lty=1)
points(c("2008", "2016"),tapply(hfi.paired.sum[hfi.paired.sum$form == i & hfi.paired.sum$fenced.year=="Unfenced", "hfi"], hfi.paired.sum[hfi.paired.sum$form == i & hfi.paired.sum$fenced.year=="Unfenced", "year"], mean) , pch=19, cex=0.5, type="b", lty=2)
text(2016.6, mean(hfi.paired.sum[hfi.paired.sum$form == i & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2016, "hfi"]), paste("Fenced\n", i, sep=""), cex = 0.8, col = "black",  pos = 1, offset = -0.5)
}
dev.off()


m1 <-lme(hfi  ~  year * form * fenced, random =~1|plot, data=hfi.paired.sum)  ##
#m2 <-lme(hfi  ~  year + form + x + terrace + water.deficit + dem.altitude +dem.slope, random =~year|form/plot, data=hfi.paired.sum)  ##
hfi.form.table.paired <- as.data.frame(round(coef(summary(m1))[,c(1:2,4:5)], digits=3)); 
hfi.form.table.paired[,names(hfi.form.table.paired)] <- lapply(names(hfi.form.table.paired), function(x) as.character(hfi.form.table.paired[,x]))
hfi.form.table.paired[hfi.form.table.paired=="0"] <-  "$<$0.001"
hfi.form.table.paired$factors <- c("Intercept. Fenced forbs", "Year", "Grasses", "Shrubs", "Unfenced", "Grasses x year", "Shrubs x year", "Unfenced x year", "Unfenced x grasses", "Unfenced x shrubs", "Unfenced x grasses x year", "Unfenced x shrubs x year") 
hfi.form.table.paired <- hfi.form.table.paired[,c(5,1:4)]



pdf(paste("/home/sean/Documents/Molesworth/Graphs/PairedPlotHFIannuals.pdf", sep="" )) ##i <- "Yes"
plot(hfi.paired.sum$hfi ~ hfi.paired.sum$year, ylim= c(min(0), max(80)), xlim=c(min(2008),  max(2017)), type="n", bty="l", ylab = "HFI", xlab = "Year", main="Paired plot biomass (HFI summed)")
## upper.se <- mean(hfi.paired.sum[hfi.paired.sum$annual == "Yes" & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2016, "hfi"]) + s.e.(hfi.paired.sum[hfi.paired.sum$annual == "Yes" & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2016, "hfi"])
## lower.se <- mean(hfi.paired.sum[hfi.paired.sum$annual == "Yes" & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2016, "hfi"]) - s.e.(hfi.paired.sum[hfi.paired.sum$annual == "Yes" & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2016, "hfi"])
## arrows(2016, upper.se, 2016, lower.se, length=0.025, angle=90, code=3, lwd = 0.25)
for (i in  unique(hfi.sum$annual)){
points(c("2008", "2016"), tapply(hfi.paired.sum[hfi.paired.sum$annual == i & hfi.paired.sum$fenced.year=="Fenced", "hfi"], hfi.paired.sum[hfi.paired.sum$annual == i & hfi.paired.sum$fenced.year=="Fenced", "year"], mean), pch=19, cex=0.5, type="b", lty=1)
points(c("2008", "2016"), tapply(hfi.paired.sum[hfi.paired.sum$annual == i & hfi.paired.sum$fenced.year=="Unfenced", "hfi"], hfi.paired.sum[hfi.paired.sum$annual == i & hfi.paired.sum$fenced.year=="Unfenced", "year"], mean), pch=19, cex=0.5, type="b", lty=2)
text(2009, mean(hfi.paired.sum[hfi.paired.sum$annual == "Yes" & hfi.paired.sum$fenced.year=="Fenced" & hfi.paired.sum$year==2008, "hfi"]), "Fenced annuals", cex = 0.8, col = "black",  pos = 3, offset = 0.25)

}
dev.off()

save(hfi.form.table.random, hfi.form.table.paired, file="/home/sean/Documents/Molesworth/Data/Objects/HFIform.Rdata") ###


#######################################################################################################################################################################
##############################################################               NPP                           ########################################################

pdf("/home/sean/Documents/Molesworth/Graphs/NPP.pdf")
par(mfrow=c(1,2), cex=1, mar=c(4,5,4,2)+0.1) #pty="s",  The default mar is ‘c(5, 4, 4, 2) +0.1’. ‘c(bottom, left, top, right)’
npp.random <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/RandomPlotNPP.csv", sep = ",", header=TRUE, skip = 0)
npp.random <-  merge(npp.random, hfi.dat[hfi.dat$year==2016, c("plot", "fenced.year")], by="plot", all.X=TRUE, all.Y=TRUE)
boxplot2(npp.random$npp ~ npp.random$fenced.year, top=TRUE, ylim= c(min(0), max(1200)), frame.plot = FALSE,  outline=FALSE, ylab =  expression("Annual NPP (g m"^~2 ~")")) #
mtext("a) Random plots", side = 3, line=2, adj = 0, cex=1)
abline(h=308); text(x= 1.5, y= 320, labels= "Global study mean")
npp.paired <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/PairedPlotNPP.csv", sep = ",", header=TRUE, skip = 0)
boxplot2(npp.paired$npp ~ npp.paired$fenced, top=TRUE, ylim= c(min(0), max(1200)), frame.plot = FALSE,  outline=FALSE, yaxt="n", names = c("Fenced", "Unfenced")) #
mtext("b) Paired plots", side = 3, line=2, adj = 0, cex=1); abline(h=308)
dev.off()

random.npp.summary <- round(c(mean(npp.random$npp), s.e.(npp.random$npp)), digits=0)

###Models for paired plots
## m1 <- glm(npp ~ site + fenced * dem.altitude + terrace + x + y + water.deficit, data = npp.paired, family=poisson) #confint(m1, level=0.95)
## m1 <- glm(npp ~ fenced ,  data = npp.paired, family=poisson) #confint(m1, level=0.95)
m1.random <- glm(npp ~ fenced.year * dem.altitude ,  data = npp.random, family=poisson) #confint(m1, level=0.95)
m1.x <- glm(npp ~ fenced * x ,  data = npp.paired, family=poisson) #confint(m1, level=0.95)
m1.y <- glm(npp ~ fenced * y ,  data = npp.paired, family=poisson) #confint(m1, level=0.95)
m1.def <- glm(npp ~ fenced * water.deficit  ,  data = npp.paired, family=poisson) #confint(m1, level=0.95)
npp.paired$site <-  as.character(npp.paired$site)
## m1.lme <-  lme(log(npp) ~  fenced,  random =~1|site, data = npp.paired) #confint(m1, level=0.95)
 x.alt= seq(0,2000,by=10); x.def=seq(0,200); x.alt.high= seq(1000,2000,by=5); x.x = seq(min(npp.paired$x), max(npp.paired$x), by=100); x.y = seq(min(npp.paired$y), max(npp.paired$y), by=100)
npp.table.random <- as.data.frame(round(coef(summary(m1.random))[,c(1:2,4)], digits=3)); 
npp.table.random[,names(npp.table.random)] <- lapply(names(npp.table.random), function(x) as.character(npp.table.random[,x]))
npp.table.random[npp.table.random=="0"] <-  "$<$0.001"


rownames(npp.table.random) <- NULL; names(npp.table.random) <- c("Random plot NPP", "sem", "P"); 
npp.table.random$Effect <- c("Intercept - Fenced", "Unfenced", "Altitude (m.a.s.l)", "Year$\\times$ Unfenced")    
table <- xtable(npp.table.random[,c(4,1:3)], label= "RandomNPP", caption = "Results from a GLM on NPP from randomly located plots established in 2007 and remeasured in 2016.") 
align(table) <-   c("l", "l","r @{$\\pm$}", "l","l")
print(table, floating = TRUE, include.rownames=FALSE, sanitize.text.function = function(x){x}, scalebox = 0.85)


m1.paired <-  glm(npp ~  fenced * dem.altitude ,  data = npp.paired, family=poisson) #confint(m1, level=0.95)
npp.table.paired <- as.data.frame(round(coef(summary(m1.paired))[,c(1:2,4)], digits=3)); 
npp.table.paired[,names(npp.table.paired)] <- lapply(names(npp.table.paired), function(x) as.character(npp.table.paired[,x]))
npp.table.paired[npp.table.paired=="0"] <-  "$<$0.001"


rownames(npp.table.paired) <- NULL; names(npp.table.paired) <- c("Paired plot NPP", "sem", "P"); 
npp.table.paired$Effect <- c("Intercept - Fenced", "Unfenced", "Altitude (m.a.s.l)", "Year$\\times$ Unfenced")    
table <- xtable(npp.table.paired[,c(4,1:3)], label= "PairedNPP", caption = "Results from a GLM on NPP from paired plots established in 2008 and remeasured in 2016.") 
align(table) <-   c("l", "l","r @{$\\pm$}", "l","l")
print(table, floating = TRUE, include.rownames=FALSE, sanitize.text.function = function(x){x}, scalebox = 0.85)




m1.alt <-  glm(npp ~ fenced + dem.altitude, data = npp.paired, family=poisson)

## m1 <- glmer(npp ~ fenced * poly(dem.altitude, 2) + x + y + terrace + water.deficit + (1|site), data = npp.paired, family=poisson) #confint(m1, level=0.95)
## predict(m1, type = "response"))


pdf("/home/sean/Documents/Molesworth2020/Analysis/Graphs/PairedPlotNPPvsAlt.pdf")
plot(npp.paired$npp ~ npp.paired$dem.altitude, ylim= c(min(0), max(2000)), xlim=c(min(800),  max(1400)), type="n", bty="l", ylab = expression("Annual NPP (g m"^~2 ~")"), xlab = "Altitude (m a.s.l.)", main="")
points(npp.paired[npp.paired$fenced=="Fe", "npp"] ~ npp.paired[npp.paired$fenced=="Fe", "dem.altitude"], pch=19, cex=0.75, col="green") 
lines(x.alt, predict(m1.alt, list(dem.altitude = x.alt, fenced=rep("Fe", length(x.alt))), type="response", se.fit = FALSE), lwd=2, col="green",  lty=1)
points(npp.paired[npp.paired$fenced=="Un", "npp"] ~ npp.paired[npp.paired$fenced=="Un", "dem.altitude"], pch="O", cex=0.75, col="brown") 
lines(x.alt, predict(m1.alt, list(dem.altitude = x.alt, fenced=rep("Un", length(x.alt))), type="response", se.fit = FALSE), lwd=2, col="brown",  lty=2)
legend("topleft",c("Fenced","Unfenced"), col=c("green","brown"), lwd=3, lty=c(1,2))
dev.off()
pdf("/home/sean/Documents/Molesworth2020/Analysis/Graphs/PairedPlotNPPvsWaterDeficit.pdf")
plot(npp.paired$npp ~ npp.paired$water.deficit, ylim= c(min(0), max(2000)), xlim=c(min(0),  max(100)), type="n", bty="l", ylab = expression("Annual NPP (g m"^~2 ~")"), xlab = "Water deficit (mm)", main="")
points(npp.paired[npp.paired$fenced=="Fe", "npp"] ~ npp.paired[npp.paired$fenced=="Fe", "water.deficit"], pch=19, cex=0.75, col="green") 
lines(x.def, predict(m1.def, list(water.deficit = x.def, fenced=rep("Fe", length(x.def))), type="response", se.fit = FALSE), lwd=2, col="green",  lty=1)
points(npp.paired[npp.paired$fenced=="Un", "npp"] ~ npp.paired[npp.paired$fenced=="Un", "water.deficit"], pch="O", cex=0.75, col="brown") 
lines(x.def, predict(m1.def, list(water.deficit = x.def, fenced=rep("Un", length(x.def))), type="response", se.fit = FALSE), lwd=2, col="brown",  lty=2)
legend("topleft",c("Fenced","Unfenced"), col=c("green","brown"), lwd=3, lty=c(1,2))
dev.off()
pdf("/home/sean/Documents/Molesworth2020/Analysis/Graphs/PairedPlotNPPvsEast.pdf")
plot(npp.paired$npp ~ npp.paired$x, ylim= c(min(0), max(2000)), xlim=c(min(min(npp.paired$x)),  max(max(npp.paired$x))), type="n", bty="l", ylab = expression("Annual NPP (g m"^~2 ~")"), xlab = "Easting (NZTM 2000)", main="")
points(npp.paired[npp.paired$fenced=="Fe", "npp"] ~ npp.paired[npp.paired$fenced=="Fe", "x"], pch=19, cex=0.75, col="green") 
lines(x.x, predict(m1.x, list(x = x.x, fenced=rep("Fe", length(x.x))), type="response", se.fit = FALSE), lwd=2, col="green",  lty=1)
points(npp.paired[npp.paired$fenced=="Un", "npp"] ~ npp.paired[npp.paired$fenced=="Un", "x"], pch="O", cex=0.75, col="brown") 
lines(x.x, predict(m1.x, list(x = x.x, fenced=rep("Un", length(x.x))), type="response", se.fit = FALSE), lwd=2, col="brown",  lty=2)
legend("topleft",c("Fenced","Unfenced"), col=c("green","brown"), lwd=3, lty=c(1,2))
dev.off()
pdf("/home/sean/Documents/Molesworth2020/Analysis/Graphs/PairedPlotNPPvsNorth.pdf")
plot(npp.paired$npp ~ npp.paired$y, ylim= c(min(0), max(2000)), xlim=c(min(min(npp.paired$y)),  max(max(npp.paired$y))), type="n", bty="l", ylab = expression("Annual NPP (g m"^~2 ~")"), xlab = "Northing (NZTM 2000)", main="")
points(npp.paired[npp.paired$fenced=="Fe", "npp"] ~ npp.paired[npp.paired$fenced=="Fe", "y"], pch=19, cex=0.75, col="green") 
lines(x.y, predict(m1.y, list(y = x.y, fenced=rep("Fe", length(x.y))), type="response", se.fit = FALSE), lwd=2, col="green",  lty=1)
points(npp.paired[npp.paired$fenced=="Un", "npp"] ~ npp.paired[npp.paired$fenced=="Un", "y"], pch="O", cex=0.75, col="brown") 
lines(x.y, predict(m1.y, list(y = x.y, fenced=rep("Un", length(x.y))), type="response", se.fit = FALSE), lwd=2, col="brown",  lty=2)
legend("topleft",c("Fenced","Unfenced"), col=c("green","brown"), lwd=3, lty=c(1,2))
dev.off()


##Models for random plots
global.excl <- read.csv("/home/sean/Documents/Molesworth2020/Analysis/Data/GlobalExclosureNPPtable.csv", sep = ",", header=TRUE, skip = 5)
# Stepwise Regression library(MASS)
m1 <-lm(npp  ~  fenced.year + x + terrace + water.deficit + dem.altitude +dem.slope, data=npp.random)  ##
step.m1 <- stepAIC(m1, direction="both")
step.m1$anova # display results 
m2 <-lm(npp  ~  fenced.year +  water.deficit, data=npp.random)  ##

npp.random <- npp.random[npp.random$npp < 2000,] ##remove two plots with high npp R1 and I2 - both in damp spots

m1 <- glm(npp ~ fenced.year * poly(dem.altitude, 2) + x + y + terrace + water.deficit, data = npp.random, family=poisson) #confint(m1, level=0.95)
step(m1,  direction = c("forward"))

m1 <- glm(npp ~  poly(dem.altitude, 2) * fenced.year, data = npp.random, family=poisson)
#m1.a <- glm(npp ~  I(dem.altitude^2), data = npp.random, family=poisson)


pdf("/home/sean/Documents/Molesworth2020/Analysis/Graphs/RandomPlotNPPvsAlt.pdf") ##i <- "Yes"
plot(npp.random$npp ~ npp.random$dem.altitude, ylim= c(min(0), max(2000)), xlim=c(min(500),  max(1500)), type="n", bty="l", ylab = expression("Annual NPP (g m"^~2 ~")"), xlab = "Altitude (m)", main="")
points(npp.random[npp.random$fenced.year=="Fenced", "npp"] ~ npp.random[npp.random$fenced.year=="Fenced", "dem.altitude"], pch=19, cex=0.75, col="green") 
lines(x.alt.high, predict(m1, list(dem.altitude = x.alt.high, fenced.year=rep("Fenced", length(x.alt.high))), type="response", se.fit = FALSE), lwd=2, col="green",  lty=1)
points(npp.random[npp.random$fenced.year=="Unfenced", "npp"] ~ npp.random[npp.random$fenced.year=="Unfenced", "dem.altitude"], pch="O", cex=0.75, col="brown") 
lines(x.alt, predict(m1, list(dem.altitude = x.alt, fenced.year=rep("Unfenced", length(x.alt))), type="response", se.fit = FALSE), lwd=2, col="brown",  lty=2)
legend("topleft",c("Fenced","Unfenced"), col=c("green","brown"), lwd=3, lty=c(1,2))
dev.off()

m1 <- gls(log(npp) ~ fenced.year + x + terrace + water.deficit + dem.altitude, data = npp.random)
vario <- Variogram(m1, form = ~x + y, resType = "pearson")
plot(vario, smooth = TRUE)#, ylim = c(0, 1.2))
m2 <- gls(npp ~ fenced.year + x + terrace + water.deficit + dem.altitude, correlation = corExp(form = ~x + y, nugget = TRUE),data = npp.random)## try also: c(corExp, corGaus, corSpher, corLin, corRatio
vario.fitted <- Variogram(m2, form = ~x + y, resType = "pearson")
plot(vario.fitted)


## The units for water deficit are in mm, higher values are areas that have a larger deficit. 
m1 <-lme(log(npp)  ~  fenced + x + terrace + water.deficit + dem.altitude +dem.slope, random=~1|site, data=npp.paired)  ##
anova.lme(m1, type="sequential", adjustSigma = TRUE)
m2 <- aov(npp  ~   site + fenced %in% site, data=npp.paired)  ##fenced nested within site
m3 <- aov(npp  ~   fenced + Error(site/fenced), data=npp.paired)  ##fenced fixed effect, site random effect. Need balanced model to use aov
m4 <- glm(npp  ~   site + fenced, data=npp.paired, family=poisson )  ##randomised block with 
m4 <- glm(npp  ~   site + fenced + dem.altitude, data=npp.paired, family=poisson )  ##randomised block with 
exp(coef(m4)) 

m4 <- glm(npp  ~   dem.altitude, data=npp.paired, family=poisson )  ##randomised block with 
x = seq(500, 2000) 
plot(npp.random$npp ~ npp.random$dem.altitude)
lines(predict(m4,list(dem.altitude = x), type="response", se.fit = FALSE),lwd=2,col="orange")


m5 <- glm(npp  ~   water.deficit, data=npp.paired, family=poisson )  ##randomised block with 
plot(npp.random$npp ~ npp.random$water.deficit)
x = seq(0,200) 
lines(predict(m5,list(water.deficit = x),type="response"),lwd=2,col="orange")


## m5 <- glm(npp  ~   fenced + site + terrace + water.deficit + dem.altitude, data=npp.paired)  ## ANCOVA
## vario <- Variogram(m5, form = ~x + y, resType = "pearson")
## plot(vario, smooth = TRUE)#, ylim = c(0, 1.2))
m6 <- gls(npp  ~   fenced + site + terrace + water.deficit + dem.altitude, correlation = corExp(form = ~x + y, nugget = TRUE), data=npp.paired)  ## ANCOVA
vario.fitted <- Variogram(m6, form = ~x + y, resType = "pearson")
plot(vario.fitted)


## coef() or summary() are on the scale of the link function, not the untransformed predictor variables.  Use the inverse of the link function to get parameter values back on the scale of x, or use the function predict with the type="response"
## attribute functions (coef, confint, resid, fitted, etc)






##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
                               ##### ##### ##### ##### ##### DCA including all species   ##### ##### ##### ##### #####
for(i in c(survey.list)) {
temp.plot.list <- plot.list[plot.list[,1]==i,2]
  
comp <- hfi.dat[which(hfi.dat$plot %in% temp.plot.list), c(28:338)]
  
##paste(toupper(substr(survey.list,1,1)), substr(survey.list, 2, regexpr(".plot", survey.list)-1), sep="")
title <- paste(toupper(substr(i,1,1)), substr(i, 2, regexpr(".plot", i)[1]-1), " Plots", sep="")
  
    
ord <- decorana(comp) ## DCA on all species              
dca.site <- as.data.frame(ord$rproj[, c(1:2)])
dca.site$plot.yr <- rownames(dca.site)
dca.site$plot <- substr(rownames(dca.site), 1, nchar(rownames(dca.site))-5)
dca.site$year <- substr(rownames(dca.site), nchar(rownames(dca.site))-3, nchar(rownames(dca.site)))
dca.spp <-  as.data.frame(ord$cproj[, c(1:2)])
dca.spp$species <- rownames(dca.spp)

### Plots Saxton1, Saxton2, Saxton4, R1, R2 may have missed out on being identified as on a terrace because of coarse pixels near the Molesworth boundary.

dca.site <- merge(dca.site, hfi.dat[,c("plot.yr","fenced.year","x","y","soil.class","rain.class","oversow.year","terrace","water.deficit","temperature","lcdb","lenz.soil","dem.altitude","dem.slope")], by="plot.yr", all.X=TRUE, all.Y=TRUE)

dca1.site.mean <- as.data.frame.table(tapply(dca.site$DCA1, list(dca.site$fenced.year, dca.site$year), mean))
dca1.site.sem <- as.data.frame.table(tapply(dca.site$DCA1, list(dca.site$fenced.year, dca.site$year), s.e.))
dca1 <- (cbind(dca1.site.mean, dca1.site.sem))[,c(1:3,6)]
names(dca1) <- c("fenced", "year", "DCA1.mean", "DCA1.sem")
dca2.site.mean <- as.data.frame.table(tapply(dca.site$DCA2, list(dca.site$fenced.year, dca.site$year), mean))
dca2.site.sem <- as.data.frame.table(tapply(dca.site$DCA2, list(dca.site$fenced.year, dca.site$year), s.e.))
dca2 <- (cbind(dca2.site.mean, dca2.site.sem))[,c(1:3,6)]
names(dca2) <- c("fenced", "year", "DCA2.mean", "DCA2.sem")
dca.site <- cbind(dca1, dca2)[,-c(5,6)]

pdf(paste("/home/sean/Documents/Molesworth2020/Analysis/Graphs/DCA", sub(" ","", title),  ".pdf", sep="" ))
    
plot(dca.site$DCA2.mean ~ dca.site$DCA1.mean, ylim=c(min(round_any(min(dca.site$DCA2.mean), 0.5, floor)), max(round_any(max(c(dca.site$DCA2.mean)), 0.5, ceiling))),
      xlim=c(min(round_any(min(dca.site$DCA1.mean), 0.5, floor)), max(round_any(max(c(dca.site$DCA1.mean)), 0.5, ceiling))),
      type="n", bty="l", ylab = "Mean DCA 2", xlab = "Mean DCA1", main=paste(title, "\nDecorana on HFI - all species", sep=""))
points(dca.site[dca.site$fenced=="Fenced", "DCA2.mean"] ~ dca.site[dca.site$fenced=="Fenced", "DCA1.mean"], pch=19, cex=0.5, type="b")
points(dca.site[dca.site$fenced=="Unfenced", "DCA2.mean"] ~ dca.site[dca.site$fenced=="Unfenced", "DCA1.mean"], pch="O", cex=0.5, type="b", lty=2)
text(dca.site[dca.site$fenced=="Unfenced", "DCA2.mean"] ~ dca.site[dca.site$fenced=="Unfenced", "DCA1.mean"], labels  = dca.site[dca.site$fenced=="Unfenced", "year"], cex = 0.8, col = "black",  pos = 1, offset = 0.5)
##text(dca.site[dca.site$fenced=="Unfenced", "DCA2.mean"] ~ dca.site[dca.site$fenced=="Fenced", "DCA1.mean"], labels  = dca.site$year, cex = 0.8, col = "black",  pos = 3, offset = 0.5)
arrows(dca.site$DCA1.mean, (dca.site$DCA2.mean + dca.site$DCA2.sem), dca.site$DCA1.mean, (dca.site$DCA2.mean - dca.site$DCA2.sem), length=0.025, angle=90, code=3, lwd = 0.25)
arrows((dca.site$DCA1.mean + dca.site$DCA1.sem), dca.site$DCA2.mean,  (dca.site$DCA1.mean - dca.site$DCA1.sem), length=0.025, angle=90, code=3, lwd = 0.25)

dev.off()

}

##i <- "haphazard.plot.list"
##### DCA only including common species

##i <- "haphazard.plot.list"
##### DCA only including common species

comp.common <- hfi.dat[hfi.dat$year==2016, c(28:338)]
for(i in names(comp.common)){comp.common[comp.common[names(comp.common)==i] != 0, names(comp.common)==i] <- 1} ## produce a matrix with 1 for occurance
col.sums.comp.common <- as.data.frame(colSums(comp.common)); col.sums.comp.common$spp <- rownames(col.sums.comp.common)
common.spp <- col.sums.comp.common[(col.sums.comp.common[1] > length(rownames(comp.common))/10)==TRUE, "spp"]  ##Produce a list of common species. 10 for species occuring in more than 10 plots in 2016 random plot survey
very.common.spp <- col.sums.comp.common[(col.sums.comp.common[1] > 80)==TRUE, "spp"]  ##Produce a list of common species. 5 for species occuring in at least half of  2016 random plots

##i <- "haphazard.plot.list"
for(i in survey.list) {
title <- paste(toupper(substr(i,1,1)), substr(i, 2, regexpr(".plot", i)[1]-1), " Plots", sep="")
temp.plot.list <- plot.list[plot.list[,1]==i,2]
  
common.spp.hfi <- hfi.dat[which(hfi.dat$plot %in% temp.plot.list), common.spp];
    
####Decorana ordination and graph common species
ord <- decorana(common.spp.hfi) ## DCA on common species              
dca.site <- as.data.frame(ord$rproj[, c(1:2)])
dca.site$plot.yr <- rownames(dca.site)
dca.site$plot <- substr(rownames(dca.site), 1, nchar(rownames(dca.site))-5)
dca.site$year <- substr(rownames(dca.site), nchar(rownames(dca.site))-3, nchar(rownames(dca.site)))
dca.spp <-  as.data.frame(ord$cproj[, c(1:2)])
dca.spp$species <- rownames(dca.spp)

dca.site <- merge(dca.site, hfi.dat[,c("plot.yr","fenced.year","x","y","soil.class","rain.class","oversow.year","terrace","water.deficit","temperature","lcdb","lenz.soil","dem.altitude","dem.slope")], by="plot.yr", all.X=TRUE, all.Y=TRUE)
dca.site$year <- as.integer(dca.site$year)
m.dca1 <- lme(DCA1 ~ year + fenced.year, random= ~1|plot, data=dca.site)
table <- as.data.frame(round(coef(summary(m.dca1))[,c(1:2,4:5)], digits=3)); 
table[,names(table)] <- lapply(names(table), function(x) as.character(table[,x]))
table[table=="0"] <-  "$<$0.001"
assign(paste("dca1.table.", substr(i, 1, nchar(i)-10), sep=""), table)

m.dca2 <- lme(DCA2 ~ year + fenced.year, random= ~1|plot, data=dca.site)
table <- as.data.frame(round(coef(summary(m.dca1))[,c(1:2,4:5)], digits=3)); 
table[,names(table)] <- lapply(names(table), function(x) as.character(table[,x]))
table[table=="0"] <-  "$<$0.001"
assign(paste("dca2.table.", substr(i, 1, nchar(i)-10), sep=""), table)
    
dca1.site.mean <- as.data.frame.table(tapply(dca.site$DCA1, list(dca.site$fenced.year, dca.site$year), mean))
dca1.site.sem <- as.data.frame.table(tapply(dca.site$DCA1, list(dca.site$fenced.year, dca.site$year), s.e.))
dca1 <- (cbind(dca1.site.mean, dca1.site.sem))[,c(1:3,6)]
names(dca1) <- c("fenced", "year", "DCA1.mean", "DCA1.sem")
dca2.site.mean <- as.data.frame.table(tapply(dca.site$DCA2, list(dca.site$fenced.year, dca.site$year), mean))
dca2.site.sem <- as.data.frame.table(tapply(dca.site$DCA2, list(dca.site$fenced.year, dca.site$year), s.e.))
dca2 <- (cbind(dca2.site.mean, dca2.site.sem))[,c(1:3,6)]
names(dca2) <- c("fenced", "year", "DCA2.mean", "DCA2.sem")
dca.site <- cbind(dca1, dca2)[,-c(5,6)]

pdf(paste("/home/sean/Documents/Molesworth2020/Analysis/Graphs/DCAcommonSpecies", sub(" ","", title),  ".pdf", sep="" ))
    
plot(dca.site$DCA2.mean ~ dca.site$DCA1.mean, ylim=c(min(round_any(min(dca.site$DCA2.mean), 0.5, floor)), max(round_any(max(c(dca.site$DCA2.mean)), 0.5, ceiling))),
      xlim=c(min(round_any(min(dca.site$DCA1.mean), 0.5, floor)), max(round_any(max(c(dca.site$DCA1.mean)), 0.5, ceiling))),
      type="n", bty="l", ylab = "Mean DCA 2", xlab = "Mean DCA1", main=paste(title, "\nDecorana on HFI - 2016 common species only", sep=""))
points(dca.site[dca.site$fenced=="Fenced", "DCA2.mean"] ~ dca.site[dca.site$fenced=="Fenced", "DCA1.mean"], pch=19, cex=0.5, type="b")
points(dca.site[dca.site$fenced=="Unfenced", "DCA2.mean"] ~ dca.site[dca.site$fenced=="Unfenced", "DCA1.mean"], pch="O", cex=0.5, type="b", lty=2)
text(dca.site[dca.site$fenced=="Unfenced", "DCA2.mean"] ~ dca.site[dca.site$fenced=="Unfenced", "DCA1.mean"], labels  = dca.site[dca.site$fenced=="Unfenced", "year"], cex = 0.8, col = "black",  pos = 1, offset = 0.5)
##text(dca.site[dca.site$fenced=="Unfenced", "DCA2.mean"] ~ dca.site[dca.site$fenced=="Fenced", "DCA1.mean"], labels  = dca.site$year, cex = 0.8, col = "black",  pos = 3, offset = 0.5)
arrows(dca.site$DCA1.mean, (dca.site$DCA2.mean + dca.site$DCA2.sem), dca.site$DCA1.mean, (dca.site$DCA2.mean - dca.site$DCA2.sem), length=0.025, angle=90, code=3, lwd = 0.25)
arrows((dca.site$DCA1.mean + dca.site$DCA1.sem), dca.site$DCA2.mean,  (dca.site$DCA1.mean - dca.site$DCA1.sem), length=0.025, angle=90, code=3, lwd = 0.25)
dev.off()

}
save(dca1.table.haphazard, dca1.table.random, dca1.table.paired, dca2.table.haphazard, dca2.table.random, dca2.table.paired, file = "/home/sean/Documents/Molesworth/Data/Objects/HFIdca.Table.Rdata") ##




