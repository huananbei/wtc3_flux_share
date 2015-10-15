#- required libraries
library(HIEv)
library(mgcv)
library(scales)
library(gplots)
library(plotBy)
library(plantecophys)
library(magicaxis)
library(lubridate)
library(doBy)
library(Hmisc)
library(zoo)
library(hexbin)
library(nlme)
library(lsmeans)
library(car)

#- load the analysis and plotting functions that do all of the actual work
source("R/functions.R")

#- export flag. Set to "T" to create pdfs of figures in "output/", or "F" to suppress output.
export=T





#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- plot R/A with a leaf-scale model (Figure 1)
plotCUE_conceptual_fig(toexport=export,Tdew=10,Ca=400,Vcmax=100,Jmax=125,Tleaf=10:42,PPFD=1500)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Plot respiration vs. temperature curves for leaves and whole-tree chambers (Figure 2)

#- get and process the whole-canopy R vs. T datasets
fits.list <- return_Rcanopy_closed()
fits.mass <- fits.list[[1]]     #- tree-level data
fits.trt <- fits.list[[2]]      #- treatment averages

#- plot R vs. T (Figure 2)
plotRvsT_figure2(fits.mass=fits.mass,fits.trt=fits.trt,export=export)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- plot Rbranch and Rleaf measured at 15 degrees C (Figure 3)
plotRleafRbranch(export=export)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Read and process the hourly flux dataset. Plot Figures 4 and 5.

#- read in the hourly flux data
dat.hr <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130910-20140530_L2_V1.csv")
dat.hr$DateTime <- as.POSIXct(dat.hr$DateTime,format="%Y-%m-%d %H:%M:%S",tz="GMT")
dat.hr$Date <- as.Date(dat.hr$DateTime)

#- partition the net fluxes into GPP and Ra components. Specify the actitation energy and the number
#-  of prior nights to include for the estimate of the basal respiration rate (lagdates)
dat.hr.p <- partitionHourlyFluxCUE_arr(dat.hr.gf=dat.hr,Ea=57.69,lagdates=3)

#- plot an example week of partitioned fluxes
plotPartitionedFluxes(dat.hr.gf3=dat.hr.p,ch_toplot="C07",startDate="2014-3-22",endDate="2014-3-27",write=F)

#- get daily sums
cue.list <- returnCUE.day(dat=dat.hr.p) # get daily sums from hourly data
cue.day <- cue.list[[1]]                # extract chamber values on each day
cue.day.trt <- cue.list[[2]]            # extract treatment averages

#- plot met, Ra, GPP, and Ra/GPP data over time (Figure 4)
plotPAR_AirT_CUE_GPP_Ra(cue.day.trt=cue.day.trt,export=export,lwidth=2.75)

#- plot PAR and Temperaure dependence of GPP, Ra, and Ra/GPP (Figure 5)
plotGPP_Ra_CUE_metdrivers(cue.day=cue.day,export=export,shading=0.7)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Plot the temperature and PAR dependence of GPP per unit leaf area (Figure 6).
#- Note that this produces four separate graphs (panels a, b, and c, plus the legend).
#   These panels were manually combined to create Figure 6.
plotGPP_hex(dat=dat.hr.p,export=export,shading=0.7)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Plot the 5 diurnal observations of leaf-level photosynthesis and stomatal conductance.
#    Set printANOVAs to "T" to print ANOVAs for each date
plotAnet_met_diurnals(export=export,lsize=2,printANOVAs=F)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
















#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Statistical analysis of Ra, GPP, and Ra/GPP (Table 1).

dat <- subset(cue.day,select=c("Date","T_treatment","Water_treatment","chamber","CUE","RtoA","GPP_la","Ra_la","PAR","leafArea"))
dat2 <- subset(dat,Water_treatment=="control")
dat2$T_treatment <- as.factor(dat2$T_treatment)
dat2$MonthFac <- as.factor(month(dat2$Date))
dat2$DateFac <- as.factor(dat2$Date)

#####
#- Ra

#- find the right transformation
sp.Ra.simple <- lm(Ra_la~T_treatment*DateFac,data=dat2)
value.r <- boxcox(object=sp.Ra.simple,lambda=seq(-1,0,1/20))
exponent.r <- value.r$x[which.max(value.r$y)]
dat2$Ra_latrans <- with(dat2,Ra_la^exponent.r)

sp.ra <- lme(Ra_latrans~T_treatment*DateFac,random=list(~1|chamber),
             #corr=corAR1(~1|chamber),
             weights=varFunc(~as.numeric(Date)),
             data=dat2)

#look at model diagnostics
plot(sp.ra,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.ra,Ra_latrans~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each species
plot(sp.ra,Ra_latrans~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.ra, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
anova(sp.ra)
lsmeans(sp.ra,"T_treatment")
####


####
#- GPP
sp.gpp <- lme(GPP_la~T_treatment*DateFac,random=list(~1|chamber),data=dat2,
              weights=varFunc(~1/leafArea))

#look at model diagnostics
plot(sp.gpp,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.gpp,GPP_la~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each species
plot(sp.gpp,GPP_la~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.gpp, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
hist(sp.gpp$residuals)
anova(sp.gpp)
lsmeans(sp.gpp,"T_treatment")
####


#####
#- RtoA, with box-cox transformation
sp.CUE.simple <- lm(RtoA~T_treatment*DateFac,data=dat2)
value <- boxcox(object=sp.CUE.simple,lambda=seq(0,1,1/20))
exponent <- value$x[which.max(value$y)]
  
dat2$RtoAtrans <- with(dat2,RtoA^exponent) # get the "best" transformation 
sp.CUE <- lme(RtoAtrans~T_treatment*DateFac,random=list(~1|chamber),data=dat2,
              #corr=corAR1(value=0.1,form=~1|chamber)
              )

#look at model diagnostics
plot(sp.CUE,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.CUE,RtoAtrans~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each species
plot(sp.CUE,RtoAtrans~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.CUE, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
anova(sp.CUE)
lsmeans(sp.CUE,"T_treatment") 
####


#- merge models together to make Table 1
table.r1 <- as.matrix(anova(sp.ra))
table.r2 <- cbind(table.r1,as.matrix(anova(sp.gpp))[1:4,3:4])
table1 <- cbind(table.r2,as.matrix(anova(sp.CUE))[1:4,3:4])[2:4,]

####
#-- analysis of just the exceedingly hot days
hotDates <- unique(dat.hr.p[which(dat.hr.p$Tair_al>40),"Date"]) # find dates with temperatures exceeding 40
dat.hr.p.hot <- subset(dat.hr.p,Date %in% hotDates)

#- daily climate metrics
hotDates_met <- summaryBy(Tair_al~T_treatment+Date,data=dat.hr.p.hot,FUN=c(mean,min,max),na.rm=T)
summaryBy(Tair_al.max~T_treatment,data=hotDates_met) # average maximum temperature on these hot dates

#- re-analyze on hot dates only
sp.CUE.hot <- lme(RtoAtrans~T_treatment*DateFac,random=list(~1|chamber),data=subset(dat2,Date %in% hotDates),
              #corr=corAR1(value=0.1,form=~1|chamber)
)

#look at model diagnostics
plot(sp.CUE.hot,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.CUE.hot,RtoAtrans~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each species
plot(sp.CUE.hot,RtoAtrans~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.CUE.hot, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
anova(sp.CUE.hot)
lsmeans(sp.CUE,"T_treatment") 
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------