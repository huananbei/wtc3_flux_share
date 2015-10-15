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
dat2$DateFac <- as.factor(dat2$Date)

#####
#- Ra
sp.ra <- lme(log(Ra_la)~T_treatment*DateFac,random=list(~1|chamber),
             #corr=corAR1(~1|chamber),
             weights=varFunc(~as.numeric(Date)),
             data=dat2)

#look at model diagnostics
plot(sp.ra,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.ra,log(Ra_la)~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each species
plot(sp.ra,log(Ra_la)~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.ra, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
anova(sp.ra)
lsmeans(sp.ra,"T_treatment")

#refit for summer months only
sp.ra.summer <- lme(Ra_la~T_treatment+DateFac,random=list(~1|chamber),
                    weights=varFunc(~as.numeric(Date)),
                    data=subset(dat2,month(Date)==1 | month(Date)==12 | month(Date)==2))
lsmeans(sp.ra.summer,"T_treatment")
#refit for periods outside of summer
sp.ra.notsummer <- lme(Ra_la~T_treatment*DateFac,random=list(~1|chamber),
                       weights=varFunc(~as.numeric(Date)),
                       data=subset(dat2,month(Date)!=1 & month(Date)!=12 & month(Date)!=2))
lsmeans(sp.ra.notsummer,"T_treatment")
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
#- RtoA
sp.CUE <- lme(log10(RtoA)~T_treatment*DateFac,random=list(~1|chamber),data=dat2,
              weights=varIdent(form=~1|T_treatment))


#look at model diagnostics
plot(sp.CUE,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.CUE,log10(RtoA)~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each species
plot(sp.CUE,log10(RtoA)~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.CUE, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
anova(sp.CUE)
lsmeans(sp.CUE,"T_treatment") 
####





#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------