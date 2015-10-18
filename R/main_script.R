#- required libraries
library(mgcv)
library(scales)
library(gplots)
library(plotBy)
library(magicaxis)
library(lubridate)
library(doBy)
library(Hmisc)
library(zoo)
library(hexbin)
library(nlme)
library(lsmeans)
library(car)
library(data.table)

#- the following libraries aren't on CRAN, but can be installed from github or bitbucket with devtools
library(devtools)
#install_bitbucket("remkoduursma/plantecophys")
#install_bitbucket("remkoduursma/HIEv")
#install_github("jslefche/piecewiseSEM")
library(piecewiseSEM) # for estimating r2 value in mixed-effects models
library(plantecophys) # for modeling leaf-level gas exchange
library(HIEv)         # is this needed?

#- load the analysis and plotting functions that do all of the actual work
source("R/functions.R")

#- export flag. Set to "T" to create pdfs of figures in "output/", or "F" to suppress output.
#- This flag is passed to many of the plotting functions below.
export=F





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