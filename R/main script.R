
source("R/functions.R")



#-------------------------------------------------------------------------------------------------------------------
#- plot R/A with a leaf-scale model (Figure 1)
plotCUE_conceptual_fig(toexport=F,Tdew=10,Ca=400,Vcmax=100,Jmax=125,Tleaf=10:42,PPFD=1500)
#-------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------
#- Plot respiration vs. temperature curves for leaves and whole-tree chambers (Figure 2)

#- get the whole-canopy R vs. T datasets
fits.list <- return_Rcanopy_closed()
fits.mass <- fits.list[[1]]     #- tree-level data
fits.trt <- fits.list[[2]]      #- treatment averages

#- get the leaf-scale RvsT data
rvt.c <- getRTWTC3()

#- plot R vs. T (Figure 2)
plotRvsT_figure2(fits.mass=fits.mass,fits.trt=fits.trt,export=F)
#-------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------
#- plot Rbranch and Rleaf measured at 15 degrees C (Figure 3)
plotRleafRbranch(export=F)
#-------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------
#- Read and process the hourly flux dataset

#- read in the hourly flux data
dat.hr <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130910-20140530_L2_V1.csv")
dat.hr$DateTime <- as.POSIXct(dat.hr$DateTime,format="%Y-%m-%d %H:%M:%S",tz="GMT")
dat.hr$Date <- as.Date(dat.hr$DateTime)

#- partition the net fluxes into GPP and Ra components. Specify the actitation energy and the number
#-  of prior nights to include for the estimate of the basal respiration rate (lagdates)
dat.hr.p <- partitionHourlyFluxCUE_arr(dat.hr.gf=dat.hr,Ea=57.69,lagdates=3)

#- plot an example week of partitioned fluxes
plotPartitionedFluxes(dat.hr.gf3=dat.hr.p,ch_toplot="C07",startDate="2014-3-22",endDate="2014-3-27",write=F)

#- plot met, Ra, GPP, and Ra/GPP data over time (Figure 4)

#-------------------------------------------------------------------------------------------------------------------


#- get the whole-tree flux data
downloadHIEv(searchHIEv("WTC_TEMP_CM_WTCFLUX_20130910-20140530_L2_V1"),topath="data/")
