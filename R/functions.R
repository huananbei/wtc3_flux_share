#----------------------------------------------------------------------------------------------------------------
#- This is a collection of many functions that do the actual work of data analysis and plotting. These functions
#    are called by just a few lines of code in "main script.R" to recreate the analyses and figures.
#----------------------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------------
# Demonstration of the short-term temperature effects on photosynthesis and respiration.

plotCUE_conceptual_fig <- function(toexport=T,Tdew=10,Ca=400,Vcmax=100,Jmax=125,Tleaf=10:45,PPFD=1500){
  # load required libraries
  #- set up parameters for model
  #VPD <- RHtoVPD(RH=RH,TdegC=Tleaf)
  VPD <- DewtoVPD(Tdew=Tdew,TdegC=Tleaf)
  
  #- predict photosynthesis and respiration with changing T (and VPD)
  output<- Photosyn(VPD=VPD,Ca=Ca,Vcmax=Vcmax,Jmax=Jmax,Tleaf=Tleaf,Tcorrect=T,Rd0=1,TrefR=20,PPFD=PPFD)
  output$AGROSS <- with(output,ALEAF+Rd)
  output.acclim <-Photosyn(VPD=VPD,Ca=Ca,Vcmax=Vcmax,Jmax=Jmax*0.95,Tleaf=Tleaf,Tcorrect=T,Rd0=0.5,TrefR=20,delsJ=600,PPFD=PPFD)
  output.acclim$AGROSS <- with(output.acclim,ALEAF+Rd)
  
  #- calculate CUE
  CUE <- 1-(output$Rd/output$ALEAF)
  CUEa <- 1-(output.acclim$Rd/output.acclim$ALEAF)
  RtoA <- (output$Rd/output$AGROSS)
  RtoA2 <- (output.acclim$Rd/output.acclim$AGROSS)
  
  #- plot
  windows(20,40);par(mfrow=c(2,1),mar=c(2,7,1,1),oma=c(3,0,2,0),xpd=F,las=1)
  plot(AGROSS~Tleaf,data=output,type="l",ylab="",cex.lab=1.4,ylim=c(0,23),
       col="forestgreen",lwd=2,axes=F)
  magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=3)
  title(ylab=expression(atop(Leaf~CO[2]~exchange,
                             (mu*mol~CO[2]~m^-2~s^-1))),cex.lab=1.3)
  lines(AGROSS~Tleaf,data=output.acclim,type="l",col="forestgreen",lwd=2,lty=2)
  lines(Rd~Tleaf,data=output,col="red",lwd=2)
  lines(Rd~Tleaf,data=output.acclim,col="red",lwd=2,lty=2)
  legend(x=12,y=27,xpd=NA,legend=c("A","R"),lwd=2,col=c("forestgreen","red"),ncol=2,bty="n")
  legend("topleft","a",cex=1.1,bty="n",inset=-0.05)
  
  plot(RtoA~output$Tleaf,type="l",ylab="",ylim=c(0,0.6),cex.lab=1.3,lwd=2,col="black",axes=F)
  magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=3)
  title(ylab=expression(R/A),cex.lab=1.3)
  
  lines(RtoA2~output.acclim$Tleaf,type="l",ylab="CUE (1-R/A)",ylim=c(0.2,1),cex.lab=1.3,lwd=2,lty=2,col="black")
  legend("topleft","b",cex=1.1,bty="n",inset=-0.05)
  
  par(xpd=NA)
  text(x=29,y=-0.15,expression(T[leaf]~(degree*C)),cex=1.4)
  
  if(toexport==T) dev.copy2pdf(file="./output/Figure1.pdf")
}
#----------------------------------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------------------
#- function to return leak parameters of WTC fluxes while in a closed state
return.leaks <- function(plotson=0){

  # I did leak tests for WTC chambers #1-4 on the night of 27 May 2014, and all 12 chambers on 28 May 2014.

  #read in all the minutely data
  dat1.1 <- read.csv(file="data/Min140527.csv")
  dat1.2 <- read.csv(file="data/Min140528.csv")
  dat1.3 <- read.csv(file="data/Min140529.csv")
  dat <- subset(rbind(dat1.1,dat1.2,dat1.3),chamber<13)
  dat$datetime <- as.POSIXct(paste(dat$date,dat$time,sep=" "),format="%Y-%m-%d %H:%M:%S",tz="GMT")
  
  #set up the times for start and end times for the two nights
  starttime1 <- as.POSIXct("2014-05-27 22:00:10",format="%Y-%m-%d %H:%M:%S",tz="GMT")
  endtime1 <- as.POSIXct("2014-05-28 07:30:00",format="%Y-%m-%d %H:%M:%S",tz="GMT")
  starttime2 <- as.POSIXct("2014-05-28 20:50:00",format="%Y-%m-%d %H:%M:%S",tz="GMT")
  endtime2 <- as.POSIXct("2014-05-29 07:30:10",format="%Y-%m-%d %H:%M:%S",tz="GMT")
  
  
  print("Processing data...")

    #------------------------------------------------------------------------------------------------------------------
  #work out the reference [CO2] data. Note that this is chamber "13", and the data has a different structure.
  ref <- subset(rbind(dat1.1,dat1.2,dat1.3),chamber==13)
  ref$datetime <- as.POSIXct(paste(ref$date,ref$time,sep=" "),format="%Y-%m-%d %H:%M:%S",tz="GMT")
  
  # just grab the [CO2] data, interpolate it
  ref2 <- ref[,c(13,19)]
  names(ref2)[1] <- "CO2ref"
  times <- data.frame(datetime =seq.POSIXt(from=starttime1,to=endtime2,by="sec"))
  
  ref.i <- merge(times,ref2,all.x=T)
  ref.i <- ref.i[with(ref.i,order(datetime)),]
  
  ref.zoo <- zoo(ref.i$CO2ref)
  ref.zoo <- na.approx(ref.zoo) #gapfill between observaitons
  
  ref.i$CO2ref <- ref.zoo #replace missing values with gapfilled numbers
  #------------------------------------------------------------------------------------------------------------------
  
  
  
  #------------------------------------------------------------------------------------------------------------------
  #subset the data for the first night
  dat1 <- subset(dat,datetime>starttime1&datetime<endtime1)
  dat1 <- merge(ref.i,dat1,by="datetime")
  dat1 <- subset(dat1,chamber<=4) #only the first four chambers were leak tested on the first day
  
  #subset the data for the second night
  dat2 <- subset(dat,datetime>starttime2&datetime<endtime2)
  dat2 <- merge(ref.i,dat2,by="datetime")
  
  
  
  #model the leak, return the squared residual sum for minimiation by optimise(), which is much faster than DEoptim().
  leak.mod <- function(theta, V, Ca, Ci,fit=1){
    
    resid <- rep(NA,length(Ca))
    resid[1] <- 0
    pred <- rep(NA,length(Ca))
    pred[1] <- Ci[1]
    
    for (i in 2:length(Ca)){
      dCO2 <- -theta*(Ci[i-1]-Ca[i-1])/V
      pred[i] <- pred[i-1] + dCO2
      
      resid[i] <- (Ci[i] - pred[i])^2
    }
    resid.sum <- sum(resid)
    if (fit==1) return(resid.sum)
    if (fit==0) return(data.frame(Ci=Ci,Ca=Ca,pred=pred))
  }
  
  
  print("Fitting the model to estimate leak parameters for the first day of measurements...")
  #fit the first day
  dat1.list <- split(dat1,dat1$chamber)
  out1 <- list()
  pred <- list()
  pred[[1]] <- data.frame()
  for (i in 1:length(dat1.list)){
    dat <- dat1.list[[i]]
    
    out1[i] <- optimise(f=leak.mod,interval=c(0,0.5),V=60,Ca=dat$CO2ref,Ci=dat$CO2L,fit=1) # volume was 30, changed to 60
    #pred[i] <- leak.mod(theta=out1[[i]][1],V=30,Ca=dat$CO2ref,Ci=dat$CO2L,fit=0)
  }
  
  print("Fitting the model to estimate leak parameters for the second day of measurements...")
  
  #fit the second day
  dat2.list <- split(dat2,dat2$chamber)
  out2 <- list()
  # pred <- list()
  # pred[[1]] <- data.frame()
  for (i in 1:length(dat2.list)){
    dat <- dat2.list[[i]]
    
    out2[i] <- optimise(f=leak.mod,interval=c(0,0.5),V=30,Ca=dat$CO2ref,Ci=dat$CO2L,fit=1)
    #pred[i] <- leak.mod(theta=out[[i]][1],V=30,Ca=dat$CO2ref,Ci=dat$CO2L,fit=0)
  }
  
  leaks_day1 <- data.frame(theta=do.call(rbind,out1))
  leaks_day1$chamber <- 1:4
  
  
  #compile leak estimates
  leaks <- data.frame(theta=do.call(rbind,out2))
  leaks$chamber <- 1:12
  #leaks$theta2 <- c(do.call(rbind,out1),rep(NA,8))
  leaks <- leaks[,c(2,1)]
  #leaks$theta <- rowMeans(leaks[,2:3],na.rm=T)
  
  
  #-- what was the variation of theta measured across the two nights?
  dtheta <- rbind(leaks_day1,leaks[1:4,])
  leaksd <- summaryBy(theta~chamber,data=dtheta,FUN=sd)
  leakmean <- summaryBy(theta~chamber,data=dtheta,FUN=mean)
  cv <- merge(leaksd,leakmean,by="chamber")
  cv$cv <- with(cv,theta.sd/theta.mean)
  #----------------------------------------------------------------------------------------------
  #plot fits, if plotson==1
  
  if(plotson==1){
    print("Plotting")
    windows(12,10)
    par(mfrow=c(4,3),mar=c(2,4,2,2))
    
    #fit the second day
    dat2.list <- split(dat2,dat2$chamber)
    out2 <- list()
    pred <- data.frame(Ci=NA,Ca=NA,pred=NA)
    # pred[[1]] <- data.frame()
    for (i in 1:length(dat2.list)){
      dat <- dat2.list[[i]]
      
      #out2[i] <- optimise(f=leak.mod,interval=c(0,0.5),V=30,Ca=dat$CO2ref,Ci=dat$CO2L,fit=1)
      pred <- leak.mod(theta=leaks$theta[i],V=30,Ca=dat$CO2ref,Ci=dat$CO2L,fit=0)
      
      plot(pred$Ci,col="black",lty=1,lwd=2,type="l",ylab="[CO2]",ylim=c(350,1100))
      lines(pred$Ca,col="grey",lty=1,lwd=2)
      lines(pred$pred,col="red",lty=1,lwd=2)
      title(main=paste("Chamber",i,sep=" "))
      if (i==3){
        legend("topright",c("data","pred","Ca"),lty=1,lwd=2,col=c("black","red","grey"),cex=1.5,bty="n")
      }
    }
  }
  
  
  print("Done with leak calculations.")
  return(leaks)
}
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#- code to return WTC flux estimates from closed chambers (whole canopy T-response curves)

return_Rcanopy_closed <- function(){

  #source("//ad.uws.edu.au/dfshare/HomesHWK$/30035219/My Documents/Work/HFE/WTC3/gx_wtc3/R/loadLibrariesWTC3.R")
  
  
  # get the leak parameter for each chamber. This will be used later to estimate the rate of respiration
  leaks <- return.leaks(plotson=0)
  
  #------------------------------------------------------
  #------------------------------------------------------
  #------------------------------------------------------
  #read in the minutely IRGA data, including a local IRGA (CO2L) and a central LI7000 (CO2CConc)
  # I manually removed the reference data from these files, sot that they could be read into R
  # These data were measured from 7pm (CPU time) on 11 Feb 14 until ~4:30am on 12 Feb 14
  #Min1 <- read.csv("./data/Rtree/Min140212_noref_day2.csv")
  #Min2 <- read.csv("./data/Rtree/Min140213_noref.csv")
  
  Min1 <- read.csv("data/Min140212.csv")
  Min2 <- read.csv("data/Min140213.csv")
  
  
  Mindat1 <- rbind(Min1,Min2)
  Mindat1$DateTime <- paste(Mindat1$date,Mindat1$time,sep=" ")
  Mindat1$DateTime <- as.POSIXct(Mindat1$DateTime,tz="GMT")
  
  
  
  #5 sets of measurements
  #start and stop times are different for this computer relative to the tdl computer
  starttime1 <- as.POSIXct("2014-02-12 19:16:09 GMT",format="%Y-%m-%d %T",tz="GMT")
  stoptime1 <- as.POSIXct("2014-02-12 20:45:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  starttime2 <- as.POSIXct("2014-02-12 21:05:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  stoptime2 <- as.POSIXct("2014-02-12 22:40:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  starttime3 <- as.POSIXct("2014-02-12 23:00:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  stoptime3 <- as.POSIXct("2014-02-13 00:30:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  starttime4 <- as.POSIXct("2014-02-13 01:00:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  stoptime4 <- as.POSIXct("2014-02-13 02:30:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  starttime5 <- as.POSIXct("2014-02-13 03:40:00 GMT",format="%Y-%m-%d %T",tz="GMT") 
  stoptime5 <- as.POSIXct("2014-02-13 04:29:09 GMT",format="%Y-%m-%d %T",tz="GMT")
  
  #Mindat2 <- subset(Mindat,DateTime>starttime1 & DateTime<= stoptime5)
  
  
  #calcualte time elapsed
  dat1.min <- subset(Mindat1,DateTime >= starttime1 & DateTime <= stoptime1);dat1.min$Measurement <- 1
  dat2.min <- subset(Mindat1,DateTime >= starttime2 & DateTime <= stoptime2);dat2.min$Measurement <- 2
  dat3.min <- subset(Mindat1,DateTime >= starttime3 & DateTime <= stoptime3);dat3.min$Measurement <- 3
  dat4.min <- subset(Mindat1,DateTime >= starttime4 & DateTime <= stoptime4);dat4.min$Measurement <- 4
  dat5.min <- subset(Mindat1,DateTime >= starttime5 & DateTime <= stoptime5);dat5.min$Measurement <- 5
  
  Mindat <- rbind(dat1.min,dat2.min,dat3.min,dat4.min,dat5.min)
  Mindat$T_treatment <- ifelse(Mindat$chamber %% 2 ==1,"ambient","elevated")
  
  #------------------------------------------------------------------------------------------------------------------
  #work out the reference [CO2] data. Note that this is chamber "13", and the data has a different structure.
  #note that the actual measured [CO2] by the central irga of the 
  ref <- subset(Mindat,chamber==13)
  
  # just grab the [CO2] data, interpolate it
  ref2 <- ref[,c(13,19)]
  names(ref2)[1] <- "CO2ref"
  times <- data.frame(DateTime =seq.POSIXt(from=starttime1,to=stoptime5,by="sec",tz="GMT"))
  
  ref.i <- merge(times,ref2,all=T)
  ref.i <- ref.i[with(ref.i,order(DateTime)),]
  
  ref.zoo <- zoo(ref.i)
  ref.zoo$CO2ref <- na.approx(ref.zoo$CO2ref) #gapfill between observaitons
  ref.df <- fortify.zoo(ref.zoo)
  ref.df$CO2ref <- as.numeric(as.character(ref.df$CO2ref))
  ref.df$DateTime <- as.POSIXct(ref.df$DateTime,tz="GMT")
  

  #------------------------------------------------------------------------------------------------------------------
  
  
  
  #------------------------------------------------------------------------------------------------------------------
  #merge the chamber data with the interpolated reference data
  dat1 <- merge(ref.df,subset(Mindat,chamber<13),by="DateTime")
  
  #-- so that worked great, except for chamber 9, where the local irga died. Ugh.
  # I"ll just fit the first three temperatures for chamber 9, which have good data
  dat9 <- subset(dat1,chamber==9)
  dat92 <- subset(dat9,DateTime <as.POSIXct("2014-02-13 00:35:00",tz="GMT"))
  
  dat1.1 <- subset(dat1,chamber !=9)
  
  dat2 <- rbind(dat1.1,dat92)
  #------------------------------------------------------------------------------------------------------------------
  
  
  
  #----------------------------------------------------------------------------------------------
  # get the chamber temperatures by downloading the inside met data from HIEv
  #downloadHIEv(searchHIEv("WTC_TEMP_CM_WTCMET_20140201-20140228_L1_v1.csv"),topath="data/from HIEv")
  met <- read.csv("C:/Repos/wtc3_flux/data/from HIEv/WTC_TEMP_CM_WTCMET_20140201-20140228_L1_v1.csv")
  met$DateTime <- as.POSIXct(met$DateTime,tz="GMT")
  #met2 <- subset(met,DateTime>starttime1 & DateTime < stoptime5)
  
  #set up the data for merging with the met
  dat2$chamber <- as.factor(paste0("C",sprintf("%02.0f",dat2$chamber)))
  dat2$DateTime1 <- nearestTimeStep(dat2$DateTime,nminutes=15,align="ceiling")
  
  dat3 <- merge(dat2,met,by.x=c("chamber","DateTime1"),by.y=c("chamber","DateTime"))
  #----------------------------------------------------------------------------------------------
  
  
  #------------------------------------------------------------------------------------------------------------------
  # create a factor for measurement number
  dat3$id <- as.factor(paste(dat3$chamber,dat3$Measurement,sep="-"))
  dat3.l <- split(dat3,dat3$id)
  #------------------------------------------------------------------------------------------------------------------
  
  
  
  
  
  #----------------------------------------------------------------------------------------------
  #model the flux, given the leak rate, return the squared residual sum for minimiation by optimise(), which is much faster than DEoptim().
  closedWTC.mod <- function(R,theta, V, Ca, Ci,fit=1){
    
    resid <- rep(NA,length(Ca))
    resid[1] <- 0
    pred <- rep(NA,length(Ca))
    pred[1] <- Ci[1]
    
    for (i in 2:length(Ca)){
      dCO2 <- (-theta*(Ci[i-1]-Ca[i-1])+R)/V
      pred[i] <- pred[i-1] + dCO2
      
      resid[i] <- (Ci[i] - pred[i])^2
    }
    resid.sum <- sum(resid)
    if (fit==1) return(resid.sum)
    if (fit==0) return(data.frame(Ci=Ci,Ca=Ca,pred=pred))
  }
  #----------------------------------------------------------------------------------------------
  
  
  
  #----------------------------------------------------------------------------------------------
  out1 <- c()
  fits1 <- c()
  #fit the respiration model using the leak estimate measured with an empty chamber
  for (i in 1:length(dat3.l)){
    dat <- dat3.l[[i]]
    theta <- leaks[dat$chamber[1],2]
    out1[[i]] <- optimise(f=closedWTC.mod,interval=c(0,200),theta=theta,V=60,Ca=dat$CO2ref,Ci=dat$CO2L,fit=1)[[1]] # V was 30
    chamber <- dat$chamber[1]
    id <- dat$id[1]
    T_treatment <- dat$T_treatment[1]
    Tair <- median(dat$Tair_al)
    RH <- median(dat$RH_al)
    Measurement <- dat$Measurement[1]
    
    fits1[[i]] <- data.frame(chamber,id,T_treatment,Measurement,theta,Tair,RH,Rcanopy=out1[[i]])
  }
  fits <- do.call(rbind,fits1)
  fits$mol_m3 <- 44.6 #(1*1000)/(0.0821*(273.15+fits$Tair))# calculate the mols of gas per m3 via pv=nrt. Remember 1m3 = 1000L. Assuming P = 1 atm
  # chamber 7 at 22 is a problem, as is ch1 at 25
  
  
  
  
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  #- Okay, now we have raw respiration data. Let's merge it with tree size to normalize the flux by mass (or leaf area)
  
  #- get an estimate of leaf area for each day of the experiment
  #treeMass <- returnTreeMass()
  treeMass <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130910-20140530_L2_V1.csv")
  treeMass$DateTime <- as.POSIXct(treeMass$DateTime,format="%Y-%m-%d %T",tz="GMT")
  treeMass2 <- subset(treeMass,DateTime==as.POSIXct("2014-02-12 01:00:00",format="%Y-%m-%d %T",tz="GMT"),
                      select=c("chamber","leafArea"))
  
  
  #- merge the respiration data with the tree mass data
  fits.mass <- merge(fits,treeMass2,by=c("chamber"))
  
  #- leaf area of ambient and warmed treatments
  summaryBy(leafArea~T_treatment,data=fits.mass,FUN=mean)
  
  #- unit conversions. Recall that Rcanopy has units of umol CO2 m^3 mol-1 s-1
  fits.mass$Rcanopy_umol <- fits.mass$Rcanopy*fits.mass$mol_m3/60
  
  fits.mass$R_la <- with(fits.mass,Rcanopy_umol/leafArea) # convert from umol CO2 s-1 to umol CO2 m-2 s-1
  fits.mass$mol_m3 <- NULL
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  
  
  #average by treatment 
  fits.mass$Tbin <- ifelse(fits.mass$Tair<20,18,
                           ifelse(fits.mass$Tair>20 & fits.mass$Tair < 24,22.5,
                                  ifelse(fits.mass$Tair>24 & fits.mass$Tair < 26, 25,
                                         ifelse(fits.mass$Tair>26 & fits.mass$Tair<30,28,32))))
  fits.trt <- summaryBy(Tair+Rcanopy_umol+R_la~T_treatment+Tbin,data=fits.mass,FUN=c(mean,standard.error),keep.names=F)
  
  #----------------------------------------------------------------------------------------------
  
  return(list(fits.mass,fits.trt))
  
}
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#- plot figure 2 (R vs. T curves at leaf and whole canopy scales)
plotRvsT_figure2 <- function(rvt.c=rvt.c,fits.mass=fits.mass,fits.trt=fits.trt,export=T){
  palette(c("black","red"))
  
  #--------------------------------------------------------------------------------------------------------
  #- fit arrheneous models for all three types of respiraiton
  # ANCOVA. No difference in slopes, so interaction removed from model
  fits.mass$logRcanopy <- log(fits.mass$Rcanopy_umol)
  fits.mass$logRarea <- log(fits.mass$R_la)
  
  
  fits.mass$invT <- 1/(fits.mass$Tair+273.15)
  
  lmRcanopy <- lm(logRcanopy~invT+T_treatment,data=fits.mass)
  lmRarea <- lm(logRarea~invT+T_treatment,data=fits.mass)
  #- extract fits for each type of R
  Rgas <- 8.314 #J per mol per K
  
  A_amb1 <- unname(exp(coef(lmRcanopy)[1]))
  A_ele1 <- unname(exp(coef(lmRcanopy)[1]+coef(lmRcanopy)[3]))
  Ea1 <- unname(-1*coef(lmRcanopy)[2]*Rgas)
  #(A_ele1-A_amb1)/A_amb1*100
  
  A_amb2 <- unname(exp(coef(lmRarea)[1]))
  A_ele2 <- unname(exp(coef(lmRarea)[1]+coef(lmRarea)[3]))
  Ea2 <- unname(-1*coef(lmRarea)[2]*Rgas)
  #(A_ele2-A_amb2)/A_amb2*100
  
  
  Q10 <- exp(10*Ea1/(Rgas*(25+273.15)^2))
  # get model predictions of R across all temperatures
  xvals <- seq(18,32.1, length=101)
  predA1 <- A_amb1*exp(-1*Ea1/(Rgas*(xvals+273.15)))
  predE1 <- A_ele1*exp(-1*Ea1/(Rgas*(xvals+273.15)))
  predA2 <- A_amb2*exp(-1*Ea2/(Rgas*(xvals+273.15)))
  predE2 <- A_ele2*exp(-1*Ea2/(Rgas*(xvals+273.15)))
  #--------------------------------------------------------------------------------------------------------
  
  
  
  
  #--------------------------------------------------------------------------------------------------------
  #- process leaf-scale R vs. T data
  rvt <- read.csv("data/WTC_TEMP_CM_GX-RdarkVsT_20140207-20140423_L1.csv")
  rvt$Date <- as.Date(rvt$Date)
  rvt.c <- subset(rvt,Date==as.Date("2014-02-07")) #- just pull out the data measured prior to the drought
  
  rvt.45 <- subset(rvt.c, Tleaf > 18 & Tleaf<=40)
  rvt.45$lnRmass <- log(rvt.45$Rmass)
  rvt.45$lnRarea <- log(rvt.45$Rarea)
  rvt.45$invT <- 1/(rvt.45$Tleaf+273.15)
  rvt.45$Treat <- as.factor(rvt.45$T_treatment)
  rvt.45$datefac <- as.factor(rvt.45$Date)
  
  rvt.45$Tleaf_bin <- cut(rvt.45$Tleaf,breaks=seq(from=18,to=40,length=25))
  rvt.45$Tleaf_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", rvt.45$Tleaf_bin), ","), function(x)sum(as.numeric(x))/2) 
  rvt.45.ch <- summaryBy(Rmass+invT+lnRmass+Rarea+lnRarea~date+chamber+Treat+Tleaf_bin_mid,data=subset(rvt.45, Tleaf_bin_mid<45),
                         keep.names=T,FUN=c(mean))
  
  rvt.treat.bin <- summaryBy(Rmass+Rarea~date+Treat+Tleaf_bin_mid,data=rvt.45.ch,keep.names=F,FUN=c(mean,standard.error))
  rvt.treat.bin$Rmass.high <- with(rvt.treat.bin,Rmass.mean+Rmass.standard.error)
  rvt.treat.bin$Rmass.low <- with(rvt.treat.bin,Rmass.mean-Rmass.standard.error)
  
  #- fit arrhenious. 
  lmRleaf <- lm(lnRarea~invT+Treat,data=rvt.45.ch)
  Rleaf_amb1 <- unname(exp(coef(lmRleaf)[1]))
  Rleaf_ele1 <- unname(exp(coef(lmRleaf)[1]+coef(lmRleaf)[3]))
  Ea_Rleaf <- unname(-1*coef(lmRleaf)[2]*Rgas)

  
  Q10_Rleaf <- exp(10*Ea_Rleaf/(Rgas*(25+273.15)^2))
  # get model predictions of R across all temperatures
  xvals_Rleaf <- seq(18,40, length=101)
  predRleafA <- Rleaf_amb1*exp(-1*Ea_Rleaf/(Rgas*(xvals_Rleaf+273.15)))
  predRleafE <- Rleaf_ele1*exp(-1*Ea_Rleaf/(Rgas*(xvals_Rleaf+273.15)))
  
  #--------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------
  
  
  
  
  
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  #- plot R vs. T with (a) raw fluxes, (b) fluxes per unit leaf area, (c) fluxes per unit total tree mass
  windows(22,32);par(mfrow=c(3,1),cex.lab=1.5,mar=c(2,7,1,2),oma=c(4,1,0,0),las=1)
  xlims=c(15,42)
  
  # #- plot AREA BASED leaf R over T for the first date of high resolution T-response curves
  plotBy(Rarea.mean~Tleaf_bin_mid|Treat,data=rvt.treat.bin,axes=F,
         ylab=expression(atop(R[leaf],
                              (mu*mol~CO[2]~m^-2~s^-1))),col=c("black","red"),pch=15,
         xlim=xlims,ylim=c(0,4),type="p",lwd=3,cex=1.6,xlab="",cex.lab=1.6,legend=F,
         panel.first=adderrorbars(x=rvt.treat.bin$Tleaf_bin_mid,y=rvt.treat.bin$Rarea.mean,SE=rvt.treat.bin$Rarea.standard.error,direction="updown"))
  magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))
  legend("bottomright","a",bty="n",inset=0.02,cex=1.5)
  lines(x=xvals_Rleaf,y=predRleafA,col="black",lwd=2)
  lines(x=xvals_Rleaf,y=predRleafE,col="red",lwd=2)
  
    #- plot raw T-response curves
  plotBy(Rcanopy_umol.mean~Tair.mean|T_treatment,data=fits.trt,type="p",xlim=xlims,ylim=c(0,30),pch=15,cex=1.6,cex.lab=1.6,axes=F,
         ylab=expression(atop(R[canopy],
                              (mu*mol~CO[2]~s^-1))),xlab="",legend=F,
         panel.first=adderrorbars(x=fits.trt$Tair.mean,y=fits.trt$Rcanopy_umol.mean,SE=fits.trt$Rcanopy_umol.standard.error,direction="updown"))
  lines(x=xvals,y=predA1,col="black",lwd=2)
  lines(x=xvals,y=predE1,col="red",lwd=2)
  magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))
  legend("bottomright","b",bty="n",inset=0.02,cex=1.5)
  
  #- plot respiration per unit leaf area
  plotBy(R_la.mean~Tair.mean|T_treatment,data=fits.trt,type="p",pch=15,ylim=c(0,2.5),xlim=xlims,cex=1.6,cex.lab=1.6,legend=F,axes=F,
         ylab=expression(atop(R["canopy, area"],
                              (mu*mol~CO[2]~m^-2~s^-1))),xlab=expression(T[air]~(degree*C)),
         panel.first=adderrorbars(x=fits.trt$Tair.mean,y=fits.trt$R_la.mean,SE=fits.trt$R_la.standard.error,direction="updown"))
  lines(x=xvals,y=predA2,col="black",lwd=2)
  lines(x=xvals,y=predE2,col="red",lwd=2)
  legend("bottomright","c",bty="n",inset=0.02,cex=1.5)
  magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))
  title(xlab=expression(Temperature~(degree*C)),xpd=NA)
  
  if(export==T) dev.copy2pdf(file="output/Figure2.pdf")
}
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#- get and plot the Rleaf and Branch data near the end of the experiment
plotRleafRbranch <- function(export=T){

  Rbranch <- read.csv("data/WTC_TEMP_CM_GX-RBRANCH_20140513-20140522_L1_v1.csv")
  Rbranch$date <- as.Date(Rbranch$date)
  Rbranch$campaign <- ifelse(Rbranch$date == as.Date("2014-05-13"),1,2)
  R1 <- subset(Rbranch,campaign==1)
  
  #- average across treatments
  R1.ch <- summaryBy(Rleaf+Rbranch~T_treatment+chamber,data=R1,FUN=c(mean),keep.names=T)
  R1.m <- summaryBy(Rleaf+Rbranch~T_treatment,data=R1.ch,FUN=c(mean,standard.error))
  
  R1.m$Rleaf.cil <- R1.m$Rleaf.mean-R1.m$Rleaf.standard.error
  R1.m$Rleaf.ciu <- R1.m$Rleaf.mean+R1.m$Rleaf.standard.error
  R1.m$Rbranch.cil <- R1.m$Rbranch.mean-R1.m$Rbranch.standard.error
  R1.m$Rbranch.ciu <- R1.m$Rbranch.mean+R1.m$Rbranch.standard.error
  
  #-----------------------------------------------------------------------------------------------------------
  # plot acclimation of tissue-specific R
  windows(8,8)
  par(mfrow=c(1,2),mar=c(5,3,3,1),oma=c(1,4,0,0),cex.lab=2.5,las=1,cex.axis=1.8,cex.lab=2)
  
  #- plot leaves
  barplot2(height=R1.m$Rleaf.mean,plot.ci=T,ci.l=R1.m$Rleaf.cil,ci.u=R1.m$Rleaf.ciu,axes=F,ci.width=0.2,col=c("darkgrey","red"),
           names.arg=c("A","W"),ylim=c(0,4))
  legend("topright","a",bty="n",inset=-0.002,cex=1.5)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  title(main="Leaves",cex.main=1.5)
  title(xlab="Treatment",outer=F,line=3)
  
  #- plot branches
  barplot2(height=R1.m$Rbranch.mean,plot.ci=T,ci.l=R1.m$Rbranch.cil,ci.u=R1.m$Rbranch.ciu,axes=F,ci.width=0.2,col=c("darkgrey","red"),
           names.arg=c("A","W"),ylim=c(0,2))
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  legend("topright","b",bty="n",inset=-0.002,cex=1.5)
  title(main="Branches",cex.main=1.5)
  title(xlab="Treatment",outer=F,line=3)
  
  
  title(ylab=expression(Respiration~(nmol~g^-1~s^-1)),outer=T,line=0)
  
  #   
  #   #- statistical analyses
  #   lm.leaf <- lme(Rleaf~T_treatment,random=~1|chamber,data=R1)
  #   
  #   
  #   #look at model diagnostics
  #   plot(lm.leaf,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
  #   plot(lm.leaf,Rleaf~fitted(.),abline=c(0,1))         #predicted vs. fitted for each species
  #   qqnorm(lm.leaf, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
  #   hist(lm.leaf$residuals)
  #   anova(lm.leaf)
  #   lsmeans(lm.leaf,"T_treatment")
  #   
  #   #- statistical analyses
  #   lm.br <- lme(Rbranch~T_treatment,random=~1|chamber,data=R1)
  #   
  #   
  #   #look at model diagnostics
  #   plot(lm.br,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
  #   plot(lm.br,Rbranch~fitted(.),abline=c(0,1))         #predicted vs. fitted for each species
  #   qqnorm(lm.br, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
  #   hist(lm.br$residuals)
  #   anova(lm.br)
  #   lsmeans(lm.br,"T_treatment")
  #   
  
  if (export==T) dev.copy2pdf(file="output/Rtree_tissue_specific_Rates.pdf")
}
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#' Adds error bars to a plot
#' 
#' @description Yet another function that adds error bars. The user must specify the length of the error bars.
#' @param x The x coordinates of start of the error bar
#' @param y The y coordinates of start of the error bar
#' @param SE The length of the error bar
#' @param direction One of 'up', 'down', 'right', 'left', 'updown' or 'rightleft'.
#' @param barlen The length of the cross-bar at the end of the error bar.
#' @param \ldots Additional parameters passed to \code{\link{arrows}}, such as the colour (\code{col}).
#' #' @details Simple wrapper for \code{\link{arrows}}, where \code{angle=90} and \code{code=3}. The \code{barlen} argument corresponds to \code{length} in \code{arrows}.
#' @examples
#' # A simple example. Also note that we can specify the colour of the error bars, or other parameters
#' # that arrows() recognizes.
#' x <- rnorm(20)
#' y <- x + rnorm(20)
#' se <- runif(20, 0.2,0.4)
#' plot(x,y,pch=21,bg="white",panel.first=adderrorbars(x,y,se,direction="updown", col="darkgrey"))
#' @export
adderrorbars <- function(x,y,SE,direction,barlen=0.04,...){
  
  if(length(direction)>1)stop("direction must be of length one.")
  if(direction == "updown")
    direction <- c("up","down")
  else if(direction == "rightleft" | direction == "leftright")direction <- c("left","right")
  
  if("up" %in% direction)
    arrows(x0=x, x1=x, y0=y, y1=y+SE, code=3, angle=90, length=barlen,...)
  if("down" %in% direction) 
    arrows(x0=x, x1=x, y0=y, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("left" %in% direction) 
    arrows(x0=x, x1=x-SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)
  if("right" %in% direction)
    arrows(x0=x, x1=x+SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)  
  
}
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
standard.error <- function(dat,na.rm=F,...){
  if(na.rm==T){
    dat <- subset(dat,is.na(dat)==F)
  }
  std <- sd(dat)
  n <- length(dat)
  se <- std/sqrt(n)
  return(se)
}
#----------------------------------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------------------------------
partitionHourlyFluxCUE_arr <- function(dat.hr.gf=dat.hr.gf,Ea=57.69,lagdates){
  rvalue = 8.134
  require(data.table)
  
  #-- convert mmolCO2 s-1 to gC hr-1
  dat.hr.gf$FluxCO2_g <- with(dat.hr.gf,FluxCO2*60*60/1000*12.0107)
  dat.hr.gf$period <- ifelse(dat.hr.gf$PAR>2,"Day","Night")
  
  #-- partition day-time net C exchange into GPP and Ra, similar to how it is done in eddy-covariance.
  #-- create a series of dates
  date.vec <- seq.Date(from=min(dat.hr.gf$Date),to=max(dat.hr.gf$Date),by="day")
  
  #-- estimate R-Tref and Tref for each date for each chamber
  #lagDates <- 3 # establish how many prior days to include
  RTdat <- expand.grid(Date=date.vec,chamber=levels(dat.hr.gf$chamber))
  RTdat$Tref <- RTdat$R_Tref <- NA
  
  
  print("Partitioning Net CO2 fluxes into GPP and Ra")
  #- set up progress bar to track that this is working
  pb <- txtProgressBar(min = 0, max = nrow(RTdat), style = 3)
  
  for (i in 1:nrow(RTdat)){
    #- trial a data.table alternative to speed this up. The filter thing was actually slower.
    
    #dat <- dplyr::filter(dat.hr.gf,chamber==RTdat$chamber[i],Date <= RTdat$Date[i], Date >= (RTdat$Date[i]-lagDates),period =="Night")
    #RTdat$Tref[i] <- mean(dat$Tair_al,na.rm=T)
    #RTdat$R_Tref[i] <- mean(dat$FluxCO2_g,na.rm=T)
    
    inds <- which(dat.hr.gf$chamber==RTdat$chamber[i] & dat.hr.gf$Date <= RTdat$Date[i] & dat.hr.gf$Date >= (RTdat$Date[i]-lagdates) & dat.hr.gf$period =="Night" )
    RTdat$Tref[i] <- mean(dat.hr.gf$Tair_al[inds],na.rm=T)
    RTdat$R_Tref[i] <- mean(dat.hr.gf$FluxCO2_g[inds],na.rm=T)
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  
  RTdat$Tref_K <- with(RTdat,Tref+273.15)
  
  #-- merge these reference data into the gap-filled flux dataframe to estimate Ra during the daytime, and hence GPP
  dat.hr.gf3 <- merge(dat.hr.gf,RTdat,by=c("Date","chamber"))
  #dat.hr.gf3$Ra_est <- with(dat.hr.gf3,R_Tref*Q10^((Tair_al-Tref)/10)) # estimate respiration rate. This is a negative number.
  dat.hr.gf3$Ra_est <- with(dat.hr.gf3,R_Tref*exp((Ea*1000/(rvalue*Tref_K))*(1-Tref_K/(Tair_al+273.15)))) # estimate respiration rate. This is a negative number.
  

  dat.hr.gf3$GPP <- ifelse(dat.hr.gf3$period=="Night",0,dat.hr.gf3$FluxCO2_g-dat.hr.gf3$Ra_est)
  dat.hr.gf3$Ra <- ifelse(dat.hr.gf3$period=="Night",dat.hr.gf3$FluxCO2_g,dat.hr.gf3$Ra_est)
  
  return(dat.hr.gf3)
}
#----------------------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------------
#-- Plots the net CO2 flux and its partitioning into GPP and Ra. Defaults to an example week for chamber 7
plotPartitionedFluxes <- function(dat.hr.gf3=dat.hr.gf3,ch_toplot="C07",startDate="2014-3-22",endDate="2014-3-27",write=F){
  
  startDate <- as.Date(startDate)
  endDate <- as.Date(endDate)
  
  #-- plot an example week, with PAR, Tair, and Cfluxes
  windows(30,15);par(cex.lab=1.5,mar=c(0,0,0,0),oma=c(7,7,2,2),cex.axis=1.5,las=1)
  layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), 
         widths=c(1,1,1), heights=c(1,1,3))
  plotBy(PAR~DateTime,axes=F,data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=startDate & 
                                              as.Date(DateTime)<=endDate),ylab="PAR",lty=1,ylim=c(0,2000),type="l",legend=F,col="black")
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  mtext(text="PAR",side=2,outer=F,line=4,cex=1.5,las=0)
  plotBy(Tair_al~DateTime,axes=F,data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=as.Date("2014-3-22") & 
                                                  as.Date(DateTime)<=as.Date("2014-3-27")),ylab="Tair",lty=1,ylim=c(11,34),type="l",legend=F,col="black")
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  mtext(text="Air T",side=2,outer=F,line=4,cex=1.5,las=0)
  
  
  plotBy(GPP~DateTime,axes=F,data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=startDate & 
                                              as.Date(DateTime)<=endDate),ylab="CO2 flux (g hr-1)",pch=16,cex=2,ylim=c(-2,10),type="b",legend=F,col="green")
  plotBy(Ra_est~DateTime,axes=F,data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=startDate & 
                                                 as.Date(DateTime)<=endDate),pch=16,col="red",type="b",cex=2,add=T,legend=F)
  plotBy(FluxCO2_g~DateTime,axes=F,data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=startDate & 
                                                    as.Date(DateTime)<=endDate),type="b",pch=16,col="black",cex=2,add=T,legend=F)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  
  mtext(text=expression(CO[2]~flux~(gC~hr^-1)),side=2,outer=F,line=4,cex=1.5,las=0)
  legend("topright",c("GPP","Measured Net Flux","Ra"),pch=16,col=c("green","black","red"),cex=2)
  abline(0,0)
  axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct(startDate),to=as.POSIXct(endDate+1),by="day"),format="%F",cex.axis=2,line=1,tick=F)
  axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct(startDate),to=as.POSIXct(endDate+1),by="day"),format="%F",cex.axis=2,labels=F,tick=T)
  
  mtext(text="Date",side=1,outer=F,line=4,cex=1.5)
  
  if(write==T){dev.copy2pdf(file="output/FigureS1_Flux_partioning_example.pdf")}
}
#----------------------------------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------------------------------
#- function to return daily sums of GPP and Ra from hourly data
returnCUE.day <- function(dat=dat.hr.p){
  
  #- calculate leaf-area specific rates of GPP, convert to umol CO2 m-2 s-1
  dat$GPP_la <- with(dat,GPP/leafArea)
  dat$GPP_la_umol <- with(dat,GPP_la/12*1*10^6/60/60)
  dat$PAR_mol <- dat$PAR*60*60*1*10^-6

  
  
  #- create a date variable that moves the window for "night" observations, to sum observations for the night period following a day period
  dat$Date2 <- as.Date(dat$DateTime)
  dat$hour <- hour(dat$DateTime)
  earlynights <- which(dat$hour < 12 & dat$period =="Night")
  dat$Date2[earlynights] <- dat$Date2[earlynights]-1
  
  #- create daily sums 
  dat.day <- dplyr::summarize(group_by(dat,Date2,chamber,T_treatment,period),
                              GPP = sum(GPP,na.rm=T),
                              Ra = sum(Ra,na.rm=T),
                              FluxCO2_g=sum(FluxCO2_g,na.rm=T),
                              Tair_al=mean(Tair_al,na.rm=T),
                              VPDair=max(VPDair,na.rm=T),
                              PAR=sum(PAR_mol,na.rm=T),
                              leafArea=mean(leafArea,na.rm=T))
  dat.day <- as.data.frame(dat.day)
  dat.day <- dat.day[with(dat.day,order(Date2,chamber)),]
  names(dat.day)[1] <- "Date"
  
  Tair_day <- dplyr::summarize(group_by(dat,Date2,chamber,T_treatment),
                               Tair_al=mean(Tair_al,na.rm=T))
  names(Tair_day)[4] <- "Tair_24hrs"                   
  
  dat.day_tair <- as.data.frame(Tair_day)
  dat.day <- as.data.frame(dat.day)
  
  #- merge night and day data
  dat.day2d <- subset(dat.day,period=="Day")[,c("Date","chamber","T_treatment","GPP","Ra","FluxCO2_g","PAR","Tair_al","VPDair","leafArea")]
  names(dat.day2d)[4:8] <- c("GPP","Raday","Cgain","PAR","T_day")
  dat.day2n <- subset(dat.day,period=='Night')[,c("Date","chamber","T_treatment","Ra","FluxCO2_g","Tair_al","VPDair")]
  names(dat.day2n)[4:7] <- c("Ranight","Closs","T_night","VPDair_night")
  names(dat.day2n)[1] <- c("Date")
  
  #- merge data such that the nightly data are one day ahead of the daily data (i.e., does tonight's respiration depend on today's photosynthesis?)
  cue.day1 <- subset(merge(dat.day2d,dat.day2n,by=c("Date","chamber","T_treatment")))
  cue.day <- merge(cue.day1,dat.day_tair,by.x=c("Date","chamber","T_treatment"),by.y=c("Date2","chamber","T_treatment"))
  cue.day$Ra <- with(cue.day,-1*(Raday+Ranight))
  cue.day$RtoA <- with(cue.day,Ra/GPP)
  cue.day$CUE <- with(cue.day,(1-(Ra/GPP)))
  cue.day$GPP_la <- with(cue.day,GPP/leafArea)
  cue.day$Ra_la <- with(cue.day,Ra/(leafArea)) #g C m-2 day-1
  
  cue.day <- cue.day[-which(cue.day$Ra<0),] # remove two days from the beginning of the dataset in C09. Problem with flux data.
  
  
  
  
  #- merge in the key
  key <- read.csv("data/WTC_TEMP_CM_TREATKEY_20121212-20140528_L1_v1.csv")
  cue.day$T_treatment <- ifelse(cue.day$T_treatment=="ambient","ambient","elevated")
  cue.day2 <- merge(cue.day,key,by=c("chamber","Date","T_treatment"))
  
  cue.day.trt <- summaryBy(GPP+GPP_la+Ra+Ra_la+RtoA+PAR+T_day+T_night+Tair_24hrs+VPDair+VPDair_night~Date+T_treatment+Water_treatment,data=cue.day2,FUN=c(mean,standard.error))
  
  return(list(cue.day2,cue.day.trt))
}
#----------------------------------------------------------------------------------------------------------------





#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- plot environment, CUE, GPP and Ra (normalized by leaf area) over time
plotPAR_AirT_CUE_GPP_Ra <- function(cue.day.trt=cue.day.trt,export=F,lwidth=2.5){
  
  windows(25,20);par(mfrow=c(4,1),mar=c(0,2,0,2),oma=c(6,8,0,6),cex.axis=2,las=1)
  palette(c("black","red"))
  
  #plot temperature
  plotBy(Tair_24hrs.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",
         ylim=c(12,45),
         ylab="",cex.lab=1.5,lwd=3,axes=F)
  mtext(expression(atop(T[air],
                        (degree*C))),side=2,las=0,cex=1.2,line=3)
  legend("topleft",c("A","W"),lty=c(1,1),lwd=3,col=c("black","red"),ncol=2,bty="n",seg.len=3,cex=1.2)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "months"),
            labels = F, tcl = 0.5)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "weeks"),
            labels = F, tcl = 0.2)
  magaxis(side=c(2),labels=c(F),frame.plot=T)
  axis(side=2,at=c(20,30,40),tick=F)
  #-overlay par
  par(new=T)
  plot(PAR.mean~Date,type="l",data=subset(cue.day.trt,Water_treatment=="control" & T_treatment=="ambient"),
       ylim=c(-40,70),col="grey",lwd=2,axes=FALSE,add=T)
  magaxis(side=4,labels=F)
  axis(side=4,at=c(0,20,40,60),tick=F)
  mtext(expression(atop(PAR,
                        (mol~d^-1))),side=4,las=0,cex=1.2,line=6.5)
  legend("topright",c("PAR"),lty=c(1),lwd=3,col=c("grey"),ncol=1,bty="n",seg.len=3,cex=1.2)
  legend("bottomright","a",bty="n",inset=-0.002,cex=1.5)
  
  
  #-- plot Ra per unit leaf area
  plotBy(Ra_la.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0.25,2.8),
         ylab="Ra (gC d-1)",cex.lab=1.5,lwd=3,axes=F)
  a <- subset(cue.day.trt,T_treatment=="ambient")
  e <- subset(cue.day.trt,T_treatment=="elevated")
  adderrorbars(x=a$Date-0.25,y=a$Ra_la.mean,SE=a$Ra_la.standard.error,direction="updown",barlen=0,col="black")
  adderrorbars(x=e$Date+0.25,y=e$Ra_la.mean,SE=e$Ra_la.standard.error,direction="updown",barlen=0,col="red")
  plotBy(Ra_la.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0,90),
         cex.lab=1.5,lwd=lwidth,add=T)
  mtext(expression(atop(R[a],
                        (gC~m^-2~d^-1))),side=2,las=0,cex=1.2,line=3)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "months"),
            labels = F, tcl = 0.5,las=2)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "weeks"),
            labels = F, tcl = 0.2)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  legend("bottomright","b",bty="n",inset=-0.002,cex=1.5)
  
  #-- plot GPP per unit leaf area
  plotBy(GPP_la.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0.5,7.7),
         ylab="",cex.lab=1.5,lwd=3,axes=F)
  a <- subset(cue.day.trt,T_treatment=="ambient")
  e <- subset(cue.day.trt,T_treatment=="elevated")
  adderrorbars(x=a$Date-.25,y=a$GPP_la.mean,SE=a$GPP_la.standard.error,direction="updown",barlen=0,col="black")
  adderrorbars(x=e$Date+0.25,y=e$GPP_la.mean,SE=e$GPP_la.standard.error,direction="updown",barlen=0,col="red")
  plotBy(GPP_la.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0,90),
         lwd=lwidth,add=T)
  mtext(expression(atop(GPP,
                        (gC~m^-2~d^-1))),side=2,las=0,cex=1.2,line=3)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "months"),
            labels = F, tcl = 0.5)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "weeks"),
            labels = F, tcl = 0.2)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  legend("bottomright","c",bty="n",inset=-0.002,cex=1.5)
  

  #-- plot RtoA
  plotBy(RtoA.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0,1),
         ylab="",cex.lab=1.5,lwd=1,axes=F)
  a <- subset(cue.day.trt,T_treatment=="ambient")
  e <- subset(cue.day.trt,T_treatment=="elevated")
  adderrorbars(x=a$Date-0.25,y=a$RtoA.mean,SE=a$RtoA.standard.error,direction="updown",barlen=0,col="black")
  adderrorbars(x=e$Date+0.25,y=e$RtoA.mean,SE=e$RtoA.standard.error,direction="updown",barlen=0,col="red")
  plotBy(RtoA.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0,90),
         cex.lab=1.5,lwd=lwidth,add=T)
  mtext(expression(R[a]/GPP),side=2,las=0,cex=1.2,line=5)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "months"),
            labels = T, tcl = 0.5,las=2)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "weeks"),
            labels = F, tcl = 0.2)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  legend("bottomright","d",bty="n",inset=-0.002,cex=1.5)
  

  if(export==T) dev.copy2pdf(file="Output/Figure4.pdf")
}
#--------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------------
# plot climate drivers of GPP, Ra, and Ra/GPP (Figure 5)

plotGPP_Ra_CUE_metdrivers <- function(cue.day=cue.day,export=T,shading=0.5,parcut=35){
  
  
  palette(alpha(c("black","red"),0.5))
  
  windows(20,30);par(mfrow=c(3,2),mar=c(0,0,0,0),oma=c(7,7,1,3),las=1,cex.lab=1.8,cex.axis=1.2)
  
  #- plot GPP
  smoothplot(PAR, GPP_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0.5,9),
             data=cue.day, kgam=5, axes=FALSE)
  title(ylab=expression(GPP~(gC~m^-2~d^-1)),outer=T,adj=0.95)
  magaxis(side=1:4,labels=c(0,1,0,0))
  legend("topright","a",bty="n",inset=-0.01,cex=1.2)
  
  smoothplot(Tair_24hrs, GPP_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0.5,9),
             data=subset(cue.day,PAR>parcut), kgam=5, axes=FALSE)
  magaxis(side=1:4,labels=c(0,0,0,0))
  legend("topright","b",bty="n",inset=0.01,cex=1.2)
  title(ylab=expression(GPP~(gC~m^-2~d^-1)),outer=T,adj=0.95)
  magaxis(side=1:4,labels=c(0,0,0,1))
  
  
  #- plot Ra
  smoothplot(PAR, Ra_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0.2,4),
             data=cue.day, kgam=5, axes=FALSE)
  title(ylab=expression(R[a]~(gC~m^-2~d^-1)),outer=T,adj=0.5)
  magaxis(side=1:4,labels=c(0,1,0,0))
  legend("topright","c",bty="n",inset=-0.01,cex=1.2)
  
  smoothplot(Tair_24hrs, Ra_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0.2,4),
             data=subset(cue.day,PAR>parcut), kgam=5, axes=FALSE)
  magaxis(side=1:4,labels=c(0,0,0,0))
  legend("topright","d",bty="n",inset=0.01,cex=1.2)
  title(ylab=expression(R[a]~(gC~m^-2~d^-1)),outer=T,adj=0.5)
  magaxis(side=1:4,labels=c(0,0,0,1))
  
  #- plot CUE
  
  smoothplot(PAR, RtoA, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,1),
             data=cue.day, kgam=5, axes=FALSE)
  title(ylab=expression(R[a]/GPP),outer=T,adj=0.15)
  magaxis(side=1:4,labels=c(1,1,0,0))
  legend("topright","e",bty="n",inset=-0.01,cex=1.2)
  
  smoothplot(Tair_24hrs, RtoA, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,1),
             data=subset(cue.day,PAR>parcut), kgam=5, axes=FALSE)
  magaxis(side=1:4,labels=c(1,0,0,0))
  legend("topright","f",bty="n",inset=0.01,cex=1.2)
  title(ylab=expression(R[a]/GPP),outer=T,adj=0.15)
  magaxis(side=1:4,labels=c(1,0,0,1))
  
  title(xlab=expression(PAR~(mol~d^-1)),outer=T,adj=0.08)
  
  title(xlab=expression(T[air]~(degree*C*", 24-"*h~mean)),outer=T,adj=0.95)
  
  if(export==T) dev.copy2pdf(file="output/Figure5.pdf")
}
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
#' Function for smoothplots of GAMs. 
fitgam <- function(X,Y,dfr, k=-1, R=NULL){
  dfr$Y <- dfr[,Y]
  dfr$X <- dfr[,X]
  if(!is.null(R)){
    dfr$R <- dfr[,R]
    model <- 2
  } else model <- 1
  dfr <- droplevels(dfr)
  
  
  if(model ==1){
    g <- gam(Y ~ s(X, k=k), data=dfr)
  }
  if(model ==2){
    g <- gamm(Y ~ s(X, k=k), random = list(R=~1), data=dfr)
  }
  
  return(g)
}


#' Plot a generalized additive model
#' @param x Variable for X axis (unquoted)
#' @param y Variable for Y axis (unquoted)
#' @param data Dataframe containing x and y
#' @param kgam the \code{k} parameter for smooth terms in gam.
#' @param random An optional random effect (quoted)
#' @param log Whether to add log axes for x or y (but no transformations are done).
#' @param fitoneline Whether to fit only 
smoothplot <- function(x,y,g=NULL,data,
                       fittype=c("gam","lm"),
                       kgam=4,
                       random=NULL,
                       randommethod=c("lmer","aggregate"),
                       log="",
                       fitoneline=FALSE,
                       pointcols=NULL,
                       linecols=NULL, 
                       xlab=NULL, ylab=NULL,
                       polycolor=alpha("lightgrey",0.7),
                       axes=TRUE,
                       ...){
  
  fittype <- match.arg(fittype)
  randommethod <- match.arg(randommethod)
  
  if(!is.null(substitute(g))){
    data$G <- as.factor(eval(substitute(g),data))
  } else {
    fitoneline <- TRUE
    data$G <- 1
  }
  data$X <- eval(substitute(x),data)
  data$Y <- eval(substitute(y),data)
  data <- droplevels(data)
  
  data <- data[!is.na(data$X) & !is.na(data$Y) & !is.na(data$G),]
  
  if(is.null(pointcols))pointcols <- palette()
  if(is.null(linecols))linecols <- palette()
  
  if(is.null(xlab))xlab <- substitute(x)
  if(is.null(ylab))ylab <- substitute(y)
  
  # If randommethod = aggregate, average by group and fit simple gam.
  if(!is.null(random) && randommethod == "aggregate"){
    data$R <- data[,random]
    
    data <- summaryBy(. ~ R, FUN=mean, na.rm=TRUE, keep.names=TRUE, data=data,
                      id=~G)
    R <- NULL
  }
  
  
  if(!fitoneline){
    
    d <- split(data, data$G)
    
    if(fittype == "gam"){
      fits <- lapply(d, function(x)try(fitgam("X","Y",x, k=kgam, R=random)))
      if(!is.null(random))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- lapply(d, function(x)lm(Y ~ X, data=x))
    }
    hran <- lapply(d, function(x)range(x$X, na.rm=TRUE))
  } else {
    if(fittype == "gam"){
      fits <- list(fitgam("X","Y",data, k=kgam, R=random))
      if(!is.null(random))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- list(lm(Y ~ X, data=data))
    }
    hran <- list(range(data$X, na.rm=TRUE))
    
  }
  
  with(data, plot(X, Y, axes=FALSE, pch=16, col=pointcols[G],
                  xlab=xlab, ylab=ylab, ...))
  
  if(axes){
    if(log=="xy")magaxis(side=1:2, unlog=1:2)
    if(log=="x"){
      magaxis(side=1, unlog=1)
      axis(2)
      box()
    }
    if(log=="y"){
      magaxis(side=2, unlog=2)
      axis(1)
      box()
    }
    if(log==""){
      axis(1)
      axis(2)
      box()
    }
  }
  
  for(i in 1:length(fits)){
    
    if(fittype == "gam"){
      nd <- data.frame(X=seq(hran[[i]][1], hran[[i]][2], length=101))
      if(!inherits(fits[[i]], "try-error")){
        p <- predict(fits[[i]],nd,se.fit=TRUE)
        addpoly(nd$X, p$fit-2*p$se.fit, p$fit+2*p$se.fit, col=polycolor[i])
        lines(nd$X, p$fit, col=linecols[i], lwd=2)
      }
    }
    if(fittype == "lm"){
      pval <- summary(fits[[i]])$coefficients[2,4]
      LTY <- if(pval < 0.05)1 else 5
      predline(fits[[i]], col=linecols[i], lwd=2, lty=LTY)
    }
  }
  
  return(invisible(fits))
}


addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.7),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}

predline <- function(fit, from=NULL, to=NULL, col=alpha("lightgrey",0.7), ...){
  
  if(is.null(from))from <- min(fit$model[,2], na.rm=TRUE)
  if(is.null(to))to <- max(fit$model[,2], na.rm=TRUE)
  
  newdat <- data.frame(X = seq(from,to, length=101))
  names(newdat)[1] <- names(coef(fit))[2]
  
  pred <- as.data.frame(predict(fit, newdat, se.fit=TRUE, interval="confidence")$fit)
  
  addpoly(newdat[[1]], pred$lwr, pred$upr, col=col)
  
  #ablinepiece(fit, from=from, to=to, ...)
  lines(pred$fit~newdat[,1])
}

#'@title Add a line to a plot
#'@description As \code{abline}, but with \code{from} and \code{to} arguments. 
#'If a fitted linear regression model is used as asn argument, it uses the min and max values of the data used to fit the model.
#'@param a Intercept (optional)
#'@param b Slope (optional)
#'@param reg A fitted linear regression model (output of \code{\link{lm}}).
#'@param from Draw from this X value
#'@param to Draw to this x value
#'@param \dots Further parameters passed to \code{\link{segments}}
#'@export
ablinepiece <- function(a=NULL,b=NULL,reg=NULL,from=NULL,to=NULL,...){
  
  # Borrowed from abline
  if (!is.null(reg)) a <- reg
  
  if (!is.null(a) && is.list(a)) {
    temp <- as.vector(coefficients(a))
    from <- min(a$model[,2], na.rm=TRUE)
    to <- max(a$model[,2], na.rm=TRUE)
    
    if (length(temp) == 1) {
      a <- 0
      b <- temp
    }
    else {
      a <- temp[1]
      b <- temp[2]
    }
  }
  
  segments(x0=from,x1=to,
           y0=a+from*b,y1=a+to*b,...)
  
}

#----------------------------------------------------------------------------------------------------------------







#----------------------------------------------------------------------------------------------------------------
#- plot heat-map style hexagons for GPP relative to air temperaure and PAR (i.e., Figure 6)
plotGPP_hex <- function(dat=dat.hr.p,export=F,shading=0.7){
  
  #- Calculate GPP per unit leaf area, convert to umol CO2 m-2 s-1
  dat$GPP_la <- with(dat,GPP/leafArea)
  dat$GPP_la_umol <- with(dat,GPP_la/12*1*10^6/60/60)

  #creates a scale of colors
  myColorRamp <- function(colors, values) {
    #v <- (values - min(values))/diff(range(values))
    v <- (values - -1)/21 
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
  }
  a <- subset(dat,PAR>10 & Water_treatment=="control" & T_treatment=="ambient")[,c("Tair_al","PAR","GPP","GPP_la","GPP_la_umol")]
  e <- subset(dat,PAR>10 & Water_treatment=="control" & T_treatment=="elevated")[,c("Tair_al","PAR","GPP","GPP_la","GPP_la_umol")]
  all <- subset(dat,PAR>10 & Water_treatment=="control")[,c("T_treatment","chamber","Date","Tair_al","PAR","GPP","GPP_la","GPP_la_umol")]
  
  
  #- create hex bins for a, e and all
  ha <- hexbin(a$Tair_al, a$PAR,xbins=30,IDs=TRUE)
  he <- hexbin(e$Tair_al, e$PAR,xbins=30,IDs=TRUE)
  hall <- hexbin(all$Tair_al, all$PAR,xbins=30,IDs=TRUE)

  #average values of points inside hexbins 
  meanHexBina<-data.frame(mean=hexTapply(ha, a$GPP_la_umol, mean)) 
  meanHexBine<-data.frame(mean=hexTapply(he, e$GPP_la_umol, mean)) 
  meanHexBinall<-data.frame(mean=hexTapply(hall, all$GPP_la_umol, mean)) 

  
  colsa <- myColorRamp(c("blue","green","yellow", "red"), meanHexBina$mean)
  colse <- myColorRamp(c("blue","green","yellow", "red"), meanHexBine$mean)

  #-- plot Ambient
  windows()
  
  ## setup coordinate system of the plot
  par(cex.lab=2,cex.axis=2)
  Pa <- plot(hall, type="n",legend=FALSE,xlab="",ylab="",main="")
  
  ##add hexagons (in the proper viewport):
  pushHexport(Pa$plot.vp)
  
  #plots hexbins based on colors of third column
  grid.hexagons(ha, style= "lattice", border = gray(.9), pen = colsa,  minarea = 1, maxarea = 1)
  if(export==T) dev.copy2pdf(file="output/Figure6a.pdf")
  
  #-- plot elevated
  windows()
  
  ## setup coordinate system of the plot
  Pe <- plot(hall, type="n",legend=FALSE,xlab="",ylab="",main="")
  
  ##add hexagons (in the proper viewport):
  pushHexport(Pe$plot.vp)
  
  #plots hexbins based on colors of third column
  grid.hexagons(he, style= "lattice", border = gray(.9), pen = colse,  minarea = 1, maxarea = 1)
  if(export==T) dev.copy2pdf(file="output/Figure6b.pdf")
  
  
  
  #- get five hexagons, add them as a legend
  windows()
  values <- seq(0,20,length=6)
  colors <- myColorRamp(c("blue","green","yellow", "red"), values)
  xloc <- rep(1,6)
  yloc <- c(1,2,3,4,5,6)
  plot(yloc~xloc,pch=18,cex=5,col=colors,ylim=c(0,10),axes=F,xlab="",ylab="")
  text(xloc+0.1,yloc,labels=values,cex=1.5)
  #dev.copy2pdf(file="/output/Figure6_legend.pdf")
  
  
  
  #-- subset data to a par range, plot Tresponse
  toplot <- subset(all, PAR >1200 & PAR <1500 & GPP_la_umol>0)
  
  #- plot smoothplots
  windows();par(las=1)
  smoothplot(Tair_al, GPP_la_umol, T_treatment,pointcols=c(alpha("black",0.3),alpha("red",0.3)),
             linecol=c("black","red"),polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",cex=1,main="",
             xlim=c(15,45),ylim=c(0,25),xlab="",ylab="",
             data=toplot, kgam=4, axes=FALSE)
  legend("topright",pch=16,legend=c("ambient","warmed"),col=c(alpha("black",0.3),alpha("red",0.3)))
  box();axis(side=1,labels=T);axis(side=2,labels=T);axis(side=4,labels=F)
  if(export==T) dev.copy2pdf(file="output/Figure6c.pdf")
  
}
