
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Statistical analysis of Ra, GPP, and Ra/GPP (Table 1).
#- Note that many of these analyses take considerable time to run. The entire script may take 10-20 minutes on a typical machine.



#- Here is a manual calculation of the expected degrees of freedom.

# T = number of timepoints. w= number of warming treatments. c= number of replicate chambers within each treatment

# Source                df calculation               df  
# Warming                  w-1                        1
# chamber(Warming)         w(c-1)                    10   # this is the random whole-plot error term to test warming main effect
# Time                    (T-1)                     254   # in this case I have 255 timepoints
# Time*Warming             (T-1)*(w-1)              254
# Time*chamber[warming]  (T-1)*(c-1)*w             2540  # this is the random sub-plot error term. it's the same as the residual
#                                                               , as there is no sub-replication here.. 254*5*2

#- note the random sub-plot error term will actually have fewer df than this, as I am excluding the drought data



dat <- subset(cue.day,select=c("Date","T_treatment","Water_treatment","chamber","CUE","RtoA","GPP_la","Ra_la","PAR","leafArea"))
dat2 <- subset(dat,Water_treatment=="control")
dat2$T_treatment <- as.factor(dat2$T_treatment)
dat2$DateFac <- as.factor(dat2$Date)


#- create "Explicitly nested" random factors. These will become the error terms in the ANOVA
dat2$plotEN <- with(dat2,interaction(T_treatment,chamber))
dat2$dateEN <- with(dat2,interaction(DateFac,plotEN)) #- this doesn't make a difference. It's the same as the residual

#####
#- Ra

#- find the "right" Box-Cox transformation for Ra_la
sp.Ra.simple <- lm(Ra_la~T_treatment*DateFac,data=dat2)
value.r <- boxcox(object=sp.Ra.simple,lambda=seq(-1,0,1/20))
exponent.r <- value.r$x[which.max(value.r$y)]
dat2$Ra_latrans <- with(dat2,Ra_la^exponent.r)

#- compare models with and without autocorrelation
sp.ra <- lme(Ra_latrans~T_treatment*DateFac,random=list(~1|chamber),
             weights=varIdent(form=~1|T_treatment),data=dat2,method="ML")
sp.ra.ar1 <- update(sp.ra,correlation=corAR1(0.7,form=~1|chamber),method="ML")
AIC(sp.ra,sp.ra.ar1)
anova(sp.ra,sp.ra.ar1) # model with autocorrelation is immensely better!

#- refit best model with REML
sp.ra.ar1.reml <- update(sp.ra.ar1,method="REML")

#look at model diagnostics
plot(sp.ra.ar1,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.ra.ar1,Ra_latrans~fitted(.)|chamber,abline=c(0,1))           #predicted vs. fitted for each chamber
plot(sp.ra.ar1,Ra_latrans~fitted(.),abline=c(0,1))                   #predicted vs. fitted
qqnorm(sp.ra.ar1, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot. Departure at high values.
anova(sp.ra.ar1)

#- exclude outliers and refit
outliers <- unname(which(abs(residuals(sp.ra.ar1,type="normalized"))>qnorm(0.999)))
noout <- dat2[-outliers,]
sp.ra.ar1.noout <- lme(Ra_latrans~T_treatment*DateFac,random=list(~1|chamber),
                       correlation=corAR1(0.7,form=~1|chamber),method="REML",
                       weights=varIdent(form=~1|T_treatment),
                       data=noout)
plot(sp.ra.ar1.noout,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.ra.ar1.noout,Ra_latrans~fitted(.)|chamber,abline=c(0,1))           #predicted vs. fitted for each chamber
plot(sp.ra.ar1.noout,Ra_latrans~fitted(.),abline=c(0,1))                   #predicted vs. fitted
qqnorm(sp.ra.ar1.noout, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot. Better.
anova(sp.ra.ar1.noout)                                     #inference doesn't change when excluding outliers

#- compare models with and without interaction terms
sp.ra.ar1.noout.noint <- update(sp.ra.ar1.noout,.~.-T_treatment:DateFac)
AIC(sp.ra.ar1.noout,sp.ra.ar1.noout.noint) # dropping the interaction results in a more parsimonious model (huge reduction in df)

#- estimate explainatory power (r2 values) of models with and without the warming by time interaction
sem.model.fits(list(sp.ra.ar1.noout,sp.ra.ar1.noout.noint))                        

#- Warming by date interaction is "significant" for Ra_la, but excluding it drops the explanatory power of the model
#    very modestly. Marginal r2 value declined from 0.70 to 0.68.
#- Therefore the warming by date interaction is not actually that important in a quantatiative sense.
####






####
#- GPP
sp.gpp <- lme(GPP_la~T_treatment*DateFac,random=list(~1|chamber),
              weights=varIdent(form=~1|T_treatment),data=dat2,method="ML")
sp.gpp.ar1 <- update(sp.gpp,correlation=corAR1(0.7,form=~1|chamber),method="ML")
AIC(sp.gpp,sp.gpp.ar1)
anova(sp.gpp,sp.gpp.ar1) # model with autocorrelation is immensely better!

#- refit best model with REML
sp.gpp.ar1.reml <- update(sp.gpp.ar1,method="REML")

#look at model diagnostics
plot(sp.gpp.ar1.reml,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.gpp.ar1.reml,GPP_la~fitted(.)|chamber,abline=c(0,1))               #predicted vs. fitted for each chamber
plot(sp.gpp.ar1.reml,GPP_la~fitted(.),abline=c(0,1))                       #predicted vs. fitted
qqnorm(sp.gpp.ar1.reml, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot
hist(sp.gpp.ar1.reml$residuals)
anova(sp.gpp.ar1.reml)
lsmeans(sp.gpp.ar1.reml,"T_treatment")

#- compare models with and without interaction terms
sp.gpp.ar1.noint <- update(sp.gpp.ar1,.~.-T_treatment:DateFac)
anova(sp.gpp.ar1,sp.gpp.ar1.noint) # dropping the interaction results in a poorer model

#- estimate explainatory power (r2 values) of models with and without the warming by time interaction
sem.model.fits(list(sp.gpp.ar1,sp.gpp.ar1.noint))   

#- so the interaction between T_treatment and date is "significant" for GPP, but it's not very important.
#- Marginal r2 drops from 0.87 to 0.83 when excluding the interaction
####







#####
#- Ra/GPP

#- find the "right" box-cox transformation
sp.CUE.simple <- lm(RtoA~T_treatment*DateFac,data=dat2)
value <- boxcox(object=sp.CUE.simple,lambda=seq(0,1,1/20))
exponent.cue <- value$x[which.max(value$y)]

dat2$RtoAtrans <- with(dat2,RtoA^exponent.cue) # get the "best" transformation 

#- compare models with and without autocorrelation
sp.cue <- lme(RtoAtrans~T_treatment*DateFac,random=list(~1|chamber),
              weights=varIdent(form=~1|T_treatment),data=dat2,method="ML")
sp.cue.ar1 <- update(sp.cue,correlation=corAR1(0.7,form=~1|chamber),method="ML")
AIC(sp.cue,sp.cue.ar1)
anova(sp.cue,sp.cue.ar1) # model with autocorrelation is immensely better!


#- refit best model with REML
sp.cue.ar1.reml <- update(sp.cue.ar1,method="REML")

#look at model diagnostics
plot(sp.cue.ar1.reml,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.cue.ar1.reml,RtoAtrans~fitted(.)|chamber,abline=c(0,1))            #predicted vs. fitted for each chamber
plot(sp.cue.ar1.reml,RtoAtrans~fitted(.),abline=c(0,1))                    #predicted vs. fitted
qqnorm(sp.cue.ar1.reml, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot. Departure at both ends
hist(sp.cue.ar1.reml$residuals)


#- exclude outliers and refit
outliers <- unname(which(abs(residuals(sp.cue.ar1.reml,type="normalized"))>qnorm(0.99)))
noout <- dat2[-outliers,]
sp.cue.ar1.reml.noout <- lme(RtoAtrans~T_treatment*DateFac,random=list(~1|chamber),
                             correlation=corAR1(0.7,form=~1|chamber),method="REML",
                             weights=varIdent(form=~1|T_treatment),
                             data=noout)
plot(sp.cue.ar1.reml.noout,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.cue.ar1.reml.noout,RtoAtrans~fitted(.)|chamber,abline=c(0,1))           #predicted vs. fitted for each chamber
plot(sp.cue.ar1.reml.noout,RtoAtrans~fitted(.),abline=c(0,1))                   #predicted vs. fitted
qqnorm(sp.cue.ar1.reml.noout, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot. Better.
anova(sp.cue.ar1.reml.noout)                                     #inference doesn't change when excluding outliers


anova(sp.cue.ar1.reml)
lsmeans(sp.cue.ar1.reml,"T_treatment") 

#- compare models with and without interaction terms
sp.cue.ar1.noint <- update(sp.cue.ar1,.~.-T_treatment:DateFac)
anova(sp.cue.ar1,sp.cue.ar1.noint) # dropping the interaction results in a poorer model

#- estimate explainatory power (r2 values) of models with and without the warming by time interaction
sem.model.fits(list(sp.cue.ar1,sp.cue.ar1.noint))   

#- so the warming by date interaction is "significant" and leads to a better model, but it is quantitatively of little 
#-  importance, as the marginal r2 drops from 0.83 to 0.82 when excluding the interaction.
####









#- merge models together to make Table 1
table.r1 <- as.matrix(anova(sp.cue.ar1.reml))
table.r2 <- cbind(table.r1,as.matrix(anova(sp.gpp.ar1.reml))[1:4,3:4])
table1 <- cbind(table.r2,as.matrix(anova(sp.cue.ar1.reml))[1:4,3:4])[2:4,]

table1 


####
#-- analysis of just the exceedingly hot days
hotDates <- unique(dat.hr.p[which(dat.hr.p$Tair_al>40),"Date"]) # find dates with temperatures exceeding 40 deg C
dat.hr.p.hot <- subset(dat.hr.p,Date %in% hotDates)

#- daily climate metrics
hotDates_met <- summaryBy(Tair_al~T_treatment+Date,data=dat.hr.p.hot,FUN=c(mean,min,max),na.rm=T)
summaryBy(Tair_al.max~T_treatment,data=hotDates_met) # average maximum temperature on these hot dates

#- re-analyze on hot dates only. Capture heteroskedasticity across dates
sp.CUE.hot <- lme(RtoA~T_treatment*DateFac,random=list(~1|chamber),
                  weights=varIdent(form=~1|DateFac),data=subset(dat2,Date %in% hotDates))

#look at model diagnostics
plot(sp.CUE.hot,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.CUE.hot,RtoA~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each species
plot(sp.CUE.hot,RtoA~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.CUE.hot, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
anova(sp.CUE.hot)
lsmeans(sp.CUE.hot,"T_treatment") 
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------