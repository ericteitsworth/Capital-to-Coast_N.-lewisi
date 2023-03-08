#Eric Teitsworth
#12/22/2020

#This script is a cleaned version of all the code leading up to this "final" analysis.
#Occupancy model using detection data from 2018-2019 and 2019-2020 sample
#seasons. Occupancy covariates are land uses, binned distances to potential pollution
#sources, benthic rating, and metrics associated with channel 
#stability. Detection covariates are daily temperature and precipitation.

rm(list=ls())

#SetWD
setwd("D:/NCSU/R/19_20Season Analysis")

library(ggplot2)
library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(corrplot)

#test for correlations in covariates
source("http://www.sthda.com/upload/rquery_cormat.r")

# Det.data can be 0 or 1 and from 8 Visits
#HAiFLS = Hydroactive inverse Flow to Stream. This is given 
# as a proportion of the watershed, clipped to exclude everything upstream of dams. 
det.data<-read.csv("cap_hist19_20.csv",header=TRUE)
surv.cov2<-read.csv("samp_cov3.csv",header=TRUE) #long format # values on day prior to capture
site.cov<-read.csv('site_covs19_20.csv',header=TRUE)

##Covariates##
#Site
HAiFLSd<-site.cov[, (5)] #developed LULC
HAiFLSf<-site.cov[, (6)] #forested LULC
HAiFLSg<-site.cov[, (7)] #grass/pasture LULC
HAiFLSc<-site.cov[, (8)] #crop LULC
HAiFLSw<-site.cov[, (9)] #wetland LULC
benthic<-site.cov[, (11)] #DEQ benthic rating
UpDamBin<-site.cov[, (12)] #binned distance to upstream dam
DownDamBin<-site.cov[, (13)] #binned distance to downstream dam
CSA<-site.cov[, (16)] #channel stability assessment
NPDES<-site.cov[, (17)] #number of NPDES sites in watershed/area=proportional number
stage<-site.cov[, (22)] #score 1-6 as measure of modification
stage2<-stage^2
substrate<- site.cov[, (23)] #categorical: dominant substrate at location
TSS<- site.cov[, (24)] #predicted total suspended solids from USGS SPARROW
TSS2<- TSS^2
ATSS<-site.cov[, (25)]
phos<-site.cov$Phosphorus #predicted total phosphorus loads from USGS SPARROW
nitro<-site.cov$Nitrogen #predicted total nitrogen loads from USGS SPARROW
phosw<-site.cov$weightedphos #phosphorus loads weighted by watershed size
nitrow<-site.cov$weightednitro #nitrogen loads weighted by watershed size
LULC<-site.cov$LULC.Sum #Total LULC percentages

#Survey
MaxAir2<-surv.cov2[, (1)] #Daily high prior to trap night
MinAir2<-surv.cov2[, (2)] #Daily low prior to trap night
precip2<-surv.cov2[, (3)] #total daily precip prior to trap night
MinAirquad2<-MinAir2^2
precipquad2<-precip2^2
count<-surv.cov2[, (4)] #raw count of nele caught by site
s1<-surv.cov2[, (6)]
s2<-surv.cov2[, (7)]

##vif testing
sitevif<- site.cov[, c(5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16)]
vif(sitevif)


##Creating the unmarked dataframe
nele.umf2<-unmarkedFrameOccu(y      =det.data,
                             siteCovs = site.cov,
                             obsCovs = surv.cov2)

##Creating the global occu model
# This model is not testing substrate, TSS, phosphorus, or nitrogen loads
# This model includes oneway interactions between all landcover types, excluding forest (overdispersed)
global.model<- occu(~MinAir2 + MinAirquad2 + precip2 + precipquad2 + s1 + s2
                     ~UpDamBin + NPDES + HAiFLSf
                     + HAiFLSd*HAiFLSw + HAiFLSc*HAiFLSw 
                     + HAiFLSg*HAiFLSw + DownDamBin
                     + stage + stage2 + substrate + benthic
                     + HAiFLSd*HAiFLSc + HAiFLSd*HAiFLSg + HAiFLSc*HAiFLSg
                     , nele.umf2, control=list(maxit=1000))
#removed MaxAir, benthic to allow for dummy variables and not be singular

#Test GOF of global model
gof<- mb.gof.test(global.model, nsim=999, parallel = FALSE)
#SS=605.58, p=0.335, c-hat=1.01

##Before Dredging
options(na.ation = "na.fail")
##Setting back to default
options(na.action = "na.omit")

## The below lines limits the dredge of the global
## to a max of 1 interaction per model, only allow interactions when both
## main effects are present, only allow quadratic effects when
## the main effect is present, and do not allow developed and forested 
## land cover types in the same model.
x<- getAllTerms(global.model)
x<- x[grep(":", x)] ##this doesn't work if there are ":" in the vector name

ss<- substitute((dc("psi(HAiFLSd)"+"psi(HAiFLSw)", "psi(HAiFLSd*HAiFLSw)")) #requires main effects with interaction
                & (dc("psi(HAiFLSc)"+"psi(HAiFLSw)", "psi(HAiFLSc*HAiFLSw)")) #requires main effects with interaction
                & (dc("psi(HAiFLSg)"+"psi(HAiFLSw)", "psi(HAiFLSg*HAiFLSw)")) #requires main effects with interaction
                & (dc("psi(HAiFLSd)"+"psi(HAiFLSc)", "psi(HAiFLSd*HAiFLSc)")) #requires main effects with interaction
                & (dc("psi(HAiFLSd)"+"psi(HAiFLSg)", "psi(HAiFLSd*HAiFLSg)")) #requires main effects with interaction
                & (dc("psi(HAiFLSc)"+"psi(HAiFLSg)", "psi(HAiFLSc*HAiFLSg)")) #requires main effects with interaction
                & (dc("p(precip2)", "p(precipquad2)")) #requires main effect with quadratic
                & (dc("p(MinAir2)", "p(MinAirquad2)")) #requires main effect with quadratic
                & (dc("p(s1)", "p(s2)"))
                & (dc("p(s2)", "p(s1)"))
                & (dc("psi(stage)", "psi(stage2)")) #requires main effect with quadratic
                & !("psi(HAiFLSd)" && "psi(HAiFLSf)") #does not allow forest and developed LULC in same model
                & !("p(s1)" && "p(MinAir2)")
                & !("p(s1)" && "p(precip2)")
                & !("p(s2)" && "p(MinAir2)")
                & !("p(s2)" && "p(precip2)")
                & (sum(x) < N), list(N = 2, x = as.call(lapply(c("c", x), as.symbol)))) #only allows 1 interaction per fit model

#This is the dredge that works. More than 8 terms per model (plus 
# intercepts) would produce NaNs with the number of base terms I have 
# chosen to include in the global model. Interactions of crop and 
# pasture with forest were removed because they always produced NaNs.
# I think forest cover is too overdispersed. 
#This produced a working dredge of 23303 models (31 covs), the top 6 of which are
# within delta 2.0, top 16 within delta 4.0
dd<- dredge(global.model, rank = "AIC", m.lim= c(0,8), subset = ss) #NOT USED, using AICc instead

ddc<- dredge(global.model, rank = "AICc", m.lim= c(0,8), subset = ss)

#This averages the all models within delta 4.0. Recommended to use
#full estimate (vs. conditional) because it is more conservative
modavg.dd<- model.avg(ddc, subset = delta < 4, fit=TRUE)
summary(modavg.dd)

modavg.2<- model.avg(dd, subset = delta <2, fit=TRUE)
summary(modavg.2)

#predict model average result
pmodavg<-predict(modavg.dd, se.fit=TRUE, type="state", backTransform=TRUE)
predict.det<-predict(modavg.dd, type="det", se.fit=TRUE, backTransform=TRUE)
write.csv(pmodavg, file = "D:/NCSU/R/19_20Season Analysis/pmodavg_ch1.csv")
write.csv(predict.det, file = "D:/NCSU/R/19_20Season Analysis/predictdet_ch1.csv")

####Predicting Effect of Pasture########
newdata<-with(site.cov,
              data.frame(HAiFLSg = seq(from = min(site.cov$HAiFLSg), #xaxis variable
                                       to = max(site.cov$HAiFLSg), length.out=116),
                         HAiFLSd = mean(HAiFLSd),
                         HAiFLSc = mean(HAiFLSc),
                         HAiFLSw = mean(HAiFLSw),
                         HAiFLSf = mean(HAiFLSf),
                         NPDES = mean(NPDES), 
                         stage = median(stage),
                         stage2 = median(stage2),
                         substrate = median(substrate),
                         benthic = median(benthic),
                         UpDamBin = median(UpDamBin),
                         DownDamBin = median(DownDamBin),
                         "HAiFLSc:HAiFLSw" = mean(HAiFLSw*HAiFLSc),
                         "HAiFLSg:HAiFLSw" = mean(HAiFLSg*HAiFLSw)))


#make fit predictions based on newdata model
pred.g<-predict(modavg.dd, newdata, full=TRUE, se.fit=FALSE, type="state", backTransform=TRUE)
newdata$fit <- pred.g$fit
newdata$SE <- pred.g$se
newdata$upr=newdata$fit+1.96*newdata$SE
newdata$lwr=newdata$fit-1.96*newdata$SE

#Plot?
ggplot(data=newdata, aes(HAiFLSg, y=fit)) +
  geom_line() +
  labs(x = "% Pasture", y = "psi")+
  ylim(-0.2,1.1) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold")) +
  theme_classic(base_size = 12) +
  geom_ribbon(aes(ymin=upr, ymax=lwr), alpha=0.25)

####Predicting Effect of Substrate########
newdata<-with(site.cov,
              data.frame(substrate = rep(seq(from = 0, #xaxis variable
                         to = 3, length.out=116))),
                         HAiFLSd = mean(HAiFLSd),
                         HAiFLSc = mean(HAiFLSc),
                         HAiFLSw = mean(HAiFLSw),
                         HAiFLSf = mean(HAiFLSf),
                         NPDES = mean(NPDES), 
                         stage = median(stage),
                         stage2 = median(stage2),
                         HAiFLSg = mean(HAiFLSg),
                         benthic = median(benthic),
                         UpDamBin = median(UpDamBin),
                         DownDamBin = median(DownDamBin),
                         "HAiFLSc:HAiFLSw" = mean(HAiFLSw*HAiFLSc),
                         "HAiFLSg:HAiFLSw" = mean(HAiFLSg*HAiFLSw))


#make fit predictions based on newdata model
pred.sub<-predict(modavg.dd, newdata, full=TRUE, se.fit=FALSE, type="state", backTransform=TRUE)
newdata$fit <- pred.sub$fit
newdata$SE <- pred.sub$se
newdata$upr=newdata$fit+1.96*newdata$SE
newdata$lwr=newdata$fit-1.96*newdata$SE

#Plot?
ggplot(data=newdata, aes(substrate, y=fit)) +
  geom_smooth(method = "lm", level=0.95, colour = "black") +
  labs(x = "Substrate Size Class", y = "psi")+
  ylim(-0.2,1.1) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold")) +
  theme_classic(base_size = 12)

####Predicting Effect of Developed########
newdata<-with(site.cov,
              data.frame(HAiFLSd = seq(from = min(site.cov$HAiFLSd), #xaxis variable
                                       to = max(site.cov$HAiFLSd), length.out=116),
                         HAiFLSg = mean(HAiFLSg),
                         HAiFLSc = mean(HAiFLSc),
                         HAiFLSw = mean(HAiFLSw),
                         HAiFLSf = mean(HAiFLSf),
                         NPDES = mean(NPDES), 
                         stage = median(stage),
                         stage2 = median(stage2),
                         substrate = median(substrate),
                         benthic = median(benthic),
                         UpDamBin = median(UpDamBin),
                         DownDamBin = median(DownDamBin),
                         "HAiFLSc:HAiFLSw" = mean(HAiFLSw*HAiFLSc),
                         "HAiFLSg:HAiFLSw" = mean(HAiFLSg*HAiFLSw)))


#make fit predictions based on newdata model
pred.d<-predict(modavg.dd, newdata, full=TRUE, se.fit=FALSE, type="state", backTransform=TRUE)
newdata$fit <- pred.d$fit
newdata$SE <- pred.d$se
newdata$upr=newdata$fit+1.96*newdata$SE
newdata$lwr=newdata$fit-1.96*newdata$SE

#Plot?
ggplot(data=newdata, aes(HAiFLSd, y=fit)) +
  geom_line() +
  labs(x = "% Developed", y = "psi")+
  ylim(-0.2,1.1) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold")) +
  theme_classic(base_size = 12) +
  geom_ribbon(aes(ymin=upr, ymax=lwr), alpha=0.25)

####Predicting Effect of Wetland########
newdata<-with(site.cov,
              data.frame(HAiFLSw = seq(from = min(site.cov$HAiFLSw), #xaxis variable
                                       to = max(site.cov$HAiFLSw), length.out=116),
                         HAiFLSg = mean(HAiFLSg),
                         HAiFLSc = mean(HAiFLSc),
                         HAiFLSd = mean(HAiFLSd),
                         HAiFLSf = mean(HAiFLSf),
                         NPDES = mean(NPDES), 
                         stage = median(stage),
                         stage2 = median(stage2),
                         substrate = median(substrate),
                         benthic = median(benthic),
                         UpDamBin = median(UpDamBin),
                         DownDamBin = median(DownDamBin),
                         "HAiFLSc:HAiFLSw" = mean(HAiFLSw*HAiFLSc),
                         "HAiFLSg:HAiFLSw" = mean(HAiFLSg*HAiFLSw)))


#make fit predictions based on newdata model
pred.w<-predict(modavg.dd, newdata, full=TRUE, se.fit=FALSE, type="state", backTransform=TRUE)
newdata$fit <- pred.w$fit
newdata$SE <- pred.w$se
newdata$upr=newdata$fit+1.96*newdata$SE
newdata$lwr=newdata$fit-1.96*newdata$SE

#Plot?
ggplot(data=newdata, aes(HAiFLSw, y=fit)) +
  geom_line() +
  labs(x = "% Wetland", y = "psi")+
  ylim(0,1.5) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold")) +
  theme_classic(base_size = 12) +
  geom_ribbon(aes(ymin=upr, ymax=lwr), alpha=0.25)


####### Scenario Predictions ######
#Scenario 1 "Best Case" - 50% reduction in pasture, gained by forest
scenario1<- with(site.cov,
                 data.frame(HAiFLSg = HAiFLSg*0.5,
                             HAiFLSf = HAiFLSf + 0.5*HAiFLSg,
                             HAiFLSd = HAiFLSd,
                             HAiFLSc = HAiFLSc,
                             HAiFLSw = HAiFLSw,
                             stage = stage,
                             stage2 = stage2,
                             substrate = substrate,
                             NPDES = NPDES,
                             UpDamBin = UpDamBin,
                             benthic = benthic,
                             DownDamBin = DownDamBin))

  # Make fit predictions based on Scenario 1 dataframe
  scenario1.pred<-predict(modavg.dd, scenario1, full=TRUE, se.fit=FALSE, type="state", backTransform=TRUE)
  scenario1$fit <- scenario1.pred$fit
  scenario1$SE <- scenario1.pred$se
  scenario1$upr=scenario1$fit+1.96*scenario1$SE
  scenario1$lwr=scenario1$fit-1.96*scenario1$SE

#Scenario 2 "Worst Case" - 100% increase in Developed, 50% increase in Pasture, offset losses in crop and forest
scenario2<- with(site.cov,
                 data.frame(HAiFLSg = HAiFLSg + HAiFLSg*0.5,
                            HAiFLSf = ifelse(HAiFLSf-0.5*HAiFLSd-0.25*HAiFLSg < 0, 0, HAiFLSf-0.5*HAiFLSd-0.25*HAiFLSg), 
                            HAiFLSd = ifelse(HAiFLSd*2 > site.cov$LULC.Sum, site.cov$LULC.Sum, HAiFLSd*2), #developed LULC cannot exceed total sum % possible
                            HAiFLSc = ifelse(HAiFLSc-0.5*HAiFLSd-0.25*HAiFLSg < 0, 0, HAiFLSc-0.5*HAiFLSd-0.25*HAiFLSg),
                            HAiFLSw = HAiFLSw,
                            stage = stage,
                            stage2 = stage2,
                            substrate = substrate,
                            NPDES = NPDES,
                            UpDamBin = UpDamBin,
                            benthic = benthic,
                            DownDamBin = DownDamBin))
    #with all of the changes to percentages found above, we do not want to 
    #increase the total possible % from truth to scenario. This statement says
    #"if the sum total % of LULC types now exceeds the true sum total %, 
    #subtract the difference in sum totals from forest LULC" because in 
    #reality this is where it would probably come from 
  scenario2$HAiFLSf<- ifelse(scenario2$HAiFLSg + scenario2$HAiFLSf + scenario2$HAiFLSd + scenario2$HAiFLSc + scenario2$HAiFLSw > site.cov$LULC.Sum, 
                           scenario2$HAiFLSf - ((scenario2$HAiFLSg + scenario2$HAiFLSf + scenario2$HAiFLSd + scenario2$HAiFLSc + scenario2$HAiFLSw) - site.cov$LULC.Sum),
                           scenario2$HAiFLSf)  
        #sometimes this makes forest cover negative, which is not possible. In this case, 
        #make it 0 and subtract that negative value from developed for the sake of being realistic
        scenario2$HAiFLSd<- ifelse(scenario2$HAiFLSf<0, scenario2$HAiFLSd + scenario2$HAiFLSf, scenario2$HAiFLSd)
        scenario2$HAiFLSf<- ifelse(scenario2$HAiFLSf<0, 0, scenario2$HAiFLSf)
    #series saying "if developed LULC equals the total sum % possible in 
    #this scenario, set this LULC to 0" ##### Not pretty, but it works...
  scenario2$HAiFLSc<- ifelse(scenario2$HAiFLSd == site.cov$LULC.Sum, 0, scenario2$HAiFLSc)
  scenario2$HAiFLSw<- ifelse(scenario2$HAiFLSd == site.cov$LULC.Sum, 0, scenario2$HAiFLSw)
  scenario2$HAiFLSg<- ifelse(scenario2$HAiFLSd == site.cov$LULC.Sum, 0, scenario2$HAiFLSg)
  scenario2$HAiFLSf<- ifelse(scenario2$HAiFLSd == site.cov$LULC.Sum, 0, scenario2$HAiFLSf)
   
     # Make fit predictions based on Scenario 2 dataframe
  scenario2.pred<-predict(modavg.dd, scenario2, full=TRUE, se.fit=FALSE, type="state", backTransform=TRUE)
  scenario2$fit <- scenario2.pred$fit
  scenario2$SE <- scenario2.pred$se
  scenario2$upr=scenario2$fit+1.96*scenario2$SE
  scenario2$lwr=scenario2$fit-1.96*scenario2$SE

#Scenario3 "Expected Case" - 100% increase in Developed LULC, but 50% decrease in Pasture LULC.
    #Forest LULC first gains the loss in Pasture, then Forest and Crop split a reduction that 
    #occurs from the gain in Developed.
scenario3<- with(site.cov,
            data.frame(HAiFLSg = HAiFLSg*0.5,
             HAiFLSf = (HAiFLSf + HAiFLSg*0.5) - HAiFLSd*0.5, 
             HAiFLSd = ifelse(HAiFLSd*2 > site.cov$LULC.Sum, site.cov$LULC.Sum, HAiFLSd*2), #Developed LULC not allowed to exceed maximum possible LULC % 
             HAiFLSc = HAiFLSc - HAiFLSd*0.5,
             HAiFLSw = HAiFLSw,
             stage = stage,
             stage2 = stage2,
             substrate = substrate,
             NPDES = NPDES,
             UpDamBin = UpDamBin,
             benthic = benthic,
             DownDamBin = DownDamBin))
   
     #with all of the changes to percentages found above, we do not want to 
    #increase the total possible % from truth to scenario. This statement says
    #"if the sum total % of LULC types now exceeds the true sum total %, 
    #subtract the difference in sum totals from forest LULC" because in 
    #reality this is where it would probably come from 
    scenario3$HAiFLSf<- ifelse(scenario3$HAiFLSg + scenario3$HAiFLSf + scenario3$HAiFLSd + scenario3$HAiFLSc + scenario3$HAiFLSw > site.cov$LULC.Sum, 
                           scenario3$HAiFLSf - ((scenario3$HAiFLSg + scenario3$HAiFLSf + scenario3$HAiFLSd + scenario3$HAiFLSc + scenario3$HAiFLSw) - site.cov$LULC.Sum),
                           scenario3$HAiFLSf)
    #sometimes this makes forest and crop LULC negative, which is not possible. In this case, 
    #make it 0 and subtract that negative value from developed for the sake of being realistic
    #Forest fix 
    scenario3$HAiFLSd<- ifelse(scenario3$HAiFLSf<0, scenario3$HAiFLSd + scenario3$HAiFLSf, scenario3$HAiFLSd)
    scenario3$HAiFLSf<- ifelse(scenario3$HAiFLSf<0, 0, scenario3$HAiFLSf)
    #Crop fix
    scenario3$HAiFLSd<- ifelse(scenario3$HAiFLSc<0, scenario3$HAiFLSd + scenario3$HAiFLSc, scenario3$HAiFLSd)
    scenario3$HAiFLSc<- ifelse(scenario3$HAiFLSc<0, 0, scenario3$HAiFLSc)

    # Make fit predictions based on Scenario 2 dataframe
    scenario3.pred<-predict(modavg.dd, scenario3, full=TRUE, se.fit=FALSE, type="state", backTransform=TRUE)
    scenario3$fit <- scenario3.pred$fit
    scenario3$SE <- scenario3.pred$se
    scenario3$upr=scenario3$fit+1.96*scenario3$SE
    scenario3$lwr=scenario3$fit-1.96*scenario3$SE
    
### Throw all 3 scenarios into csv file for visualization in ArcMAP
    write.csv(cbind(scenario1.pred$fit, scenario1.pred$upr, scenario1.pred$lwr,
                      scenario2.pred$fit, scenario2.pred$upr, scenario2.pred$lwr,
                      scenario3.pred$fit, scenario3.pred$upr, scenario3.pred$lwr),
                      file = "C:/Users/ewteitsw/Documents/R/19_20Season Analysis/3scenario_preds.csv")
    