###########################################################################################################################################################
#Codes for "Non-optimum ambient temperature and cause-specific hospital admissions for cardiovascular diseases in the context of particulate air pollution"
#Authors for codes: Huimeng Liu, Jian Lei, Yunxing Jiang, Lijun Bai et.al.
#Correspondence to Shaowei Wu, Xi Chen.
##########################################################################################################################################################

######################################################################################################################################
#Numbers and fractions of hospital admissions for total major cardiovascular events attributable to non-optimum temperatures
########################################################################################################
remove(list = ls())
library(tidyverse)
library(dlnm) 
library(mvmeta) 
library(splines)
library(survival)
library(mgcv)
########mgcv###############STEP 1#########################
#IMPORT DATA
data_use <- read.csv("E:/OneDrive/11_Doctor research/17.PM2.5组分健康效应/住院/温度code/Github code1126/Temperature-and-modification-by-PM2.5/test_data_fin.csv") %>% 
  mutate(date = as.Date(date))
regions <- unique(data_use$district)


#ARRANGE THE DATA AS A LIST OF DATA SETS
data <- lapply(regions,function(x) data_use[data_use$district==x,])
names(data) <- regions

#SET EMPTY MATRIX TO STORE THE RESULTSyall<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
yall<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
Sall <- vector("list",length(regions))
names(Sall) <- regions


for(i in regions) {
  # PRINT
  cat(i,"")
  
  #GENERATRE COVARIATES
  sub$year <- year(sub$date)
  sub$date1 <- 1:length(sub$date)
  sub$rhu <- sub$l0rhu
  sub$dow    <- as.factor(weekdays(sub$date))
  #######################################################################################
  #TRIM THE OUTLIERS 
  
  #DEFINE THE CROSS-BASIS
  lag <- 28
  
  arglag <- list(knots=logknots(lag,3))
  argvar <- list(fun="ns",df=3) 
  cb <- crossbasis(sub$l0avg_temp, lag = lag, argvar=argvar, arglag=arglag)

    tryCatch({
      mfirst <- gam(PSN_NO~cb+ ns(rhu, df=3) +ns(date1,df=8*length(unique(year)))
                    +as.factor(Holiday)
                    +as.factor(dow),
                    family=quasipoisson,sub,na.action="na.exclude")
  
    pred <- crosspred(cb, mfirst,by = 0.1,cumul = TRUE);
    mmt <- pred$predvar[which.min(pred$allRRfit)];
    
    crall <- crossreduce(cb,mfirst,cen=19,by = 0.1)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    
  }, error=function(e){cat("Error",conditionMessage(e), "\n")})    
}
###########STEP 2#############
#META
method <- "reml"
#OVERALL
mvall <- mvmeta(yall~1,Sall,method=method) 
summary(mvall)
################################################################################
# OBTAIN BLUPS

blup <- blup(mvall,vcov=T)

# RE-CENTERING VALUE AS MINIMUM MORBIDITY TEMPERATURE IN EACH CITY
# GENERATE THE MATRIX FOR STORING THE RESULTS
mintempcity <- rep(NA,length(regions))
names(mintempcity) <- regions
q25city <- rep(NA,length(regions))
  names(q25city) <- regions
 q975city <- rep(NA,length(regions))
names(q75city) <- regions
for (i in regions) {
  sub <- data[[i]]
  tperc <- quantile(sub$l0avg_temp,0.5:0.995/100,na.rm=T)
  
  
  argvar <- list(x = tperc,fun="ns",df=3) 
  
  bvar <- do.call(onebasis,argvar)
  mintempcity[i] <-tperc[which.min((bvar%*%blup[[i]]$blup))]
  q25city[i] <- quantile(sub$l0avg_temp,0.025,na.rm = T)
  q975city[i] <- quantile(sub$l0avg_temp,0.975,na.rm = T)
}

##########STEP 3###############################################
#COMPUTE THE AN IN EACH CITY WITH EMPIRICAL CI
##############################################################
source("E:\\OneDrive\\11_Doctor research\\17.PM2.5组分健康效应\\住院\\温度code\\attrdl.R")
# CREATE THE MATRIX TO STORE THE ATTRIBUTABLE DEATHS
totan <- rep(NA,length(regions))
names(totan) <- regions
# NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
nsim <- 1000

# CREATE THE ARRAY TO STORE THE CI OF ATTRIBUTABLE NUMBER OF HOSPITAL ADMISSIONS
matsim<-array(NA,dim=c(length(regions),5,1),dimnames = list(regions,c("nationwide","cold","heat","excold","exheat")))
af.arraysim<-array(NA,dim=c(length(regions),5,nsim),dimnames = list(regions,c("nationwide","cold","heat","excold","exheat")))
an.arraysim<-array(NA,dim=c(length(regions),5,nsim),dimnames = list(regions,c("nationwide","cold","heat","excold","exheat")))
city.sim<-array(NA,dim=c(length(regions),5,8),dimnames = list(regions,
                                                              c("nationwide","cold","heat","excold","exheat"),
                                                              c("an.sim.mean","an.sim.se",'an.low','an.high','af.sim.mean','af.sim.se',
                                                                'af.low','af.high')))
#RUN THE LOOP FOR EACH CITY
for(i in regions){
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  sub <- data[[i]]
  sub <- sub [order(sub$date),] 
  #DEFINE PARAMETER FOR DLNM
  lag <- 28
  
  arglag <- list(knots=logknots(lag,3)) 
  argvar <- list(fun="ns",df=3) 
  # DERIVE THE CROSS-BASIS
  cb <- crossbasis(sub$l0avg_temp, lag = lag, argvar=argvar, arglag=arglag)
   cen <- mintempcity[i]
  #DEFINE CENTERING VALUE
  cen = mintempcity[i]
  matsim[i,"nationwide",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                                    vcov=blup[[i]]$vcov,type="an",dir="forw",cen=cen)
  
  matsim[i,"cold",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                              vcov=blup[[i]]$vcov,type="an",dir="forw",cen=cen,
                              range=c(-100,cen))
  matsim[i,"heat",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                              vcov=blup[[i]]$vcov,type="an",dir="forw",cen=cen,
                              range=c(cen,100))
  
  matsim[i,"excold",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                                vcov=blup[[i]]$vcov,type="an",dir="forw",cen=cen,
                                range=c(-100,q25city[i]))
  
  matsim[i,"exheat",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                                vcov=blup[[i]]$vcov,type="an",dir="forw",cen= mintempcity[i],
                                range=c(q975city[i],100))
  
  #DERIVE CONFIDENCE INTERVALS
  an.arraysim[i,"nationwide",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                                         vcov=blup[[i]]$vcov,type="an",dir="forw",cen=cen,sim=T,nsim=nsim)
  
  an.arraysim[i,"cold",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                                   vcov=blup[[i]]$vcov,type="an",dir="forw",cen=cen,
                                   range=c(-100,cen),sim=T,nsim=nsim)
  
  an.arraysim[i,"heat",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                                   vcov=blup[[i]]$vcov,type="an",dir="forw",cen=cen,
                                   range=c(cen,100),sim=T,nsim=nsim)
  
  an.arraysim[i,"excold",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                                     vcov=blup[[i]]$vcov,type="an",dir="forw",cen=cen,
                                     range=c(-100,q25city[i]),sim=T,nsim=nsim)
  
  
  an.arraysim[i,"exheat",] <- attrdl(sub$l0avg_temp,cb,sub$PSN_NO,coef=blup[[i]]$blup,
                                     vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
                                     range=c(q975city[i],100),sim=T,nsim=nsim)
  
  # STORE THE DENOMINATOR OF ATTRIBUTABLE DEATHS, I.E. TOTAL OBSERVED MORTALITY
  # CORRECT DENOMINATOR TO COMPUTE THE ATTRIBUTABLE FRACTION LATER, AS IN attrdl
  totan[i] <- sum(sub$PSN_NO,na.rm=T)
}


# TOTAL
# NB: FIRST SUM THROUGH CITIES
antot <- colSums(matsim,na.rm = T)
antotlow <- apply(apply(an.arraysim,c(2,3),sum,na.rm = T),1,quantile,0.025)

antothigh <- apply(apply(an.arraysim,c(2,3),sum,na.rm = T),1,quantile,0.975)
################################################################################
# TOTAL ATRRIBUTBALE NUMBER
totantot <- sum(totan)
################################################################################
# TOTAL ATTRIBUTABLE FRACTIONS

aftot <- antot/totantot*100
aftotlow <- antotlow/totantot*100 
aftothigh <- antothigh/totantot*100
AF_all <- data.frame(aftot = aftot,aftotlow= aftotlow,aftothigh=aftothigh,antot = antot,antotlow = antotlow,antothigh = antothigh, tot = totantot)
