###########################################################################################################################################################
#Codes for "Non-optimum ambient temperature and cause-specific hospital admissions for cardiovascular diseases in the context of particulate air pollution"
#Authors for codes: Huimeng Liu, Jian Lei, Yunxing Jiang, Lijun Bai et.al.
#Correspondence to Shaowei Wu, Xi Chen.
##########################################################################################################################################################

######################################################################################################################################
#Overall exposure-response relationships between non-optimum temperature and hospital admissions (HAs) for major CVDs
##########################################################################################################################################
library(tidyverse)
library(lubridate)
library(dlnm) 
library(mvmeta) 
library(splines)
library(mgcv)

###############STEP 1#########################
#IMPORT DATA
data_use <- read.csv("E:/OneDrive/11_Doctor research/17.PM2.5组分健康效应/住院/温度code/Github code1126/test_data_fin.csv") %>% 
  mutate(date = as.Date(date))
regions <- unique(data_use$district)


#ARRANGE THE DATA AS A LIST OF DATA SETS
data <- lapply(regions,function(x) data_use[data_use$district==x,])
names(data) <- regions

# TEMPERATURE RANGES
ranges <- t(sapply(data, function(x) 
  range(x$l0avg_temp,na.rm=T)))
bound <- colMeans(ranges)
value_cold <- -3.1
value_hot<-  27.9

#SET EMPTY MATRIX TO STORE THE RESULTS
yall<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
yhot <- matrix(NA,length(regions),5,dimnames=list(regions,paste("b",seq(5),sep="")))
ycold <- matrix(NA,length(regions),5,dimnames=list(regions,paste("b",seq(5),sep="")))
  
Sall <- vector("list",length(regions))
names(Sall) <- regions

Shot <- c()
Scold <- c()


  
#LOOP FOR CITIES
for(i in regions) {
    # PRINT
    cat(i,"")
    
    # LOAD CITY-SPECIFIC DATA
    sub <- data[[i]]
    
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
    
    #RUN THE MODELS 
    tryCatch({
      mfirst <- gam(PSN_NO~cb+ ns(rhu, df=3) +ns(date1,df=8*length(unique(year)))
                    +as.factor(Holiday)
                    +as.factor(dow),
                    family=quasipoisson,sub,na.action="na.exclude")
      
      # PREDICTION AND REDUCTION TO OVERALL CUMULATIVE EXPOSURE-RESPONSE
      pred <- crosspred(cb, mfirst,by = 0.1,cumul = TRUE);
    
      
      crall <- crossreduce(cb,mfirst,cen=20,by = 0.1) #
      crhot <- crossreduce(cb,mfirst,type="var",value=value_hot,cen=20)
      crcold <- crossreduce(cb,mfirst,type="var",value=value_cold,cen=20)
      
      # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
      yall[i,] <- coef(crall)
      Sall[[i]] <- vcov(crall)
      yhot[i,] <- coef(crhot)
      Shot[[i]] <- vcov(crhot)
      ycold[i,] <- coef(crcold)
      Scold[[i]] <- vcov(crcold)
      
    
      
      
    }, error=function(e){cat("Error",conditionMessage(e), "\n")})    
  }
  
 ###########STEP 2#############
  
  #META
  method <- "reml"
  #OVERALL
  mvall <- mvmeta(yall~1,Sall,method=method) 
  summary(mvall)
  #HOT
  mvhot <- mvmeta(yhot~1,Shot,method=method)
  summary(mvhot)
  #COLD
  mvcold <- mvmeta(ycold~1,Scold,method=method)
  summary(mvcold)
  
  # OBTAIN BLUPS
  
  blup <- blup(mvall,vcov=T)
  
  # RE-CENTERING
  # DEFINE RELATED AVERAGE TEMPERATURES
 
  tperc <- seq(bound[1],bound[2],by = 0.1)
  
  cb.tm1 <- crossbasis(tperc,lag=lag,argvar=argvar,arglag=arglag)
  
  # DEFINE INDICATOR FOR CENTERING PERCENTILE FROM AVERAGE ASSOCIATION
  bvar <- do.call("onebasis",c(list(x=tperc),attr(cb.tm1,"argvar")))
  xlag <- 0:lag
  blag <- do.call("onebasis",c(list(x=xlag),attr(cb.tm1,"arglag")))
  
  #PREDICT THE CITY-SPECIFIC CUMULATIVE ASSOCIATIONS
  regall <- lapply(seq(nrow(yall)),function(i) 
    crosspred(bvar,coef=yall[i,],
              vcov=Sall[[i]],model.link="log",cen=23))
  reghot <- lapply(seq(nrow(yhot)),function(i) 
    crosspred(blag,coef=yhot[i,],
              vcov=Shot[[i]],model.link="log",cen=23))
  regcold <- lapply(seq(nrow(ycold)),function(i) 
    crosspred(blag,coef=ycold[i,],
              vcov=Scold[[i]],model.link="log",cen=23))
  
  
  #PREDICT THE OVERALL CUMULATIVE ASSOCIATIONS
  cpall.total <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
                           model.link="log",by=0.1,from=bound[1],to=bound[2],cen=mmt)
 
  #MINIUM MORBIDITY TEMPERATURE
  mmt <- cpall.total$predvar[which.min(cpall.total$allRRfit)]
  #REPREDICT AT MINIMUM MORBIDITY TEMPERATURE
  cpall.total <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
                           model.link="log",by=0.1,from=bound[1],to=bound[2],cen=mmt)
  #PREDICT THE LAG PATTERNS OF EXTREME HOT TEMPERATURES
  cphot.total <- crosspred(blag,coef=coef(mvhot),vcov=vcov(mvhot),
                           model.link="log",at=0:lag)
  #PREDICT THE LAG PATTERNS OF EXTREME COLD TEMPERATURES
  cpcold.total <- crosspred(blag,coef=coef(mvcold),vcov=vcov(mvcold),
                            model.link="log",at=0:lag)
  
 ##############STEP #EEEEEEEEEE
#PLOT OVERALL EXPOSURE-RESPONSE CURVES
xlab <- paste("Daily Mean Temperature")
plot(cpall.total,type="n",ylab="Relative Risk",
     ci.arg = list(col = c("#D0E7Ed")),cex.axis = 1.3, cex.lab = 1.3)

lines(cpall.total,col=1,lwd=1)
lines(c(mmt,mmt),c(0,10),lty = 2,lwd = 2,col = "red")
mtext("Example",cex=1.2,side = 3,line = 0.5)
mtext(paste0(mmt),at = c(mmt+4),line = -0.8,cex=1)

#PLOT LAG PATTERN FOR EXTREME HOT TEMPERATURE
plot(cphot.total,type="n",ylab="Relative Risk",cex.axis = 1.3, cex.lab = 1.5,
     ci.arg = list(col = c("#D0E7Ed")),ylim = c(0.95,1.04))
lines(cphot.total,col=1,lwd=1)
mtext("Lag pattern for hot",cex=1,line = 0)
#PLOT LAG PATTERN FOR EXTREME HOT TEMPERATURE

plot(cpcold.total,type="n",ylab="Relative Risk",cex.axis = 1.3, cex.lab = 1.5,
     ci.arg = list(col = c("#D0E7Ed")),ylim = c(0.85,1.10))
lines(cpcold.total,col=1,lwd=1)
mtext("Lag pattern for cold",cex=1,line = 0)


p2_5 <- -3.1
p975 <- 27.9


#OUT PUT ESTIMATION AT SPECIFIC TEMPERATURE
results <- round(cbind(cpall.total$allRRfit,cpall.total$allRRlow,cpall.total$allRRhigh)
                 [c(as.character(p2_5),as.character(p975),
                    as.character(mmt)),],digits=3)

