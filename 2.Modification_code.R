###########################################################################################################################################################
#Codes for "Non-optimum ambient temperature and cause-specific hospital admissions for cardiovascular diseases in the context of particulate air pollution"
#Authors for codes: Huimeng Liu, Jian Lei, Yunxing Jiang, Lijun Bai et.al.
#Correspondence to Shaowei Wu, Xi Chen.
##########################################################################################################################################################

######################################################################################################################################
#MODIFICATION BY PM2.5 AND BC ON THE ASSOCIATIONS BETWEEN NON-OPTIMUM TEMPERATURES AND HOSPITAL ADMISSIONS 
##########################################################################################################################################

library(tidyverse);library(dlnm);library(splines);
library(mgcv);library(mixmeta);
###############STEP 1#########################
#IMPORT DATA
data_use <- read.csv("test_data_fin.csv") %>% 
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


#LOOP FOR PM2.5 AND TIS MAJOR CONSTITUENTS AT LAG01
pollag <- c("PM2.5lag01","BClag01")

for(pol_i in pollag) {
  
  #DEFINE MATRIX TO STORE THE RESULTS

  coef<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  coefm1<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  coefm2<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  coefm3<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  coefm4<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  
  vcovl <- vector("list",length(regions))
  names(vcovl) <- regions
  vcovl1 <- vector("list",length(regions))
  names(vcovl1) <- regions
  vcovl2 <- vector("list",length(regions))
  names(vcovl2) <- regions
  vcovl3 <- vector("list",length(regions))
  names(vcovl3) <- regions
  vcovl4 <- vector("list",length(regions))
  names(vcovl4) <- regions
  
  
  MMT1 <- vector()
  MMT2 <- vector()
  MMT3 <- vector()
  MMT4 <- vector()
  #LOOP FOR CITIES
  for(i in regions) {
    # PRINT
    cat(i,"")
    
    
    #LOAD THE DATA
    sub <- data_use %>% 
      filter(district == i) 
    
    #GENERATRE COVARIATES
    sub$year <- year(sub$date)
    sub$date1 <- 1:length(sub$date)
    sub$rhu <- sub$l0rhu
    sub$dow    <- as.factor(weekdays(sub$date))
   
    
    #DEFINE QUARTILE GROUPS OF PM2.5 AND TIS MAJOR CONSTITUENTS
    qk<-quantile(sub[,pol_i],na.rm=T) 
    sub$qgroup[sub[,pol_i]<=qk[2]]<-"q1"
    sub$qgroup[sub[,pol_i]>qk[2] & sub[,pol_i]<=qk[3]]<-"q2"
    sub$qgroup[sub[,pol_i]>qk[3] & sub[,pol_i]<=qk[4]]<-"q3"
    sub$qgroup[sub[,pol_i]>qk[4] & sub[,pol_i]<=qk[5]]<-"q4"
    
    
    #DEFINE THE CROSS-BASIS
    lag <- 28
    arglag <- list(knots=logknots(lag,3))
    argvar <- list(fun="ns",df=3)
    cb <- crossbasis(sub$l0avg_temp, lag = lag, argvar=argvar, arglag=arglag)
    sub["pol_adj"] <- sub[,pol_i]
    cen <- quantile(sub$l0avg_temp,probs = 0.75,na.rm = T)
    #RUN THE MODEL
    tryCatch({
      
      fit <- gam(PSN_NO~cb:qgroup+ pol_adj+ns(rhu, df=3) + ns(date1,df=8*length(unique(year)))
                 +as.factor(Holiday)
                 +as.factor(dow),
                 control=list(maxit=100),
                 family=quasipoisson,sub,na.action="na.exclude")
      
   
      
      #EXTRACT COEFFICIENT OF EACH GROUP
      ## group 1
      indexcoef1<-grep("qgroupq1",names(coef(fit)))
      coef1<-coef(fit)[indexcoef1]
      
      indexvcov1<-grep("qgroupq1",rownames(vcov(fit)))
      vcov1<-vcov(fit)[indexvcov1,indexvcov1]
      
      pred1<-crossreduce(cb, coef=coef1, vcov=vcov1, type="overall",model.link="log",cen= cen)
      coefm1[i,]<-coef(pred1)
      vcovl1[[i]]<-vcov(pred1)
      mmt1 <- pred1$predvar[which.min(pred1$RRfit)]
    
      ## group 2
      indexcoef2<-grep("qgroupq2",names(coef(fit)))
      coef2<-coef(fit)[indexcoef2]
      
      indexvcov2<-grep("qgroupq2",rownames(vcov(fit)))
      vcov2<-vcov(fit)[indexvcov2,indexvcov2]
      
      pred2<-crossreduce(cb, coef=coef2, vcov=vcov2, type="overall",model.link="log",cen=cen)
      coefm2[i,]<-coef(pred2)
      vcovl2[[i]]<-vcov(pred2)
      mmt2 <- pred2$predvar[which.min(pred2$RRfit)];
      
      
      ## group 3
      indexcoef3<-grep("qgroupq3",names(coef(fit)))
      coef3<-coef(fit)[indexcoef3]
      
      indexvcov3<-grep("qgroupq3",rownames(vcov(fit)))
      vcov3<-vcov(fit)[indexvcov3,indexvcov3]
      
      pred3<-crossreduce(cb, coef=coef3, vcov=vcov3, type="overall",model.link="log",cen=cen)
      coefm3[i,]<-coef(pred3)
      vcovl3[[i]]<-vcov(pred3)
      mmt3 <- pred3$predvar[which.min(pred3$RRfit)];
      
      
      ## group 4
      indexcoef4<-grep("qgroupq4",names(coef(fit)))
      coef4<-coef(fit)[indexcoef4]
      
      indexvcov4<-grep("qgroupq4",rownames(vcov(fit)))
      vcov4<-vcov(fit)[indexvcov4,indexvcov4]
      
      pred4<-crossreduce(cb, coef=coef4, vcov=vcov4, type="overall",model.link="log",cen=cen)
      coefm4[i,]<-coef(pred4)
      vcovl4[[i]]<-vcov(pred4)
      mmt4 <- pred4$predvar[which.min(pred4$RRfit)];
      
      MMT1[i] <- mmt1
      MMT2[i] <- mmt2
      MMT3[i] <- mmt3
      MMT4[i] <- mmt4
     
    }, error=function(e){cat("Error",conditionMessage(e), "\n")})    
  }
 
  
  ### META FOR NATIONAL RESULT
  
  meta_nest1<-mixmeta(coefm1~1,vcovl1,method = "reml")
  meta_nest2<-mixmeta(coefm2~1,vcovl2,method = "reml")
  meta_nest3<-mixmeta(coefm3~1,vcovl3,method = "reml")
  meta_nest4<-mixmeta(coefm4~1,vcovl4,method = "reml")

  # PREDICT THE POOLED COEFFICIENTS
  avgtmean <- sapply(data,function(x) mean(x$l0avg_temp,na.rm=T))
  rangetmean <- sapply(data,function(x) diff(range(x$l0avg_temp,na.rm=T)))
  
  datanew <- data.frame(avgtmean=mean(tapply(avgtmean,regions,mean)),
                        rangetmean=mean(tapply(rangetmean,regions,mean))) 
  mvpred1 <- predict(meta_nest1,datanew,vcov=T,format="list")
  mvpred2 <- predict(meta_nest2,datanew,vcov=T,format="list")
  mvpred3 <- predict(meta_nest3,datanew,vcov=T,format="list")
  mvpred4 <- predict(meta_nest4,datanew,vcov=T,format="list") 
 
  #RE-CENTERING
  ranges <- t(sapply(data, function(x) 
    range(x$l0avg_temp,na.rm=T)))
  bound <- round(colMeans(ranges,1),1)
  
  tperc <- seq(bound[1],bound[2],length=50)
  argvar1 <-list(fun="ns",df=3)
  cb.tm1 <- crossbasis(tperc,maxlag=lag,argvar=argvar1,arglag=arglag)
  
  # DEFINE INDICATOR FOR CENTERING PERCENTILE FROM AVERAGE ASSOCIATIONF
  bvar <- do.call("onebasis",c(list(x=tperc),attr(cb.tm1,"argvar")))
  xlag <- 0:lag
  
  #PREDICTION FOR EACH POLLUTANTS' GROUP
  cp1 <- crosspred(bvar,coef=mvpred1$fit,vcov=mvpred1$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt1)
  cp2 <- crosspred(bvar,coef=mvpred2$fit,vcov=mvpred2$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt2)
  cp3 <- crosspred(bvar,coef=mvpred3$fit,vcov=mvpred3$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt3)
  cp4 <- crosspred(bvar,coef=mvpred4$fit,vcov=mvpred4$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt4)
  
  #MINIMUM MORBIDITY TEMPERATURE
  mmt1 <- cp1$predvar[which.min(cp1$allRRfit)]
  mmt2 <- cp2$predvar[which.min(cp2$allRRfit)]
  mmt3 <- cp3$predvar[which.min(cp3$allRRfit)]
  mmt4 <- cp4$predvar[which.min(cp4$allRRfit)]
  #REPREDICT AT MINIMUM MORBIDITY TEMPERATURE
  cp1 <- crosspred(bvar,coef=mvpred1$fit,vcov=mvpred1$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt1)
  cp2 <- crosspred(bvar,coef=mvpred2$fit,vcov=mvpred2$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt2)
  cp3 <- crosspred(bvar,coef=mvpred3$fit,vcov=mvpred3$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt3)
  cp4 <- crosspred(bvar,coef=mvpred4$fit,vcov=mvpred4$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt4)
}

#PLOT
plot(cp1,ylab="Relative risk",cex.axis = 1.3, cex.lab = 1.3,xlab=" ",
     lwd=2,col=(alpha("#362FBB",0.8)),ci.arg=list(density = 20,col=alpha("#E6C6FF",0.8)))
lines(cp2,ci="area",lwd=2,col=alpha("#50BCB9",0.8),ci.arg=list(density = 20,col=alpha("#BEFCFF",0.8)))

lines(cp3,ci="area",lwd=2,col=alpha("#FFB845",0.8),ci.arg=list(densit = 20,col=alpha("#F4D44E",0.8)))
lines(cp4,ci="area",lwd=2,col=alpha("#F97698",0.8),ci.arg=list(density = 20,col=alpha("#F4B998",0.8)))

