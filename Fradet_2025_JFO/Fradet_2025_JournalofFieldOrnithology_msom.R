#### Journal of Field Ornithology, Fradet_2025
#### Multi-species Occupancy Modeling
#### Danielle Fradet, MWFB, dtfradet@ucdavis.edu
#### Brett Furnas, CDFW, brett.furnas@wildlife.ca.gov


#### Set up (run functions at bottom first!!!)
library(jagsUI)
library(chron)
library(rjags)




#### Get data
data<-read.csv("Fradet_2025_JFO_data.csv",as.is=T)




#### Format the data
history<-get.history(data)
y<-history[["y"]] # detection history
pas<-history[["pas"]] # passerines = 1
yday<-history[["yday"]] # survey date
range(yday,na.rm=T)
mn.yday<-round(mean(yday,na.rm=T))
sd.yday<-round(sd(yday,na.rm=T))
standards<-c(mn.yday,sd.yday) # standardize dates for modeling
yday<-(yday-mn.yday)/sd.yday
yday[is.na(yday)]<-0
covars.site<-history[["covars.site"]] # site level covariates




### Inspect Data
dim(y)
dimnames(y)
sum(pas)
naives<-get.naive(y) # Naives estimates of occupancy
naives[order(naives$naive,decreasing=T),]

index.naives<-which(naives$naive>0.25) # Common species
naives.25<-naives[which(naives$naive>0.25),]
nrow(naives.25)




#### Run model (set nc to match parallel processing capacity)
msom.1<-run.msom(y,pas,yday,covars.site,standards,ni=10000,nb=5000,nc=3,nt=3,parallel=T)



#### Evaluate model results (set to 90% credible intervals)
results.main<-get.model.estimates(msom.1,"mu.A.rice","sd.B.qjd",alpha=0.1)


################################################
################################################
################################################
###### FUNCTIONS (Run these first!!!)

##### Generate detection history from species interpretation results 
get.history<-function(data){
  
  sites<-sort(unique(data$Concat))
  n.sites<-length(sites)
  surveys<-array(c(rep(1:3,each=3),rep(1:3,3)),dim=c(9,2))
  colnames(surveys)<-c("day","period")
  n.surveys<-nrow(surveys)
  species<-sort(unique(data$Event))
  species<-species[!species%in%c("")]
  n.species<-length(species)
  
  pas<-rep(0,length(species))
  
  
  pas[which(species%in%c("AMCR", "AMGO", "AMRO", "ATFL", "BEVI","BEWR","BGGN", "BHCO", "BHGR", "BLGR", "BLPH", "BRBL","BRCR","BTYW","BUOR", 
                         "BUSH", "CALT","CASJ", "CEDW", "CHSP","CLSW","CORA","COYE","DEJU","EUST","GRKI","GTGR", "HOFI","HOSP", "HOWR", "LASP",
                         "LAZB","LEGO", "MAWR","NAWA","NOMO", "NOPA", "NRWS", "OATI","PUFI","RWBL","SAVS","SOSP","SPTO", "SWTH", "TRES", "WAVI",
                         "WBNU", "WEBL", "WEKI", "WEME","WETA", "WEWP", "YBCU","YEWA"))]<-1 
  
  rice<-c("G1","G2","H2","I3","J2","J3","K2","K3", "K5","L1","L4","M2","N1","N2","N3", "O1","O3","O4","P2","P3", "P4","Q1","Q3", "R1","R2", "R4","S1","S2","T1", "T2", "U1","U2")
  marsh<-c("A6","B6","D5", "F5","G4", "H3", "H4","J4","J5", "M5","M7","M8","N200P1","N200P2","N280P3","N280P5", "N5","N8","O6", "P5","POPSP1","POPSP5","Q6", "R6", "S3","S4","S5","S6","T6","U5", "U6")
  riparian<-c("A5","A7","CP1","D7","E6","F6","G5","G6","H5","H6", "I5", "I6","J6","K6", "L9","M9","MP1","MP5","N60P3","N60P4","N280P1", "O7","O72", "O8","O82", "P7","R7","Q7","S7","T71","T72", "U71","U72")
  
  y<-array(NA,dim=c(n.sites,n.surveys,n.species))
  dimnames(y)<-list(sites,1:n.surveys,species)
  covars.site<-array(NA,dim=c(n.sites,4))
  dimnames(covars.site)<-list(sites,c("year","rice","marsh","riparian"))
  
  yday<-array(NA,dim=c(n.sites,n.surveys))
  dimnames(yday)<-list(sites,1:n.surveys)
  
  for(j in 1:n.sites){
    covars.site[j,"year"]<-unique(subset(data,Concat==sites[j])$Year)
    covars.site[j,"rice"]<-ifelse(substr(sites[j],5,7)%in%rice,1,0)
    covars.site[j,"marsh"]<-ifelse(substr(sites[j],5,7)%in%marsh,1,0)
    covars.site[j,"riparian"]<-ifelse(substr(sites[j],5,7)%in%riparian,1,0)
    
    dates<-sort(unique(subset(data,Concat==sites[j])$Julian))
    start<-dates[1]
    for(k in 1:n.surveys){
      
      temp.data<-subset(data,Concat==sites[j] & 
                          Julian==as.numeric(start+surveys[k,"day"]-1) &
                          Period==as.numeric(surveys[k,"period"]))
      if(nrow(temp.data)>0){
        yday[j,k]<-as.numeric(start+surveys[k,"day"]-1)
        for(i in 1:n.species){
          temp.spp<-subset(temp.data,Event==species[i])
          y[j,k,i]<-ifelse(nrow(temp.spp)>0,1,0)
        } # i
      } # if
    } # k
  } # j
  
  out<-list(y,yday,pas,covars.site)
  names(out)<-c("y","yday","pas","covars.site")
  
  return(out)
}



#### Multi-species Occupancy Model, no data augmentation
run.msom<-function(y,pas,yday,covars.site,standards,ni,nb,nc,nt,parallel){
  
  naives<-get.naive(y) # for restricting species-level reporting
  index.naives<-which(naives$naive>0.25)
  
  
  rice<-as.numeric(covars.site[,"rice"]) # site-level covariates
  marsh<-as.numeric(covars.site[,"marsh"])
  riparian<-as.numeric(covars.site[,"riparian"])
  w.rice<-0.412
  w.marsh<-0.515
  w.riparian<-0.073
  
  
  p1<-c(1,0,0,1,0,0,1,0,0) # survey times
  p2<-c(0,1,0,0,1,0,0,1,0)
  p3<-c(0,0,1,0,0,1,0,0,1)
  
  jd.fit<-140:200 # for phenology 
  n.jd.fit<-length(jd.fit)
  jd.fit.c<-(jd.fit-standards[1])/standards[2]
  
  ## Specify model in BUGS language
  sink("model.txt")
  cat("
      model {
      
      # Flat priors
      mu.A.rice ~ dunif(-20,20)
      tau.A.rice <- 1/(sd.A.rice*sd.A.rice)
      sd.A.rice ~ dunif(0,20)
      
      mu.A.marsh ~ dunif(-20,20)
      tau.A.marsh <- 1/(sd.A.marsh*sd.A.marsh)
      sd.A.marsh ~ dunif(0,20)
      
      mu.A.riparian ~ dunif(-20,20)
      tau.A.riparian <- 1/(sd.A.riparian*sd.A.riparian)
      sd.A.riparian ~ dunif(0,20)
      
      mu.B.1 ~ dunif(-20,20)
      tau.B.1 <- 1/(sd.B.1*sd.B.1)
      sd.B.1 ~ dunif(0,20)
      
      mu.B.2 ~ dunif(-20,20)
      tau.B.2 <- 1/(sd.B.2*sd.B.2)
      sd.B.2 ~ dunif(0,20)
      
      mu.B.3 ~ dunif(-20,20)
      tau.B.3 <- 1/(sd.B.3*sd.B.3)
      sd.B.3 ~ dunif(0,20)
      
      mu.B.jd ~ dunif(-20,20)
      tau.B.jd <- 1/(sd.B.jd*sd.B.jd)
      sd.B.jd ~ dunif(0,20)
      
      mu.B.qjd ~ dunif(-20,20)
      tau.B.qjd <- 1/(sd.B.qjd*sd.B.qjd)
      sd.B.qjd ~ dunif(0,20)
      
      # Likelihoods
      for (i in 1:n.species) { # hyper parameters by species
      A.rice[i] ~ dnorm(mu.A.rice, tau.A.rice)# site-level covariates
      A.marsh[i] ~ dnorm(mu.A.marsh, tau.A.marsh)
      A.riparian[i] ~ dnorm(mu.A.riparian, tau.A.riparian) 
      B.1[i] ~ dnorm(mu.B.1,tau.B.1) # 30 min before sunrise  (survey-level covariates)
      B.2[i] ~ dnorm(mu.B.2,tau.B.2) # sunrise
      B.3[i] ~ dnorm(mu.B.3,tau.B.3) # 30 min after sunrise 
      B.jd[i] ~ dnorm(mu.B.jd,tau.B.jd) # control for phenology
      B.qjd[i] ~ dnorm(mu.B.qjd,tau.B.qjd)
      
      for (j in 1:n.sites) { # Ecological model for true occurrence
      z[j,i] ~ dbern(psi[j,i])
      logit(psi[j,i]) <- lpsi.lim[j,i] 
      lpsi.lim[j,i] <- min(999, max(-999, lpsi[j,i])) # logit stabilization	
      lpsi[j,i] <- A.rice[i]*rice[j]+ A.marsh[i]*marsh[j]
      + A.riparian[i]*riparian[j]
      
      # Observation model for replicated detection/nondetection observations
      for (k in 1:n.surveys) {
      y[j,k,i] ~ dbern(p.eff[j,k,i])
      p.eff[j,k,i] <- z[j,i] * p[j,k,i]
      logit(p[j,k,i]) <- lp.lim[j,k,i] 
      lp.lim[j,k,i] <- min(999, max(-999, lp[j,k,i])) # logit stabilization
      lp[j,k,i] <- B.1[i]*p1[k] + B.2[i]*p2[k]+ B.3[i]*p3[k] 
      + B.jd[i]*jd[j,k] + B.qjd[i]*jd[j,k]*jd[j,k]
      } # k
      } # j
      } # i
      
      # Derived quantities
      for(i in 1:n.index.naives){ # Overall occupancies and detection probs for common species
      psi.spp[i] <- sum(psi[,index.naives[i]])/n.sites
      psi.spp.weight[i]<-w.rice*psi.rice[i]+w.marsh*psi.marsh[i]+w.riparian*psi.riparian[i]
      
      p1.spp[i] <- sum(p[,1,index.naives[i]]+p[,4,index.naives[i]]+p[,7,index.naives[i]])/(n.sites*3)
      p2.spp[i] <- sum(p[,2,index.naives[i]]+p[,5,index.naives[i]]+p[,8,index.naives[i]])/(n.sites*3)
      p3.spp[i] <- sum(p[,3,index.naives[i]]+p[,6,index.naives[i]]+p[,9,index.naives[i]])/(n.sites*3)
      A.rice.spp[i] <- A.rice[index.naives[i]]
      A.marsh.spp[i] <- A.marsh[index.naives[i]]
      A.riparian.spp[i] <- A.riparian[index.naives[i]]
      B.jd.spp[i] <- B.jd[index.naives[i]]
      B.qjd.spp[i] <- B.qjd[index.naives[i]]
      } # i
      
      for(i in 1:n.index.naives){ # occupancies by habitat type for common species
        psi.rice[i]<-sum(z.rice[,i])/sum(rice[])
        psi.marsh[i]<-sum(z.marsh[,i])/sum(marsh[])
        psi.riparian[i]<-sum(z.riparian[,i])/sum(riparian[])
        for (j in 1:n.sites) {
          z.rice[j,i] <- psi[j,index.naives[i]]*rice[j]
          z.marsh[j,i] <- psi[j,index.naives[i]]*marsh[j]
          z.riparian[j,i] <- psi[j,index.naives[i]]*riparian[j]
        } # j 
      } # i
      
      
      for(i in 1:n.species){ # occupancies by habitat type for ALL species 
        psi.rice.all[i]<-sum(z.rice.all[,i])/sum(rice[])
        psi.marsh.all[i]<-sum(z.marsh.all[,i])/sum(marsh[])
        psi.riparian.all[i]<-sum(z.riparian.all[,i])/sum(riparian[])
        for (j in 1:n.sites) {
          z.rice.all[j,i] <- psi[j,i]*rice[j]
          z.marsh.all[j,i] <- psi[j,i]*marsh[j]
          z.riparian.all[j,i] <- psi[j,i]*riparian[j]
        } # j 
      } # i
      
      
      for (j in 1:n.sites) { # richness by taxa
      for(i in 1:n.species){
      z.rice.pas[j,i] <- z[j,i]*rice[j]*pas[i]
      z.rice.notpas[j,i] <- z[j,i]*rice[j]*(1-pas[i])
      z.marsh.pas[j,i] <- z[j,i]*marsh[j]*pas[i]
      z.marsh.notpas[j,i] <- z[j,i]*marsh[j]*(1-pas[i])
      z.riparian.pas[j,i] <- z[j,i]*riparian[j]*pas[i]
      z.riparian.notpas[j,i] <- z[j,i]*riparian[j]*(1-pas[i])
      } # i3
      temp.rich.rice.pas[j] <- sum(z.rice.pas[j,])
      temp.rich.rice.notpas[j] <- sum(z.rice.notpas[j,])
      temp.rich.marsh.pas[j] <- sum(z.marsh.pas[j,])
      temp.rich.marsh.notpas[j] <- sum(z.marsh.notpas[j,])
      temp.rich.riparian.pas[j] <- sum(z.riparian.pas[j,])
      temp.rich.riparian.notpas[j] <- sum(z.riparian.notpas[j,])
      } # j2
      
      for(i in 1:n.index.naives){    #loop to calculate differences in average occupancy between habitat types for each indivdiual species
      delta.rice.marsh[i] <- psi.rice[i] -psi.marsh[i]
      delta.rice.riparian[i] <- psi.rice[i] -psi.riparian[i] 
      delta.riparian.marsh[i] <- psi.riparian[i] -psi.marsh[i]
      }

      
      
      rich.rice.pas <- sum(temp.rich.rice.pas[])/ sum(rice[])
      rich.rice.notpas <- sum(temp.rich.rice.notpas[])/ sum(rice[])
      
      rich.marsh.pas <- sum(temp.rich.marsh.pas[])/ sum(marsh[])
      rich.marsh.notpas <- sum(temp.rich.marsh.notpas[])/ sum(marsh[])
      
      rich.riparian.pas <- sum(temp.rich.riparian.pas[])/ sum(riparian[])
      rich.riparian.notpas <- sum(temp.rich.riparian.notpas[])/ sum(riparian[])
      
    ##tests for species richness passerines and nonpasserines
      
      test.rich.notpas.rip.marsh<-rich.riparian.notpas-rich.marsh.notpas
      test.rich.notpas.rip.rice<-rich.riparian.notpas-rich.rice.notpas
      test.rich.notpas.marsh.rice<-rich.marsh.notpas-rich.rice.notpas
      
      test.rich.pas.rip.marsh<-rich.riparian.pas-rich.marsh.pas
      test.rich.pas.rip.rice<-rich.riparian.pas-rich.rice.pas
      test.rich.pas.marsh.rice<-rich.marsh.pas-rich.rice.pas
      
  
      
    ##loop for simpson's for each habitat type
    for(i in 1:n.species){
      p.psi.spp.rice.all[i]<-psi.rice.all[i]/sum(psi.rice.all[])
      p.psi.spp.marsh.all[i]<-psi.marsh.all[i]/sum(psi.marsh.all[])
      p.psi.spp.riparian.all[i]<-psi.riparian.all[i]/sum(psi.riparian.all[])
      } # i
    D.rice<-sum(p.psi.spp.rice.all[]*p.psi.spp.rice.all[])
    D.marsh<-sum(p.psi.spp.marsh.all[]*p.psi.spp.marsh.all[])
    D.riparian<-sum(p.psi.spp.riparian.all[]*p.psi.spp.riparian.all[])
    index.simpson.rice<-(1/D.rice)/n.species
    index.simpson.marsh<-(1/D.marsh)/n.species
    index.simpson.riparian<-(1/D.riparian)/n.species
    even.simpson.rice<-index.simpson.rice*n.species 
    even.simpson.marsh<-index.simpson.marsh*n.species
    even.simpson.riparian<-index.simpson.riparian*n.species
    
    ##tests for even comparison
    test.even.marsh.rice<-even.simpson.marsh-even.simpson.rice
    test.even.rip.rice<-even.simpson.riparian-even.simpson.rice
    test.even.rip.marsh<-even.simpson.riparian-even.simpson.marsh
    
      
  ##loop for passerines (from the pass vector)
      for(i in 1:n.species){ # phenology
      optimum.date.spp[i]<-sum(max.date.spp[,i]) 
      temp.optimum.pas[i]<-optimum.date.spp[i]*pas[i]
      for(j in 1:n.jd.fit){ 
      max.date.spp[j,i]<-equals(p.date[j,i],max.p.date[i])*jd.fit[j]
      p.date[j,i] <- B.jd[i]*jd.fit.c[j] + B.qjd[i]*jd.fit.c[j]*jd.fit.c[j]
      } # j
      max.p.date[i]<-max(p.date[,i])
      } # i
      optimum.date.pas<-sum(temp.optimum.pas[])/sum(pas[])
      } # end model
      
      ",fill = TRUE)
  sink()
  
  
  ## Bundle data
  win.data <- list(y = y, 
                   pas=pas,
                   rice=rice,
                   marsh=marsh,
                   riparian=riparian,
                   jd=yday,
                   p1=p1,
                   p2=p2,
                   p3=p3,
                   w.rice=w.rice,
                   w.marsh=w.marsh,
                   w.riparian=w.riparian,
                   n.species = length(y[1,1,]),
                   n.sites = length(y[,1,1]), 
                   n.surveys = length(y[1,,1]),
                   index.naives=index.naives,
                   n.index.naives=length(index.naives),
                   jd.fit=jd.fit,
                   jd.fit.c=jd.fit.c,
                   n.jd.fit=length(jd.fit)
  )
  
  ## Initial values
  zst <- sapply(1:length(y[1,1,]),f<-function(i){apply(y[,,i], 1, max,na.rm=TRUE)})# Observed occurrence as starting values for z
  inits <- function() list(z = zst,
                           mu.A.rice= -2,sd.A.rice=0.5,
                           mu.A.marsh=0,sd.A.marsh=0.5,
                           mu.A.riparian=0,sd.A.riparian=0.5,
                           mu.B.1= -1, sd.B.1=0.5,
                           mu.B.2=0, sd.B.2=0.5,
                           mu.B.3=0, sd.B.3=0.5,
                           mu.B.jd=0, sd.B.jd=0.5,
                           mu.B.qjd=0, sd.B.qjd=0.5
  )
  
  ## Parameters monitored
  params <- c("mu.A.rice","sd.A.rice",
              "mu.A.marsh","sd.A.marsh",
              "mu.A.riparian","sd.A.riparian",
              "mu.B.1","sd.B.1",
              "mu.B.2","sd.B.2",
              "mu.B.3","sd.B.3",
              "mu.B.jd","sd.B.jd",
              "mu.B.qjd","sd.B.qjd",
              "psi.spp","p1.spp","p2.spp","p3.spp",
              "A.rice.spp","A.marsh.spp","A.riparian.spp",
              "B.jd.spp", "B.qjd.spp",
              "rich.rice.pas","rich.marsh.pas","rich.riparian.pas",
              "rich.rice.notpas","rich.marsh.notpas","rich.riparian.notpas",
              "psi.rice","psi.marsh","psi.riparian",
              "psi.rice.all","psi.marsh.all","psi.riparian.all",
              "optimum.date.spp","optimum.date.pas",
              "psi.spp.weight",
              "even.simpson.rice", "even.simpson.marsh","even.simpson.riparian",
              "index.simpson.rice", "index.simpson.marsh","index.simpson.riparian",
              "test.even.marsh.rice", "test.even.rip.rice", "test.even.rip.marsh",
              "test.rich.notpas.rip.marsh","test.rich.notpas.rip.rice","test.rich.notpas.marsh.rice",
              "test.rich.pas.rip.marsh","test.rich.pas.rip.rice","test.rich.pas.marsh.rice",
              "delta.rice.marsh", "delta.rice.riparian", "delta.riparian.marsh" 
              
              
  )
  
  ## Call JAGS from R
  out <- jags(win.data, inits, params, "model.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb, parallel=parallel)
  
  
  ## Summarize posteriors
  return(out)}




####
get.naive<-function(y,alpha){
  out<-data.frame(species=dimnames(y)[[3]],naive=NA)
  for(k in 1:nrow(out)){
    temp.1<-y[,,k]
    out$naive[k]<-mean(apply(temp.1,1,max,na.rm=T))
  }
  return(out)
}




#### Model parameters with CI
get.model.estimates<-function(model,index.1,index.2,alpha){
  index.1<-which(rownames(model$summary)==index.1)
  index.2<-which(rownames(model$summary)==index.2)
  names<-rownames(model$summary)[index.1:index.2]
  out<-data.frame(Parameter=names,
                  MEAN=as.numeric(model$summary[index.1:index.2,"mean"]),
                  SD=as.numeric(model$summary[index.1:index.2,"sd"]),
                  CV=NA,
                  CI.LO=NA,CI.UP=NA,
                  Rhat=as.numeric(model$summary[index.1:index.2,"Rhat"]),
                  p.val=NA)
  out$CV<-abs(out$SD/out$MEAN)
  L<-model$sims.list
  temp.1<-array(unlist(L), dim = c(length(unlist(L))/(dim(model$summary)[1]), 
                                   dim(model$summary)[1]))
  colnames(temp.1)<-rownames(model$summary)
  for(i in 1:nrow(out)){
    temp.2<-temp.1[,(index.1+i-1)]
    out$CI.LO[i]<-as.numeric(quantile(temp.2,alpha/2)) 
    out$CI.UP[i]<-as.numeric(quantile(temp.2,1-alpha/2))
    out$p.val[i]<-ifelse(out$CI.UP[i]<0 | out$CI.LO[i]>0, "***","")
  }
  return(out)
}


### END
