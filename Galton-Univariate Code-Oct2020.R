############SIMULATIONS#####################################################
rm(list=ls())

#########################################################
#### EXACT MLE###########################################
###LOG L FUNCTION###
Logll <- function(TL,freq,theta){
  m<- length(TL)
  
  if( (pnorm(theta[2]*TL[2]-theta[1])) < 1e-16  ){
    a <- -1e+6}else{
      a <- freq[1]*log(pnorm(theta[2]*TL[2]- theta[1]))}
  #print(a)
  if( (1-pnorm(theta[2]*TL[m]-theta[1])) < 1e-16  ) {
    b <- -1e+6 }else{
      b <- freq[m]*log(1-pnorm(theta[2]*TL[m]-theta[1]))}
  #print(b)
  c<-0
  for(i in 2:(m-1)){
    #print(i)
    #print(freq[i])
    if ( (pnorm(theta[2]*TL[i+1]-theta[1]) - pnorm(theta[2]*TL[i]-theta[1])) < 1e-16 ){
      c <- c -1e+6 }else{
        c <- c + freq[i]*(log( pnorm(theta[2]*TL[i+1]-theta[1]) - pnorm(theta[2]*TL[i]-theta[1] ))) }
    #print(c)
    #print(TL[i])
  }
  L <- -(a+b+c)
  return(L)
}


#######################################################################
##########################################################
##########################################################
### EM ALGORITHM########
#########################################
### Calculation of Updates mu###
###E-Step For Mean Estimate###
Mest<- function(theta,bl,bu,Freq){
  Aj<- rep(0,length(bl))
  astar<- rep(0,length(bl))
  bstar<- rep(0,length(bl))
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2]
    astar[i]<- (bl[i]-theta[1])/theta[2]
    
  }
  
  for(i in 1:length(bl)){
    
    Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    
  }
  M<- sum(Aj*Freq)/sum(Freq) 
  #return(Aj)
  return(M)
  
}

############################################
#############################################
### Calculation of Updates Sigma###
###E-Step For Variance Estimate ###
SSest<- function(theta,bl,bu,muupdate,Freq){
  Bj<- rep(0,length(bl))
  bstar<- rep(0,length(bl))
  astar<- rep(0,length(bl))
  
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2]
    astar[i]<- (bl[i]-theta[1])/theta[2]
    
  }
  
  astar[1]<- -1000
  bstar[length(bl)]<- 1000
  
  
  for(i in 1:length(bl)){
    Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    +(muupdate-theta[1])^2+
      (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i]))))
    
  }
  SS<- sum(Bj*Freq)/sum(Freq)
  #S<- sqrt(SS)
  return(SS)
}

#####################################################
##########################################################
####MAXIMAZATION STEP:######
EM<- function(bl,bu,Freq,theta_init,maxit=1000,tol1=1e-3,tol2=1e-4){
  flag<- 0
  Mu_cur<- theta_init[1]
  S_cur<- theta_init[2]
  
  for (i in 1:maxit){
    #print(paste("Iteration number=", i))
    cur<- c(Mu_cur,S_cur)
    
    Munew<- Mest(theta=cur,bl,bu,Freq)
    
    SSnew<- SSest(theta=cur,bl,bu,
                  muupdate=Mest(theta=cur,bl,bu,Freq) ,Freq)
    
    
    Mu_new<- Munew
    S_new<- sqrt(SSnew)
    new_step<- c(Mu_new,S_new)
    
    if(abs(cur[1]-new_step[1])<tol1 & abs(cur[2]-new_step[2])<tol2){flag<-1 ;break}
    
    
    Mu_cur<- Mu_new
    S_cur<- S_new
  }
  if(!flag) warning("Didn't Converge \n")
  #list(paste("updated Mu of the EM=" , Mu_cur),paste("Updated Sigma of the EM=" ,S_cur),
  #    paste("iteration number=",i))
  updateres<- c(Mu_cur,(S_cur)^2)
  return(updateres)
}

###########################################################
###########################################################

######################################################################
###Galton data###
Galton<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton.csv",header=TRUE)
X<- Galton$Parent 
MP<- mean(X)
SSP<- ((1/length(X))*sum((X-mean(X))^2))


Galton<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton.csv",header=TRUE)
Y<- Galton$Children
MeC<- mean(Y)
SSC<- ((1/length(Y))*sum((Y-mean(Y))^2))

mle1<- c(MP,SSP,MeC,SSC)
###For Galton Data:Parent###
Galton2 <- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton2.csv",header=TRUE)
Galton2[1,1]<- -Inf
Galton2[11,2]<- Inf
#Galton2

Galton3 <- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton3.csv",header=TRUE)

Galton3[1,1]<- -Inf
Galton3[14,2]<- Inf
#Galton3

thetaintP<- c(MP,SSP)
blP<- Galton2$TL
buP<- Galton2$TU
FrP<- Galton2$Freq


thetaintC<- c(MeC,SSC)
blC<- Galton3$TL
buC<- Galton3$TU
FrC<- Galton3$Freq
###################################################
outParent<- EM(bl=blP,bu=buP,Freq=FrP,theta_init=thetaintP,maxit=1000,tol1=1e-3,tol2=1e-4)
outParent

outChildren<- EM(bl=blC,bu=buC,Freq=FrC,theta_init=thetaintC,maxit=1000,tol1=1e-3,tol2=1e-4)
outChildren

Emres<- c(outParent,outChildren) 
###################################################################################################

res <- optim(c((68/2),(1/2)),fn=Logll,TL=blP,freq=FrP,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))
estimate<- res$par
muP1<- round(res$par[1]/res$par[2],4)
SP1<- round(1/res$par[2],4)

#print(paste("Mle estimate of mu for exact likelihood =" , round(muP1,4)))
#print(paste("Mle estimate of Sigma for exact likelihood=" , round(SP1,4)))


res1 <- optim(c((68/2),(1/2)),fn=Logll,TL=blC,freq=FrC,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))
estimate1<- res1$par
muC1<- res1$par[1]/res1$par[2]
SC1<- 1/res1$par[2]

Mleexc<- c(muP1,(SP1)^2,muC1,(SC1)^2)
#####################################################################################################

ls()
rm("TL","freq","i","j")
rm("Fr","Mest","Fr","i","SSest","EM","Logll")
####################################################################################################
###################################################################################
library(truncnorm)
###1. Update Both Mu & Sigma based on the Formulas###

###Calculation of updates for Mu###
##E-Step for Mu estimation: Simulating Z's###
ZMCEM<- function(theta,data){
  k<- 1000
  sim<- matrix(rep(0,k*nrow(data)),ncol=nrow(data))
  #sim<- matrix(rep(0,k*3),ncol=3)
  
  for(i in 1 :nrow(data))
    sim[,i]<- rtruncnorm(k,a=data[i,1],b=data[i,2],mean=theta[1],sd=theta[2])
  
  return(sim)
}  
### E-Step for Mu & Sigma###
MuMCEM<- function(data,simZ){
  n<- sum(data[,3])
  Z<- colMeans(simZ)
  numerator<- rep(0,nrow(data))
  for (i in 1:nrow(data))
    numerator[i]<- data[i,3]*Z[i]
  sum(numerator)
  MuN<- (1/n)*sum(numerator)
  return(MuN)
}

################################################
sigmaMCEM<- function(data,simZZ,mupd){
  n<- sum(data[,3])
  ZZ<- simZZ
  #print(ZZ)
  NewZ<- (ZZ-mupd)^2
  #print(NewZ)
  SZNEW<- colMeans(NewZ) 
  #print(SZNEW)
  numerator<- rep(0,nrow(data))
  for (i in 1:nrow(data))
    numerator[i]<- data[i,3]*SZNEW[i]
  
  sigmaNN<- (1/n)*sum(numerator)
  sig<- sqrt(sigmaNN)
  
  # Z<- colMeans(ZZ)
  #print(Z)
  
  
  #zsqr<- (Z-mupd)^2
  
  #numerator<- rep(0,nrow(data))
  #for (i in 1:nrow(data))
  # numerator[i]<- data[i,3]*zsqr[i]
  
  #sigmaNN<- (1/n)*sum(numerator)
  #sig<- sqrt(sigmaNN)
  return(sig)
}

##########################################
###MONTE CARLO EM###
##########################################################
MCEM<- function(data,theta_init,maxit=1000,tol1=1e-2,tol2=1e-3){
  flag<- 0
  Mu_cur<- theta_init[1]
  S_cur<- theta_init[2]
  #theta_cur<- c(Mu_cur,S_cur)
  iter<- rep(0,maxit)
  Svec<- rep(0,maxit)
  Mvec<- rep(0,maxit)
  for (i in 1:maxit){
    cur<- c(Mu_cur,S_cur)
    Munew<- MuMCEM(data=mydat,simZ=ZMCEM(theta=cur,data=mydat))
    
    Snew<- sigmaMCEM(data=mydat,simZZ=ZMCEM(theta=cur,data=mydat),
                     mupd=MuMCEM(data=mydat,simZ=ZMCEM(theta=cur,data=mydat)))
    
    
    
    Mu_new<- Munew
    S_new<- Snew
    
    new_step<- c(Mu_new,S_new)
    
    
    if(abs(cur[1]-new_step[1])<tol1 & abs(cur[2]-new_step[2])<tol2){flag<-1 ;break}
    
    
    Mu_cur<- Mu_new
    
    S_cur<- S_new
    
    iter[i]<- i
    Svec[i]<- S_new
    Mvec[i]<- Mu_new 
    
  }
  
  if(!flag) warning("Didn't Converge \n")
  
  list(paste("updated Mu of the EM=" , Mu_cur),paste("Updated Sigma of the EM=" ,S_cur),
       paste("iteration number=",i))
  update=c(Mu_cur,(S_cur)^2)
  return(update)
}

##################################################################

##GALTON DATA: Parent##
###For Galton Data###
Galton2 <- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton2.csv",header=TRUE)
#Galton2
BL <- Galton2$TL
BL[1] <- -Inf
BU <- Galton2$TU
BU[11] <- Inf

Freq <- Galton2$Freq
GaltonDat<- cbind(BL,BU,Freq)
GD<- as.data.frame(GaltonDat)

mydat<- GD
#theta_init<- c(Me,SD)
MCEMParent<- MCEM(data=mydat,theta_init=thetaintP,maxit = 1000,tol1=1e-3,tol2=1e-4)
MCEMParent

#############################################################
### GALTON DATA: Children###
Galton3 <- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton3.csv",header=TRUE)
#Galton3
#dim(Galton3)
BL1 <- Galton3$TL
BL1[1] <- -Inf
BU1 <- Galton3$TU
BU1[14] <- Inf

Freqc <- Galton3$Freq
GaltonDatC<- cbind(BL1,BU1,Freqc)
GDC<- as.data.frame(GaltonDatC)

mydat<- GDC
MCEMChildren<- MCEM(data=mydat,theta_init=thetaintC,maxit = 1000,tol1=1e-3,tol2=1e-4)
MCEMChildren
######################################################################################
mcemres<- c(MCEMParent,MCEMChildren)
############################################################################################
GaltonRes<- matrix(rep(0,4*4),ncol=4)
colnames(GaltonRes)<- c("MuParent","VParent","MuChildren","VChildren")
row.names(GaltonRes)<- c("MLEungrp","MLEExact","EM","MCEM")
GaltonRes
GaltonRes[1,]<- mle1
GaltonRes[2,]<- Mleexc
GaltonRes[3,]<- Emres
GaltonRes[4,]<- mcemres
GaltonRes
library(xtable)
xtable(GaltonRes,digits = 5,"Parameters estimation for Univariate Galton Data")


