###Galton data###
Galton <- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Galton Data/Galton.csv",header=TRUE)  
X<- Galton$Parent
MP<- mean(X)
SSP<- ((1/length(X))*sum((X-mean(X))^2))


Galton <- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Galton Data/Galton.csv",header=TRUE)  
Y<- Galton$Children
MeC<- mean(Y)
SSC<- ((1/length(Y))*sum((Y-mean(Y))^2))

mle1<- c(MP,SSP,MeC,SSC)
###For Galton Data:Parent###
Galton2 <- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Galton Data/Galton2.csv",header=TRUE)

Galton2[1,1]<- -Inf
Galton2[11,2]<- Inf
#Galton2

Galton3 <- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Galton Data/Galton3.csv",header=TRUE)

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
##########################################################################################
#### EXACT MLE###########################################
###LOG L FUNCTION###
Logll <- function(TL,freq,theta){ 
  m<- length(TL) 
  
  if( (pnorm(theta[2]*TL[2]-theta[1])) < 1e-16  ){ 
    a <- -1e+6} 
  else{
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

############################################################################################################
res <- optim(c((68/2.5),(1/2.5)),fn=Logll,TL=blP,freq=FrP,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))
estimate<- res$par
muP1<- round(res$par[1]/res$par[2],4)
SP1<- round(1/res$par[2],4)

#print(paste("Mle estimate of mu for exact likelihood =" , round(muP1,4)))
#print(paste("Mle estimate of Sigma for exact likelihood=" , round(SP1,4)))


res1 <- optim(c((68/2.5),(1/2.5)),fn=Logll,TL=blC,freq=FrC,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))
estimate1<- res1$par
muC1<- res1$par[1]/res1$par[2]
SC1<- 1/res1$par[2]

Mleexc<- c(muP1,muC1,(SP1)^2,(SC1)^2)
Mleexc
############################################################################################################
### EM ALGORITHM For simulated data and Galton data###
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
  dinom<- NULL
  for(i in 1:length(bl)){ 
    #print(i)
    #print(pnorm(bstar[i]))
    #print(pnorm(astar[i]))
    dinom[i]<- (pnorm(bstar[i])-pnorm(astar[i]))
    #print(dinom[i])
    if(dinom[i]==0) { Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001)}
    else {Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i])}
  }
  #Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i])     
  #}
  M<- sum(Aj*Freq)/sum(Freq) 
  return(M) 
  
}


############################################)
#############################################
### Calculation of Updates Sigma###
###E-Step For Variance Estimate ###
SSest<- function(theta,bl,bu,muupdate,Freq){ 
  Bj<- rep(0,length(bl)) 
  bstar<- rep(0,length(bl)) 
  astar<- rep(0,length(bl)) 
  
  dinom<- NULL
  
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2] 
    astar[i]<- (bl[i]-theta[1])/theta[2] 
    
  }
  
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  
  
  for(i in 1:length(bl)){ 
    dinom[i]<- (pnorm(bstar[i])-pnorm(astar[i])) 
    if(dinom[i]==0) {Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/0.0001)
    +(muupdate-theta[1])^2+
      (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001))}
    
    else{Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/dinom[i])
    +(muupdate-theta[1])^2+
      (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i]))}
    
    
    # Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    #+(muupdate-theta[1])^2+
    # (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i]))))
    
  }
  SS<- sum(Bj*Freq)/sum(Freq) 
  #S<- sqrt(SS)
  return(SS) 
}


#####################################################
##########################################################
####Iterating MAXIMAZATION STEP:######
### This is the maximization step of the EM algorithm (M-step) which I have defined it using the function EM
EM<- function(bl,bu,Freq,theta_init,maxit=1000,tol1=1e-3,tol2=1e-4){ 
  flag<- 0 
  Mu_cur<- theta_init[1] 
  S_cur<- theta_init[2]  
  
  for (i in 1:maxit){ 
    cur<- c(Mu_cur,S_cur) 
    
    Munew<- Mest(theta=cur,bl,bu,Freq) 
    
    SSnew<- SSest(theta=cur,bl,bu,
                  muupdate=Mest(theta=cur,bl,bu,Freq) ,Freq) 
    
    Mu_new<- Munew  
    #print(Mu_new)
    S_new<- sqrt(SSnew)  
    #print(S_new)
    new_step<- c(Mu_new,S_new) 
    #print(abs(cur[1]-new_step[1])<tol1)
    #print(abs(cur[2]-new_step[2])<tol2)
    
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
#########################################################################################################
outParent<- EM(bl=blP,bu=buP,Freq=FrP,theta_init=thetaintP,maxit=1000,tol1=1e-3,tol2=1e-4)
outParent

outChildren<- EM(bl=blC,bu=buC,Freq=FrC,theta_init=thetaintC,maxit=1000,tol1=1e-3,tol2=1e-4)
outChildren

Emres<- c(outParent[1],outChildren[1],outParent[2],outChildren[2])
Emres
#########################################################################################################
### MCEM ALGORITHM For simulated data and Galton data###

library(truncnorm) 
###Calculation of updates for Mu###
##E-Step for Mu estimation: Simulating Z's###
ZMCEM<- function(theta,data){ 
  k<- 1000 
  sim<- matrix(rep(0,k*nrow(data)),ncol=nrow(data)) 
  #sim<- matrix(rep(0,k*3),ncol=3)
  
  for(i in 1 :nrow(data)) {
    sim[,i]<- rtruncnorm(k,a=data[i,1],b=data[i,2],mean=theta[1],sd=theta[2]) }
  return(sim) 
}  
### E-Step for Mu & Sigma###
MuMCEM<- function(data,simZ){ 
  n<- sum(data[,3])
  Z<- colMeans(simZ)
  numerator<- rep(0,nrow(data)) 
  for (i in 1:nrow(data)) {
    numerator[i]<- data[i,3]*Z[i]}
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
  for (i in 1:nrow(data)){
    numerator[i]<- data[i,3]*SZNEW[i] }
  
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
  return(sig) ### return the sig (updated sigma from the E-step)
}

##########################################
###MONTE CARLO EM###
##########################################################
### This is the Iteration of maximization step of the MCEM algorithm (M-step) which I have defined it using the function MCEM

MCEM<- function(data,theta_init,maxit=1000,tol1=1e-2,tol2=1e-3){ 
  flag<- 0 
  Mu_cur<- theta_init[1] 
  S_cur<- theta_init[2]  
  #theta_cur<- c(Mu_cur,S_cur)
  iter<- rep(0,maxit)  
  Svec<- rep(0,maxit)
  Mvec<- rep(0,maxit)
  for (i in 1:maxit){ 
    #print(paste("Iteration number=", i))
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
  
  #update=list(paste("updated Mu of the EM=" , Mu_cur),paste("Updated Sigma of the EM=" ,S_cur),
  #    paste("iteration number=",i)) 
  update=c(Mu_cur,(S_cur)^2) 
  return(update) 
}
#############################################################################################################
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

###########################################
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

mcemres<- c(MCEMParent[1],MCEMChildren[1],MCEMParent[2],MCEMChildren[2])
mcemres
######################################################################################################
##################################################################################################
#####Standard Errors of EM estimates-Univariate ##################
Muvarstd<- function(thetaupd,bl,bu,Freq,Data){ 
  Wj<- rep(0,length(bl)) 
  WjN<- rep(0,length(bl)) 
  
  astar<- rep(0,length(bl)) 
  bstar<- rep(0,length(bl)) 
  Mstdj<- rep(0,length(bl))
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-thetaupd[1])/thetaupd[2] 
    astar[i]<- (bl[i]-thetaupd[1])/thetaupd[2] 
    
  }
  
  
  dinom1<- NULL
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  for(i in 1:length(bl)){ 
    dinom1[i]<- (pnorm(bstar[i])-pnorm(astar[i]))
    if(dinom1[i]==0) { WjN[i]<- (-thetaupd[2])*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001)}
    else {WjN[i]<- (-thetaupd[2])*((dnorm(bstar[i])-dnorm(astar[i]))/dinom1[i])}
  }
  
  Ej2<- rep(0,length(bl)) 
  
  dinom<- NULL
  
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  
  for(i in 1:length(bl)){ 
    dinom[i]<- (pnorm(bstar[i])-pnorm(astar[i])) 
    if(dinom[i]==0) {Ej2[i]<- (bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/0.0001}
    
    
    else{Ej2[i]<- (bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/dinom[i]}
  } 
  SMU<- WjN/(thetaupd[2]^2)
  SVAR<- (-1/(2*(thetaupd[2]^2)))*Ej2
  Sj1<- rbind(SMU,SVAR)
 
  IME<- array(rep(0,2*2*(nrow(Data))),c(2,2,nrow(Data))) 
  for(i in 1:nrow(Data)){
    IME[,,i]<- Sj1[,i]%*%t(Sj1[,i])*Data[i,3]
  }
  
  InfME<- apply(IME,c(1,2),sum)
  
  var_EM_est<- solve(InfME)
  std_EM_est<- sqrt(diag(var_EM_est))
  

  return(std_EM_est) 
  
}
##################################################################################################################
c(outParent[1],sqrt(outParent[2]))

#c(outParent[1],sqrt(outParent[2]))
#[1] 68.300256  1.801338
c(outChildren[1],sqrt(outChildren[2]))
#[1] 68.098344  2.551413


se_EM_Parent<- Muvarstd(thetaupd=c(outParent[1],sqrt(outParent[2])),
                        bl=Galton2$TL,bu=Galton2$TU,Freq=Galton2$Freq,Data=Galton2)
se_EM_Parent

#[1] 0.03818231 0.06807437


se_EM_Children<- Muvarstd(thetaupd=c(outChildren[1],sqrt(outChildren[2])),
                          bl=Galton3$TL,bu=Galton3$TU,Freq=Galton3$Freq,Data=Galton3)
se_EM_Children

#[1] 0.05232407 0.14653837


lower_MU_P<- outParent[1]-qnorm(0.975)*se_EM_Parent[1]
Upper_MU_P<- outParent[1]+qnorm(0.975)*se_EM_Parent[1]
CI_Mu_P<- c(lower_MU_P,Upper_MU_P)
CI_Mu_P

lower_Sigma2_P<- outParent[2]-qnorm(0.975)*se_EM_Parent[2]
Upper_Sigma2_P<- outParent[2]+qnorm(0.975)*se_EM_Parent[2]
CI_Sigma2_P<- c(lower_Sigma2_P,Upper_Sigma2_P)
CI_Sigma2_P
sqrt(CI_Sigma2_P)

lower_MU_C<- outChildren[1]-qnorm(0.975)*se_EM_Children[1]
Upper_MU_C<- outChildren[1]+qnorm(0.975)*se_EM_Children[1]
CI_Mu_C<- c(lower_MU_C,Upper_MU_C)
CI_Mu_C

lower_Sigma2_C<- outChildren[2]-qnorm(0.975)*se_EM_Children[2]
Upper_Sigma2_C<- outChildren[2]+qnorm(0.975)*se_EM_Children[2]
CI_Sigma2_C<- c(lower_Sigma2_C,Upper_Sigma2_C)
CI_Sigma2_C
sqrt(CI_Sigma2_C)


outParent[1]
CI_Mu_P
outParent[2]
CI_Sigma2_P

outChildren[1]
outChildren[2]
CI_Mu_C
CI_Sigma2_C
############################################################################################################
### Simulate data from MCEM Estimates###
library(truncnorm) 

ZsimMCEM<- function(theta,data){ 
  k<- 1000 
  sim<- matrix(rep(0,k*nrow(data)),ncol=nrow(data)) 
  #sim<- matrix(rep(0,k*3),ncol=3)
  
  for(i in 1 :nrow(data)) {
    sim[,i]<- rtruncnorm(k,a=data[i,1],b=data[i,2],mean=theta[1],sd=theta[2]) }
  
  wj<- matrix(rep(0,k*nrow(data)),ncol=nrow(data))
  wj2<- matrix(rep(0,k*nrow(data)),ncol=nrow(data))
  
  for(i in 1:nrow(data)){
    wj[,i]<- sim[,i]-theta[1]
    wj2[,i]<- (sim[,i]-theta[1])^2
  }
  #print(wj)
  #print(wj2)
  
  d2<- array(rep(0,2*2*nrow(data)),c(2,2,nrow(data)))
  for (i in 1: nrow(data)){
    
    d2[1,1,i]<- 1/theta[2]^2
    d2[1,2,i]<- d2[2,1,i]<- (1/theta[2]^4)*mean(wj[,i])
    d2[2,2,i]<- (-1/(2*theta[2]^4))+((1/theta[2]^6)*mean(wj2[,i]))
  }
  
  
  #print(d2)
  
  d1mu<- wj/(theta[2]^2)  
  #print(d1mu)
  d1sigma2<- (-1/(2*theta[2]^2))+(wj2/(2*theta[2]^4))  
  #print(d1sigma2)
  
  dd1<- array(rep(0,2*k*nrow(data)),c(k,2,nrow(data)))
  for(i in 1:nrow(data)){
    dd1[,1,i]<- d1mu[,i]
    dd1[,2,i]<- d1sigma2[,i]  
  }
  #print(dd1)
  
  
  d1<- array(rep(0,2*1*nrow(data)),c(2,1,nrow(data)))
  for(i in 1:nrow(data)){
    d1[1,1,i]<- ((mean(wj[,i]))/(theta[2]^2))
    
    d1[2,1,i]<- (-1/(2*theta[2]^2))+((mean(wj2[,i]))/(2*theta[2]^4))  
  }
  #print(d1)
  
  diff1<- array(rep(0,2*k*nrow(data)),c(k,2,nrow(data)))
  for(i in 1:nrow(data)){
    diff1[,1,i]<- dd1[,1,i]-d1[1,1,i] 
    diff1[,2,i]<- dd1[,2,i]-d1[2,1,i]  
  }
  #print(diff1)
  
  diff1_2<- array(rep(0,2*2*k*nrow(data)),c(2,2,k,nrow(data)))
  for(i in 1:nrow(data)){
    for(j in 1:k){
      diff1_2[,,j,i]<- diff1[j,,i]%*%t(diff1[j,,i])
      
    }
  }
  #print(diff1_2)
  
  mydif<- apply(diff1_2,c(1,2,4),mean)
  #print(mydif)
  
  mydifmat<- array(rep(0,2*2*nrow(data)),c(2,2,nrow(data)))
  for(i in 1:nrow(data)){
    mydifmat[,,i]<- (d2[,,i]-mydif[,,i])*data[i,ncol(data)]
  }
  #print(mydifmat)
  
  myinf_init<- apply(mydifmat,c(1,2),sum)
  #print(myinf_init)
  myinf_final<- solve(myinf_init)
  #print(myinf_final)
  
  #se_MCEM_est<- c(sqrt(myinf_final[1,1]),sqrt(myinf_final[2,2]))
  se_MCEM_est<- sqrt(diag(myinf_final))
  names(se_MCEM_est)<- c("std_mu","std_sigma2")
  return(se_MCEM_est) 
}  
##########################################################################################################
MCEMParent

theta1<- c(MCEMParent[1],sqrt(MCEMParent[2]))

theta2<- c(MCEMChildren[1],sqrt(MCEMChildren[2]))

se_MCEM_Parent<- ZsimMCEM(theta=theta1,data=GD) 
se_MCEM_Parent

se_EM_Parent
#mydat
se_MCEM_Children<- ZsimMCEM(theta=theta2,data=mydat) 
se_MCEM_Children

se_EM_Children
######################################################################################################
######################################################################################################
GaltonRes<- matrix(rep(0,6*3),ncol=6)
colnames(GaltonRes)<- c("MuParent","std_mu_P","MuChildren","std_mu_C","VParent","VChildren")
row.names(GaltonRes)<- c("MLEExact","EM","MCEM")
GaltonRes

#Mleexc
#GaltonRes[1,]<- mle1
GaltonRes[,1]<- c(Mleexc[1],outParent[1],MCEMParent[1])
GaltonRes[,2]<- c(0,se_EM_Parent[1],se_MCEM_Parent[1])
GaltonRes[,3]<- c(Mleexc[2],outChildren[1],MCEMChildren[1]) 
GaltonRes[,4]<- c(0,se_EM_Children[1],se_MCEM_Children[1])
GaltonRes[,5]<- c(Mleexc[3],outParent[2],MCEMParent[2])
#GaltonRes[,2]<- c(0,se_EM_Parent[2],se_MCEM_Parent[2])
GaltonRes[,6]<- c(Mleexc[4],outChildren[2],MCEMChildren[2])
#GaltonRes[,4]<- c(0,se_EM_Children[2],se_MCEM_Children[2])


GaltonRes

library(xtable)
xtable(GaltonRes,digits = 5,"Parameters estimation for Univariate Galton Data")
