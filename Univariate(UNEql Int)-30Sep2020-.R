############SIMULATIONS#####################################################
rm(list=ls())

#####Univariate Simulations#####
### n=50 & 10 Classes###
set.seed(5555)

sim2<- matrix(rep(0,50*30),ncol=30)
#ncol(sim1)

for(i in 1:ncol(sim2)){
  sim2[,i]<- rnorm(50,68,1.80)
}


#dim(sim2)
Fr<- matrix(rep(0,10*ncol(sim2)),ncol=ncol(sim2))
#Fr<- matrix(rep(0,8*ncol(sim)),ncol=ncol(sim))

for(i in 1:ncol(sim2)){
  Fr[,i]<- table(cut(sim2[,i],breaks=c(-Inf,62,65,67.5,68,70,71,73.5,74,77,Inf)))
  
  }
#print(Fr)

simdata2<- array(rep(0,10*3*30),c(10,3,30))
med2<- array(rep(0,10*2*30),c(10,2,30))
for(i in 1:ncol(sim2)){
  simdata2[,1,i]<- c(-Inf,62,65,67.5,68,70,71,73.5,74,77)
  simdata2[,2,i]<- c(62,65,67.5,68,70,71,73.5,74,77,Inf)
  simdata2[,3,i]<- Fr[,i]
  med2[,1,i]<- (simdata2[,1,i]+simdata2[,2,i])/2
  med2[1,1,i]<- 60.5
  med2[10,1,i]<- 78.5
  med2[,2,i]<- Fr[,i]
}
#simdata2
#simdata2
#med2


Mynewarray2<- array(rep(0*10*2*30),c(10,2,30))
for (i in 1:30){
  Mynewarray2[,1,i]<- (simdata2[,1,i]+simdata2[,2,i])/2
  Mynewarray2[,2,i]<- simdata2[,3,i]
  Mynewarray2[1,1,i]<- 60.5
  Mynewarray2[10,1,i]<- 78.5
}
#Mynewarray

newdata2<- matrix(rep(0*30*50),ncol = 30,nrow=50)
for (i in 1:30){
  newdata2[,i]<- rep(Mynewarray2[,1,i],Mynewarray2[,2,i])
}
#newdata2[,1]

Me2<- colMeans(newdata2)
SS2<- rep(0,30)
SDD2<- rep(0,30)
for(i in 1:30){
  SS2[i]<- (sum((newdata2[,i]-Me2[i])^2))/nrow(newdata2)
  SDD2[i]<- sqrt(SS2[i])
}
#SS2
#SDD2


#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
#######################################################################
#####Univariate Simulations#####
### n=150 & 10 Classes###

sim4<- matrix(rep(0,150*30),ncol=30)
#ncol(sim1)

#set.seed(789)
for(i in 1:ncol(sim4)){
  sim4[,i]<- rnorm(150,68,1.80)
}


Fr<- matrix(rep(0,10*ncol(sim4)),ncol=ncol(sim4))
#Fr<- matrix(rep(0,8*ncol(sim)),ncol=ncol(sim))

for(i in 1:ncol(sim4)){
  Fr[,i]<- table(cut(sim4[,i],breaks=c(-Inf,62,65,67.5,68,70,71,73.5,74,77,Inf)))
  
}
#print(Fr)

simdata4<- array(rep(0,10*3*30),c(10,3,30))
med4<- array(rep(0,10*2*30),c(10,2,30))
for(i in 1:ncol(sim4)){
  simdata4[,1,i]<- c(-Inf,62,65,67.5,68,70,71,73.5,74,77)
  simdata4[,2,i]<- c(62,65,67.5,68,70,71,73.5,74,77,Inf)
  simdata4[,3,i]<- Fr[,i]
  med4[,1,i]<- (simdata4[,1,i]+simdata4[,2,i])/2
  med4[1,1,i]<- 60.5
  med4[10,1,i]<- 78.5
  med4[,2,i]<- Fr[,i]
}
#simdata4


Mynewarray4<- array(rep(0*10*2*30),c(10,2,30))
for (i in 1:30){
  Mynewarray4[,1,i]<- (simdata4[,1,i]+simdata4[,2,i])/2
  Mynewarray4[,2,i]<- simdata4[,3,i]
  Mynewarray4[1,1,i]<- 60.5
  Mynewarray4[10,1,i]<- 78.5
}
#Mynewarray

newdata4<- matrix(rep(0*30*150),ncol = 30,nrow=150)
for (i in 1:30){
  newdata4[,i]<- rep(Mynewarray4[,1,i],Mynewarray4[,2,i])
}
#newdata4

Me4<- colMeans(newdata4)
SS4<- rep(0,30)
SDD4<- rep(0,30)
for(i in 1:30){
  SS4[i]<- (sum((newdata4[,i]-Me4[i])^2))/nrow(newdata4)
  SDD4[i]<- sqrt(SS4[i])
}
#SS4
#SDD4

#Me4<- Memle4
#SDD4<- stdmle4

############################################################################
#########################################################
###n=1000, 10 Claases
sim6<- matrix(rep(0,1000*30),ncol=30)

#set.seed(666)
for(i in 1:ncol(sim6)){
  sim6[,i]<- rnorm(1000,68,1.80)
}


Fr<- matrix(rep(0,10*ncol(sim6)),ncol=ncol(sim6))

for(i in 1:ncol(sim6)){
  #sim1.cut<- cut(sim[,1],breaks=c(-Inf,65,66,67,68,69,70,71,Inf))
  Fr[,i]<- table(cut(sim6[,i],breaks=c(-Inf,62,65,67.5,68,70,71,73.5,74,77,Inf)))
  #print(Fr)
}


simdata6<- array(rep(0,10*3*30),c(10,3,30))
med6<- array(rep(0,10*2*30),c(10,2,30))

for(i in 1:ncol(sim6)){
  simdata6[,1,i]<- c(-Inf,62,65,67.5,68,70,71,73.5,74,77)
  simdata6[,2,i]<- c(62,65,67.5,68,70,71,73.5,74,77,Inf)
  simdata6[,3,i]<- Fr[,i]
  med6[,1,i]<- (simdata6[,1,i]+simdata6[,2,i])/2
  med6[1,1,i]<- 60.5
  med6[10,1,i]<- 78.5
  med6[,2,i]<- Fr[,i]
}
#simdata6


Mynewarray6<- array(rep(0*10*2*30),c(10,2,30))
for (i in 1:30){
  Mynewarray6[,1,i]<- (simdata6[,1,i]+simdata6[,2,i])/2
  Mynewarray6[,2,i]<- simdata6[,3,i]
  Mynewarray6[1,1,i]<- 60.5
  Mynewarray6[10,1,i]<- 78.5
}
#Mynewarray

newdata6<- matrix(rep(0*30*1000),ncol = 30,nrow=1000)
for (i in 1:30){
  newdata6[,i]<- rep(Mynewarray6[,1,i],Mynewarray6[,2,i])
}
#newdata6

Me6<- colMeans(newdata6)
SS6<- rep(0,30)
SDD6<- rep(0,30)
for(i in 1:30){
  SS6[i]<- (sum((newdata6[,i]-Me6[i])^2))/nrow(newdata6)
  SDD6[i]<- sqrt(SS6[i])
}
#SS6
#SDD6

#Me6<- Memle6
#SDD6<- stdmle6
############################################################################
###n=500, 10 Claases
sim7<- matrix(rep(0,500*30),ncol=30)

#set.seed(666)
for(i in 1:ncol(sim7)){
  sim7[,i]<- rnorm(500,68,1.80)
}


Fr<- matrix(rep(0,10*ncol(sim7)),ncol=ncol(sim7))

for(i in 1:ncol(sim7)){
  #sim1.cut<- cut(sim[,1],breaks=c(-Inf,65,66,67,68,69,70,71,Inf))
  Fr[,i]<- table(cut(sim7[,i],breaks=c(-Inf,62,65,67.5,68,70,71,73.5,74,77,Inf)))
  #print(Fr)
}


simdata7<- array(rep(0,10*3*30),c(10,3,30))
med7<- array(rep(0,10*2*30),c(10,2,30))

for(i in 1:ncol(sim7)){
  simdata7[,1,i]<- c(-Inf,62,65,67.5,68,70,71,73.5,74,77)
  simdata7[,2,i]<- c(62,65,67.5,68,70,71,73.5,74,77,Inf)
  simdata7[,3,i]<- Fr[,i]
  med7[,1,i]<- (simdata7[,1,i]+simdata7[,2,i])/2
  med7[1,1,i]<- 60.5
  med7[10,1,i]<- 78.5
  med7[,2,i]<- Fr[,i]
}
#simdata6


Mynewarray7<- array(rep(0*10*2*30),c(10,2,30))
for (i in 1:30){
  Mynewarray7[,1,i]<- (simdata7[,1,i]+simdata7[,2,i])/2
  Mynewarray7[,2,i]<- simdata7[,3,i]
  Mynewarray7[1,1,i]<- 60.5
  Mynewarray7[10,1,i]<- 78.5
}
#Mynewarray

newdata7<- matrix(rep(0*30*500),ncol = 30,nrow=500)
for (i in 1:30){
  newdata7[,i]<- rep(Mynewarray7[,1,i],Mynewarray7[,2,i])
}
#newdata6

Me7<- colMeans(newdata7)
SS7<- rep(0,30)
SDD7<- rep(0,30)
for(i in 1:30){
  SS7[i]<- (sum((newdata7[,i]-Me7[i])^2))/nrow(newdata7)
  SDD7[i]<- sqrt(SS7[i])
}
#SS6
#SDD6

#Me6<- Memle6
#SDD6<- stdmle6
##################################################################
#############################################################################################################


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

#####FOR SIM DATA 2############################################
MUMult2<- rep(0,30)
SDMult2<- rep(0,30)

for(j in 1:30){
  TL<- simdata2[,1,j]
  freq<- simdata2[,3,j]
  res2 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))
  #estimsim<- ressim$par
  MUMult2[j]<- res2$par[1]/res2$par[2]
  SDMult2[j]<- 1/res2$par[2]
}

#MUMult2
#SDMult2
########################################################
rm("TL","freq","i","j")

#####FOR SIM DATA 4 ############################################
MUMult4<- rep(0,30)
SDMult4<- rep(0,30)


for(j in 1:30){
  TL<- simdata4[,1,j]
  freq<- simdata4[,3,j]
  res4 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))
  #estimsim<- ressim$par
  MUMult4[j]<- res4$par[1]/res4$par[2]
  SDMult4[j]<- 1/res4$par[2]
}

#MUMult4
#SDMult4
########################################################
rm("TL","freq","i","j")

#####FOR SIM DATA 6############################################
MUMult6<- rep(0,30)
SDMult6<- rep(0,30)


for(j in 1:30){
  TL<- simdata6[,1,j]
  freq<- simdata6[,3,j]
  res6 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))
  #estimsim<- ressim$par
  MUMult6[j]<- res6$par[1]/res6$par[2]
  SDMult6[j]<- 1/res6$par[2]
}

#MUMult6
#SDMult6
########################################################
rm("TL","freq","i","j")
#####FOR SIM DATA 7############################################33
MUMult7<- rep(0,30)
SDMult7<- rep(0,30)


for(j in 1:30){
  TL<- simdata7[,1,j]
  freq<- simdata7[,3,j]
  res7 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))
  #estimsim<- ressim$par
  MUMult7[j]<- res7$par[1]/res7$par[2]
  SDMult7[j]<- 1/res7$par[2]
}

#MUMult7
#SDMult7
########################################################
#ls()
rm("TL","freq","i","j")
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
########################################################
########################################################
### 30 SIMULATION RESULT, n=50 & 10 Classes####
output2<- matrix(rep(0,2*30),ncol=2)
colnames(output2)<- c("mean","std")


for(i in 1:30){
  output2[i,]<- EM(bl=simdata2[,1,i],bu=simdata2[,2,i],Freq=simdata2[,3,i],theta_init=c(67,2),
                   maxit = 1000,tol1=1e-3,tol2=1e-4)
  
}
#print (output2)
#simdata1    
####################################################################

########################################################
########################################################
### 30 SIMULATION RESULT, n=150 & 10 Classes####
output4<- matrix(rep(0,2*30),ncol=2)
colnames(output4)<- c("mean","std")


for(i in 1:30){
  output4[i,]<- EM(bl=simdata4[,1,i],bu=simdata4[,2,i],Freq=simdata4[,3,i],theta_init=c(67,2)
                   ,maxit = 1000,tol1=1e-3,tol2=1e-4)
  
}
#print (output4)
#simdata1    
####################################################################
### 30 SIMULATION RESULT, n=100 & 7 Classes####

output5<- matrix(rep(0,2*30),ncol=2)
colnames(output5)<- c("mean","std")

for(i in 1:30){
  output5[i,]<- EM(bl=simdata5[,1,i],bu=simdata5[,2,i],Freq=simdata5[,3,i],theta_init=c(67,2),
                   maxit = 1000,tol1=1e-3,tol2=1e-4)
  
}

#print (output5)
##########################################################
### 30 SIMULATION RESULT, n=1000 & 10 Classes####
output6<- matrix(rep(0,2*30),ncol=2)
colnames(output6)<- c("mean","std")


for(i in 1:30){
  output6[i,]<- EM(bl=simdata6[,1,i],bu=simdata6[,2,i],Freq=simdata6[,3,i],theta_init=c(67,2),
                   maxit = 1000,tol1=1e-3,tol2=1e-4)
  
}

#print(output6)
##########################################################
### 30 SIMULATION RESULT, n=500 & 10 Classes####
output7<- matrix(rep(0,2*30),ncol=2)
colnames(output7)<- c("mean","std")


for(i in 1:30){
  output7[i,]<- EM(bl=simdata7[,1,i],bu=simdata7[,2,i],Freq=simdata7[,3,i],theta_init=c(67,2),
                   maxit = 1000,tol1=1e-3,tol2=1e-4)
  
}


######################################################
######################################################
ls()
rm("Fr","Mest","Fr","i","SSest","EM","Logll")
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
##### Simulated data, n=50 & 10 Classes
outputMCEM2<- matrix(rep(0,2*30),ncol=2)
colnames(outputMCEM2)<- c("mean","std")

for(i in 1:30){
  #i<- 1 
  mydat<- simdata2[,,i]
  outputMCEM2[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3)
  
}
#print (outputMCEM2)

######################################################################

##################################################################
##### Simulated data, n=150 & 10 Classes
outputMCEM4<- matrix(rep(0,2*30),ncol=2)
colnames(outputMCEM4)<- c("mean","std")

for(i in 1:30){
  #i<- 1 
  mydat<- simdata4[,,i]
  outputMCEM4[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3)
  
}
#print (outputMCEM4)


#########################################################################
##### Simulated data, n=1000 & 10 Classes

outputMCEM6<- matrix(rep(0,2*30),ncol=2)
colnames(outputMCEM6)<- c("mean","std")

for(i in 1:30){
  #i<- 1 
  mydat<- simdata6[,,i]
  outputMCEM6[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3)
  
}

#print (outputMCEM6)
#######################################################################
#######Simulated data, n=500 & 10 Classes

outputMCEM7<- matrix(rep(0,2*30),ncol=2)
colnames(outputMCEM7)<- c("mean","std")

for(i in 1:30){
  #i<- 1 
  mydat<- simdata7[,,i]
  outputMCEM7[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3)
  
}

##########################################################################################
###########COMPARISONS, n=50 & 10 classes#########################
###MEAN###
MCompare2<- matrix(rep(0,30*4),ncol=4)
MCompare2[,1]<- Me2
MCompare2[,2]<- output2[,1]
MCompare2[,3]<- outputMCEM2[,1]
MCompare2[,4]<- MUMult2
#MUMult2
#SDMult2
#MCompare
colnames(MCompare2)<- c("MLEUngrp","EM","MCEM","MLEExact")
#MCompare2


MmuMLE2<-mean(MCompare2[,1])
MmuEM2<-mean(MCompare2[,2])
MmuMCEM2<-mean(MCompare2[,3])

MSDMLE2<-sd(MCompare2[,1])
MSDEM2<-sd(MCompare2[,2])
MSDMCEM2<-sd(MCompare2[,3])

BMLm2<- mean(MCompare2[,1]-68)
BEMm2<- mean(MCompare2[,2]-68)
BMCEMm2<- mean(MCompare2[,3]-68)

MSEMLm2<-(sum((MCompare2[,1]-68)^2))/(length(MCompare2[,1]))
MSEEMm2<-(sum((MCompare2[,2]-68)^2))/(length(MCompare2[,2]))
MSEMCEMm2<-(sum((MCompare2[,3]-68)^2))/(length(MCompare2[,3]))


MmuMLEX2<-mean(MCompare2[,4])
MSDMLEX2<-sd(MCompare2[,4])
BMLEXm2<- mean(MCompare2[,4]-68)
MSEMLEXm2<-(sum((MCompare2[,4]-68)^2))/(length(MCompare2[,4]))


#############################################################################################
###STD###
SDCompare2<- matrix(rep(0,30*4),ncol=4)
SDCompare2[,1]<- (SS2)
SDCompare2[,2]<- (output2[,2])
SDCompare2[,3]<- outputMCEM2[,2]
SDCompare2[,4]<- (SDMult2)^2

colnames(SDCompare2)<- c("MLEUngrp","EM","MCEM","MLEExact")
#SDCompare2
#mean(SS2)
MSSDMLE2<-mean(SDCompare2[,1])
MSSDEM2<-mean(SDCompare2[,2])
MSSDMCEM2<-mean(SDCompare2[,3])

SSDMLE2<-sd(SDCompare2[,1])
SSDEM2<-sd(SDCompare2[,2])
SSDMCEM2<-sd(SDCompare2[,3])

BMLs2<- mean(SDCompare2[,1]-3.24)
BEMs2<- mean(SDCompare2[,2]-3.24)
BMCEMs2<- mean(SDCompare2[,3]-3.24)



MSEMLs2<-(sum((SDCompare2[,1]-3.24)^2))/(length(SDCompare2[,1]))
MSEEMs2<-(sum((SDCompare2[,2]-3.24)^2))/(length(SDCompare2[,2]))
MSEMCEMs2<-(sum((SDCompare2[,3]-3.24)^2))/(length(SDCompare2[,3]))

#SDCompare2[,4]<- SDMult2
MSSDMLEX2<-mean(SDCompare2[,4])
SSDMLEX2<-sd(SDCompare2[,4])
BMLEXs2<- mean(SDCompare2[,4]-3.24)
MSEMLEXs2<-(sum((SDCompare2[,4]-3.24)^2))/(length(SDCompare2[,4]))



##################################################################################
#####################################################################
###########COMPARISONS, n=150 & 10 classes#########################
MCompare4<- matrix(rep(0,30*4),ncol=4)
MCompare4[,1]<- Me4
MCompare4[,2]<- output4[,1]
MCompare4[,3]<- outputMCEM4[,1]
MCompare4[,4]<- MUMult4
#MCompare
colnames(MCompare4)<- c("MLEUngrp","EM","MCEM","MLEExact")
#MCompare4

MmuMLE4<-mean(MCompare4[,1])
MmuEM4<-mean(MCompare4[,2])
MmuMCEM4<-mean(MCompare4[,3])

MSDMLE4<-sd(MCompare4[,1])
MSDEM4<-sd(MCompare4[,2])
MSDMCEM4<-sd(MCompare4[,3])


BMLm4<- mean(MCompare4[,1]-68)
BEMm4<- mean(MCompare4[,2]-68)
BMCEMm4<- mean(MCompare4[,3]-68)

MSEMLm4<-(sum((MCompare4[,1]-68)^2))/(length(MCompare4[,1]))
MSEEMm4<-(sum((MCompare4[,2]-68)^2))/(length(MCompare4[,2]))
MSEMCEMm4<-(sum((MCompare4[,3]-68)^2))/(length(MCompare4[,3]))


MmuMLEX4<-mean(MCompare4[,4])
MSDMLEX4<-sd(MCompare4[,4])
BMLEXm4<- mean(MCompare4[,4]-68)
MSEMLEXm4<-(sum((MCompare4[,4]-68)^2))/(length(MCompare4[,4]))
###############################################
###STD###
SDCompare4<- matrix(rep(0,30*4),ncol=4)
SDCompare4[,1]<- (SS4)
SDCompare4[,2]<- output4[,2]
SDCompare4[,3]<- outputMCEM4[,2]
SDCompare4[,4]<- (SDMult4)^2
colnames(SDCompare4)<- c("MLEUngrp","EM","MCEM","MLEExact")
#SDCompare2

MSSDMLE4<-mean(SDCompare4[,1])
MSSDEM4<-mean(SDCompare4[,2])
MSSDMCEM4<-mean(SDCompare4[,3])

SSDMLE4<-sd(SDCompare4[,1])
SSDEM4<-sd(SDCompare4[,2])
SSDMCEM4<-sd(SDCompare4[,3])


BMLs4<- mean(SDCompare4[,1]-3.24)
BEMs4<- mean(SDCompare4[,2]-3.24)
BMCEMs4<- mean(SDCompare4[,3]-3.24)

MSEMLs4<-(sum((SDCompare4[,1]-3.24)^2))/(length(SDCompare4[,1]))
MSEEMs4<-(sum((SDCompare4[,2]-3.24)^2))/(length(SDCompare4[,2]))
MSEMCEMs4<-(sum((SDCompare4[,3]-3.24)^2))/(length(SDCompare4[,3]))


#SDCompare4[,4]<- SDMult4
MSSDMLEX4<-mean(SDCompare4[,4])
SSDMLEX4<-sd(SDCompare4[,4])
BMLEXs4<- mean(SDCompare4[,4]-3.24)
MSEMLEXs4<-(sum((SDCompare4[,4]-3.24)^2))/(length(SDCompare4[,4]))

##################################################################################
#################################################################
##########COMPARISONS, n=1000 & 10 classes#########################
MCompare6<- matrix(rep(0,30*4),ncol=4)
MCompare6[,1]<- Me6
MCompare6[,2]<- output6[,1]
MCompare6[,3]<- outputMCEM6[,1]
MCompare6[,4]<- MUMult6

#MCompare
colnames(MCompare6)<- c("MLEUngrp","EM","MCEM","MLEExact")
#MCompare5


MmuMLE6<-mean(MCompare6[,1])
MmuEM6<-mean(MCompare6[,2])
MmuMCEM6<-mean(MCompare6[,3])

MSDMLE6<-sd(MCompare6[,1])
MSDEM6<-sd(MCompare6[,2])
MSDMCEM6<-sd(MCompare6[,3])


BMLm6<- mean(MCompare6[,1]-68)
BEMm6<- mean(MCompare6[,2]-68)
BMCEMm6<- mean(MCompare6[,3]-68)

MSEMLm6<-(sum((MCompare6[,1]-68)^2))/(length(MCompare6[,1]))
MSEEMm6<-(sum((MCompare6[,2]-68)^2))/(length(MCompare6[,2]))
MSEMCEMm6<-(sum((MCompare6[,3]-68)^2))/(length(MCompare6[,3]))


MmuMLEX6<-mean(MCompare6[,4])
MSDMLEX6<-sd(MCompare6[,4])
BMLEXm6<- mean(MCompare6[,4]-68)
MSEMLEXm6<-(sum((MCompare6[,4]-68)^2))/(length(MCompare6[,4]))

###########################################
###STD###

SDCompare6<- matrix(rep(0,30*4),ncol=4)
SDCompare6[,1]<- (SS6)
SDCompare6[,2]<- output6[,2]
SDCompare6[,3]<- outputMCEM6[,2]
SDCompare6[,4]<- (SDMult6)^2
colnames(SDCompare6)<- c("MLEUngrp","EM","MCEM","MLEExact")
#SDCompare5

MSSDMLE6<-mean(SDCompare6[,1])
MSSDEM6<-mean(SDCompare6[,2])
MSSDMCEM6<-mean(SDCompare6[,3])


SSDMLE6<-sd(SDCompare6[,1])
SSDEM6<-sd(SDCompare6[,2])
SSDMCEM6<-sd(SDCompare6[,3])


BMLs6<- mean(SDCompare6[,1]-3.24)
BEMs6<- mean(SDCompare6[,2]-3.24)
BMCEMs6<- mean(SDCompare6[,3]-3.24)

MSEMLs6<-(sum((SDCompare6[,1]-3.24)^2))/(length(SDCompare6[,1]))
MSEEMs6<-(sum((SDCompare6[,2]-3.24)^2))/(length(SDCompare6[,2]))
MSEMCEMs6<-(sum((SDCompare6[,3]-3.24)^2))/(length(SDCompare6[,3]))


#SDCompare6[,4]<- SDMult6
MSSDMLEX6<-mean(SDCompare6[,4])
SSDMLEX6<-sd(SDCompare6[,4])
BMLEXs6<- mean(SDCompare6[,4]-3.24)
MSEMLEXs6<-(sum((SDCompare6[,4]-3.24)^2))/(length(SDCompare6[,4]))
####################################################################


##########COMPARISONS, n=500 & 10 classes#########################
MCompare7<- matrix(rep(0,30*4),ncol=4)
MCompare7[,1]<- Me7
MCompare7[,2]<- output7[,1]
MCompare7[,3]<- outputMCEM7[,1]
MCompare7[,4]<- MUMult7

#MCompare
colnames(MCompare7)<- c("MLEUngrp","EM","MCEM","MLEExact")
#MCompare7


MmuMLE7<-mean(MCompare7[,1])
MmuEM7<-mean(MCompare7[,2])
MmuMCEM7<-mean(MCompare7[,3])

MSDMLE7<-sd(MCompare7[,1])
MSDEM7<-sd(MCompare7[,2])
MSDMCEM7<-sd(MCompare7[,3])


BMLm7<- mean(MCompare7[,1]-68)
BEMm7<- mean(MCompare7[,2]-68)
BMCEMm7<- mean(MCompare7[,3]-68)

MSEMLm7<-(sum((MCompare7[,1]-68)^2))/(length(MCompare7[,1]))
MSEEMm7<-(sum((MCompare7[,2]-68)^2))/(length(MCompare7[,2]))
MSEMCEMm7<-(sum((MCompare7[,3]-68)^2))/(length(MCompare7[,3]))


MmuMLEX7<-mean(MCompare7[,4])
MSDMLEX7<-sd(MCompare7[,4])
BMLEXm7<- mean(MCompare7[,4]-68)
MSEMLEXm7<-(sum((MCompare7[,4]-68)^2))/(length(MCompare7[,4]))

###########################################
###STD###

SDCompare7<- matrix(rep(0,30*4),ncol=4)
SDCompare7[,1]<- (SS7)
SDCompare7[,2]<- output7[,2]
SDCompare7[,3]<- outputMCEM7[,2]
SDCompare7[,4]<- (SDMult7)^2
colnames(SDCompare7)<- c("MLEUngrp","EM","MCEM","MLEExact")
#SDCompare5

MSSDMLE7<-mean(SDCompare7[,1])
MSSDEM7<-mean(SDCompare7[,2])
MSSDMCEM7<-mean(SDCompare7[,3])


SSDMLE7<-sd(SDCompare7[,1])
SSDEM7<-sd(SDCompare7[,2])
SSDMCEM7<-sd(SDCompare7[,3])


BMLs7<- mean(SDCompare7[,1]-3.24)
BEMs7<- mean(SDCompare7[,2]-3.24)
BMCEMs7<- mean(SDCompare7[,3]-3.24)

MSEMLs7<-(sum((SDCompare7[,1]-3.24)^2))/(length(SDCompare7[,1]))
MSEEMs7<-(sum((SDCompare7[,2]-3.24)^2))/(length(SDCompare7[,2]))
MSEMCEMs7<-(sum((SDCompare7[,3]-3.24)^2))/(length(SDCompare7[,3]))


#SDCompare6[,4]<- SDMult6
MSSDMLEX7<-mean(SDCompare7[,4])
SSDMLEX7<-sd(SDCompare7[,4])
BMLEXs7<- mean(SDCompare7[,4]-3.24)
MSEMLEXs7<-(sum((SDCompare7[,4]-3.24)^2))/(length(SDCompare7[,4]))

####################################################################

####################################################################
### MEAN, n=50, Classes=10 #################################################################
a1<- c("MLEUngrp","MLEExact","EM","MCEM")
n1<- c(rep(50,4))
#cl1<- c(rep(10,4))
b1<- c(round(MmuMLE2,6),round(MmuMLEX2,6),round(MmuEM2,6),round(MmuMCEM2,6))
c1<- c(round(MSDMLE2,6),round(MSDMLEX2,6),round(MSDEM2,6),round(MSDMCEM2,6))
#d1<- c(round(BMLm2,5),round(BEMm2,5),round(BMCEMm2,5))
e1<- c(round(MSEMLm2,6),round(MSEMLEXm2,6),round(MSEEMm2,6),round(MSEMCEMm2,6))

#BMLEXm2
#BMLEXs2
#rm("m1")
m1<- data.frame(cbind(a1,n1,b1,c1,e1))
colnames(m1)<- c("Method","n","mean","std","MSE")
#print(m1)
####### STD ,n=50, Classes=10 ##########################
a2<- c("MLEUngrp","MLEExact","EM","MCEM")
n2<- c(rep(50,4))
#cl2<- c(rep(10,4))
b2<- c(round(MSSDMLE2,6),round(MSSDMLEX2,6),round(MSSDEM2,6),round(MSSDMCEM2,6))
c2<- c(round(SSDMLE2,6),round(SSDMLEX2,6),round(SSDEM2,6),round(SSDMCEM2,6))
#d2<- c(round(BMLs2,5),round(BEMs2,5),round(BMCEMs2,5))
e2<- c(round(MSEMLs2,6),round(MSEMLEXs2,6),round(MSEEMs2,6),round(MSEMCEMs2,6))

m2<- data.frame(cbind(a2,n2,b2,c2,e2))
colnames(m2)<- c("Method","n","mean","std","MSE")
#m2<- as.data.frame(m2)
#print(m2)
#rbind(m1,m2)
########################################
#########################################
####### Mean ,n=150, Classes=10####################################
a4<- c("MLEUngrp","MLEExact","EM","MCEM")
n4<- c(rep(100,4))
#cl4<- c(rep(10,4))
b4<- c(round(MmuMLE4,6),round(MmuMLEX4,6),round(MmuEM4,6),round(MmuMCEM4,6))
c4<- c(round(MSDMLE4,6),round(MSDMLEX4,6),round(MSDEM4,6),round(MSDMCEM4,6))
#d4<- c(round(BMLm4,5),round(BEMm4,5),round(BMCEMm4,5))
e4<- c(round(MSEMLm4,6),round(MSEMLEXm4,6),round(MSEEMm4,6),round(MSEMCEMm4,6))


m4<- data.frame(cbind(a4,n4,b4,c4,e4))
colnames(m4)<- c("Method","n","mean","std","MSE")
#m4<- as.data.frame(m4)
#print(m4)
###### STD ,n=150, Classes=10###########################
a44<- c("MLEUngrp","MLEExact","EM","MCEM")
n44<- c(rep(100,4))
#cl44<- c(rep(10,4))
b44<- c(round(MSSDMLE4,6),round(MSSDMLEX4,6),round(MSSDEM4,6),round(MSSDMCEM4,6))
c44<- c(round(SSDMLE4,6),round(SSDMLEX4,6),round(SSDEM4,6),round(SSDMCEM4,6))
#d44<- c(round(BMLs4,5),round(BEMs4,5),round(BMCEMs4,5))
e44<- c(round(MSEMLs4,6),round(MSEMLEXs4,6),round(MSEEMs4,6),round(MSEMCEMs4,6))

m44<- data.frame(cbind(a44,n44,b44,c44,e44))
colnames(m44)<- c("Method","n","mean","std","MSE")

#m44<- as.data.frame(m44)
#print(m44)
###########################################
######################################################
##### Mean ,n=1000, Classes=10##################################################

a6<- c("MLEUngrp","MLEExact","EM","MCEM")
n6<- c(rep(1000,4))
#cl6<- c(rep(10,4))
b6<- c(round(MmuMLE6,6),round(MmuMLEX6,6),round(MmuEM6,6),round(MmuMCEM6,6))
c6<- c(round(MSDMLE6,6),round(MSDMLEX6,6),round(MSDEM6,6),round(MSDMCEM6,6))
#d6<- c(round(BMLm6,5),round(BEMm6,5),round(BMCEMm6,5))
e6<- c(round(MSEMLm6,6),round(MSEMLEXm6,6),round(MSEEMm6,6),round(MSEMCEMm6,6))

m6<- data.frame(cbind(a6,n6,b6,c6,e6))
colnames(m6)<- c("Method","n","mean","std","MSE")

#m6<- as.data.frame(m6)
#print(m5)
#### STD ,n=1000, Classes=10#############################
a66<- c("MLEUngrp","MLEExact","EM","MCEM")
n66<- c(rep(1000,4))
#cl66<- c(rep(10,4))
b66<- c(round(MSSDMLE6,6),round(MSSDMLEX6,6),round(MSSDEM6,6),round(MSSDMCEM6,6))
c66<- c(round(SSDMLE6,6),round(SSDMLEX6,6),round(SSDEM6,6),round(SSDMCEM6,6))
#d66<- c(round(BMLs6,5),round(BEMs6,5),round(BMCEMs6,5))
e66<- c(round(MSEMLs6,6),round(MSEMLEXs6,6),round(MSEEMs6,6),round(MSEMCEMs6,6))


m66<- data.frame(cbind(a66,n66,b66,c66,e66))
colnames(m66)<- c("Method","n","mean","std","MSE")

#m66<- as.data.frame(m66)
#print(m55)

######################################################
##### Mean ,n=500, Classes=10##################################################

a7<- c("MLEUngrp","MLEExact","EM","MCEM")
n7<- c(rep(500,4))
#cl7<- c(rep(10,4))
b7<- c(round(MmuMLE7,6),round(MmuMLEX7,6),round(MmuEM7,6),round(MmuMCEM7,6))
c7<- c(round(MSDMLE7,6),round(MSDMLEX7,6),round(MSDEM7,6),round(MSDMCEM7,6))
#d7<- c(round(BMLm7,5),round(BEMm7,5),round(BMCEMm7,5))
e7<- c(round(MSEMLm7,6),round(MSEMLEXm7,6),round(MSEEMm7,6),round(MSEMCEMm7,6))

m7<- data.frame(cbind(a7,n7,b7,c7,e7))
colnames(m7)<- c("Method","n","mean","std","MSE")

#m7<- as.data.frame(m7)
#print(m5)
#### STD ,n=500, Classes=10#############################
a77<- c("MLEUngrp","MLEExact","EM","MCEM")
n77<- c(rep(500,4))
#cl77<- c(rep(10,4))
b77<- c(round(MSSDMLE7,6),round(MSSDMLEX7,6),round(MSSDEM7,6),round(MSSDMCEM7,6))
c77<- c(round(SSDMLE7,6),round(SSDMLEX7,6),round(SSDEM7,6),round(SSDMCEM7,6))
#d77<- c(round(BMLs7,5),round(BEMs7,5),round(BMCEMs7,5))
e77<- c(round(MSEMLs7,6),round(MSEMLEXs7,6),round(MSEEMs7,6),round(MSEMCEMs7,6))


m77<- data.frame(cbind(a77,n77,b77,c77,e77))
colnames(m77)<- c("Method","n","mean","std","MSE")

#m77<- as.data.frame(m77)
#print(m55)

###################################################################
####################################################################
###MEAN###
#MUEST7<- rbind(m3,m5,m7)
MUEST10<- rbind(m1,m4,m7,m6)
###STD###
#STDEST7<- rbind(m33,m55,m77)
STDEST10<- rbind(m2,m44,m77,m66)
#mean(SS6)
#output3
library(xtable)
#xtable(MUEST8,caption = "Estimate of MU in 4 methods  for sample sizes: 50,100,1000 and 7 classes")
xtable(MUEST10,caption = "Estimate of MU along with its standard deviation (std) and mean squared error (MSE)
       using 4 methods for sample sizes: 50,100,500,1000 over Unequal Intervals (number of classes=10) ")
#xtable(STDEST8,"Estimate of STD in 4 methods  for sample sizes: 50,100,1000 and 7 classes")
xtable(STDEST10,"Estimate of variance along with its standard deviation (std) and mean squared error (MSE)
       using 4 methods for sample sizes: 50,100,500,1000 over Unequal Intervals (number of classes=10)")
##################################################################
pdf("MU Unequal.pdf")
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))
boxplot(MCompare2,col="orchid",ylab="Mean",ylim=c(67,68.7),las=2,cex.axis=0.70)
abline(h=68)
title("n=50")

boxplot(MCompare4,col="orchid",ylab="Mean",ylim=c(67.3,68.7),las=2,cex.axis=0.70)
abline(h=68)
title("n=150")

boxplot(MCompare7,col="orchid",ylab="Mean",ylim=c(67.3,68.7),las=2,cex.axis=0.70)
abline(h=68)
title("n=500")

boxplot(MCompare6,col="orchid",ylab="Mean",ylim=c(67.3,68.7),las=2,cex.axis=0.70)
abline(h=68)
title("n=1000")
dev.off()
###################################################################
pdf("Var UneQUAL.pdf")
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))
boxplot(SDCompare2,col="yellow",ylab="Var",ylim=c(1.6,6),las=2,cex.axis=0.70)
abline(h=3.24)
title("n=50")

boxplot(SDCompare4,col="yellow",ylab="var",ylim=c(1.6,4.8),las=2,cex.axis=0.70)
abline(h=3.24)
title("n=150")

boxplot(SDCompare7,col="yellow",ylab="var",ylim=c(1.6,4.8),las=2,cex.axis=0.70)
abline(h=3.24)
title("n=500")

boxplot(SDCompare6,col="yellow",ylab="var",ylim=c(1.6,4.8),las=2,cex.axis=0.70)
abline(h=3.24)
title("n=1000")

dev.off()
#########################################################
########################################################################
