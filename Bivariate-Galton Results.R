rm(list=ls())
library(tmvtnorm)
library(MASS)
library(xtable)
#################################################################
#Galton<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton.csv",header=TRUE)
Galton<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Galton.csv",header=TRUE)

X<- Galton$Parent
MeP<- mean(X)
SDP<- sqrt((1/length(X))*sum((X-mean(X))^2))
VP<- SDP^2
#MeP
#SDP
#VP
Y<- Galton$Children
MeC<- mean(Y)
SDC<- sqrt((1/length(Y))*sum((Y-mean(Y))^2))
VC<- SDC^2
#MeC
#SDC
#VC
#cor(Galton)
#cov(Galton)

rho<- sum((X-MeP)*(Y-MeC))/(SDC*SDP*length(X))
#rho
MLEUngrouped<- c(MeP,MeC,VP,VC,rho)
MLEUngrouped
#mean(Galton$Parent)
#mean(Galton$Children)

#muint<- c(MeP,MeC)
#sigmaint<- cov(Galton)

#Galton
#summary(Galton)
#cov(Galton)
#cor(Galton)
####MLE EXACT#############################
###DATA###
#DataBiv<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton-Bivariate.csv")
DataBiv<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Galton-Bivariate.csv")


DataBiv[1:14,1]<- -Inf
DataBiv[141:154,2]<- Inf
DataBiv[c(1,15,29,43,57,71,85,99,113,127,141),3]<- -Inf
DataBiv[c(14,28,42,56,70,84,98,112,126,140,154),4]<- Inf
#DataBiv
#nrow(DataBiv)
###############################################################################################################################
library(tmvtnorm)
library(MASS)

#############################################################################################################################
Billik<- function(Data,theta){
  #mu<- c(theta[1],theta[2])
  Sigma<- matrix(c(theta[3],theta[5]*sqrt(theta[3])*sqrt(theta[4]),theta[5]*sqrt(theta[3])*sqrt(theta[4])
                   , theta[4]),byrow=T,ncol=2)
  #print(Sigma)
  #Sigma<- matrix(c(theta[3],theta[5],theta[5],theta[4]),byrow=T,ncol=2)
  b<- 0
  for(i in 1:nrow(Data)){
    a<- pmvnorm(lower=c(Data[i,1],Data[i,3]),upper=c(Data[i,2],Data[i,4]),mean=c(theta[1],theta[2]),sigma=Sigma)
    if (a[1]==0) {a<- 1e-03}
    else {a<- a[1]}
    b<- b+Data[i,5]*log(a)
    #print(b)
  }
  ##The joint likelihood
  #print(-b)
  return(-b)
  
}

#theta<- c(65,65,3,6,2.12132)

#Billik(DataBiv,theta)
################################
###Parameter Estimation Using EM Algorithm####
#######################################################################
###E-STEP:####
###Moments###
MEXI<- function(Data,mu,sigma){ 
  n<- nrow(Data) 
  d<- (ncol(Data)-1)/2  
  
  ind<- ncol(Data)-1 
  lowerb<- NULL 
  upperb<- NULL 
  for (i in 1:ind){ 
    if (i%%2!=0) {lowerb<- append(lowerb,i)}
    else {upperb<- append(upperb,i)} 
    
  }
  Data1<- as.matrix(Data) 
  
  
  EX<- matrix(rep(0,n*d),ncol=d,nrow=n) 
  for(i in 1:nrow(Data)){ 
    mnts<- mtmvnorm(mean=mu, sigma=sigma, 
                    lower=Data1[i,c(lowerb)],
                    upper=Data1[i,c(upperb)]) 
    EX[i,]<- mnts$tmean 
  }
  return(EX) 
  
  
}


#mu1 <- c(67,67)
#sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)



#############################################################
McovXI<- function(Data,mu,sigma){ 
  n<- nrow(Data) 
  d<- (ncol(Data)-1)/2 
  
  ind<- ncol(Data)-1
  lowerb<- NULL
  upperb<- NULL
  for (i in 1:ind){ 
    if (i%%2!=0) {lowerb<- append(lowerb,i)} 
    else {upperb<- append(upperb,i)} 
    
  }
  Data1<- as.matrix(Data)
  
  s1<- array(rep(0,n*d*d),c(d,d,n)) 
  for(i in 1:nrow(Data)){ 
    mnts<- mtmvnorm(mean=mu, sigma=sigma,
                    lower=Data1[i,c(lowerb)],
                    upper=Data1[i,c(upperb)]) 
    s1[,,i]<- mnts$tvar 
  }
  return(s1) 
  
}


#############################################################
###MU:E-step###
MEM<- function(Data,EXi){ 
  F<- ncol(Data) 
  Freq<- Data[,F] 
  A<- EXi*Data[,F] 
  
  mupd<- apply(A,2,sum)/sum(Data[,F]) 
  return(mupd) 
  
}

################################################################################
###E(XX) estimate###
EXXest<- function(Data,EXi,ss1){ 
  n<- nrow(Data)
  d<- (ncol(Data)-1)/2 
  ExxNew<- array(rep(0,n*d*d),c(d,d,n)) 
  
  for(i in 1:n){ 
    ExxNew[,,i]<- ss1[,,i]+EXi[i,]%*%t(EXi[i,])
    
    
  }
  return(ExxNew)  
  
}




######################################################################
Mysig<- function(Data,EXX,EX,UPDMU){ 
  n<- nrow(Data) 
  d<- (ncol(Data)-1)/2 
  s2<- array(rep(0,n*d*d),c(d,d,n)) 
  
  for(i in 1:n){ 
    
    s2[,,i]<- (EXX[,,i]-UPDMU%*%t(EX[i,])-EX[i,]%*%t(UPDMU)+(UPDMU%*%t(UPDMU)))*Data[i,ncol(Data)] 
  }
  
  SigmaNew<- apply(s2,c(1,2),sum)/sum(Data[,ncol(Data)]) 
  
  if (ncol(SigmaNew)>2){ 
    for (i in 1:ncol(SigmaNew)){ 
      if (SigmaNew[upper.tri(SigmaNew)][i]!=SigmaNew[lower.tri(SigmaNew)][i]){ 
        
        SigmaNew[lower.tri(SigmaNew)][i]<- SigmaNew[upper.tri(SigmaNew)][i] 
        
      }
    }
    
  }
  
  
  return(SigmaNew) 
}




########################################################################
########################################################################
###M-STEP:###
### M-Step###
EMBivG<- function(Data,mu_init,sigma_init,maxit=1000,tol1=1e-4,tol2=1e-3){ 
  flag<- 0 
  M_cur<- mu_init 
  S_cur<- sigma_init 
  #iter<- rep(0,maxit) 
  
  for (i in 1:maxit){ 
    
    
    MU_cur<- M_cur 
    SS_cur<- S_cur 
    
    
    MuNew<- MEM(Data=mydat,EXi=MEXI(Data=mydat,mu= MU_cur,sigma=SS_cur)) 
    
    
    SigmaNew<- Mysig(Data=mydat,
                     EXX=EXXest(Data=mydat,
                                EXi=MEXI(Data=mydat,mu=MU_cur,sigma=SS_cur),
                                ss1=McovXI(Data=mydat,mu=MU_cur,sigma=SS_cur)),
                     EX=MEXI(Data=mydat,mu=MU_cur,sigma=SS_cur),
                     UPDMU=MEM(Data=mydat,EXi=MEXI(Data=mydat,mu=MU_cur,sigma=SS_cur)))  
    
    
    
    Mu_new<- MuNew 
    Sigma_new<- SigmaNew 
    
    
    diff1<- MU_cur-Mu_new 
    diff2<- SS_cur-Sigma_new 
    #print(diff1)
    #print(diff2)
    
    D1<- abs(mean(diff1)) 
    D2<- abs(mean(diff2)) 
    
    if(D1<tol1 & D2<tol2){flag<- 1 ;break}  
    
    M_cur<- Mu_new 
    
    S_cur<- Sigma_new 
    
    #iter[i]<- i
    
  }
  
  if(!flag) warning("Didn't Converge \n") ## the condition if there is no convergence in the stopping rule
  
  #list(paste("updated Mu of the EM=" , M_cur),paste("Updated Sigma of the EM=" ,S_cur),
  #    paste("iteration number=",i))
  #print(M_cur,S_cur)
  #print(S_cur)
  update=list(M_cur,(S_cur)) ## the list of the updated estimates of the parameters
  return(update) ## return the list of updated estimates of the parameters
  
}
######################################################################################################################
set.seed(4321)
library(tmvtnorm) 
library(MASS) 

#####E-STEP####
###SIMULATE Z###
#############################################
ZMCEM1<- function(Data,mu,sigma){ 
  #k=5
  k<- 5000  
  n<- nrow(Data) 
  d<- (ncol(Data)-1)/2
  
  ind<- ncol(Data)-1 
  lowerb<- NULL 
  upperb<- NULL 
  for (i in 1:ind){ 
    if (i%%2!=0) {lowerb<- append(lowerb,i)} 
    else {upperb<- append(upperb,i)} 
    
  }
  Data1<- as.matrix(Data) 
  
  m<- array(rep(0,k*d*n),c(k,d,n))  
  #print(m)
  
  for(i in 1:n){ 
    m[,,i]<- rtmvnorm(k, mean=mu,
                      sigma=sigma, lower=Data1[i,c(lowerb)], upper=Data1[i,c(upperb)],
                      algorithm="gibbs") 
    #  print(m)
  }
  #print(m)
  
  #mm<- matrix(rep(0,2*n),ncol=2)
  #print(mm)
  #for(i in 1:n){
  #mm[i,]<- apply(m[,,i],2,mean)
  #}
  #print(mm)
  return(m)
  
}


###E-step: Mean Calculation###
######################################################
Mu<- function(Data,Y){ 
  n<- nrow(Data) 
  F<- ncol(Data) 
  T<- sum(Data[,F]) 
  #print(T)
  d<- (ncol(Data)-1)/2 
  mm<- matrix(rep(0,d*n),ncol=d)  
  #print(mm)
  for(i in 1:n){ 
    mm[i,]<- apply(Y[,,i],2,mean)  
  }
  #print(mm)
  Ym<- mm*Data[,F] 
  #print(Ym)
  sm<- rep(0,d) 
  for(i in 1:d){ 
    sm[i]<- sum(Ym[,i]) 
  }
  #print(sm)
  mu<- (sm/T) 
  return(mu) 
}


### E-Step: Sigma Calculation###
#######################################################
#######################################################################
### SIGMA###
Sigma1<- function(Data,Y,M){## Using the above two functions, in this function to find the covariance matrix, 
  ### The arguments of this function are: Data= the data set we are working with, Y= the generated (simulated) samples 
  ### from the function ZMCEM1, M= the mean vector from the Mu function above
  n<- nrow(Data) ## set n equals number of row of data (equals number of rectangles/surfaces) 
  F<- ncol(Data) ## number of column of dataset
  total<- sum(Data[,F]) ## Total sum of observations (Data)
  #print(Y)
  #k=5
  k=5000 ## number of generated samples
  d<- (ncol(Data)-1)/2 ## To find the number of dimension of data 
  myarray<- array(rep(0,d*d*k*n),c(d,d,k,n)) ## create and initialize an array of zero's with the dims: d,d,k,n to create 
  ### d*d covariance matrices over the generated samples and the rectangles/surfaces
  for(i in 1:n){ ## use for loop to implement over all rows (rectangles/surfaces)
    for(j in 1:k){ ## use another for loop to implement over all generated samples
      myarray[,,j,i]<- (Y[j,,i]-M)%*%t(Y[j,,i]-M) ## over the array of generated samples, for each row, and each generated 
      ### samples, we have to find the difference between the generated samples and the mean vectors and square it by 
      ### multiply it to its transpose
    }
    
  }
  #print(myarray)
  #m<- apply(myarray[,,,1],c(1,2),mean)
  #print(m)
  Newarray<- array(rep(0,d*d*n),c(d,d,n)) ## create and initialize an array of zero's for the covariances over the 
  ### rectangles/surfaces
  #print(Newarray)
  for(i in 1:n){ ## use for loop to do the function over all rectangles/surfaces
    Newarray[,,i]<- apply(myarray[,,,i],c(1,2),mean) ## take mean over the first and second dimensions of array from last
    ### step (myarray) and assign the results for each rectangles of Newarray  
  }
  #print(Newarray)
  
  NNarray<- array(rep(0,d*d*n),c(d,d,n))  ## create and initialize an array of zero's for the covariances over the 
  ### rectangles/surfaces
  for(i in 1:n){ ## use for loop to do the function over all rectangles/surfaces
    NNarray[,,i]<- Newarray[,,i]*Data[i,F] ## Multiply the covariance matrices over the rectangles/surfaces
    ### from above (Newaaray) by the frequency of each surface
  }
  #print(NNarray)
  S1<- apply(NNarray,c(1,2),sum) ## sum over all the rectangles/surfaces
  #print(S1)
  mcov<- S1/total ## devide the result by total number of observations to find the covariance matrix
  #print(mcov)
  #print(myarray)
  #m1<- array(rep(0,2*2*n),c(2,2,n))
  #m2<- array(rep(0,2*2*n),c(2,2,n))
  #for(i in 1:n){
  # m1[,,i]<- (Y[i,]-M)%*%t(Y[i,]-M)
  #m2[,,i]<- m1[,,i]*Data[i,5]
  
  #  }
  #print(m1)
  #print(m2)
  # m3<- apply(m2,1:2,sum)
  #print(m3)
  #mcov<- m3/total
  #print(mcov)
  return(mcov) ## return the covariance matrix
}




### M-Step###
MCEMBiv<- function(data,mu_init,sigma_init,maxit=1000,tol1=1e-3,tol2=1e-3){
  flag<- 0 
  M_cur<- mu_init
  S_cur<- sigma_init
  #theta_cur<- c(Mu_cur,S_cur)
  iter<- rep(0,maxit) 
  #Svec<- rep(0,maxit)
  #Mvec<- rep(0,maxit)
  for (i in 1:maxit){ 
    #print(paste("Iteration number=", i))
    #cur<- c(Mu_cur,S_cur)
    MU_cur<- M_cur
    SS_cur<- S_cur 
    
    MuNew<- Mu(Data=mydat,Y=ZMCEM1(Data=mydat,mu=MU_cur,sigma=SS_cur))
    
    SigmaNew<- Sigma1(Data=mydat,Y=ZMCEM1(Data=mydat,mu=MU_cur,sigma=SS_cur),
                      M=Mu(Data=mydat,Y=ZMCEM1(Data=mydat,mu=MU_cur,sigma=SS_cur) ) ) 
    
    Mu_new<- MuNew 
    Sigma_new<- SigmaNew 
    
    diff1<- MU_cur-Mu_new 
    diff2<- SS_cur-Sigma_new
    
    D1<- abs(mean(diff1)) 
    D2<- abs(mean(diff2))  
    
    
    #print(D1)
    #print(D2)
    if(D1<tol1 & D2<tol2){flag<- 1 ;break} 
    M_cur<- Mu_new 
    
    S_cur<- Sigma_new 
  }
  
  if(!flag) warning("Didn't Converge \n") 
  #updatedest=list(paste("updated Mu of the EM=" , M_cur),paste("Updated Sigma of the EM=" ,S_cur),
  #  paste("iteration number=",i))
  
  
  updatedest=list(M_cur,(S_cur)) 
  return(updatedest) 
  
}
#####################################################################################################################
#####################################################################################################################
### Standard Erroe of EM Estimates####
Std_Mu_EM<- function(Data,mu,sigma){ 
  n<- nrow(Data) 
  d<- (ncol(Data)-1)/2  
  Freq<- Data[,ncol(Data)]
  
  ind<- ncol(Data)-1 
  lowerb<- NULL 
  upperb<- NULL 
  for (i in 1:ind){ 
    if (i%%2!=0) {lowerb<- append(lowerb,i)}
    else {upperb<- append(upperb,i)} 
    
  }
  Data1<- as.matrix(Data) 
  
  EX<- matrix(rep(0,n*d),ncol=d,nrow=n) 
  CovX<- array(rep(0,n*d*d),c(d,d,n))
  for(i in 1:nrow(Data)){ 
    mnts<- mtmvnorm(mean=mu, sigma=sigma, 
                    lower=Data1[i,c(lowerb)],
                    upper=Data1[i,c(upperb)]) 
    EX[i,]<- mnts$tmean 
    CovX[,,i]<- mnts$tvar
  }
  EX_minus_MU <- matrix(rep(0,n*d),ncol=d,nrow=n) 
  EX_minus_MU<- EX-mu

  dS_mu <- matrix(rep(0,n*d),ncol=d,nrow=n) 
  
  for(i in 1:nrow(dS_mu)){
    dS_mu[i,]<- (EX_minus_MU[i,])%*%solve(sigma)
    
  }
  
  S1_S1T<- array(rep(0,d*d*n),c(d,d,n))  
  for (i in 1:n){
    S1_S1T[,,i]<- dS_mu[i,]%*%t(dS_mu[i,])
  }
  
  Inf_j_m<- array(rep(0,d*d*n),c(d,d,n))
  for(i in 1:n){
    Inf_j_m[,,i]<- S1_S1T[,,i]*Data[i,5]
  }
  Inf_mat<- apply(Inf_j_m,c(1,2),sum)
  cov_Mu_EM<- solve(Inf_mat)
  Se_mu_EM<- c(sqrt(diag(cov_Mu_EM)))
  names(Se_mu_EM)<- c("std_mu_x1","std_mu_x2")
  
  

  return(Se_mu_EM)
  
  
}
#######################################################################################################################
#### Standard Errors of MCEM Estimates###
set.seed(9876)
library(tmvtnorm) 
library(MASS) 
#######################################################################
STD_MU_MCEM_Biv<- function(Data,mu,sigma){
  k<- 5000  
  n<- nrow(Data) 
  d<- (ncol(Data)-1)/2
  
  ind<- ncol(Data)-1 
  lowerb<- NULL 
  upperb<- NULL 
  for (i in 1:ind){ 
    if (i%%2!=0) {lowerb<- append(lowerb,i)} 
    else {upperb<- append(upperb,i)} 
    
  }
  Data1<- as.matrix(Data) 
  
  zsim<- array(rep(0,k*d*n),c(k,d,n))  
  #print(m)
  
  for(i in 1:n){ 
    zsim[,,i]<- rtmvnorm(k, mean=mu,
                         sigma=sigma, lower=Data1[i,c(lowerb)], upper=Data1[i,c(upperb)],
                         algorithm="gibbs") 
    #  print(m)
  }
  d2<- solve(sigma)
  d1<- array(rep(0,k*d*n),c(k,d,n)) 
  for(i in 1:n){
    d1[,,i]<- (zsim[,,i]-mu)%*%solve(sigma)
  }

  
  md1<- matrix(rep(0,d*n),ncol=d)
  for(i in 1:n){
    md1[i,]<- apply(d1[,,i],2,mean)
  }
  diff1<- array(rep(0,k*d*n),c(k,d,n)) 
  for(i in 1:n){
    diff1[,,i]<- d1[,,i]-md1[i,] 
  }
  der1<- matrix(rep(0,d*n),ncol=d)
  for(i in 1:n){
    der1[i,]<- apply(diff1[,,i],2,mean)
  }

  diff2<- array(rep(0,d*d*n),c(d,d,n))
  for(i in 1:n){
    diff2[,,i]<- der1[i,]%*%t(der1[i,])
  }
  sj<- array(rep(0,d*d*n),c(d,d,n))
  for(i in 1:n){
    sj[,,i]<- (d2+diff2[,,i])*Data[i,ncol(Data)]
  }
  IM<- apply(sj,c(1,2),sum)
  Inf_mat<- solve(IM)
  se_MCEM_MU<- diag(sqrt(Inf_mat))
  names(se_MCEM_MU)<- c("se_mu_x1","se_mu_x2")
  return(se_MCEM_MU)
}


#######################################################################################################################

MLEGalton<- nlm(Billik,Data=DataBiv,theta<- c(67,67,3.2,6.2,0.5),hessian=TRUE)$estimate
  
MLEGalton


#t1<- Null
#t1<- system.time(nlm(Billik,Data=DataBiv,theta<- c(67,67,3.2,6.2,0.5),hessian=TRUE)$estimate)
#t1

MLEGalton
#[1] 68.300475 68.098651  3.243895  6.513746  0.470162
########################################################################################################################
mu1 <- c(MeP,MeC)
mu1
sigma1 <- matrix(c(VP, rho*SDP*SDC, rho*SDP*SDC, VC), 2, 2)
sigma1

mydat<- DataBiv

EM_Galton<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-3))
EM_Galton

mu_em<- c(EM_Galton[1],EM_Galton[2])
mu_em
sigma_em<- matrix(c(EM_Galton[3],EM_Galton[4],EM_Galton[5],EM_Galton[6]),ncol=2)
sigma_em
se_EM_Galton<- Std_Mu_EM(Data=mydat,mu=mu_em,sigma=sigma_em)   
se_EM_Galton

#std_mu_x1  std_mu_x2 
#0.05965676 0.08425964
########################################################################################################################
MCEM_Galton<- unlist(MCEMBiv(data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-3,tol2=1e-3))
MCEM_Galton

mu_mcem<- c(MCEM_Galton[1],MCEM_Galton[2])
mu_mcem
sigma_mcem<- matrix(c(MCEM_Galton[3],MCEM_Galton[4],MCEM_Galton[5],MCEM_Galton[6]),ncol=2)
sigma_mcem

se_MCEM<- STD_MU_MCEM_Biv(Data=mydat,mu=mu_mcem,sigma=sigma_mcem)   
se_MCEM

#se_mu_x1   se_mu_x2 
#0.05807384 0.07091751
