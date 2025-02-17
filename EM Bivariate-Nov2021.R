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
#############################################################################
##########################################################################
### Check if the function works on data ###

###INITIAL VALUES for dim=2 & dim =3####
## d=2##
#sqrt(3.2*6.2)*0.5
mu1 <- c(67,67)
sigma1 <- matrix(c(3.2, 2.227106, 2.227106, 6.2), 2, 2)
## d=3##
#mu2<- c(67,67,67)
#sigma2<- matrix(c(3.1,2.2,1.6,2.2,4.3,2.45,1.6,2.45,4.8),3,3)
#sigma2
################################################################
EM50b10<- matrix(rep(0*6*500),ncol=6,nrow=500) 
colnames(EM50b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")
for(i in 1:nrow(EM50b10)){ 
  mydat<- simulateddata50b10[,,i]
  EM50b10[i,]<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-3))
  
}


EM50b10
#######################################################################################################
EM100b10<- matrix(rep(0*6*500),ncol=6,nrow=500) 
colnames(EM100b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")
for(i in 1:nrow(EM100b10)){ 
  mydat<- simulateddata100b10[,,i]
  EM100b10[i,]<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-3))
  
}

EM100b10
#########################################################################################################
#read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/Simulations1/simulateddata300b10.csv")


EM300b10<- matrix(rep(0*6*500),ncol=6,nrow=500) 
colnames(EM300b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")
for(i in 1:nrow(EM300b10)){ 
  mydat<- simulateddata300b10[,,i]
  EM300b10[i,]<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-3))
  
}

EM300b10
###########################################################################################################
EM600b10<- matrix(rep(0*6*500),ncol=6,nrow=500) 
colnames(EM600b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")
for(i in 1:nrow(EM600b10)){ 
  mydat<- simulateddata600b10[,,i]
  EM600b10[i,]<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-3))
  
}

EM600b10

###########################################################################################################################
EM1000b10<- matrix(rep(0*6*500),ncol=6,nrow=500) 
colnames(EM1000b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")
for(i in 1:nrow(EM1000b10)){ 
  mydat<- simulateddata1000b10[,,i]
  EM1000b10[i,]<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-3))
  
}

EM1000b10

######################################################################################################################
###########################################################################################################################
EM1500b10<- matrix(rep(0*6*500),ncol=6,nrow=500) 
colnames(EM1500b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")
for(i in 1:nrow(EM1500b10)){ 
  mydat<- simulateddata1500b10[,,i]
  EM1500b10[i,]<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-3))
  
}

EM1500b10

dataEM_n1500_bin10<- write.csv(EM1500b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM1500b10.csv")

##############################################################################################################################

#EM1000b10<- matrix(rep(0*6*10),ncol=6,nrow=10) 
#colnames(EM1000b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")
#for(i in 1:nrow(EM1000b10)){ 
 # mydat<- simulateddata1000b10[,,i]
  #EM1000b10[i,]<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-5))
  
#}

#EM1000b10

#rm("EM1000b10")
#rm("i")
#EM1000b10
#dim(simulateddata1000b10)

#EM1000b10<- matrix(rep(0*6*500),ncol=6,nrow=500) 

#i=1
#mydat<- simulateddata1000b10[,,i]

#EM1000b10[i,]<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-3))

#EM1000b10

#EM1000b10<- matrix(rep(0*5*500),ncol=6,nrow=500) 
#colnames(EM1000b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")

#tt1<- NULL
#for(i in 1:nrow(EM1000b10)){ 
 # mydat<- simulateddata1000b10[,,i]
  #tt1<- system.time(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-5))
  
#}

#tt1
#tt1
#user  system elapsed 
#5.56    0.01    5.58 

#tt
#68+(0.008*1.96)
#68-(0.008*1.96)

#########################################################################################################################
dataEM_n50_bin10<- write.csv(EM50b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM50b10.csv")
dataEM_n100_bin10<- write.csv(EM100b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM100b10.csv")
dataEM_n300_bin10<- write.csv(EM300b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM300b10.csv")
dataEM_n600_bin10<- write.csv(EM600b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM600b10.csv")
dataEM_n1000_bin10<- write.csv(EM1000b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM1000b10.csv")







