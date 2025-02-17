#thetaupd<- c(output1000m15EM[1,1],sqrt(output1000m15EM[1,2]))
#thetaupd
#output1000m15EM
#######################################################################################
Mustd<- function(thetaupd,bl,bu,Freq,Data){ 
  Wj<- rep(0,length(bl)) 
  WjN<- rep(0,length(bl)) 
  
  astar<- rep(0,length(bl)) 
  bstar<- rep(0,length(bl)) 
  Mstdj<- rep(0,length(bl))
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-thetaupd[1])/thetaupd[2] 
    astar[i]<- (bl[i]-thetaupd[1])/thetaupd[2] 
    
  }
  dinom<- NULL
  for(i in 1:length(bl)){ 
    #print(i)
    #print(pnorm(bstar[i]))
    #print(pnorm(astar[i]))
    dinom[i]<- (pnorm(bstar[i])-pnorm(astar[i]))
    #print(dinom[i])
    if(dinom[i]==0) { Wj[i]<- thetaupd[1]-thetaupd[2]*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001)}
    else {Wj[i]<- thetaupd[1]-thetaupd[2]*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i])}
  }
  dinom1<- NULL
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  for(i in 1:length(bl)){ 
    dinom1[i]<- (pnorm(bstar[i])-pnorm(astar[i]))
    #print(dinom[i])
    if(dinom1[i]==0) { WjN[i]<- (-thetaupd[2])*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001)}
    else {WjN[i]<- (-thetaupd[2])*((dnorm(bstar[i])-dnorm(astar[i]))/dinom1[i])}
  }
  #Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i])     
  #}
  #print(Wj)
  #print(WjN)
  Sj_mu<- WjN/(thetaupd[2]^2)
  IMU<- rep(0,nrow(Data))
  for(i in 1:nrow(Data)){
    IMU[i]<- Sj_mu[i]*Sj_mu[i]*Data[i,3]
  }
  
  var_mu<- 1/(sum(IMU))
  std_mu<- sqrt(var_mu)
  #print(std_mu)
  #print((Wj-thetaupd[1]))
  #Mstdj<- Freq*((Wj-thetaupd[1])^2)
  #print(WjN%*%t(WjN))
  #print(WjN^2)
  #print(diag(WjN%*%t(WjN)))
  
  
  Mstdj<- Freq*(diag(WjN%*%t(WjN)))
  #print(Mstdj)
  #print(Mstdj1)
  #print(Mstdj)
  Mstd<- (sum(Mstdj))/(thetaupd[2]^4)
  #results<- list(paste("std for Mu=" , 1/Mstd),Wj,WjN)
  #return(1/sqrt(Mstd)) 
  return(std_mu)
  
}
#i=1
#Mustd(thetaupd=c(output1000m15EM[i,1],sqrt(output1000m15EM[i,2])),
 #     bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i],Data=simdata1000m15[,,i])
#simdata1000m15[,,1]
#sqrt(0.005668293)

#############################################################################################################
#############################################################################################################
### N=1000###
se_mu<- rep(0,nrow(output1000m15EM))


for(i in 1:nrow(output1000m15EM)){
  se_mu[i]<- Mustd(thetaupd=c(output1000m15EM[i,1],sqrt(output1000m15EM[i,2])),
                   bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i],Data=simdata1000m15[,,i])
  
 # se_sigma2[i]<- VarstdNew(thetaupd=c(output1000m15EM[i,1],sqrt(output1000m15EM[i,2])),
  #                         bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i])
  
  
}



CI_mu<- matrix(rep(0,2*nrow(output1000m15EM)),nrow=nrow(output1000m15EM),ncol=2)
colnames(CI_mu)<- c("Lower bound","Upper bound")

for(i in 1:nrow(output1000m15EM)){
  CI_mu[i,1]<- output1000m15EM[i,1]-(se_mu[i]*qnorm(0.975))
  CI_mu[i,2]<- output1000m15EM[i,1]+(se_mu[i]*qnorm(0.975))
  
}



true_mu<- 68


Empirical_conf_mu<- 0

for(i in 1:nrow(output1000m15EM)){
  if(true_mu>=CI_mu[i,1] & true_mu<=CI_mu[i,2]) Empirical_conf_mu<- Empirical_conf_mu+1
  
 # if(true_sigma2>=CI_sigma2[i,1] & true_sigma2<=CI_sigma2[i,2]) Empirical_conf_sigma2<- Empirical_conf_sigma2+1
}


Empirical_conf_mu/500
#Empirical_conf_sigma2/500
################################################################################################################
################################################################################################################
### N=50###
se_mu50<- rep(0,nrow(output50m15EM))
for(i in 1:nrow(output50m15EM)){
  se_mu50[i]<- Mustd(thetaupd=c(output50m15EM[i,1],sqrt(output50m15EM[i,2])),
                     bl=simdata50m15[,1,i],bu=simdata50m15[,2,i],Freq=simdata50m15[,3,i],Data=simdata50m15[,,i])
  
  #se_sigma2_50[i]<- VarstdNew(thetaupd=c(output50m15EM[i,1],sqrt(output50m15EM[i,2])),
   #                           bl=simdata50m15[,1,i],bu=simdata50m15[,2,i],Freq=simdata50m15[,3,i])
  
  
}

CI_mu50<- matrix(rep(0,2*nrow(output50m15EM)),nrow=nrow(output50m15EM),ncol=2)
colnames(CI_mu50)<- c("Lower bound","Upper bound")

for(i in 1:nrow(output50m15EM)){
  CI_mu50[i,1]<- output50m15EM[i,1]-(se_mu50[i]*qnorm(0.975))
  CI_mu50[i,2]<- output50m15EM[i,1]+(se_mu50[i]*qnorm(0.975))
  
}

true_mu<- 68
Empirical_conf_mu50<- 0

for(i in 1:nrow(output50m15EM)){
  if(true_mu>=CI_mu50[i,1] & true_mu<=CI_mu50[i,2]) Empirical_conf_mu50<- Empirical_conf_mu50+1
  
  #if(true_sigma2>=CI_sigma2_50[i,1] & true_sigma2<=CI_sigma2_50[i,2]) Empirical_conf_sigma2_50<- Empirical_conf_sigma2_50+1
}

Empirical_conf_mu50/500
###########################################################################################################################
###########################################################################################################################
### N=100###
se_mu100<- rep(0,nrow(output100m15EM))


for(i in 1:nrow(output100m15EM)){
  se_mu100[i]<- Mustd(thetaupd=c(output100m15EM[i,1],sqrt(output100m15EM[i,2])),
                      bl=simdata100m15[,1,i],bu=simdata100m15[,2,i],Freq=simdata100m15[,3,i],Data=simdata100m15[,,i])
  
  #se_sigma2_100[i]<- VarstdNew(thetaupd=c(output100m15EM[i,1],sqrt(output100m15EM[i,2])),
   #                            bl=simdata100m15[,1,i],bu=simdata100m15[,2,i],Freq=simdata100m15[,3,i])
  
}


CI_mu100<- matrix(rep(0,2*nrow(output100m15EM)),nrow=nrow(output100m15EM),ncol=2)
colnames(CI_mu100)<- c("Lower bound","Upper bound")

for(i in 1:nrow(output100m15EM)){
  CI_mu100[i,1]<- output100m15EM[i,1]-(se_mu100[i]*qnorm(0.975))
  CI_mu100[i,2]<- output100m15EM[i,1]+(se_mu100[i]*qnorm(0.975))
  
}



true_mu<- 68


Empirical_conf_mu100<- 0

for(i in 1:nrow(output100m15EM)){
  if(true_mu>=CI_mu100[i,1] & true_mu<=CI_mu100[i,2]) Empirical_conf_mu100<- Empirical_conf_mu100+1
  
  #if(true_sigma2>=CI_sigma2_100[i,1] & true_sigma2<=CI_sigma2_100[i,2]) Empirical_conf_sigma2_100<- Empirical_conf_sigma2_100+1
}

Empirical_conf_mu100/500
########################################################################################################################
#######################################################################################################################
### N=300###
se_mu300<- rep(0,nrow(output300m15EM))

for(i in 1:nrow(output300m15EM)){
  se_mu300[i]<- Mustd(thetaupd=c(output300m15EM[i,1],sqrt(output300m15EM[i,2])),
                      bl=simdata300m15[,1,i],bu=simdata300m15[,2,i],Freq=simdata300m15[,3,i],Data=simdata300m15[,,i])
  
  #se_sigma2_300[i]<- VarstdNew(thetaupd=c(output300m15EM[i,1],sqrt(output300m15EM[i,2])),
   #                            bl=simdata300m15[,1,i],bu=simdata300m15[,2,i],Freq=simdata300m15[,3,i])
  
}


CI_mu300<- matrix(rep(0,2*nrow(output300m15EM)),nrow=nrow(output300m15EM),ncol=2)
colnames(CI_mu300)<- c("Lower bound","Upper bound")

for(i in 1:nrow(output300m15EM)){
  CI_mu300[i,1]<- output300m15EM[i,1]-(se_mu300[i]*qnorm(0.975))
  CI_mu300[i,2]<- output300m15EM[i,1]+(se_mu300[i]*qnorm(0.975))
  
}




CI_sigma2_300<- matrix(rep(0,2*nrow(output300m15EM)),nrow=nrow(output300m15EM),ncol=2)
colnames(CI_sigma2_300)<- c("Lower bound","Upper bound")


true_mu<- 68


Empirical_conf_mu300<- 0

for(i in 1:nrow(output300m15EM)){
  if(true_mu>=CI_mu300[i,1] & true_mu<=CI_mu300[i,2]) Empirical_conf_mu300<- Empirical_conf_mu300+1
  
#  if(true_sigma2>=CI_sigma2_300[i,1] & true_sigma2<=CI_sigma2_300[i,2]) Empirical_conf_sigma2_300<- Empirical_conf_sigma2_300+1
}

Empirical_conf_mu300/500
########################################################################################################################
#######################################################################################################################
### N=600###
se_mu600<- rep(0,nrow(output600m15EM))


for(i in 1:nrow(output600m15EM)){
  se_mu600[i]<- Mustd(thetaupd=c(output600m15EM[i,1],sqrt(output600m15EM[i,2])),
                      bl=simdata600m15[,1,i],bu=simdata600m15[,2,i],Freq=simdata600m15[,3,i],Data=simdata600m15[,,i])
  
  #se_sigma2_600[i]<- VarstdNew(thetaupd=c(output600m15EM[i,1],sqrt(output600m15EM[i,2])),
   #                            bl=simdata600m15[,1,i],bu=simdata600m15[,2,i],Freq=simdata600m15[,3,i])
  
}




CI_mu600<- matrix(rep(0,2*nrow(output600m15EM)),nrow=nrow(output600m15EM),ncol=2)
colnames(CI_mu600)<- c("Lower bound","Upper bound")


for(i in 1:nrow(output600m15EM)){
  CI_mu600[i,1]<- output600m15EM[i,1]-(se_mu600[i]*qnorm(0.975))
  CI_mu600[i,2]<- output600m15EM[i,1]+(se_mu600[i]*qnorm(0.975))
  
}



true_mu<- 68


Empirical_conf_mu600<- 0

for(i in 1:nrow(output600m15EM)){
  if(true_mu>=CI_mu600[i,1] & true_mu<=CI_mu600[i,2]) Empirical_conf_mu600<- Empirical_conf_mu600+1
  
  #if(true_sigma2>=CI_sigma2_600[i,1] & true_sigma2<=CI_sigma2_600[i,2]) Empirical_conf_sigma2_600<- Empirical_conf_sigma2_600+1
}

Empirical_conf_mu600/500

EMP_Conf_Mu_prop<- c(Empirical_conf_mu50/500,Empirical_conf_mu100/500,Empirical_conf_mu300/500,
                     Empirical_conf_mu600/500,Empirical_conf_mu/500)



print(EMP_Conf_Mu_prop*100)

#n<- c(50,100,300,600,1000)


#output1000m15EM[,1]  
sd(output1000m15EM[,1])
#se_mu
mean(se_mu)
#########################
#output50m15EM[,1]  
sd(output50m15EM[,1])
#se_mu50
mean(se_mu50)
############################
#output100m15EM[,1]  
sd(output100m15EM[,1])
#se_mu100
mean(se_mu100)
################################
#output300m15EM[,1]  
sd(output300m15EM[,1])
#se_mu300
mean(se_mu300)
##################################
#output600m15EM[,1]  
sd(output600m15EM[,1])
#se_mu600
mean(se_mu600)
##################################
##################################
n<- c(50,100,300,600,1000)
EMP_CI_MU_Prop<- EMP_Conf_Mu_prop*100
EMP_CI_MU_Prop
sd_MU<- c(sd(output50m15EM[,1]),sd(output100m15EM[,1]),sd(output300m15EM[,1]),
          sd(output600m15EM[,1]),sd(output1000m15EM[,1]))

Se_MU_hat<- c(mean(se_mu50),mean(se_mu100),mean(se_mu300),mean(se_mu600),mean(se_mu))
Se_MU_hat
sd_MU
##########################################################################################################################################
#######################################################################################################################
std_Mu_estimates<- cbind(n,Se_MU_hat,sd_MU,EMP_CI_MU_Prop)
std_Mu_estimates
library(xtable)
xtable(std_Mu_estimates,digits=6)
se_mu[1]

############################################################################################################################
##############################################################################################################################
##############################################################################################################################
c(outParent[1],sqrt(outParent[2]))
c(outChildren[1],sqrt(outChildren[2]))

outParent

se_mu_Parent<- Mustd(thetaupd=c(outParent[1],sqrt(outParent[2])),
                     bl=Galton2$TL,bu=Galton2$TU,Freq=Galton2$Freq,Data=Galton2)
se_mu_Parent

se_EM_Parent[1]


se_mu_Children<- Mustd(thetaupd=c(outChildren[1],sqrt(outChildren[2])),
                       bl=Galton3$TL,bu=Galton3$TU,Freq=Galton3$Freq,Galton3)
se_mu_Children
se_EM_Children[1]


lower_MU_P<- outParent[1]-qnorm(0.975)*se_mu_Parent
Upper_MU_P<- outParent[1]+qnorm(0.975)*se_mu_Parent
CI_Mu_P<- c(lower_MU_P,Upper_MU_P)
CI_Mu_P


lower_MU_C<- outChildren[1]-qnorm(0.975)*se_mu_Children
Upper_MU_C<- outChildren[1]+qnorm(0.975)*se_mu_Children
CI_Mu_C<- c(lower_MU_C,Upper_MU_C)
CI_Mu_C



outParent[1]
CI_Mu_P


outChildren[1]
CI_Mu_C

