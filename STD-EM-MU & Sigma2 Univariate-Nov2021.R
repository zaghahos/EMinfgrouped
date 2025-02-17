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
  
  
  #dinom<- NULL
  #for(i in 1:length(bl)){ 
   # dinom[i]<- (pnorm(bstar[i])-pnorm(astar[i]))
#    if(dinom[i]==0) { Wj[i]<- thetaupd[1]-thetaupd[2]*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001)}
 #   else {Wj[i]<- thetaupd[1]-thetaupd[2]*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i])}
  #}
  dinom1<- NULL
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  for(i in 1:length(bl)){ 
    dinom1[i]<- (pnorm(bstar[i])-pnorm(astar[i]))
    #print(dinom[i])
    if(dinom1[i]==0) { WjN[i]<- (-thetaupd[2])*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001)}
    else {WjN[i]<- (-thetaupd[2])*((dnorm(bstar[i])-dnorm(astar[i]))/dinom1[i])}
  }
 
  Ej2<- rep(0,length(bl)) 
  
  #astar<- rep(0,length(bl)) 
  #bstar<- rep(0,length(bl)) 
  #for(i in 1:length(bl)){
   # bstar[i]<- (bu[i]-thetaupd[1])/thetaupd[2] 
    #astar[i]<- (bl[i]-thetaupd[1])/thetaupd[2] 
    
  #}
  dinom<- NULL
  
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  
  for(i in 1:length(bl)){ 
    dinom[i]<- (pnorm(bstar[i])-pnorm(astar[i])) 
    if(dinom[i]==0) {Ej2[i]<- (bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/0.0001}
    
    
    else{Ej2[i]<- (bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/dinom[i]}
  } 
  #print(WjN)
  #print(Ej2)
  SMU<- WjN/(thetaupd[2]^2)
  SVAR<- (-1/(2*(thetaupd[2]^2)))*Ej2
  Sj1<- rbind(SMU,SVAR)
  #print(SMU)
  #print(SVAR)
  #print(Sj1)
  #Sj<- rbind(WjN/(thetaupd[2]^2),Ej2/(4*thetaupd[2]^4))
  IME<- array(rep(0,2*2*(nrow(Data))),c(2,2,nrow(Data))) 
  for(i in 1:nrow(Data)){
    IME[,,i]<- Sj1[,i]%*%t(Sj1[,i])*Data[i,3]
  }
  
  #print(IME)
  InfME<- apply(IME,c(1,2),sum)
  
  var_EM_est<- solve(InfME)
  #print(var_EM_est)
  std_EM_est<- sqrt(diag(var_EM_est))
  #print(std_EM_est)
  #print(dim(IM))
  
  
  #Mstdj<- Freq*(diag(WjN%*%t(WjN)))
  #Mstdj1<- (1/thetaupd[2]^4)*(diag(WjN%*%t(WjN)))
  #Mstdj2<- Freq*Mstdj1
  #Mstdj3<- sqrt(1/sum(Mstdj2))
  #print(Mstdj3)
  #print(Mstdj)
  #Mstd<- (sum(Mstdj))/(thetaupd[2]^4)
  #results<- list(paste("std for Mu=" , 1/Mstd),Wj,WjN)
  
  
  return(std_EM_est) 
  
}
#i=1
#Muvarstd(thetaupd=c(output1000m15EM[i,1],sqrt(output1000m15EM[i,2])),
 #     bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i], Data=simdata1000m15[,,i])

#sqrt(0.005668293)
#simdata1000m15[,,1]

#sqrt(1/2)
#1/(sqrt(2))
################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
### N=1000###
se_EM_1000<- matrix(rep(0,2*nrow(output1000m15EM)),ncol=2)

for(i in 1:nrow(output1000m15EM)){
  se_EM_1000[i,]<- Muvarstd(thetaupd=c(output1000m15EM[i,1],sqrt(output1000m15EM[i,2])),
                   bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i],Data=simdata1000m15[,,i])
  
}

se_EM_1000
#se_EM_1000[1:5,]

#se_MCEM1000m15[1:5,]

CI_mu_1000<- matrix(rep(0,2*nrow(output1000m15EM)),nrow=nrow(output1000m15EM),ncol=2)
colnames(CI_mu_1000)<- c("Lower bound","Upper bound")

for(i in 1:nrow(output1000m15EM)){
  CI_mu_1000[i,1]<- output1000m15EM[i,1]-(se_EM_1000[i,1]*qnorm(0.975))
  CI_mu_1000[i,2]<- output1000m15EM[i,1]+(se_EM_1000[i,1]*qnorm(0.975))
  
}

#CI_mu_1000

CI_sigma2_1000<- matrix(rep(0,2*nrow(output1000m15EM)),nrow=nrow(output1000m15EM),ncol=2)
colnames(CI_sigma2_1000)<- c("Lower bound","Upper bound")


for(i in 1:nrow(output1000m15EM)){
  CI_sigma2_1000[i,1]<- output1000m15EM[i,2]-(se_EM_1000[i,2]*qnorm(0.975))
  CI_sigma2_1000[i,2]<- output1000m15EM[i,2]+(se_EM_1000[i,2]*qnorm(0.975))
  
}
#CI_sigma2_1000

true_mu<- 68
true_sigma2<- (2.5)^2
#true_sigma2


Empirical_conf_mu_1000<- 0
Empirical_conf_sigma2_1000<- 0

for(i in 1:nrow(output1000m15EM)){
  if(true_mu>=CI_mu_1000[i,1] & true_mu<=CI_mu_1000[i,2]) Empirical_conf_mu_1000<- Empirical_conf_mu_1000+1
  
  if(true_sigma2>=CI_sigma2_1000[i,1] & true_sigma2<=CI_sigma2_1000[i,2]) Empirical_conf_sigma2_1000<- Empirical_conf_sigma2_1000+1
}


Empirical_conf_mu_1000/500
Empirical_conf_sigma2_1000/500
############################################################################################################################
### N=50###
se_EM_50<- matrix(rep(0,2*nrow(output50m15EM)),ncol=2)

for(i in 1:nrow(output50m15EM)){
  se_EM_50[i,]<- Muvarstd(thetaupd=c(output50m15EM[i,1],sqrt(output50m15EM[i,2])),
                            bl=simdata50m15[,1,i],bu=simdata50m15[,2,i],Freq=simdata50m15[,3,i],Data=simdata50m15[,,i])
  
}
se_EM_50

CI_mu_50<- matrix(rep(0,2*nrow(output50m15EM)),nrow=nrow(output50m15EM),ncol=2)
colnames(CI_mu_50)<- c("Lower bound","Upper bound")

for(i in 1:nrow(output50m15EM)){
  CI_mu_50[i,1]<- output50m15EM[i,1]-(se_EM_50[i,1]*qnorm(0.975))
  CI_mu_50[i,2]<- output50m15EM[i,1]+(se_EM_50[i,1]*qnorm(0.975))
  
}

#CI_mu_100

CI_sigma2_50<- matrix(rep(0,2*nrow(output50m15EM)),nrow=nrow(output50m15EM),ncol=2)
colnames(CI_sigma2_50)<- c("Lower bound","Upper bound")


for(i in 1:nrow(output50m15EM)){
  CI_sigma2_50[i,1]<- output50m15EM[i,2]-(se_EM_50[i,2]*qnorm(0.975))
  CI_sigma2_50[i,2]<- output50m15EM[i,2]+(se_EM_50[i,2]*qnorm(0.975))
  
}
#CI_sigma2_100

true_mu<- 68
true_sigma2<- (2.5)^2
#true_sigma2


Empirical_conf_mu_50<- 0
Empirical_conf_sigma2_50<- 0

for(i in 1:nrow(output50m15EM)){
  if(true_mu>=CI_mu_50[i,1] & true_mu<=CI_mu_50[i,2]) Empirical_conf_mu_50<- Empirical_conf_mu_50+1
  
  if(true_sigma2>=CI_sigma2_50[i,1] & true_sigma2<=CI_sigma2_50[i,2]) 
    Empirical_conf_sigma2_50<- Empirical_conf_sigma2_50+1
}


Empirical_conf_mu_50/500
Empirical_conf_sigma2_50/500
########################################################################################################################
### N=100##
se_EM_100<- matrix(rep(0,2*nrow(output100m15EM)),ncol=2)

for(i in 1:nrow(output100m15EM)){
  se_EM_100[i,]<- Muvarstd(thetaupd=c(output100m15EM[i,1],sqrt(output100m15EM[i,2])),
                            bl=simdata100m15[,1,i],bu=simdata100m15[,2,i],Freq=simdata100m15[,3,i],Data=simdata100m15[,,i])
  
}
se_EM_100

CI_mu_100<- matrix(rep(0,2*nrow(output100m15EM)),nrow=nrow(output100m15EM),ncol=2)
colnames(CI_mu_100)<- c("Lower bound","Upper bound")

for(i in 1:nrow(output100m15EM)){
  CI_mu_100[i,1]<- output100m15EM[i,1]-(se_EM_100[i,1]*qnorm(0.975))
  CI_mu_100[i,2]<- output100m15EM[i,1]+(se_EM_100[i,1]*qnorm(0.975))
  
}

#CI_mu_100

CI_sigma2_100<- matrix(rep(0,2*nrow(output100m15EM)),nrow=nrow(output100m15EM),ncol=2)
colnames(CI_sigma2_100)<- c("Lower bound","Upper bound")


for(i in 1:nrow(output100m15EM)){
  CI_sigma2_100[i,1]<- output100m15EM[i,2]-(se_EM_100[i,2]*qnorm(0.975))
  CI_sigma2_100[i,2]<- output100m15EM[i,2]+(se_EM_100[i,2]*qnorm(0.975))
  
}
#CI_sigma2_100

true_mu<- 68
true_sigma2<- (2.5)^2
#true_sigma2


Empirical_conf_mu_100<- 0
Empirical_conf_sigma2_100<- 0

for(i in 1:nrow(output100m15EM)){
  if(true_mu>=CI_mu_100[i,1] & true_mu<=CI_mu_100[i,2]) Empirical_conf_mu_100<- Empirical_conf_mu_100+1
  
  if(true_sigma2>=CI_sigma2_100[i,1] & true_sigma2<=CI_sigma2_100[i,2]) 
    Empirical_conf_sigma2_100<- Empirical_conf_sigma2_100+1
}


Empirical_conf_mu_100/500
Empirical_conf_sigma2_100/500
###################################################################################################################
### N=300###
se_EM_300<- matrix(rep(0,2*nrow(output300m15EM)),ncol=2)

for(i in 1:nrow(output300m15EM)){
  se_EM_300[i,]<- Muvarstd(thetaupd=c(output300m15EM[i,1],sqrt(output300m15EM[i,2])),
                           bl=simdata300m15[,1,i],bu=simdata300m15[,2,i],Freq=simdata300m15[,3,i],Data=simdata300m15[,,i])
  
}
se_EM_300

CI_mu_300<- matrix(rep(0,2*nrow(output300m15EM)),nrow=nrow(output300m15EM),ncol=2)
colnames(CI_mu_300)<- c("Lower bound","Upper bound")

for(i in 1:nrow(output300m15EM)){
  CI_mu_300[i,1]<- output300m15EM[i,1]-(se_EM_300[i,1]*qnorm(0.975))
  CI_mu_300[i,2]<- output300m15EM[i,1]+(se_EM_300[i,1]*qnorm(0.975))
  
}

#CI_mu_300

CI_sigma2_300<- matrix(rep(0,2*nrow(output300m15EM)),nrow=nrow(output300m15EM),ncol=2)
colnames(CI_sigma2_300)<- c("Lower bound","Upper bound")


for(i in 1:nrow(output300m15EM)){
  CI_sigma2_300[i,1]<- output300m15EM[i,2]-(se_EM_300[i,2]*qnorm(0.975))
  CI_sigma2_300[i,2]<- output300m15EM[i,2]+(se_EM_300[i,2]*qnorm(0.975))
  
}
#CI_sigma2_100

true_mu<- 68
true_sigma2<- (2.5)^2
#true_sigma2


Empirical_conf_mu_300<- 0
Empirical_conf_sigma2_300<- 0

for(i in 1:nrow(output300m15EM)){
  if(true_mu>=CI_mu_300[i,1] & true_mu<=CI_mu_300[i,2]) Empirical_conf_mu_300<- Empirical_conf_mu_300+1
  
  if(true_sigma2>=CI_sigma2_300[i,1] & true_sigma2<=CI_sigma2_300[i,2]) 
    Empirical_conf_sigma2_300<- Empirical_conf_sigma2_300+1
}


Empirical_conf_mu_300/500
Empirical_conf_sigma2_300/500
##################################################################################################################
###N=600###
se_EM_600<- matrix(rep(0,2*nrow(output600m15EM)),ncol=2)

for(i in 1:nrow(output600m15EM)){
  se_EM_600[i,]<- Muvarstd(thetaupd=c(output600m15EM[i,1],sqrt(output600m15EM[i,2])),
                           bl=simdata600m15[,1,i],bu=simdata600m15[,2,i],Freq=simdata600m15[,3,i],Data=simdata600m15[,,i])
  
}
se_EM_600

CI_mu_600<- matrix(rep(0,2*nrow(output600m15EM)),nrow=nrow(output600m15EM),ncol=2)
colnames(CI_mu_600)<- c("Lower bound","Upper bound")

for(i in 1:nrow(output600m15EM)){
  CI_mu_600[i,1]<- output600m15EM[i,1]-(se_EM_600[i,1]*qnorm(0.975))
  CI_mu_600[i,2]<- output600m15EM[i,1]+(se_EM_600[i,1]*qnorm(0.975))
  
}

#CI_mu_600

CI_sigma2_600<- matrix(rep(0,2*nrow(output600m15EM)),nrow=nrow(output600m15EM),ncol=2)
colnames(CI_sigma2_600)<- c("Lower bound","Upper bound")


for(i in 1:nrow(output600m15EM)){
  CI_sigma2_600[i,1]<- output600m15EM[i,2]-(se_EM_600[i,2]*qnorm(0.975))
  CI_sigma2_600[i,2]<- output600m15EM[i,2]+(se_EM_600[i,2]*qnorm(0.975))
  
}
#CI_sigma2_100

true_mu<- 68
true_sigma2<- (2.5)^2
#true_sigma2


Empirical_conf_mu_600<- 0
Empirical_conf_sigma2_600<- 0

for(i in 1:nrow(output600m15EM)){
  if(true_mu>=CI_mu_600[i,1] & true_mu<=CI_mu_600[i,2]) Empirical_conf_mu_600<- Empirical_conf_mu_600+1
  
  if(true_sigma2>=CI_sigma2_600[i,1] & true_sigma2<=CI_sigma2_600[i,2]) 
    Empirical_conf_sigma2_600<- Empirical_conf_sigma2_600+1
}


Empirical_conf_mu_600/500
Empirical_conf_sigma2_600/500



################################################################################################################
################################################################################################################

EMP_Conf_Mu_prop<- c(Empirical_conf_mu_50/500,Empirical_conf_mu_100/500,Empirical_conf_mu_300/500,
                     Empirical_conf_mu_600/500,Empirical_conf_mu_1000/500)


EMP_Conf_sigma2_prop<- c(Empirical_conf_sigma2_50/500,Empirical_conf_sigma2_100/500,Empirical_conf_sigma2_300/500,
                         Empirical_conf_sigma2_600/500,Empirical_conf_sigma2_1000/500)

print(EMP_Conf_Mu_prop*100)
print(EMP_Conf_sigma2_prop*100)

##################################
n<- c(50,100,300,600,1000)
EMP_CI_MU_Prop<- EMP_Conf_Mu_prop*100
EMP_CI_Sigma2_Prop<- EMP_Conf_sigma2_prop*100
EMP_CI_MU_Prop
EMP_CI_Sigma2_Prop

sd_MU<- c(sd(output50m15EM[,1]),sd(output100m15EM[,1]),sd(output300m15EM[,1]),
          sd(output600m15EM[,1]),sd(output1000m15EM[,1]))
sd_MU
ave_MU<- c(mean(output50m15EM[,1]),mean(output100m15EM[,1]),mean(output300m15EM[,1]),
          mean(output600m15EM[,1]),mean(output1000m15EM[,1]))
ave_MU

Se_MU_hat<- c(mean(se_EM_50[,1]),mean(se_EM_100[,1]),mean(se_EM_300[,1]),mean(se_EM_600[,1]),mean(se_EM_1000[,1]))
Se_MU_hat

sd_var<- c(sd(output50m15EM[,2]),sd(output100m15EM[,2]),sd(output300m15EM[,2]),
           sd(output600m15EM[,2]),sd(output1000m15EM[,2]))
sd_var

Ave_var<- c(mean(output50m15EM[,2]),mean(output100m15EM[,2]),mean(output300m15EM[,2]),
           mean(output600m15EM[,2]),mean(output1000m15EM[,2]))
Ave_var

Se_Sigma2_hat<-  c(mean(se_EM_50[,2]),mean(se_EM_100[,2]),mean(se_EM_300[,2]),mean(se_EM_600[,2]),mean(se_EM_1000[,2]))
Se_Sigma2_hat
########################################################################################################################


##########################################################################################################################################
#######################################################################################################################
std_Mu_estimates<- cbind(n,ave_MU,sd_MU,Se_MU_hat,EMP_CI_MU_Prop)
std_Mu_estimates

std_Sigma2_estimates<- cbind(n,Ave_var,sd_var,Se_Sigma2_hat,EMP_CI_Sigma2_Prop)
std_Sigma2_estimates

library(xtable)
xtable(std_Mu_estimates,digits=6)
xtable(std_Sigma2_estimates,digits=6)



############################################################################################################################
##############################################################################################################################
##############################################################################################################################
c(outParent[1],sqrt(outParent[2]))
c(outChildren[1],sqrt(outChildren[2]))

se_EM_Parent<- Muvarstd(thetaupd=c(outParent[1],sqrt(outParent[2])),
                     bl=Galton2$TL,bu=Galton2$TU,Freq=Galton2$Freq,Data=Galton2)
se_EM_Parent




se_EM_Children<- Muvarstd(thetaupd=c(outChildren[1],sqrt(outChildren[2])),
                       bl=Galton3$TL,bu=Galton3$TU,Freq=Galton3$Freq,Data=Galton3)
se_EM_Children



se_mu_Parent
se_mu_Children

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
