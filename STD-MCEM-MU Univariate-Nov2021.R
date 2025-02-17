outputMCEM50m15<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out50m15MCEM.csv")
outputMCEM100m15<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out100m15MCEM.csv")
outputMCEM300m15<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out300m15MCEM.csv")
outputMCEM600m15<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out600m15MCEM.csv")
outputMCEM1000m15<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out1000m15MCEM.csv")



outputMCEM50m15<- outputMCEM50m15[,-1]
outputMCEM100m15<- outputMCEM100m15[,-1]
outputMCEM300m15<- outputMCEM300m15[,-1]
outputMCEM600m15<- outputMCEM600m15[,-1]
outputMCEM1000m15<- outputMCEM1000m15[,-1]
outputMCEM1000m15
######################################################################################################################

#rm("i")
#theta<- c(outputMCEM1000m15[i,1],sqrt(outputMCEM1000m15[i,2]))
#theta
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
############################################################################################################
#############################################################################################################
############################################################################################################  
outputMCEM1000m15

theta1000<-matrix(rep(0,2*nrow(output1000m15EM)),ncol=2)
for(i in 1:nrow(output1000m15EM)){
  theta1000[i,]<- c(output1000m15EM[i,1],sqrt(output1000m15EM[i,2]))
}
theta1000


se_MCEM1000m15<- matrix(rep(0,2*nrow(output1000m15EM)),ncol=2)
for(i in 1:nrow(output1000m15EM)){
  se_MCEM1000m15[i,]<-  ZsimMCEM(theta1000[i,],data=simdata1000m15[,,i])
}
se_MCEM1000m15
se_MCEM1000m15[1:5,]
se_EM_1000[1:5,]
trueval<- c(68,(2.5)^2)
#########################################################################################
CI_mu_1000_MCEM<- matrix(rep(0,2*nrow(outputMCEM1000m15)),nrow=nrow(outputMCEM1000m15),ncol=2)
colnames(CI_mu_1000_MCEM)<- c("Lower bound","Upper bound")

for(i in 1:nrow(outputMCEM1000m15)){
  CI_mu_1000_MCEM[i,1]<- outputMCEM1000m15[i,1]-(se_MCEM1000m15[i,1]*qnorm(0.975))
  CI_mu_1000_MCEM[i,2]<- outputMCEM1000m15[i,1]+(se_MCEM1000m15[i,1]*qnorm(0.975))
  
}

CI_mu_1000_MCEM

CI_sigma2_1000_MCEM<- matrix(rep(0,2*nrow(outputMCEM1000m15)),nrow=nrow(outputMCEM1000m15),ncol=2)
colnames(CI_sigma2_1000_MCEM)<- c("Lower bound","Upper bound")


for(i in 1:nrow(outputMCEM1000m15)){
  CI_sigma2_1000_MCEM[i,1]<- outputMCEM1000m15[i,2]-(se_MCEM1000m15[i,2]*qnorm(0.975))
  CI_sigma2_1000_MCEM[i,2]<- outputMCEM1000m15[i,2]+(se_MCEM1000m15[i,2]*qnorm(0.975))
  
}
CI_sigma2_1000_MCEM

true_mu<- 68
true_sigma2<- (2.5)^2
true_sigma2


Empirical_conf_mu_1000_MCEM<- 0
Empirical_conf_sigma2_1000_MCEM<- 0

for(i in 1:nrow(outputMCEM1000m15)){
  if(true_mu>=CI_mu_1000_MCEM[i,1] & true_mu<=CI_mu_1000_MCEM[i,2]) Empirical_conf_mu_1000_MCEM<- Empirical_conf_mu_1000_MCEM+1
  
  if(true_sigma2>=CI_sigma2_1000_MCEM[i,1] & true_sigma2<=CI_sigma2_1000_MCEM[i,2]) 
    Empirical_conf_sigma2_1000_MCEM<- Empirical_conf_sigma2_1000_MCEM+1
}


Empirical_conf_mu_1000_MCEM/500
Empirical_conf_sigma2_1000_MCEM/500
###############################################################################################################################
################################################################################################################################
outputMCEM50m15

theta50<-matrix(rep(0,2*nrow(output50m15EM)),ncol=2)
for(i in 1:nrow(output50m15EM)){
  theta50[i,]<- c(output50m15EM[i,1],sqrt(output50m15EM[i,2]))
}
theta50


se_MCEM50m15<- matrix(rep(0,2*nrow(output50m15EM)),ncol=2)
for(i in 1:nrow(output50m15EM)){
  se_MCEM50m15[i,]<-  ZsimMCEM(theta50[i,],data=simdata50m15[,,i])
}
se_MCEM50m15
se_MCEM50m15[1:5,]
#se_EM_1000[1:5,]
trueval<- c(68,(2.5)^2)
#########################################################################################
CI_mu_50_MCEM<- matrix(rep(0,2*nrow(outputMCEM50m15)),nrow=nrow(outputMCEM50m15),ncol=2)
colnames(CI_mu_50_MCEM)<- c("Lower bound","Upper bound")

for(i in 1:nrow(outputMCEM50m15)){
  CI_mu_50_MCEM[i,1]<- outputMCEM50m15[i,1]-(se_MCEM50m15[i,1]*qnorm(0.975))
  CI_mu_50_MCEM[i,2]<- outputMCEM50m15[i,1]+(se_MCEM50m15[i,1]*qnorm(0.975))
  
}

CI_mu_50_MCEM

CI_sigma2_50_MCEM<- matrix(rep(0,2*nrow(outputMCEM50m15)),nrow=nrow(outputMCEM50m15),ncol=2)
colnames(CI_sigma2_50_MCEM)<- c("Lower bound","Upper bound")


for(i in 1:nrow(outputMCEM50m15)){
  CI_sigma2_50_MCEM[i,1]<- outputMCEM50m15[i,2]-(se_MCEM50m15[i,2]*qnorm(0.975))
  CI_sigma2_50_MCEM[i,2]<- outputMCEM50m15[i,2]+(se_MCEM50m15[i,2]*qnorm(0.975))
  
}
CI_sigma2_50_MCEM

true_mu<- 68
true_sigma2<- (2.5)^2
true_sigma2


Empirical_conf_mu_50_MCEM<- 0
Empirical_conf_sigma2_50_MCEM<- 0

for(i in 1:nrow(outputMCEM50m15)){
  if(true_mu>=CI_mu_50_MCEM[i,1] & true_mu<=CI_mu_50_MCEM[i,2]) Empirical_conf_mu_50_MCEM<- Empirical_conf_mu_50_MCEM+1
  
  if(true_sigma2>=CI_sigma2_50_MCEM[i,1] & true_sigma2<=CI_sigma2_50_MCEM[i,2]) 
    Empirical_conf_sigma2_50_MCEM<- Empirical_conf_sigma2_50_MCEM+1
}


Empirical_conf_mu_50_MCEM/500
Empirical_conf_sigma2_50_MCEM/500
###############################################################################################################################
################################################################################################################################
outputMCEM100m15

theta100<-matrix(rep(0,2*nrow(output100m15EM)),ncol=2)
for(i in 1:nrow(output100m15EM)){
  theta100[i,]<- c(output100m15EM[i,1],sqrt(output100m15EM[i,2]))
}
theta100


se_MCEM100m15<- matrix(rep(0,2*nrow(output100m15EM)),ncol=2)
for(i in 1:nrow(output100m15EM)){
  se_MCEM100m15[i,]<-  ZsimMCEM(theta100[i,],data=simdata100m15[,,i])
}
se_MCEM100m15
se_MCEM100m15[1:5,]
#se_EM_1000[1:5,]
#trueval<- c(68,(2.5)^2)
#########################################################################################
CI_mu_100_MCEM<- matrix(rep(0,2*nrow(outputMCEM100m15)),nrow=nrow(outputMCEM100m15),ncol=2)
colnames(CI_mu_100_MCEM)<- c("Lower bound","Upper bound")

for(i in 1:nrow(outputMCEM100m15)){
  CI_mu_100_MCEM[i,1]<- outputMCEM100m15[i,1]-(se_MCEM100m15[i,1]*qnorm(0.975))
  CI_mu_100_MCEM[i,2]<- outputMCEM100m15[i,1]+(se_MCEM100m15[i,1]*qnorm(0.975))
  
}

CI_mu_100_MCEM

CI_sigma2_100_MCEM<- matrix(rep(0,2*nrow(outputMCEM100m15)),nrow=nrow(outputMCEM100m15),ncol=2)
colnames(CI_sigma2_100_MCEM)<- c("Lower bound","Upper bound")


for(i in 1:nrow(outputMCEM100m15)){
  CI_sigma2_100_MCEM[i,1]<- outputMCEM100m15[i,2]-(se_MCEM100m15[i,2]*qnorm(0.975))
  CI_sigma2_100_MCEM[i,2]<- outputMCEM100m15[i,2]+(se_MCEM100m15[i,2]*qnorm(0.975))
  
}
CI_sigma2_100_MCEM

true_mu<- 68
true_sigma2<- (2.5)^2
true_sigma2


Empirical_conf_mu_100_MCEM<- 0
Empirical_conf_sigma2_100_MCEM<- 0

for(i in 1:nrow(outputMCEM100m15)){
  if(true_mu>=CI_mu_100_MCEM[i,1] & true_mu<=CI_mu_100_MCEM[i,2]) Empirical_conf_mu_100_MCEM<- Empirical_conf_mu_100_MCEM+1
  
  if(true_sigma2>=CI_sigma2_100_MCEM[i,1] & true_sigma2<=CI_sigma2_100_MCEM[i,2]) 
    Empirical_conf_sigma2_100_MCEM<- Empirical_conf_sigma2_100_MCEM+1
}


Empirical_conf_mu_100_MCEM/500
Empirical_conf_sigma2_100_MCEM/500
###############################################################################################################################
###############################################################################################################################
################################################################################################################################
outputMCEM300m15

theta300<-matrix(rep(0,2*nrow(output300m15EM)),ncol=2)
for(i in 1:nrow(output300m15EM)){
  theta300[i,]<- c(output300m15EM[i,1],sqrt(output300m15EM[i,2]))
}
theta300


se_MCEM300m15<- matrix(rep(0,2*nrow(output300m15EM)),ncol=2)
for(i in 1:nrow(output300m15EM)){
  se_MCEM300m15[i,]<-  ZsimMCEM(theta300[i,],data=simdata300m15[,,i])
}
se_MCEM300m15
se_MCEM300m15[1:5,]
#se_EM_1000[1:5,]
#trueval<- c(68,(2.5)^2)
#########################################################################################
CI_mu_300_MCEM<- matrix(rep(0,2*nrow(outputMCEM300m15)),nrow=nrow(outputMCEM300m15),ncol=2)
colnames(CI_mu_300_MCEM)<- c("Lower bound","Upper bound")

for(i in 1:nrow(outputMCEM300m15)){
  CI_mu_300_MCEM[i,1]<- outputMCEM300m15[i,1]-(se_MCEM300m15[i,1]*qnorm(0.975))
  CI_mu_300_MCEM[i,2]<- outputMCEM300m15[i,1]+(se_MCEM300m15[i,1]*qnorm(0.975))
  
}

CI_mu_300_MCEM

CI_sigma2_300_MCEM<- matrix(rep(0,2*nrow(outputMCEM300m15)),nrow=nrow(outputMCEM300m15),ncol=2)
colnames(CI_sigma2_300_MCEM)<- c("Lower bound","Upper bound")


for(i in 1:nrow(outputMCEM300m15)){
  CI_sigma2_300_MCEM[i,1]<- outputMCEM300m15[i,2]-(se_MCEM300m15[i,2]*qnorm(0.975))
  CI_sigma2_300_MCEM[i,2]<- outputMCEM300m15[i,2]+(se_MCEM300m15[i,2]*qnorm(0.975))
  
}
CI_sigma2_300_MCEM

true_mu<- 68
true_sigma2<- (2.5)^2
true_sigma2


Empirical_conf_mu_300_MCEM<- 0
Empirical_conf_sigma2_300_MCEM<- 0

for(i in 1:nrow(outputMCEM300m15)){
  if(true_mu>=CI_mu_300_MCEM[i,1] & true_mu<=CI_mu_300_MCEM[i,2]) Empirical_conf_mu_300_MCEM<- Empirical_conf_mu_300_MCEM+1
  
  if(true_sigma2>=CI_sigma2_300_MCEM[i,1] & true_sigma2<=CI_sigma2_300_MCEM[i,2]) 
    Empirical_conf_sigma2_300_MCEM<- Empirical_conf_sigma2_300_MCEM+1
}


Empirical_conf_mu_300_MCEM/500
Empirical_conf_sigma2_300_MCEM/500
###############################################################################################################################
################################################################################################################################
outputMCEM600m15

theta600<-matrix(rep(0,2*nrow(output600m15EM)),ncol=2)
for(i in 1:nrow(output600m15EM)){
  theta600[i,]<- c(output600m15EM[i,1],sqrt(output600m15EM[i,2]))
}
theta600


se_MCEM600m15<- matrix(rep(0,2*nrow(output600m15EM)),ncol=2)
for(i in 1:nrow(output600m15EM)){
  se_MCEM600m15[i,]<-  ZsimMCEM(theta600[i,],data=simdata600m15[,,i])
}
se_MCEM600m15
se_MCEM600m15[1:5,]
#se_EM_1000[1:5,]
#trueval<- c(68,(2.5)^2)
#########################################################################################
CI_mu_600_MCEM<- matrix(rep(0,2*nrow(outputMCEM600m15)),nrow=nrow(outputMCEM600m15),ncol=2)
colnames(CI_mu_600_MCEM)<- c("Lower bound","Upper bound")

for(i in 1:nrow(outputMCEM600m15)){
  CI_mu_600_MCEM[i,1]<- outputMCEM600m15[i,1]-(se_MCEM600m15[i,1]*qnorm(0.975))
  CI_mu_600_MCEM[i,2]<- outputMCEM600m15[i,1]+(se_MCEM600m15[i,1]*qnorm(0.975))
  
}

CI_mu_600_MCEM

CI_sigma2_600_MCEM<- matrix(rep(0,2*nrow(outputMCEM600m15)),nrow=nrow(outputMCEM600m15),ncol=2)
colnames(CI_sigma2_600_MCEM)<- c("Lower bound","Upper bound")


for(i in 1:nrow(outputMCEM600m15)){
  CI_sigma2_600_MCEM[i,1]<- outputMCEM600m15[i,2]-(se_MCEM600m15[i,2]*qnorm(0.975))
  CI_sigma2_600_MCEM[i,2]<- outputMCEM600m15[i,2]+(se_MCEM600m15[i,2]*qnorm(0.975))
  
}
CI_sigma2_600_MCEM

true_mu<- 68
true_sigma2<- (2.5)^2
true_sigma2


Empirical_conf_mu_600_MCEM<- 0
Empirical_conf_sigma2_600_MCEM<- 0

for(i in 1:nrow(outputMCEM600m15)){
  if(true_mu>=CI_mu_600_MCEM[i,1] & true_mu<=CI_mu_600_MCEM[i,2]) Empirical_conf_mu_600_MCEM<- Empirical_conf_mu_600_MCEM+1
  
  if(true_sigma2>=CI_sigma2_600_MCEM[i,1] & true_sigma2<=CI_sigma2_600_MCEM[i,2]) 
    Empirical_conf_sigma2_600_MCEM<- Empirical_conf_sigma2_600_MCEM+1
}


Empirical_conf_mu_600_MCEM/500
Empirical_conf_sigma2_600_MCEM/500
######################################################################################################################
######################################################################################################################
###############################################################################################################
################################################################################################################

Conf_Mu_prop_MCEM<- c(Empirical_conf_mu_50_MCEM/500,Empirical_conf_mu_100_MCEM/500,Empirical_conf_mu_300_MCEM/500,
                     Empirical_conf_mu_600_MCEM/500,Empirical_conf_mu_1000_MCEM/500)


Conf_sigma2_prop_MCEM<- c(Empirical_conf_sigma2_50_MCEM/500,Empirical_conf_sigma2_100_MCEM/500,Empirical_conf_sigma2_300_MCEM/500,
                         Empirical_conf_sigma2_600_MCEM/500,Empirical_conf_sigma2_1000_MCEM/500)

print(Conf_Mu_prop_MCEM*100)
print(Conf_sigma2_prop_MCEM*100)

##################################
n<- c(50,100,300,600,1000)
CI_MU_Prop_MCEM<- Conf_Mu_prop_MCEM*100
CI_Sigma2_Prop_MCEM<- Conf_sigma2_prop_MCEM*100
CI_MU_Prop_MCEM
CI_Sigma2_Prop_MCEM


outputMCEM1000m15
sd_MU_MCEM<- c(sd(outputMCEM50m15[,1]),sd(outputMCEM100m15[,1]),sd(outputMCEM300m15[,1]),
          sd(outputMCEM600m15[,1]),sd(outputMCEM1000m15[,1]))
sd_MU_MCEM

AVE_MU_MCEM<- c(mean(outputMCEM50m15[,1]),mean(outputMCEM100m15[,1]),mean(outputMCEM300m15[,1]),
                mean(outputMCEM600m15[,1]),mean(outputMCEM1000m15[,1]))
AVE_MU_MCEM


Se_MU_hat_MCEM<- c(mean(se_MCEM50m15[,1]),mean(se_MCEM100m15[,1]),mean(se_MCEM300m15[,1]),
                   mean(se_MCEM600m15[,1]),mean(se_MCEM1000m15[,1]))
Se_MU_hat_MCEM

sd_var_MCEM<- c(sd(outputMCEM50m15[,2]),sd(outputMCEM100m15[,2]),sd(outputMCEM300m15[,2]),
                sd(outputMCEM600m15[,2]),sd(outputMCEM1000m15[,2]))
sd_var_MCEM

Ave_var_MCEM<- c(mean(outputMCEM50m15[,2]),mean(outputMCEM100m15[,2]),mean(outputMCEM300m15[,2]),
                 mean(outputMCEM600m15[,2]),mean(outputMCEM1000m15[,2]))
Ave_var_MCEM

Se_Sigma2_hat_MCEM<-  c(mean(se_MCEM50m15[,2]),mean(se_MCEM100m15[,2]),mean(se_MCEM300m15[,2]),
                        mean(se_MCEM600m15[,2]),mean(se_MCEM1000m15[,2]))
Se_Sigma2_hat_MCEM
##########################################################################################################################################
#######################################################################################################################
std_Mu_estimates_MCEM<- cbind(n,AVE_MU_MCEM,sd_MU_MCEM,Se_MU_hat_MCEM,CI_MU_Prop_MCEM)
std_Mu_estimates_MCEM

std_Sigma2_estimates_MCEM<- cbind(n,Ave_var_MCEM,sd_var_MCEM,Se_Sigma2_hat_MCEM,CI_Sigma2_Prop_MCEM)
std_Sigma2_estimates_MCEM

library(xtable)
xtable(std_Mu_estimates_MCEM,digits=6)
xtable(std_Sigma2_estimates_MCEM,digits=6)



############################################################################################################################
##############################################################################################################################
##############################################################################################################################
