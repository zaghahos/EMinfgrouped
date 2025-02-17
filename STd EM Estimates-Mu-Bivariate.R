dataEM_n50_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM50b10.csv")
dataEM_n100_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM100b10.csv")
dataEM_n300_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM300b10.csv")
dataEM_n600_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM600b10.csv")
dataEM_n1000_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/EM Oct2021/EM1000b10.csv")

dataEM_n50_bin10<- dataEM_n50_bin10[,-1]
dataEM_n100_bin10<- dataEM_n100_bin10[,-1]
dataEM_n300_bin10<- dataEM_n300_bin10[,-1]
dataEM_n600_bin10<- dataEM_n600_bin10[,-1]
dataEM_n1000_bin10<- dataEM_n1000_bin10[,-1]


dataEM_n50_bin10
dataEM_n100_bin10
dataEM_n300_bin10
dataEM_n600_bin10
dataEM_n1000_bin10
##########################################################################################################################
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
  #print(EX)
  #print(CovX)
  EX_minus_MU <- matrix(rep(0,n*d),ncol=d,nrow=n) 
  EX_minus_MU<- EX-mu
  #print(EX_minus_MU)
  
  #print(EX_minus_MU[1,]%*%solve(sigma))
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
  #print(Inf_j_m)
  Inf_mat<- apply(Inf_j_m,c(1,2),sum)
  #print(Inf_mat)
  cov_Mu_EM<- solve(Inf_mat)
  #Se_mu_EM<- c(sqrt(cov_Mu_EM[1,1]),sqrt(cov_Mu_EM[2,2]))
  Se_mu_EM<- c(sqrt(diag(cov_Mu_EM)))
  names(Se_mu_EM)<- c("std_mu_x1","std_mu_x2")
  
  
  #print(cov_Mu_EM)
  #print(Se_mu_EM)
  
  return(Se_mu_EM)


}

#i=1
#Result<- Std_Mu_EM(Data=simulateddata1000b10[,,i],mu=mu_EM1000b10[i,],sigma=sigma_EM1000b10[,,i])   

#Result
############################################################################################################
#### N=1000###
EM1000b10<- dataEM_n1000_bin10
EM1000b10

mu_EM1000b10<- matrix(rep(0,2*nrow(EM1000b10)),ncol=2)
for(i in 1:nrow(EM1000b10)){
  mu_EM1000b10[i,]<- c(EM1000b10[i,1],EM1000b10[i,2])
}
mu_EM1000b10


sigma_EM1000b10<- array(rep(0,2*2*nrow(EM1000b10)),c(2,2,nrow(EM1000b10)))
#sigma_EM1000b10
for(i in 1:nrow(EM1000b10)){
  sigma_EM1000b10[,,i]<- matrix(c(EM1000b10[i,3],EM1000b10[i,4],EM1000b10[i,5],EM1000b10[i,6]),ncol=2)
}
sigma_EM1000b10
###############################################################################################################
Result_std_EM1000b10_Mu<- matrix(rep(0,nrow(EM1000b10)*2),ncol=2)


for(i in 1:nrow(EM1000b10)){
  Result_std_EM1000b10_Mu[i,]<- Std_Mu_EM(Data=simulateddata1000b10[,,i],mu=mu_EM1000b10[i,],sigma=sigma_EM1000b10[,,i])   
  
}

Result_std_EM1000b10_Mu


##################################################################################################################
CI_mu_x1_1000<- matrix(rep(0,2*nrow(EM1000b10)),ncol=2)
colnames(CI_mu_x1_1000)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_1000)){
  CI_mu_x1_1000[i,]<- c(EM1000b10[i,1]-(qnorm(0.975)*Result_std_EM1000b10_Mu[i,1]),
                   EM1000b10[i,1]+(qnorm(0.975)*Result_std_EM1000b10_Mu[i,1]))
  
}

CI_mu_x1_1000

##############################################################################################
CI_mu_x2_1000<- matrix(rep(0,2*nrow(EM1000b10)),ncol=2)
colnames(CI_mu_x2_1000)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_1000)){
  CI_mu_x2_1000[i,]<- c(EM1000b10[i,2]-(qnorm(0.975)*Result_std_EM1000b10_Mu[i,2]),
                   EM1000b10[i,2]+(qnorm(0.975)*Result_std_EM1000b10_Mu[i,2]))
  
}

CI_mu_x2_1000


#Result_std_EM1000b10_sigma
###############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_1000<- 0
for(i in 1:nrow(EM1000b10)){
  if (true_mu[1]>= CI_mu_x1_1000[i,1]& true_mu[1]<= CI_mu_x1_1000[i,2]) {counter_mu_x1_1000<- counter_mu_x1_1000+1}
}

counter_mu_x1_1000
counter_mu_x1_1000_percent<- counter_mu_x1_1000/500
counter_mu_x1_1000_percent
#####################################################################################################
counter_mu_x2_1000<- 0
for(i in 1:nrow(EM1000b10)){
  if (true_mu[2]>= CI_mu_x2_1000[i,1]& true_mu[2]<= CI_mu_x2_1000[i,2]) {counter_mu_x2_1000<- counter_mu_x2_1000+1}
}

counter_mu_x2_1000
counter_mu_x2_1000_percent<- counter_mu_x2_1000/500
counter_mu_x2_1000_percent
#################################################################################################################
#################################################################################################################
### N=50###
EM50b10<- dataEM_n50_bin10
EM50b10

mu_EM50b10<- matrix(rep(0,2*nrow(EM50b10)),ncol=2)
for(i in 1:nrow(EM50b10)){
  mu_EM50b10[i,]<- c(EM50b10[i,1],EM50b10[i,2])
}
mu_EM50b10


sigma_EM50b10<- array(rep(0,2*2*nrow(EM50b10)),c(2,2,nrow(EM50b10)))
#sigma_EM1000b10
for(i in 1:nrow(EM50b10)){
  sigma_EM50b10[,,i]<- matrix(c(EM50b10[i,3],EM50b10[i,4],EM50b10[i,5],EM50b10[i,6]),ncol=2)
}
sigma_EM50b10
###############################################################################################################
Result_std_EM50b10_Mu<- matrix(rep(0,nrow(EM50b10)*2),ncol=2)


for(i in 1:nrow(EM50b10)){
  Result_std_EM50b10_Mu[i,]<- Std_Mu_EM(Data=simulateddata50b10[,,i],mu=mu_EM50b10[i,],sigma=sigma_EM50b10[,,i])   
  
}

Result_std_EM50b10_Mu


##################################################################################################################
CI_mu_x1_50<- matrix(rep(0,2*nrow(EM50b10)),ncol=2)
colnames(CI_mu_x1_50)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_50)){
  CI_mu_x1_50[i,]<- c(EM50b10[i,1]-(qnorm(0.975)*Result_std_EM50b10_Mu[i,1]),
                        EM50b10[i,1]+(qnorm(0.975)*Result_std_EM50b10_Mu[i,1]))
  
}

CI_mu_x1_50

##############################################################################################
CI_mu_x2_50<- matrix(rep(0,2*nrow(EM50b10)),ncol=2)
colnames(CI_mu_x2_50)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_50)){
  CI_mu_x2_50[i,]<- c(EM50b10[i,2]-(qnorm(0.975)*Result_std_EM50b10_Mu[i,2]),
                        EM50b10[i,2]+(qnorm(0.975)*Result_std_EM50b10_Mu[i,2]))
  
}

CI_mu_x2_50


#Result_std_EM1000b10_sigma
###############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_50<- 0
for(i in 1:nrow(EM50b10)){
  if (true_mu[1]>= CI_mu_x1_50[i,1]& true_mu[1]<= CI_mu_x1_50[i,2]) {counter_mu_x1_50<- counter_mu_x1_50+1}
}

counter_mu_x1_50
counter_mu_x1_50_percent<- counter_mu_x1_50/500
counter_mu_x1_50_percent
#####################################################################################################
counter_mu_x2_50<- 0
for(i in 1:nrow(EM50b10)){
  if (true_mu[2]>= CI_mu_x2_50[i,1]& true_mu[2]<= CI_mu_x2_50[i,2]) {counter_mu_x2_50<- counter_mu_x2_50+1}
}

counter_mu_x2_50
counter_mu_x2_50_percent<- counter_mu_x2_50/500
counter_mu_x2_50_percent
#################################################################################################################
#################################################################################################################
### N=100###
############################################################################################################
EM100b10<- dataEM_n100_bin10
EM100b10

mu_EM100b10<- matrix(rep(0,2*nrow(EM100b10)),ncol=2)
for(i in 1:nrow(EM100b10)){
  mu_EM100b10[i,]<- c(EM100b10[i,1],EM100b10[i,2])
}
mu_EM100b10


sigma_EM100b10<- array(rep(0,2*2*nrow(EM100b10)),c(2,2,nrow(EM100b10)))
#sigma_EM1000b10
for(i in 1:nrow(EM100b10)){
  sigma_EM100b10[,,i]<- matrix(c(EM100b10[i,3],EM100b10[i,4],EM100b10[i,5],EM100b10[i,6]),ncol=2)
}
sigma_EM100b10
###############################################################################################################
Result_std_EM100b10_Mu<- matrix(rep(0,nrow(EM100b10)*2),ncol=2)


for(i in 1:nrow(EM100b10)){
  Result_std_EM100b10_Mu[i,]<- Std_Mu_EM(Data=simulateddata100b10[,,i],mu=mu_EM100b10[i,],sigma=sigma_EM100b10[,,i])   
  
}

Result_std_EM100b10_Mu


##################################################################################################################
CI_mu_x1_100<- matrix(rep(0,2*nrow(EM100b10)),ncol=2)
colnames(CI_mu_x1_100)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_100)){
  CI_mu_x1_100[i,]<- c(EM100b10[i,1]-(qnorm(0.975)*Result_std_EM100b10_Mu[i,1]),
                        EM100b10[i,1]+(qnorm(0.975)*Result_std_EM100b10_Mu[i,1]))
  
}

CI_mu_x1_100

##############################################################################################
CI_mu_x2_100<- matrix(rep(0,2*nrow(EM100b10)),ncol=2)
colnames(CI_mu_x2_100)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_100)){
  CI_mu_x2_100[i,]<- c(EM100b10[i,2]-(qnorm(0.975)*Result_std_EM100b10_Mu[i,2]),
                        EM100b10[i,2]+(qnorm(0.975)*Result_std_EM100b10_Mu[i,2]))
  
}

CI_mu_x2_100


#Result_std_EM1000b10_sigma
###############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_100<- 0
for(i in 1:nrow(EM100b10)){
  if (true_mu[1]>= CI_mu_x1_100[i,1]& true_mu[1]<= CI_mu_x1_100[i,2]) {counter_mu_x1_100<- counter_mu_x1_100+1}
}

counter_mu_x1_100
counter_mu_x1_100_percent<- counter_mu_x1_100/500
counter_mu_x1_100_percent
#####################################################################################################
counter_mu_x2_100<- 0
for(i in 1:nrow(EM100b10)){
  if (true_mu[2]>= CI_mu_x2_100[i,1]& true_mu[2]<= CI_mu_x2_100[i,2]) {counter_mu_x2_100<- counter_mu_x2_100+1}
}

counter_mu_x2_100
counter_mu_x2_100_percent<- counter_mu_x2_100/500
counter_mu_x2_100_percent
#################################################################################################################
#################################################################################################################
### N=300###
############################################################################################################
EM300b10<- dataEM_n300_bin10
EM300b10

mu_EM300b10<- matrix(rep(0,2*nrow(EM300b10)),ncol=2)
for(i in 1:nrow(EM300b10)){
  mu_EM300b10[i,]<- c(EM300b10[i,1],EM300b10[i,2])
}
mu_EM300b10


sigma_EM300b10<- array(rep(0,2*2*nrow(EM300b10)),c(2,2,nrow(EM300b10)))
#sigma_EM1000b10
for(i in 1:nrow(EM300b10)){
  sigma_EM300b10[,,i]<- matrix(c(EM300b10[i,3],EM300b10[i,4],EM300b10[i,5],EM300b10[i,6]),ncol=2)
}
sigma_EM300b10
###############################################################################################################
Result_std_EM300b10_Mu<- matrix(rep(0,nrow(EM300b10)*2),ncol=2)


for(i in 1:nrow(EM300b10)){
  Result_std_EM300b10_Mu[i,]<- Std_Mu_EM(Data=simulateddata300b10[,,i],mu=mu_EM300b10[i,],sigma=sigma_EM300b10[,,i])   
  
}

Result_std_EM300b10_Mu


##################################################################################################################
CI_mu_x1_300<- matrix(rep(0,2*nrow(EM300b10)),ncol=2)
colnames(CI_mu_x1_300)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_300)){
  CI_mu_x1_300[i,]<- c(EM300b10[i,1]-(qnorm(0.975)*Result_std_EM300b10_Mu[i,1]),
                       EM300b10[i,1]+(qnorm(0.975)*Result_std_EM300b10_Mu[i,1]))
  
}

CI_mu_x1_300

##############################################################################################
CI_mu_x2_300<- matrix(rep(0,2*nrow(EM300b10)),ncol=2)
colnames(CI_mu_x2_300)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_300)){
  CI_mu_x2_300[i,]<- c(EM300b10[i,2]-(qnorm(0.975)*Result_std_EM300b10_Mu[i,2]),
                       EM300b10[i,2]+(qnorm(0.975)*Result_std_EM300b10_Mu[i,2]))
  
}

CI_mu_x2_300


#Result_std_EM1000b10_sigma
###############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_300<- 0
for(i in 1:nrow(EM300b10)){
  if (true_mu[1]>= CI_mu_x1_300[i,1]& true_mu[1]<= CI_mu_x1_300[i,2]) {counter_mu_x1_300<- counter_mu_x1_300+1}
}

counter_mu_x1_300
counter_mu_x1_300_percent<- counter_mu_x1_300/500
counter_mu_x1_300_percent
#####################################################################################################
counter_mu_x2_300<- 0
for(i in 1:nrow(EM300b10)){
  if (true_mu[2]>= CI_mu_x2_300[i,1]& true_mu[2]<= CI_mu_x2_300[i,2]) {counter_mu_x2_300<- counter_mu_x2_300+1}
}

counter_mu_x2_300
counter_mu_x2_300_percent<- counter_mu_x2_300/500
counter_mu_x2_300_percent
#################################################################################################################
#################################################################################################################
### N=600###
############################################################################################################
EM600b10<- dataEM_n600_bin10
EM600b10

mu_EM600b10<- matrix(rep(0,2*nrow(EM600b10)),ncol=2)
for(i in 1:nrow(EM600b10)){
  mu_EM600b10[i,]<- c(EM600b10[i,1],EM600b10[i,2])
}
mu_EM600b10


sigma_EM600b10<- array(rep(0,2*2*nrow(EM600b10)),c(2,2,nrow(EM600b10)))
#sigma_EM1000b10
for(i in 1:nrow(EM600b10)){
  sigma_EM600b10[,,i]<- matrix(c(EM600b10[i,3],EM600b10[i,4],EM600b10[i,5],EM600b10[i,6]),ncol=2)
}
sigma_EM600b10
###############################################################################################################
Result_std_EM600b10_Mu<- matrix(rep(0,nrow(EM600b10)*2),ncol=2)


for(i in 1:nrow(EM600b10)){
  Result_std_EM600b10_Mu[i,]<- Std_Mu_EM(Data=simulateddata600b10[,,i],mu=mu_EM600b10[i,],sigma=sigma_EM600b10[,,i])   
  
}

Result_std_EM600b10_Mu


##################################################################################################################
CI_mu_x1_600<- matrix(rep(0,2*nrow(EM600b10)),ncol=2)
colnames(CI_mu_x1_600)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_600)){
  CI_mu_x1_600[i,]<- c(EM600b10[i,1]-(qnorm(0.975)*Result_std_EM600b10_Mu[i,1]),
                       EM600b10[i,1]+(qnorm(0.975)*Result_std_EM600b10_Mu[i,1]))
  
}

CI_mu_x1_600

##############################################################################################
CI_mu_x2_600<- matrix(rep(0,2*nrow(EM600b10)),ncol=2)
colnames(CI_mu_x2_600)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_600)){
  CI_mu_x2_600[i,]<- c(EM600b10[i,2]-(qnorm(0.975)*Result_std_EM600b10_Mu[i,2]),
                       EM600b10[i,2]+(qnorm(0.975)*Result_std_EM600b10_Mu[i,2]))
  
}

CI_mu_x2_600


#Result_std_EM1000b10_sigma
###############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_600<- 0
for(i in 1:nrow(EM600b10)){
  if (true_mu[1]>= CI_mu_x1_600[i,1]& true_mu[1]<= CI_mu_x1_600[i,2]) {counter_mu_x1_600<- counter_mu_x1_600+1}
}

counter_mu_x1_600
counter_mu_x1_600_percent<- counter_mu_x1_600/500
counter_mu_x1_600_percent
#####################################################################################################
counter_mu_x2_600<- 0
for(i in 1:nrow(EM600b10)){
  if (true_mu[2]>= CI_mu_x2_600[i,1]& true_mu[2]<= CI_mu_x2_600[i,2]) {counter_mu_x2_600<- counter_mu_x2_600+1}
}

counter_mu_x2_600
counter_mu_x2_600_percent<- counter_mu_x2_600/500
counter_mu_x2_600_percent
#################################################################################################################
#################################################################################################################
### REsults for MU_x1###
sample_size<- c(50,50,100,100,300,300,600,600,1000,1000)
sample_size1<- c(50,100,300,600,1000)
method<- rep(c("EM","MCEM"),5)

Ave_mu_x1_EM<- c(mean(EM50b10[,1]),mean(EM100b10[,1]),mean(EM300b10[,1]),mean(EM600b10[,1]),mean(EM1000b10[,1]))

sd_mu_x1_EM<- c(sd(EM50b10[,1]),sd(EM100b10[,1]),sd(EM300b10[,1]),sd(EM600b10[,1]),sd(EM1000b10[,1]))

se_mu_x1_EM<- c(mean(Result_std_EM50b10_Mu[,1]),mean(Result_std_EM100b10_Mu[,1]),mean(Result_std_EM300b10_Mu[,1]),
                mean(Result_std_EM600b10_Mu[,1]),mean(Result_std_EM1000b10_Mu[,1]))

CI_mu_X1_EM<- c(counter_mu_x1_50_percent,counter_mu_x1_100_percent,counter_mu_x1_300_percent,
                counter_mu_x1_600_percent,counter_mu_x1_1000_percent)

Res1<- as.data.frame(cbind(sample_size1,Ave_mu_x1_EM,sd_mu_x1_EM,se_mu_x1_EM,CI_mu_X1_EM))
Res1


Ave_mu_x2_EM<- c(mean(EM50b10[,2]),mean(EM100b10[,2]),mean(EM300b10[,2]),mean(EM600b10[,2]),mean(EM1000b10[,2]))

sd_mu_x2_EM<- c(sd(EM50b10[,2]),sd(EM100b10[,2]),sd(EM300b10[,2]),sd(EM600b10[,2]),sd(EM1000b10[,2]))

se_mu_x2_EM<- c(mean(Result_std_EM50b10_Mu[,2]),mean(Result_std_EM100b10_Mu[,2]),mean(Result_std_EM300b10_Mu[,2]),
                mean(Result_std_EM600b10_Mu[,2]),mean(Result_std_EM1000b10_Mu[,2]))

CI_mu_X2_EM<- c(counter_mu_x2_50_percent,counter_mu_x2_100_percent,counter_mu_x2_300_percent,
                counter_mu_x2_600_percent,counter_mu_x2_1000_percent)

Res2<- as.data.frame(cbind(sample_size1,Ave_mu_x2_EM,sd_mu_x2_EM,se_mu_x2_EM,CI_mu_X2_EM))
Res2
