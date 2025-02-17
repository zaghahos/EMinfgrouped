dataMCEM_n50_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MCEM Oct2021/MCEM50b10.csv")
dataMCEM_n100_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MCEM Oct2021/MCEM100b10.csv")
dataMCEM_n300_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MCEM Oct2021/MCEM300b10.csv")
dataMCEM_n600_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MCEM Oct2021/MCEM600b10.csv")
dataMCEM_n1000_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MCEM Oct2021/MCEM1000b10.csv")


dataMCEM_n50_bin10<- dataMCEM_n50_bin10[,-1]
dataMCEM_n100_bin10<- dataMCEM_n100_bin10[,-1]
dataMCEM_n300_bin10<- dataMCEM_n300_bin10[,-1]
dataMCEM_n600_bin10<- dataMCEM_n600_bin10[,-1]
dataMCEM_n1000_bin10<- dataMCEM_n1000_bin10[,-1]


dataMCEM_n50_bin10
dataMCEM_n100_bin10
dataMCEM_n300_bin10
dataMCEM_n600_bin10
dataMCEM_n1000_bin10
#############################################################################################################
#############################################################################################################
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
  #print(zsim)
  #print(dim(zsim))
  d2<- solve(sigma)
  #print(-d2)
  d1<- array(rep(0,k*d*n),c(k,d,n)) 
  for(i in 1:n){
    d1[,,i]<- (zsim[,,i]-mu)%*%solve(sigma)
  }
  #print(d1)
  
  
  #print(d1[,,1])
  #print(apply(d1[,,1],2,mean))
  md1<- matrix(rep(0,d*n),ncol=d)
  for(i in 1:n){
    md1[i,]<- apply(d1[,,i],2,mean)
  }
  #print(md1)
  diff1<- array(rep(0,k*d*n),c(k,d,n)) 
  for(i in 1:n){
    diff1[,,i]<- d1[,,i]-md1[i,] 
  }
  #print(diff1)
  der1<- matrix(rep(0,d*n),ncol=d)
  for(i in 1:n){
    der1[i,]<- apply(diff1[,,i],2,mean)
  }
  #print(der1)
  
  diff2<- array(rep(0,d*d*n),c(d,d,n))
  for(i in 1:n){
    diff2[,,i]<- der1[i,]%*%t(der1[i,])
  }
  #print(diff2)
  sj<- array(rep(0,d*d*n),c(d,d,n))
  for(i in 1:n){
    sj[,,i]<- (d2+diff2[,,i])*Data[i,ncol(Data)]
  }
  #print(sj)
  IM<- apply(sj,c(1,2),sum)
  #print(IM)
  Inf_mat<- solve(IM)
  #print(Inf_mat)
  se_MCEM_MU<- diag(sqrt(Inf_mat))
  names(se_MCEM_MU)<- c("se_mu_x1","se_mu_x2")
  #print(se_MCEM_MU)
  return(se_MCEM_MU)
}

############################################################################################################
##### TEST####
#i=1
#Result<- STD_MU_MCEM_Biv(Data=simulateddata1000b10[,,i],mu=mu_MCEM1000b10[i,],sigma=sigma_MCEM1000b10[,,i])   
#Result
####N=1000###
MCEM1000b10<- dataMCEM_n1000_bin10
MCEM1000b10

mu_MCEM1000b10<- matrix(rep(0,2*nrow(MCEM1000b10)),ncol=2)
for(i in 1:nrow(MCEM1000b10)){
  mu_MCEM1000b10[i,]<- c(MCEM1000b10[i,1],MCEM1000b10[i,2])
}
mu_MCEM1000b10


sigma_MCEM1000b10<- array(rep(0,2*2*nrow(MCEM1000b10)),c(2,2,nrow(MCEM1000b10)))
#sigma_EM1000b10
for(i in 1:nrow(MCEM1000b10)){
  sigma_MCEM1000b10[,,i]<- matrix(c(MCEM1000b10[i,3],MCEM1000b10[i,4],MCEM1000b10[i,5],MCEM1000b10[i,6]),ncol=2)
}
sigma_MCEM1000b10
###############################################################################################################
Result_std_MCEM1000b10_Mu<- matrix(rep(0,nrow(MCEM1000b10)*2),ncol=2)


for(i in 1:nrow(MCEM1000b10)){
  Result_std_MCEM1000b10_Mu[i,]<- STD_MU_MCEM_Biv(Data=simulateddata1000b10[,,i],mu=mu_MCEM1000b10[i,],sigma=sigma_MCEM1000b10[,,i])   
  
}

Result_std_MCEM1000b10_Mu

################
CI_mu_x1_1000_MCEM<- matrix(rep(0,2*nrow(MCEM1000b10)),ncol=2)
colnames(CI_mu_x1_1000_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_1000_MCEM)){
  CI_mu_x1_1000_MCEM[i,]<- c(MCEM1000b10[i,1]-(qnorm(0.975)*Result_std_MCEM1000b10_Mu[i,1]),
                        MCEM1000b10[i,1]+(qnorm(0.975)*Result_std_MCEM1000b10_Mu[i,1]))
  
}

CI_mu_x1_1000_MCEM

CI_mu_x2_1000_MCEM<- matrix(rep(0,2*nrow(MCEM1000b10)),ncol=2)
colnames(CI_mu_x2_1000_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_1000_MCEM)){
  CI_mu_x2_1000_MCEM[i,]<- c(MCEM1000b10[i,2]-(qnorm(0.975)*Result_std_MCEM1000b10_Mu[i,2]),
                             MCEM1000b10[i,2]+(qnorm(0.975)*Result_std_MCEM1000b10_Mu[i,2]))
  
}

CI_mu_x2_1000_MCEM



###############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_1000_MCEM<- 0
for(i in 1:nrow(MCEM1000b10)){
  if (true_mu[1]>= CI_mu_x1_1000_MCEM[i,1]& true_mu[1]<= CI_mu_x1_1000_MCEM[i,2])
    {counter_mu_x1_1000_MCEM<- counter_mu_x1_1000_MCEM+1}
}

counter_mu_x1_1000_MCEM
counter_mu_x1_1000_MCEM_percent<- counter_mu_x1_1000_MCEM/500
counter_mu_x1_1000_MCEM_percent
#####################################################################################################
counter_mu_x2_1000_MCEM<- 0
for(i in 1:nrow(MCEM1000b10)){
  if (true_mu[2]>= CI_mu_x2_1000_MCEM[i,1]& true_mu[2]<= CI_mu_x2_1000_MCEM[i,2]) 
    {counter_mu_x2_1000_MCEM<- counter_mu_x2_1000_MCEM+1}
}

counter_mu_x2_1000_MCEM
counter_mu_x2_1000_MCEM_percent<- counter_mu_x2_1000_MCEM/500
counter_mu_x2_1000_MCEM_percent

###################################################################################################################
###################################################################################################################
### N=50###
MCEM50b10<- dataMCEM_n50_bin10
MCEM50b10

mu_MCEM50b10<- matrix(rep(0,2*nrow(MCEM50b10)),ncol=2)
for(i in 1:nrow(MCEM50b10)){
  mu_MCEM50b10[i,]<- c(MCEM50b10[i,1],MCEM50b10[i,2])
}
mu_MCEM50b10


sigma_MCEM50b10<- array(rep(0,2*2*nrow(MCEM50b10)),c(2,2,nrow(MCEM50b10)))
#sigma_EM1000b10
for(i in 1:nrow(MCEM50b10)){
  sigma_MCEM50b10[,,i]<- matrix(c(MCEM50b10[i,3],MCEM50b10[i,4],MCEM50b10[i,5],MCEM50b10[i,6]),ncol=2)
}
sigma_MCEM50b10
###############################################################################################################
Result_std_MCEM50b10_Mu<- matrix(rep(0,nrow(MCEM50b10)*2),ncol=2)


for(i in 1:nrow(MCEM50b10)){
  Result_std_MCEM50b10_Mu[i,]<- STD_MU_MCEM_Biv(Data=simulateddata50b10[,,i],mu=mu_MCEM50b10[i,],sigma=sigma_MCEM50b10[,,i])   
  
}

Result_std_MCEM50b10_Mu
############################
CI_mu_x1_50_MCEM<- matrix(rep(0,2*nrow(MCEM50b10)),ncol=2)
colnames(CI_mu_x1_50_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_50_MCEM)){
  CI_mu_x1_50_MCEM[i,]<- c(MCEM50b10[i,1]-(qnorm(0.975)*Result_std_MCEM50b10_Mu[i,1]),
                             MCEM50b10[i,1]+(qnorm(0.975)*Result_std_MCEM50b10_Mu[i,1]))
  
}

CI_mu_x1_50_MCEM

CI_mu_x2_50_MCEM<- matrix(rep(0,2*nrow(MCEM50b10)),ncol=2)
colnames(CI_mu_x2_50_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_50_MCEM)){
  CI_mu_x2_50_MCEM[i,]<- c(MCEM50b10[i,2]-(qnorm(0.975)*Result_std_MCEM50b10_Mu[i,2]),
                             MCEM50b10[i,2]+(qnorm(0.975)*Result_std_MCEM50b10_Mu[i,2]))
  
}

CI_mu_x2_50_MCEM

###############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_50_MCEM<- 0
for(i in 1:nrow(MCEM50b10)){
  if (true_mu[1]>= CI_mu_x1_50_MCEM[i,1]& true_mu[1]<= CI_mu_x1_50_MCEM[i,2])
  {counter_mu_x1_50_MCEM<- counter_mu_x1_50_MCEM+1}
}

counter_mu_x1_50_MCEM
counter_mu_x1_50_MCEM_percent<- counter_mu_x1_50_MCEM/500
counter_mu_x1_50_MCEM_percent
#####################################################################################################
counter_mu_x2_50_MCEM<- 0
for(i in 1:nrow(MCEM50b10)){
  if (true_mu[2]>= CI_mu_x2_50_MCEM[i,1]& true_mu[2]<= CI_mu_x2_50_MCEM[i,2]) 
  {counter_mu_x2_50_MCEM<- counter_mu_x2_50_MCEM+1}
}

counter_mu_x2_50_MCEM
counter_mu_x2_50_MCEM_percent<- counter_mu_x2_50_MCEM/500
counter_mu_x2_50_MCEM_percent

###################################################################################################################
###################################################################################################################
### N=100###
MCEM100b10<- dataMCEM_n100_bin10
MCEM100b10

mu_MCEM100b10<- matrix(rep(0,2*nrow(MCEM100b10)),ncol=2)
for(i in 1:nrow(MCEM100b10)){
  mu_MCEM100b10[i,]<- c(MCEM100b10[i,1],MCEM100b10[i,2])
}
mu_MCEM100b10


sigma_MCEM100b10<- array(rep(0,2*2*nrow(MCEM100b10)),c(2,2,nrow(MCEM100b10)))
#sigma_EM1000b10
for(i in 1:nrow(MCEM100b10)){
  sigma_MCEM100b10[,,i]<- matrix(c(MCEM100b10[i,3],MCEM100b10[i,4],MCEM100b10[i,5],MCEM100b10[i,6]),ncol=2)
}
sigma_MCEM100b10
###############################################################################################################
Result_std_MCEM100b10_Mu<- matrix(rep(0,nrow(MCEM100b10)*2),ncol=2)


for(i in 1:nrow(MCEM100b10)){
  Result_std_MCEM100b10_Mu[i,]<- STD_MU_MCEM_Biv(Data=simulateddata100b10[,,i],mu=mu_MCEM100b10[i,],sigma=sigma_MCEM100b10[,,i])   
  
}

Result_std_MCEM100b10_Mu
#################################################
CI_mu_x1_100_MCEM<- matrix(rep(0,2*nrow(MCEM100b10)),ncol=2)
colnames(CI_mu_x1_100_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_100_MCEM)){
  CI_mu_x1_100_MCEM[i,]<- c(MCEM100b10[i,1]-(qnorm(0.975)*Result_std_MCEM100b10_Mu[i,1]),
                             MCEM100b10[i,1]+(qnorm(0.975)*Result_std_MCEM100b10_Mu[i,1]))
  
}

CI_mu_x1_100_MCEM

CI_mu_x2_100_MCEM<- matrix(rep(0,2*nrow(MCEM100b10)),ncol=2)
colnames(CI_mu_x2_100_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_100_MCEM)){
  CI_mu_x2_100_MCEM[i,]<- c(MCEM100b10[i,2]-(qnorm(0.975)*Result_std_MCEM100b10_Mu[i,2]),
                             MCEM100b10[i,2]+(qnorm(0.975)*Result_std_MCEM100b10_Mu[i,2]))
  
}

CI_mu_x2_100_MCEM
##############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_100_MCEM<- 0
for(i in 1:nrow(MCEM100b10)){
  if (true_mu[1]>= CI_mu_x1_100_MCEM[i,1]& true_mu[1]<= CI_mu_x1_100_MCEM[i,2])
  {counter_mu_x1_100_MCEM<- counter_mu_x1_100_MCEM+1}
}

counter_mu_x1_100_MCEM
counter_mu_x1_100_MCEM_percent<- counter_mu_x1_100_MCEM/500
counter_mu_x1_100_MCEM_percent
#####################################################################################################
counter_mu_x2_100_MCEM<- 0
for(i in 1:nrow(MCEM100b10)){
  if (true_mu[2]>= CI_mu_x2_100_MCEM[i,1]& true_mu[2]<= CI_mu_x2_100_MCEM[i,2]) 
  {counter_mu_x2_100_MCEM<- counter_mu_x2_100_MCEM+1}
}

counter_mu_x2_100_MCEM
counter_mu_x2_100_MCEM_percent<- counter_mu_x2_100_MCEM/500
counter_mu_x2_100_MCEM_percent


###################################################################################################################
###################################################################################################################
### N=300###
MCEM300b10<- dataMCEM_n300_bin10
MCEM300b10

mu_MCEM300b10<- matrix(rep(0,2*nrow(MCEM300b10)),ncol=2)
for(i in 1:nrow(MCEM300b10)){
  mu_MCEM300b10[i,]<- c(MCEM300b10[i,1],MCEM300b10[i,2])
}
mu_MCEM300b10


sigma_MCEM300b10<- array(rep(0,2*2*nrow(MCEM300b10)),c(2,2,nrow(MCEM300b10)))
#sigma_EM1000b10
for(i in 1:nrow(MCEM300b10)){
  sigma_MCEM300b10[,,i]<- matrix(c(MCEM300b10[i,3],MCEM300b10[i,4],MCEM300b10[i,5],MCEM300b10[i,6]),ncol=2)
}
sigma_MCEM300b10
###############################################################################################################
Result_std_MCEM300b10_Mu<- matrix(rep(0,nrow(MCEM300b10)*2),ncol=2)


for(i in 1:nrow(MCEM300b10)){
  Result_std_MCEM300b10_Mu[i,]<- STD_MU_MCEM_Biv(Data=simulateddata300b10[,,i],mu=mu_MCEM300b10[i,],sigma=sigma_MCEM300b10[,,i])   
  
}

Result_std_MCEM300b10_Mu
#################################################
CI_mu_x1_300_MCEM<- matrix(rep(0,2*nrow(MCEM300b10)),ncol=2)
colnames(CI_mu_x1_300_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_300_MCEM)){
  CI_mu_x1_300_MCEM[i,]<- c(MCEM300b10[i,1]-(qnorm(0.975)*Result_std_MCEM300b10_Mu[i,1]),
                            MCEM300b10[i,1]+(qnorm(0.975)*Result_std_MCEM300b10_Mu[i,1]))
  
}

CI_mu_x1_300_MCEM

CI_mu_x2_300_MCEM<- matrix(rep(0,2*nrow(MCEM300b10)),ncol=2)
colnames(CI_mu_x2_300_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_300_MCEM)){
  CI_mu_x2_300_MCEM[i,]<- c(MCEM300b10[i,2]-(qnorm(0.975)*Result_std_MCEM300b10_Mu[i,2]),
                            MCEM300b10[i,2]+(qnorm(0.975)*Result_std_MCEM300b10_Mu[i,2]))
  
}

CI_mu_x2_300_MCEM
##############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_300_MCEM<- 0
for(i in 1:nrow(MCEM300b10)){
  if (true_mu[1]>= CI_mu_x1_300_MCEM[i,1]& true_mu[1]<= CI_mu_x1_300_MCEM[i,2])
  {counter_mu_x1_300_MCEM<- counter_mu_x1_300_MCEM+1}
}

counter_mu_x1_300_MCEM
counter_mu_x1_300_MCEM_percent<- counter_mu_x1_300_MCEM/500
counter_mu_x1_300_MCEM_percent
#####################################################################################################
counter_mu_x2_300_MCEM<- 0
for(i in 1:nrow(MCEM300b10)){
  if (true_mu[2]>= CI_mu_x2_300_MCEM[i,1]& true_mu[2]<= CI_mu_x2_300_MCEM[i,2]) 
  {counter_mu_x2_300_MCEM<- counter_mu_x2_300_MCEM+1}
}

counter_mu_x2_300_MCEM
counter_mu_x2_300_MCEM_percent<- counter_mu_x2_300_MCEM/500
counter_mu_x2_300_MCEM_percent

###################################################################################################################
###################################################################################################################
### N=600###
MCEM600b10<- dataMCEM_n600_bin10
MCEM600b10

mu_MCEM600b10<- matrix(rep(0,2*nrow(MCEM600b10)),ncol=2)
for(i in 1:nrow(MCEM600b10)){
  mu_MCEM600b10[i,]<- c(MCEM600b10[i,1],MCEM600b10[i,2])
}
mu_MCEM600b10


sigma_MCEM600b10<- array(rep(0,2*2*nrow(MCEM600b10)),c(2,2,nrow(MCEM600b10)))
#sigma_EM1000b10
for(i in 1:nrow(MCEM600b10)){
  sigma_MCEM600b10[,,i]<- matrix(c(MCEM600b10[i,3],MCEM600b10[i,4],MCEM600b10[i,5],MCEM600b10[i,6]),ncol=2)
}
sigma_MCEM600b10
###############################################################################################################
Result_std_MCEM600b10_Mu<- matrix(rep(0,nrow(MCEM600b10)*2),ncol=2)


for(i in 1:nrow(MCEM600b10)){
  Result_std_MCEM600b10_Mu[i,]<- STD_MU_MCEM_Biv(Data=simulateddata600b10[,,i],mu=mu_MCEM600b10[i,],sigma=sigma_MCEM600b10[,,i])   
  
}

Result_std_MCEM600b10_Mu

#################################################
CI_mu_x1_600_MCEM<- matrix(rep(0,2*nrow(MCEM600b10)),ncol=2)
colnames(CI_mu_x1_600_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1_600_MCEM)){
  CI_mu_x1_600_MCEM[i,]<- c(MCEM600b10[i,1]-(qnorm(0.975)*Result_std_MCEM600b10_Mu[i,1]),
                            MCEM600b10[i,1]+(qnorm(0.975)*Result_std_MCEM600b10_Mu[i,1]))
  
}

CI_mu_x1_600_MCEM

CI_mu_x2_600_MCEM<- matrix(rep(0,2*nrow(MCEM600b10)),ncol=2)
colnames(CI_mu_x2_600_MCEM)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x2_600_MCEM)){
  CI_mu_x2_600_MCEM[i,]<- c(MCEM600b10[i,2]-(qnorm(0.975)*Result_std_MCEM600b10_Mu[i,2]),
                            MCEM600b10[i,2]+(qnorm(0.975)*Result_std_MCEM600b10_Mu[i,2]))
  
}

CI_mu_x2_600_MCEM
##############################################################################################
true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
#true_cov
counter_mu_x1_600_MCEM<- 0
for(i in 1:nrow(MCEM600b10)){
  if (true_mu[1]>= CI_mu_x1_600_MCEM[i,1]& true_mu[1]<= CI_mu_x1_600_MCEM[i,2])
  {counter_mu_x1_600_MCEM<- counter_mu_x1_600_MCEM+1}
}

counter_mu_x1_600_MCEM
counter_mu_x1_600_MCEM_percent<- counter_mu_x1_600_MCEM/500
counter_mu_x1_600_MCEM_percent
#####################################################################################################
counter_mu_x2_600_MCEM<- 0
for(i in 1:nrow(MCEM600b10)){
  if (true_mu[2]>= CI_mu_x2_600_MCEM[i,1]& true_mu[2]<= CI_mu_x2_600_MCEM[i,2]) 
  {counter_mu_x2_600_MCEM<- counter_mu_x2_600_MCEM+1}
}

counter_mu_x2_600_MCEM
counter_mu_x2_600_MCEM_percent<- counter_mu_x2_600_MCEM/500
counter_mu_x2_600_MCEM_percent

###########################################################################################################
###########################################################################################################
###########################################################################################################
### REsults for MU_x1###
sample_size<- c(50,50,100,100,300,300,600,600,1000,1000)
sample_size1<- c(50,100,300,600,1000)
method<- rep(c("EM","MCEM"),5)

Ave_mu_x1_MCEM<- c(mean(MCEM50b10[,1]),mean(MCEM100b10[,1]),mean(MCEM300b10[,1]),mean(MCEM600b10[,1]),mean(MCEM1000b10[,1]))

sd_mu_x1_MCEM<- c(sd(MCEM50b10[,1]),sd(MCEM100b10[,1]),sd(MCEM300b10[,1]),sd(MCEM600b10[,1]),sd(MCEM1000b10[,1]))

se_mu_x1_MCEM<- c(mean(Result_std_MCEM50b10_Mu[,1]),mean(Result_std_MCEM100b10_Mu[,1]),mean(Result_std_MCEM300b10_Mu[,1]),
                mean(Result_std_MCEM600b10_Mu[,1]),mean(Result_std_MCEM1000b10_Mu[,1]))

CI_mu_X1_MCEM<- c(counter_mu_x1_50_MCEM_percent,counter_mu_x1_100_MCEM_percent,counter_mu_x1_300_MCEM_percent,
                counter_mu_x1_600_MCEM_percent,counter_mu_x1_1000_MCEM_percent)

Res3<- as.data.frame(cbind(sample_size1,Ave_mu_x1_MCEM,sd_mu_x1_MCEM,se_mu_x1_MCEM,CI_mu_X1_MCEM))
Res3



Ave_mu_x2_MCEM<- c(mean(MCEM50b10[,2]),mean(MCEM100b10[,2]),mean(MCEM300b10[,2]),mean(MCEM600b10[,2]),mean(MCEM1000b10[,2]))

sd_mu_x2_MCEM<- c(sd(MCEM50b10[,2]),sd(MCEM100b10[,2]),sd(MCEM300b10[,2]),sd(MCEM600b10[,2]),sd(MCEM1000b10[,2]))

se_mu_x2_MCEM<- c(mean(Result_std_MCEM50b10_Mu[,2]),mean(Result_std_MCEM100b10_Mu[,2]),mean(Result_std_MCEM300b10_Mu[,2]),
                  mean(Result_std_MCEM600b10_Mu[,2]),mean(Result_std_MCEM1000b10_Mu[,2]))

CI_mu_X2_MCEM<- c(counter_mu_x2_50_MCEM_percent,counter_mu_x2_100_MCEM_percent,counter_mu_x2_300_MCEM_percent,
                  counter_mu_x2_600_MCEM_percent,counter_mu_x2_1000_MCEM_percent)

Res4<- as.data.frame(cbind(sample_size1,Ave_mu_x2_MCEM,sd_mu_x2_MCEM,se_mu_x2_MCEM,CI_mu_X2_MCEM))
Res4
################################################################################################################

Ave_MU_X1<- c(mean(EM50b10[,1]),mean(MCEM50b10[,1]),
              mean(EM100b10[,1]),mean(MCEM100b10[,1]),
              mean(EM300b10[,1]),mean(MCEM300b10[,1]),
              mean(EM600b10[,1]),mean(MCEM600b10[,1]),
              mean(EM1000b10[,1]),mean(MCEM1000b10[,1]))


Ave_MU_X1

sd_MU_X1<- c(sd(EM50b10[,1]),sd(MCEM50b10[,1]),
              sd(EM100b10[,1]),sd(MCEM100b10[,1]),
              sd(EM300b10[,1]),sd(MCEM300b10[,1]),
              sd(EM600b10[,1]),sd(MCEM600b10[,1]),
              sd(EM1000b10[,1]),sd(MCEM1000b10[,1]))


sd_MU_X1


se_Mu_X1<- c(mean(Result_std_EM50b10_Mu[,1]),mean(Result_std_MCEM50b10_Mu[,1]),
             mean(Result_std_EM100b10_Mu[,1]),mean(Result_std_MCEM100b10_Mu[,1]),
             mean(Result_std_EM300b10_Mu[,1]),mean(Result_std_MCEM300b10_Mu[,1]),
             mean(Result_std_EM600b10_Mu[,1]),mean(Result_std_MCEM600b10_Mu[,1]),
             mean(Result_std_EM1000b10_Mu[,1]),mean(Result_std_MCEM1000b10_Mu[,1]))

se_Mu_X1
             
CI_X1<- c(counter_mu_x1_50_percent,counter_mu_x1_50_MCEM_percent,
          counter_mu_x1_100_percent,counter_mu_x1_100_MCEM_percent,
          counter_mu_x1_300_percent,counter_mu_x1_300_MCEM_percent,
          counter_mu_x1_600_percent,counter_mu_x1_600_MCEM_percent,
          counter_mu_x1_1000_percent,counter_mu_x1_1000_MCEM_percent)
          
          
          
         
CI_X1             
Res_mu_x1<- as.data.frame(cbind(sample_size,method,Ave_MU_X1,sd_MU_X1,se_Mu_X1,CI_X1))

#################################################################################################
Ave_MU_X2<- c(mean(EM50b10[,2]),mean(MCEM50b10[,2]),
              mean(EM100b10[,2]),mean(MCEM100b10[,2]),
              mean(EM300b10[,2]),mean(MCEM300b10[,2]),
              mean(EM600b10[,2]),mean(MCEM600b10[,2]),
              mean(EM1000b10[,2]),mean(MCEM1000b10[,2]))


Ave_MU_X2


sd_MU_X2<- c(sd(EM50b10[,2]),sd(MCEM50b10[,2]),
             sd(EM100b10[,2]),sd(MCEM100b10[,2]),
             sd(EM300b10[,2]),sd(MCEM300b10[,2]),
             sd(EM600b10[,2]),sd(MCEM600b10[,2]),
             sd(EM1000b10[,2]),sd(MCEM1000b10[,2]))


sd_MU_X2


se_Mu_X2<- c(mean(Result_std_EM50b10_Mu[,2]),mean(Result_std_MCEM50b10_Mu[,2]),
             mean(Result_std_EM100b10_Mu[,2]),mean(Result_std_MCEM100b10_Mu[,2]),
             mean(Result_std_EM300b10_Mu[,2]),mean(Result_std_MCEM300b10_Mu[,2]),
             mean(Result_std_EM600b10_Mu[,2]),mean(Result_std_MCEM600b10_Mu[,2]),
             mean(Result_std_EM1000b10_Mu[,2]),mean(Result_std_MCEM1000b10_Mu[,2]))

se_Mu_X2

CI_X2<- c(counter_mu_x2_50_percent,counter_mu_x2_50_MCEM_percent,
          counter_mu_x2_100_percent,counter_mu_x2_100_MCEM_percent,
          counter_mu_x2_300_percent,counter_mu_x2_300_MCEM_percent,
          counter_mu_x2_600_percent,counter_mu_x2_600_MCEM_percent,
          counter_mu_x2_1000_percent,counter_mu_x2_1000_MCEM_percent)




CI_X2
Res_mu_x2<- as.data.frame(cbind(sample_size,method,Ave_MU_X2,sd_MU_X2,se_Mu_X2,CI_X2))

Res_mu_x1
Res_mu_x2
library(xtable)
xtable(Res_mu_x1,digit=6)
xtable(Res_mu_x2,digits=6)
