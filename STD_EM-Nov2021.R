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
dataEM_n1000_bin10
####################################################################################################################

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
  for(i in 1:nrow(Data)){ 
    mnts<- mtmvnorm(mean=mu, sigma=sigma, 
                    lower=Data1[i,c(lowerb)],
                    upper=Data1[i,c(upperb)]) 
    EX[i,]<- mnts$tmean 
  }
  EX_minus_MU <- matrix(rep(0,n*d),ncol=d,nrow=n) 
  EX_minus_MU<- EX-mu
  #print(EX_minus_MU)
  
  #print(EX_minus_MU[1,]%*%solve(sigma))
  S1 <- matrix(rep(0,n*d),ncol=d,nrow=n) 
  
  for(i in 1:nrow(S1)){
    S1[i,]<- (EX_minus_MU[i,])%*%solve(sigma)
  }
  #print(S1)
  #S1<- t(EX_minus_MU)%*%solve(sigma)
  #print(S1[1,])
  #print(S1[1,]%*%t(S1[1,]))
  S1_S1T<- array(rep(0,d*d*n),c(d,d,n))  
  #print(S1_S1T)
  for (i in 1:n){
    S1_S1T[,,i]<- S1[i,]%*%t(S1[i,])
  }
  #print(S1_S1T)
  Inf_j_m<- array(rep(0,d*d*n),c(d,d,n))
  for(i in 1:n){
    Inf_j_m[,,i]<- S1_S1T[,,i]*Data[i,5]
  }
  #print(Inf_j_m)
  Inf_mat<- apply(Inf_j_m,c(1,2),sum)
  #print(Inf_mat)
  cov_Mu_EM<- solve(Inf_mat)
  Se_mu_EM<- c(sqrt(cov_Mu_EM[1,1]),sqrt(cov_Mu_EM[2,2]))
  
  #print(cov_Mu_EM)
  #print(Se_mu_EM)
  
  return(Se_mu_EM)
  
  
}

#mu1 <- c(67,67)
#sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
#EM1000b10

#mu_EM1000b10[1,]
#sigma_EM1000b10[,,1]

result1<- Std_Mu_EM(Data=simulateddata1000b10[,,1],mu=mu_EM1000b10[1,],sigma=sigma_EM1000b10[,,1])
                    
result1

simulateddata1000b10

result1[1]*1.96



a<- matrix(c(1,2,3,4),ncol=2)
a

a^2

fm1
a<- c(1,2)
a
t(a)
t(a)%*%sigma1
dim(t(a))
dim(a)
t(fm1[1,])%*%solve(sigma1)  
solve(sigma1)
am
MEXI
#############################################################
Std_Sigma_EM<- function(Data,mu,sigma){ 
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
  
  s1<- array(rep(0,n*d*d),c(d,d,n)) 
  for(i in 1:nrow(Data)){ 
    mnts<- mtmvnorm(mean=mu, sigma=sigma,
                    lower=Data1[i,c(lowerb)],
                    upper=Data1[i,c(upperb)]) 
    s1[,,i]<- mnts$tvar 
  }
  
  
  S1_cov<- array(rep(0,n*d*d),c(d,d,n))
  for(i in 1:nrow(Data)){
    S1_cov[,,i]<- sigma-s1[,,i]
    }
  #print(S1_cov)
  S1_cov_square<- S1_cov^2
  S1_cov_square1<- S1_cov_square/4
  
  S1_cov_Freq<- array(rep(0,n*d*d),c(d,d,n))
  for(i in 1:nrow(Data)){
    S1_cov_Freq[,,i]<- S1_cov_square1[,,i]*Freq[i]
  }
  
  #print(S1_cov_Freq)
  I_emp_cov<- apply(S1_cov_Freq,c(1,2),sum)
  std_emp_cov<- 1/(sqrt(I_emp_cov))
  #print(std_emp_cov)
  return(std_emp_cov) 
  
}


Result2<- Std_Sigma_EM(Data=simulateddata1000b10[,,1],mu=mu_EM,sigma=sigma_EM)  
Result2  

EM1000b10

#Result_std_EM_mu<- matrix(rep(0,nrow(EM1000b10)*2),ncol=2)
#Result_std_EM_mu

#Result_std_EM_mu<- Std_Mu_EM(Data=simulateddata1000b10[,,1],mu=mu_EM,sigma=sigma_EM)   
  
#Result_std_EM_sigma<- Std_Sigma_EM(Data=simulateddata1000b10[,,1],mu=mu_EM,sigma=sigma_EM)


###################################################################################################################
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

###################################################################################################################
Result_std_EM1000b10_mu<- matrix(rep(0,nrow(EM1000b10)*2),ncol=2)


for(i in 1:nrow(EM1000b10)){
  Result_std_EM1000b10_mu[i,]<- Std_Mu_EM(Data=simulateddata1000b10[,,1],mu=mu_EM1000b10[i,],sigma=sigma_EM1000b10[,,i])   
  
}

Result_std_EM1000b10_mu


Result_std_EM1000b10_sigma<- array(rep(0,2*2*nrow(EM1000b10)),c(2,2,nrow(EM1000b10)))
for(i in 1:nrow(EM1000b10)){
  Result_std_EM1000b10_sigma[,,i]<- Std_Sigma_EM(Data=simulateddata1000b10[,,1],mu=mu_EM1000b10[i,],sigma=sigma_EM1000b10[,,i])
  
  }
Result_std_EM1000b10_sigma
##################################################################################################################
CI_mu_x1<- matrix(rep(0,2*nrow(EM1000b10)),ncol=2)
colnames(CI_mu_x1)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1)){
  CI_mu_x1[i,]<- c(EM1000b10[i,1]-(qnorm(0.975)*Result_std_EM1000b10_mu[i,1]),
                   EM1000b10[i,1]+(qnorm(0.975)*Result_std_EM1000b10_mu[i,1]))
  
}

CI_mu_x1

##############################################################################################
CI_mu_x2<- matrix(rep(0,2*nrow(EM1000b10)),ncol=2)
colnames(CI_mu_x2)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_mu_x1)){
  CI_mu_x2[i,]<- c(EM1000b10[i,2]-(qnorm(0.975)*Result_std_EM1000b10_mu[i,2]),
                   EM1000b10[i,2]+(qnorm(0.975)*Result_std_EM1000b10_mu[i,2]))
  
}

CI_mu_x2


#Result_std_EM1000b10_sigma
###############################################################################################
CI_sigma_x1<- matrix(rep(0,2*nrow(EM1000b10)),ncol=2)
colnames(CI_sigma_x1)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_sigma_x1)){
  CI_sigma_x1[i,]<- c(EM1000b10[i,3]-(qnorm(0.975)*Result_std_EM1000b10_sigma[1,1,i]),
                   EM1000b10[i,3]+(qnorm(0.975)*Result_std_EM1000b10_sigma[1,1,i]))
  
}

CI_sigma_x1


CI_sigma_x2<- matrix(rep(0,2*nrow(EM1000b10)),ncol=2)
colnames(CI_sigma_x2)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_sigma_x2)){
  CI_sigma_x2[i,]<- c(EM1000b10[i,6]-(qnorm(0.975)*Result_std_EM1000b10_sigma[2,2,i]),
                      EM1000b10[i,6]+(qnorm(0.975)*Result_std_EM1000b10_sigma[2,2,i]))
  
}

CI_sigma_x2

CI_sigma_x1X2<- matrix(rep(0,2*nrow(EM1000b10)),ncol=2)
colnames(CI_sigma_x1X2)<- c("lower_bound","upper_bound")
for(i in 1:nrow(CI_sigma_x1X2)){
  CI_sigma_x1X2[i,]<- c(EM1000b10[i,4]-(qnorm(0.975)*Result_std_EM1000b10_sigma[1,2,i]),
                      EM1000b10[i,4]+(qnorm(0.975)*Result_std_EM1000b10_sigma[1,2,i]))
  
}

CI_sigma_x1X2
CI_mu_x1

true_mu<- c(68,68)
true_cov<- matrix(c(3,2,2,6),ncol=2)
true_cov
counter_mu_x1<- 0
for(i in 1:nrow(EM1000b10)){
  if (true_mu[1]>= CI_mu_x1[i,1]& true_mu[1]<= CI_mu_x1[i,2]) {counter_mu_x1<- counter_mu_x1+1}
}

counter_mu_x1
counter_mu_x1_percent<- counter_mu_x1/500
counter_mu_x1_percent
#####################################################################################################
counter_mu_x2<- 0
for(i in 1:nrow(EM1000b10)){
  if (true_mu[2]>= CI_mu_x2[i,1]& true_mu[2]<= CI_mu_x2[i,2]) {counter_mu_x2<- counter_mu_x2+1}
}

counter_mu_x2
counter_mu_x2_percent<- counter_mu_x2/500
counter_mu_x2_percent
######################################################################################################
counter_sigma2_x1<- 0
for(i in 1:nrow(EM1000b10)){
  if (true_cov[1,1]>= CI_sigma_x1[i,1]& true_cov[1,1]<= CI_sigma_x1[i,2]) {counter_sigma2_x1<- counter_sigma2_x1+1}
}

counter_sigma2_x1
counter_sigma2_x1_percent<- counter_sigma2_x1/500
counter_sigma2_x1_percent
###########################################################################################################
counter_sigma2_x2<- 0
for(i in 1:nrow(EM1000b10)){
  if (true_cov[2,2]>= CI_sigma_x2[i,1]& true_cov[2,2]<= CI_sigma_x2[i,2]) {counter_sigma2_x2<- counter_sigma2_x2+1}
}

counter_sigma2_x2
counter_sigma2_x2_percent<- counter_sigma2_x2/500
counter_sigma2_x2_percent
#################################################################################################################
counter_sigma_x1x2<- 0
for(i in 1:nrow(EM1000b10)){
  if (true_cov[1,2]>= CI_sigma_x1X2[i,1]& true_cov[1,2]<= CI_sigma_x1X2[i,2]) {counter_sigma_x1x2<- counter_sigma_x1x2+1}
}


counter_sigma_x1x2
counter_sigma_x1x2_percent<- counter_sigma_x1x2/500
counter_sigma_x1x2_percent
#######################################################################################################################
a<- matrix(1:6,ncol=2)
a
a^2
t(a)
t(a)%*%a
a%*%t(a)
diag(3)
1*2: score function
2*2: score funnction
2*2
5*5
3
1*3
3*3
3*3