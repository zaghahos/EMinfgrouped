### To calculate exact MLE we are using Eq (16), where F(a,b) is bivariate normal cumulative function
##To be able to use the cumulative distribution fuction(cdf) of bivariate normal we need to use the followig libraries:
library(tmvtnorm)
library(MASS)
## Then we have to define the log-likelihood function according to Eq (16), and for the log-likelihhod funcion, in order 
## to find the maximum likelihood of the fuction (MLE) we have to use the optimization function, where we use the *nlm* 
## function for the bivariate.


###Parameters Estimation Using nlm Function###
##Exact MLE Function for Bivariate Normal (Grouped data)###


Billik<- function(Data,theta){ ## defining the log-likelihood function, The function has two arguments, Data: bivariate grouped
  ##data and theta: the value of the parameters (it is a vector as the input) 
  ### Theta: in the vector of all parameters, the first d elements of it are the means
  ### the other elements from d+1: length(theta), are d: Varinces (s11,s22,...,sdd), and then the covariances from the 
  ### upper part of the covariance matrix: s12,s13,...,s1d,s2(2+1),...,s2d,s3(3+1),...,s3d,...,s(d-1)d
  ### That means it has a form of: theta=(mu1,...,mud,Var1,...,Vard,cov(1,(1+1)),cov(1,d),cov(2,(2+1)),...,cov(2,d),...,cov(d-1,d))
  
  d<- (ncol(Data)-1)/2 ##To find the number of variables in our grouped data
  
  mu<- theta[c(1:d)] ## Build the vector of mu
      ###To Build the matrix of sigma
  sig<- matrix(rep(0,d*d),ncol=d) ## Create an initial zero matrix
  
  subtheta<- theta[-c(1:d)] ## Step 1: sub select the elements for sigma from theta
  
  diag(sig)<- subtheta[1:d] ## step2: replace the diagonal element of sig matrix with the first d elements of subtheta 
  ###(subtheta is the sub vector o theta after removing its first d )
  
  subtheta1<- subtheta[-c(1:d)] ## Step 3: Create another subvector from the new (above) subselected elements to replace the 
  ###off-diagonal elements of sig matrix
  
  sig[lower.tri(sig)] <- sig[upper.tri(sig)]<- subtheta1 ## step 4: replace the lower & upper triangle of the sig matrix 
  ### with the elements of subtheta1, and now the sig matrix which is the sigma matrix is created from the vector of theta 
  ### (the parameters)
  #print(sig)

  ind<- ncol(Data)-1 ## select the number of columns of data needed to create the lower & upper bounds for the variables, to 
  ### select the index of which columns shuold be used for the lower bounds and which ones are for the upper bound
  #ind
  
  lowerb<- NULL ## The initialized lowerbound incides of the variables (the column numbers that contain the lowerbound of the variables)
  upperb<- NULL ##The initialized upperbound incides of the variables (the column numbers that contain the upperbound of the variables)
  for (i in 1:ind){ ## set the loop to select the columns for the variables lower & upper bounds
    if (i%%2!=0) {lowerb<- append(lowerb,i)} ## to select the lower bounds columns
    else {upperb<- append(upperb,i)} ## to select the upper bound columns
    
  }
  
  #print(lowerb)
  #print(upperb)
  Data1<- as.matrix(Data) ## As pmvnorm function, needed the numerics for the lower & upper bounds and the results of the previous
  ### lines of code were in the format of lists, so we have to transform data to the matrix to have the results as numeric/vectors
  ### otherwise, we will get error when running the pmvnorm function.
  b<- 0 ### ## define initial assignment for the summation in Eq (16)
  for(i in 1:nrow(Data1)){ ## define a loop to do the summation over the rectangles/surfaces
    a<- pmvnorm(lower=Data1[i,c(lowerb)],upper=Data1[i,c(upperb)],mean=mu,sigma=sig) ##Calculating
    ## P_ij(theta) [inside paranthesis of log in Eq (16)], lower: lower bound rectangle, upper: upper bound rectanle, these 
    ## should be computed over all of the rectangles
    if (a[1]==0) {a<- 1e-03} ## define this condition to avoid zero results which lead to error for log
    else {a<- a[1]}
    b<- b+Data[i,ncol(Data)]*log(a) ## multiply the log(P_ij(theta)) by the frequencies and sum over all rectangles (Eq 16)
    #print(b)
  }
  #print(-b)
  ##The joint likelihood
  return(-b) ## here we have to return minus b,as for nlm function, in order to maximize the log-likelihood, we have to 
  ## find the optimum piont of -b. (why? because nlm gives us the minimum, to get the maximum, we use -b)
}
#############################################
### Initialize theta for d=2 & d=3###
#mu1 <- c(67,67)
#sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
#2.16/(sqrt(3.1*6.05))

#mu2<- c(67,67,67)
#sigma2<- matrix(c(3.1,2.2,1.6,2.2,4.3,2.45,1.6,2.45,4.8),3,3)


#theta1<- c(mu1,as.vector(sigma1))
#theta2<- c(mu2,as.vector(sigma2))


#MLEExactBiv<- optim(theta<- c(65,65,3,6,0.5),fn=Billik,Data=simulateddata1000b10[,,1],method="Nelder-Mead") ## use optim function with the same arguments described
#MLEExactBiv$par ## to see the Exact MLE of the parameters


#MLEExactBiv<- optim(theta<- c(67,67,3.10,6.05,0.4987641),fn=Billik,Data=simulateddata1000b10[,,1],method="Nelder-Mead") ## use optim function with the same arguments described
#MLEExactBiv$par ## to see the Exact MLE of the parameters


#MLEExact<- nlm(Billik,Data=simulateddata1000b10[,,i],theta<- c(65,65,3,6,0.5),hessian=TRUE)$estimate

#MLEExact<- nlm(Billik,Data=simulateddata1000b10[,,1],theta<- c(67,67,3.10,6.05,0.4987641),hessian=TRUE)$estimate

#MLEExact


######JOHN:TO PRODuCE THE MATRIX OF ESTIMATE FOR n=1000 &Binned=10 
MLEBINNED1000b10<- matrix(rep(0,5*10),ncol=5)
colnames(MLEBINNED1000b10)<- c("mux1","mux2","Vx1","Vx2","rhox1x2")

for(i in 1:nrow(MLEBINNED1000b10)){
  MLEBINNED1000b10[i,]<- nlm(Billik,Data=simulateddata1000b10[,,i],theta<- c(67,67,3.10,6.05,0.4987641),hessian=TRUE)$estimate
  #####JOHN: NOTE THAT the argument Data should be replaced properly
}
MLEBINNED1000b10


OUTMLEBINNED1000b10<- write.csv(MLEBINNED1000b10,"C:/Users/sh_za/Desktop/Results/Bivariate/MLE Binned/Sample size1000/MLEBINNED1000b10.csv")


