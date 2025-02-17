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
  print(sig)

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
mu1 <- c(67,67)
sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)


mu2<- c(67,67,67)
sigma2<- matrix(c(3.1,2.2,1.6,2.2,4.3,2.45,1.6,2.45,4.8),3,3)


theta1<- c(mu1,as.vector(sigma1))
theta2<- c(mu2,as.vector(sigma2))

################################
####MLE EXACT#############################
## TO check the function on a bivariate grouped data: Galton data
###DATA###
DataBiv<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton-Bivariate.csv") ## read the data in the form of grouped data (bivariate)

DataBiv[1:14,1]<- -Inf ## set the first lower 14 limits of parent data equals -inf (you can see the table of 
##data to understand about the intervals and see the gropued form of data using that table )   
DataBiv[141:154,2]<- Inf ## set the last 14 limilts of children equals -inf 
DataBiv[c(1,15,29,43,57,71,85,99,113,127,141),3]<- -Inf ##set the last upper limit of the last intervals (parent data) equals inf  
DataBiv[c(14,28,42,56,70,84,98,112,126,140,154),4]<- Inf ##set the last upper limit of the last intervals (children data) equals inf
#DataBiv
#nrow(DataBiv)
#MLEExact<- nlm(Billik,Data=DataBiv,theta<- c(65,65,3,6,2.12132),hessian=TRUE) ## use nlm function with the same arguments described
##above, now for the Galton data and assign it to the MLEEact
#MLEExact$estimate ## to see the Exact MLE of the parameters 

MLEExactGalton<- optim(theta<- c(65,65,3,6,2.12132),fn=Billik,Data=DataBiv,method="Nelder-Mead") ## use optim function with the same arguments described
MLEExactGalton$par ## to see the Exact MLE of the parameters

rm(DataBiv)
##############################################################################################
###MLE EXACT For Simulted Data####
###data4: n=1000##########################
### Here is a data set simulated for n=1000 , and bivariate grouped data
### Please just run the codes to see the simulated data and then implement the Exact MLE code, like the one for the Bivariate
### Galton, on these simulated data set. How data are simulated is explained in other files.

## TO check the function on a simulation of bivariate grouped data set
### Implement the function on a smulation of Bivariate grouped data 
### Simulate the data d=2###
mm<- c(68,68)
ss<- matrix(c(3,2,2,6),2,2)

ssdata<- matrix(rep(0,1000*2),c(1000,2))

ssdata<- mvrnorm(n=1000,mm,ss)

#ssdata

x<- rep(0,1000)
y<- rep(0,1000)


x<- ssdata[,1]
y<- ssdata[,2]

Freqtable1<- array(rep(0,10*10),c(10,10))

x.cut<- cut(x,breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf))
y.cut<- cut(y,breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf))
Freqtable1<- table(x.cut,y.cut)


simulateddata2<- array(rep(0,5*100),c(100,5))
lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10)
lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10))
upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10)
upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10))
simulateddata2[,1]<- lower.x
simulateddata2[,3]<- lower.y
simulateddata2[,2]<- upper.x
simulateddata2[,4]<- upper.y
simulateddata2[,5]<- c(Freqtable1)

#simulateddata2
##########################################################################################
### run the unction###

MLEExactBiv<- optim(theta<- c(65,65,3,6,2.12132),fn=Billik,Data=simulateddata2,method="Nelder-Mead") ## use optim function with the same arguments described
MLEExactBiv$par ## to see the Exact MLE of the parameters

########################################################################################

rm(simulateddata2)
##############################################################################################
### Implement the function on a smulation of Trivariate grouped data 
### Simulate the data d=3###
#####################################################################################
mm<- c(68,68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,1.5,2,4,2.5,1.5,2.5,5),3,3) ## The covariance matrix for the simulation

ssdata<- array(rep(0,3000*3),c(3000,3)) 
ssdata<- mvrnorm(n=3000,mm,ss) 

x<- rep(0,3000) 
y<- rep(0,3000) 
z<- rep(0,3000) 

x<- ssdata[,1] 
y<- ssdata[,2] 
z<- ssdata[,3] 

Freqtable1<- array(rep(0,8*8*8),c(8,8,8)) 
x.cut<- cut(x,breaks=c(-Inf,65,66,67,68,69,70,71,Inf))  
y.cut<- cut(y,breaks=c(-Inf,63,65,67,69,71,73,75,Inf)) 
z.cut<- cut(z,breaks=c(-Inf,64,65.5,67,68.5,70,71.5,73,Inf)) 

Freqtable1<- table(x.cut,y.cut,z.cut) 


simulateddata3<- array(rep(0,7*512),c(512,7)) 
lower.x<- rep(c(-Inf,65,66,67,68,69,70,71),64) 
ly<- c(rep(-Inf,8),rep(63,8),rep(65,8),rep(67,8),rep(69,8),rep(71,8),rep(73,8),rep(75,8)) 
lower.y<- rep(ly,8) 
lower.z<- c(rep(-Inf,64),rep(64,64),rep(65.5,64),rep(67,64),rep(68.5,64),rep(70,64),rep(71.5,64),rep(73,64)) 

upper.x<- rep(c(65,66,67,68,69,70,71,Inf),64) 
uy<- c(rep(63,8),rep(65,8),rep(67,8),rep(69,8),rep(71,8),rep(73,8),rep(75,8),rep(Inf,8)) 
upper.y<- rep(uy,8) 
upper.z<- c(rep(64,64),rep(65.5,64),rep(67,64),rep(68.5,64),rep(70,64),
            rep(71.5,64),rep(73,64),rep(Inf,64)) 

simulateddata3[,1]<- lower.x 
simulateddata3[,3]<- lower.y 
simulateddata3[,5]<- lower.z  
simulateddata3[,2]<- upper.x 
simulateddata3[,4]<- upper.y 
simulateddata3[,6]<- upper.z 
simulateddata3[,7]<- c(Freqtable1) 
#simulateddata3

#############################################################################
### Run the function###
#######n=3000#####################
## TO check the function on a simulation of trivariate grouped data set
MLEExacttri<- optim(theta<- c(67,67,67,3.1,4.2,5.5,1.8,1.7,2.3),fn=Billik,Data=simulateddata[,,15],method="Nelder-Mead") ## use nlm function with the same arguments described
MLEExacttri$par
