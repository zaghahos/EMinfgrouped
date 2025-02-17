## This part is based on section 2.2.3 & 2.2.5 EM algorithm
## In the first part, we have find and define the functions for the estimate of the parameters based on fomulas 19-23,
##and the equations on page 15,
## then in the second part, we iterate over the current estimates of the parameters to find the most updated estimates based on 
##M-step for the parameters


library(tmvtnorm)
library(MASS)
############################
###Parameter Estimation Using EM Algorithm####
#######################################################################
###E-STEP:####
###Moments###
MEXI<- function(Data,mu,sigma){ ## the function to find the first moment of the truncated bivariate normal distributions for 
  ## the rectangles which will be used to calculate the estimate of the covariance matrix 
  ##it has three arguments:  Data: in the form of bivariate grouped data, mu: mean vector, sigma: covariance matrix
  n<- nrow(Data) ## to find the number of all rectangles
  EX<- matrix(rep(0,n*2),ncol=2,nrow=n) ## define an initialized matrix of zeros to put the mean vectors in it
  #s1<- array(rep(0,n*2*2),c(2,2,n))
  for(i in 1:nrow(Data)){## use for loop to do the iterations over all rectangles
    mnts<- mtmvnorm(mean=mu, sigma=sigma, 
                    lower=c(Data[i,1],Data[i,3]),
                    upper=c(Data[i,2],Data[i,4])) ## to find the first moments vectors for each rectangle, we have to use the mtmvnorm
    ## function from the library of tmvtnorm, the lower: is the lower limit of rectangle, the upper: is the upper limit of rectangle
    ## for each separated interval
    EX[i,]<- mnts$tmean ## we select the tmeans (to get only the first moment vectors results from the mnts' output) and put the 
    ## resulted mean vectors for each rectangle on each raw of the matrix
    # s1[,,i]<- mnts$tvar
  }
  return(EX) ## return the first moment vector matrix of the rectagles (the expectation part of mu on page 15 for the rectangles)
  #s1
}
#############################################################
McovXI<- function(Data,mu,sigma){ ## the function to find the covariance matrices using the current theta of the truncated 
  ##bivariate normal distributions for the rectangles which will be used to calculate the estimate of the covariance matrix 
  ##it has three arguments:  Data: in the form of bivariate grouped data, mu: mean vector, sigma: covariance matrix
  #### Attention: in this part, the covariance is only found from the mtmvnorm function for the current state of theta, 
  #### however, in the EM algorithm, to calculate the estimate of the covariance matrix, we have to use the 
  #### updated means/mean vectors to find the new estimate for covariances, so the result of this part, will be used
  #### indirectly, to calculate the estimate of the covariances (If you have any question regarding this part, I can explain it
  #### to Camila, and she will explain it simply in your language to be easily understandable for you )/
  #### also see the detail on the following part related to function with name: EXXest
  n<- nrow(Data) ## to caclulate the number of all rectangles
  #EX<- matrix(rep(0,n*2),ncol=2,nrow=n)
  s1<- array(rep(0,n*2*2),c(2,2,n)) ## define and assigned an initialized array for the covariance matrices using the current 
      ##theta for all the rectangles 
  for(i in 1:nrow(Data)){## use for loop to iterate the calcultion over the rectanlges
    mnts<- mtmvnorm(mean=mu, sigma=sigma,
                    lower=c(Data[i,1],Data[i,3]),
                    upper=c(Data[i,2],Data[i,4])) ## again to find the covariance matrices using the current theta for 
    ##each rectangle, we have to use the mtmvnorm
    ## function from the library of tmvtnorm, the lower: is the lower limit of rectangle, the upper: is the upper limit of rectangle
    ## for each separated interval
    
    s1[,,i]<- mnts$tvar ## assing the tvar output from the above results to the i element of the defined array for the rectangles 
  }
  return(s1) ## return the s1, the array of the covariance matrices (using the current theta) results (This is used in the expectation in the formula of 
  ## the estimate of the covariance matrix on page 15)
  #s
}
#############################################################
###MU:E-step###
MEM<- function(Data,EXi){## the function to calculate the  mean vectors , formula is on page 15,
  ##or the for the bivariate case each element of the vector are eaxctly according to the formulas 19,20 
  ## this function has two arguments: Data in the form of bivariate grouped data, EXi: the first moments matrix of the
  ## rectangles from above
  Freq<- Data[,5] ## The vector of frequencies over the rectangles
  A<- EXi*Data[,5] ## the numerator of the estimate of mean vector (see page 15, also for the bivariate case, the elementes
  ## of the vecor are exactly the numerator of formulas 19, 20
  
  mupd<- apply(A,2,sum)/sum(Data[,5]) ## the estimate of mean vector (see page 15, also for the bivariate case, the elementes
  ## of the vecor are exactly the estimate of the means of formulas 19, 20)
  return(mupd) ## return the mean vector
  
  #mupdPP<- sum(EXi[,1]*Freq)/sum(Freq)
  #mupdCC<- sum(EXi[,2]*Freq)/sum(Freq)
  #muupd<-c(mupdPP,mupdCC)
  #return(muupd)
}


################################################################################
###E(XX) estimate###
EXXest<- function(Data,EXi,ss1){ ## the function to calculate the second moments, using the formulas is on page 15,
  ##or the for the bivariate case each element of the vector are eaxctly according to the formulas 20,21,22 
  ## this function has three arguments: Data in the form of bivariate grouped data, EXi: the first moments matrix of the
  ## rectangles from above, ss1: the array of covariances of the rectangles from above 
  n<- nrow(Data)## the number of all rectangles
  ExxNew<- array(rep(0,n*2*2),c(2,2,n)) ## define and assigne an initialized array, where the elements of the array, are the 
  ## second moments for the rectangles
  for(i in 1:n){## use for loop to do all the computations over all rectangles
    ExxNew[,,i]<- ss1[,,i]+EXi[i,]%*%t(EXi[i,]) ## using the estimated first moments and the covariances (for the current theta)
    ## to find the second moments matrices
    #### Explanation to clarify this part: If you see the formula of the expectation for calculating the estimate of the 
    #### covariance matrix (page 16, also for bivariate: formulas 21,22,23), we have to extend the expectation exactly, 
    #### and E[(x-muupdate)^t * (x-muupdate)]=E[x^t*x]-muupdate^t*E[x]-E[x^t]*muupdate+muupdate^t*muupdate, where the
    #### expectation is done over the current state of theta
    #### in the extended form above, E[x^t*x]: is the second moment that is used in this formula (using the current state of 
    #### theta) and replace in this summation to find the required expectation in the formulas (page 15, formulas 21,22,23)
    #### in the functions defined above, from mtmvnorm output over the rectangles, we can find the first moments and the 
    #### covariances where the covariances are not updated version, however, we use those covariances to find the 
    #### second moments (on the extended form of the expectations) under the current state of theta's (parameters)
    
  }
  return(ExxNew) ## return the array where its elements are the second moment matrics for the current theta 
  
}

######################################################################
Mysig<- function(Data,EXX,EX,UPDMU){## the function to calculate the estimate of the covariances (formula on page 16, also 
  ## for bivariate case, each resulted elements are in the form of formulas 21,22, 23)
  ## this function has four arguments: Data in the form of bivariate grouped data, EX: the first moments matrix of the
  ## rectangles from above, EXX: the array of second moments matrices over the rectangles from above, UPDMU: the updated mu
  ## from the above functions
  n<- nrow(Data) ## the number of rectangles
  s2<- array(rep(0,n*2*2),c(2,2,n)) ## define and assign an array to put the estimate of the covariance matrices
  for(i in 1:n){## use for loop to do the computations over all rectangles
    s2[,,i]<- (EXX[,,i]-UPDMU%*%t(EX[i,])-EX[i,]%*%t(UPDMU)+(UPDMU%*%t(UPDMU)))*Data[i,5] ## Here, we use the extended form 
    ## of the expectations for each rectangle and multiply it by the frequency of that relevant rectangle 
    
  }
  
  #ssx1N<- sum(s2[1,1,])/sum(Data[,5])  
  #ssx2N<- sum(s2[2,2,])/sum(Data[,5])
  #covx1x2N<- sum(s2[1,2,])/sum(Data[,5])
  #rhoN<- covx1x2N/(sqrt(ssx1N)*sqrt(ssx2N))
  #ssupd<- apply(s2,c(1,2),sum)/sum(Data[,5])
  
  SigmaNew<- apply(s2,c(1,2),sum)/sum(Data[,5]) ## use the fuction apply to find the summation in the numerator of the 
  ## formulas and devide it by total frequency in the bivariate grouped data
  #return(ssupd)
  
  return(SigmaNew) ## Return the estimate of the covariance matrix
}



#DataBiv<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton-Bivariate.csv")

#DataBiv[1:14,1]<- -Inf
#DataBiv[141:154,2]<- Inf
#DataBiv[c(1,15,29,43,57,71,85,99,113,127,141),3]<- -Inf
#DataBiv[c(14,28,42,56,70,84,98,112,126,140,154),4]<- Inf
#DataBiv

#mydat<- DataBiv

#mu1 <- c(67,67)
#sigma1 <- matrix(c(3.2, 2.2, 2.2, 6.2), 2, 2)



#Mysig(Data=mydat,
 #     EXX=EXXest(Data=mydat,
  #               EXi=MEXI(Data=mydat,mu=mu1,sigma=sigma1),
   #              ss1=McovXI(Data=mydat,mu=mu1,sigma=sigma1)),
    #  EX=MEXI(Data=mydat,mu=mu1,sigma=sigma1),
     # UPDMU=MEM(Data=mydat,EXi=MEXI(Data=mydat,mu=mu1,sigma=sigma1)))

########################################################################
########################################################################
###M-STEP:###
### M-Step###
EMBivG<- function(Data,mu_init,sigma_init,maxit=1000,tol1=1e-4,tol2=1e-5){## The EM function, to find the estimate of 
  ## the mean vector/means and covariance matrix
  ## This function has 6 elements: 
  ## Data: the grouped bivariate data, mu_iniit: the initial values of mu vector, sigma_init: the initial covariance matrix,
  ## tol1: stopping criteria for means, tol2: stopping criteria for variances
  flag<- 0 ## the initial value of flag which will be used in the stopping criteria conditions
  M_cur<- mu_init ## assing the current value of mean equal to the initial point 
  S_cur<- sigma_init ## assing the current value of covariance equal to the initial point 
  #theta_cur<- c(Mu_cur,S_cur)
  #iter<- rep(0,maxit) 
  #Svec<- rep(0,maxit)
  #Mvec<- rep(0,maxit)
  for (i in 1:maxit){ ## define a for loop to do update the estimates up to the maximum of 1000 times
    #print(paste("Iteration number=", i))
    #cur<- c(Mu_cur,S_cur)
    MU_cur<- M_cur ## the current state of mean calling inside the for loop
    SS_cur<- S_cur ## the current state of covariance calling inside the for loop
    
    #   MuNew<- Mu(Data=mydat,Y=ZMCEM1(Data=mydat,mu=MU_cur,sigma=SS_cur) )
    
    #  SigmaNew<- Sigma1(Data=mydat,Y=ZMCEM1(Data=mydat,mu=MU_cur,sigma=SS_cur),
    #                   M=Mu(Data=mydat,Y=ZMCEM1(Data=mydat,mu=MU_cur,sigma=SS_cur) ) )
    
    
    MuNew<- MEM(Data=mydat,EXi=MEXI(Data=mydat,mu= MU_cur,sigma=SS_cur)) ## updating step for mean
    
    
    SigmaNew<- Mysig(Data=mydat,
                   EXX=EXXest(Data=mydat,
                              EXi=MEXI(Data=mydat,mu=MU_cur,sigma=SS_cur),
                              ss1=McovXI(Data=mydat,mu=MU_cur,sigma=SS_cur)),
                   EX=MEXI(Data=mydat,mu=MU_cur,sigma=SS_cur),
                   UPDMU=MEM(Data=mydat,EXi=MEXI(Data=mydat,mu=MU_cur,sigma=SS_cur))) ## updating step for covariance matrix
    
    #ssx1N<- sum(SigNew[1,1,])/sum(Data[,5])
    #ssx2N<- sum(SigNew[2,2,])/sum(Data[,5])
    #covx1x2N<- sum(SigNew[1,2,])/sum(Data[,5])
    #rhoN<- covx1x2N/(sqrt(ssx1N)*sqrt(ssx2N))
    
    
    #SigmaNew<- matrix(c(ssx1N,covx1x2N,covx1x2N,ssx2N),ncol=2,nrow=2)
    
    #print(MuNew)
    #print(SigmaNew) 
    
    Mu_new<- MuNew ## reassign the new estimate of mean (renaming)
    Sigma_new<- SigmaNew ## reassign the new estimate of covariance (renaming)
    
    #new_stepM<- Mu_new
    #new_stepS<- Sigma_new
    
    diff1<- MU_cur-Mu_new ## find the differences between the current mean vector and the new estimate of mean vector
    diff2<- SS_cur-Sigma_new ## find the differences between the current covariance matrix and the new estimate of covariance matrix
    #print(diff1)
    #print(diff2)
    
    D1<- abs((diff1[1]+diff1[2])/2) ## calculate the absolute value of the mean of differences in means
    #D2<- (abs(diff2[1,1])+abs(diff2[1,2])+abs(diff2[2,1])+abs(diff2[2,2]))/4
    D2<- abs((diff2[1,1]+diff2[1,2]+diff2[2,1]+diff2[2,2])/4) ## calculate the absolute value of the mean of differences 
    ##in the elements of the covariance matrix
    
    #print(D1)
    #print(D2)
    if(D1<tol1 & D2<tol2){flag<- 1 ;break} ## define the condition for the absolute values of the differences according to the 
    ## stopping criteria, continue updating until the values of the parameter estimates meet both conditions, then stop and 
    ## take the final new estimates as the updated estimates of the parameters 
    
    M_cur<- Mu_new ## update the current value of mean vector with the new estimate
    
    S_cur<- Sigma_new ## update the current value of covariance matrix with the new estimate
    
    #iter[i]<- i
    
  }
  
  if(!flag) warning("Didn't Converge \n") ## the condition if there is no convergence in the stopping rule
  #list(paste("updated Mu of the EM=" , M_cur),paste("Updated Sigma of the EM=" ,S_cur),
  #    paste("iteration number=",i))
  #print(M_cur,S_cur)
  #print(S_cur)
  #update=c(M_cur,(S_cur))
  #return(update)
  MX1<- M_cur[1] ## assign the last & final cur to the final estimate for mean X1
  MX2<- M_cur[2] ## assign the last & final cur to the final estimate for mean X2
  SX1<- (S_cur[1,1]) ## assign the last & final cur to the final estimate for var X1
  SX2<- (S_cur[2,2]) ## assign the last & final cur to the final estimate for var X2
  RhoX1X2<- S_cur[1,2]/(sqrt(SX1*SX2)) ## assign the last & final cur to the final estimate for correlation
  #print(iter)
  updatedest<- c(MX1,MX2,SX1,SX2,RhoX1X2) ## the vector of the updated estimates of the parameters
  names(updatedest)<- c("mux1","mux2","Sx1","Sx2","rhox1x2") ## name the updated estimates
  return(updatedest) ## return the updated estimates of the parameters
  
}
#############################################################################
### GALTON DATA
DataBiv<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton-Bivariate.csv") ## Read the bivariat grouped data

DataBiv[1:14,1]<- -Inf ## ## set the first lower 14 limits of parent data equals -inf (you can see the table of 
##data to understand about the intervals and see the gropued form of data using that table ) 
DataBiv[141:154,2]<- Inf ## set the last 14 limilts of children equals -inf
DataBiv[c(1,15,29,43,57,71,85,99,113,127,141),3]<- -Inf ##set the last upper limit of the last intervals (parent data) equals inf
DataBiv[c(14,28,42,56,70,84,98,112,126,140,154),4]<- Inf ##set the last upper limit of the last intervals (children data) equals inf
#DataBiv

mydat<- DataBiv ## rename the data to match with the names inside the function

mu1 <- c(67,67) ## initial value for mean vector
sigma1 <- matrix(c(3.2, 2.2, 2.2, 6.2), 2, 2) ## initial value for covariance matrices



myout1<- EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-5) ## run the EMBivG funtion to see the 
## estimate of the parameters using the EM algorithm 

myout1 ## see the output above

##########################################################################
###INITIAL VALUES####
mu1 <- c(67,67) ## initial value for mean vector
sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2) ## initial value for covariance matrices
##########################################################################
####n=50####
myoutsim50<- matrix(rep(0*5*30),ncol=5,nrow=30) ## define and assign a matrix to put the resulted estimates of 
## of the parameters using EM algorithm for 30 simulated datasets
colnames(myoutsim50)<- c("mux1","mux2","Sx1","Sx2","rhox1x2") ## name the columns of the matrix
for(i in 1:nrow(myoutsim50)){ ## use for loop to find the estimates for all 30 simulated datasets 
  mydat<- simulateddata50[,,i]## assign the new mydat (each of the 30 datasets) to mydat
  myoutsim50[i,]<- EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-5)## implement the EMBivG function
  ## and put the results on the relevant raw of the matrix
  
}
#colMeans(myoutsim50) 