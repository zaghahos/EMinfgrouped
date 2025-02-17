thetaupd<- c(output1000m15EM[1],sqrt(output1000m15EM[2]))
thetaupd

Mustd<- function(thetaupd,bl,bu,Freq){ 
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
  #print((Wj-thetaupd[1]))
  #Mstdj<- Freq*((Wj-thetaupd[1])^2)
  Mstdj<- Freq*(WjN^2)
  #print(Mstdj)
  #print(Mstdj1)
  #print(Mstdj)
  Mstd<- (sum(Mstdj))/(thetaupd[2]^4)
  #results<- list(paste("std for Mu=" , 1/Mstd),Wj,WjN)
  return(1/Mstd) 
  
}

Mustd(thetaupd=c(output1000m15EM[1],sqrt(output1000m15EM[2])),
      bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i])

################################################################################################
Varstd<- function(thetaupd,bl,bu,Freq){ 
  Wj<- rep(0,length(bl)) 
  Wj2<- rep(0,length(bl)) 
  
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
  #Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i])     
  #}
  #print(Wj)
  dinom1<- NULL
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  for(i in 1:length(bl)){ 
    dinom1[i]<- (pnorm(bstar[i])-pnorm(astar[i])) 
    if(dinom1[i]==0) {Wj2[i]<- (thetaupd[2]^2-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/0.0001)
    +(thetaupd[1]^2)-
      (2*thetaupd[2]*(thetaupd[1])*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001))}
    
    else{Wj2[i]<- (thetaupd[2]^2-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/dinom1[i])
    +(thetaupd[1])^2-
      (2*thetaupd[2]*(thetaupd[1])*((dnorm(bstar[i])-dnorm(astar[i]))/dinom1[i]))}
    
    
    # Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    #+(muupdate-theta[1])^2+
    # (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i]))))
    
  }
  #print(Wj)
  #print(Wj2)
  
  V1<- (Wj2-2*thetaupd[1]*Wj+thetaupd[1]^2)
  V2<- (1-V1/(thetaupd[2]^2))^2
  V3<- Freq*V2
  V4<- sum(V3)
  Vstd<- (1/4)*(V4/(thetaupd[2]^4))
  #print(V4)
  
  #SS<- sum(Bj*Freq)/sum(Freq) 
  #S<- sqrt(SS)
  return(1/Vstd) 
}


Varstd(thetaupd=c(output1000m15EM[1],sqrt(output1000m15EM[2])),
      bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i])
#############################################################################################################

VarstdNew<- function(thetaupd,bl,bu,Freq){ 
  #Wj<- rep(0,length(bl)) 
  Ej2<- rep(0,length(bl)) 
  
  astar<- rep(0,length(bl)) 
  bstar<- rep(0,length(bl)) 
  #Mstdj<- rep(0,length(bl))
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-thetaupd[1])/thetaupd[2] 
    astar[i]<- (bl[i]-thetaupd[1])/thetaupd[2] 
    
  }
  dinom<- NULL
  
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  
  for(i in 1:length(bl)){ 
    dinom[i]<- (pnorm(bstar[i])-pnorm(astar[i])) 
    if(dinom[i]==0) {Ej2[i]<- (bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/0.0001}

    
    else{Ej2[i]<- (bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/dinom[i]}
      

    # Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    #+(muupdate-theta[1])^2+
    # (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i]))))
    
  }
  #print(Ej2)
  V1<- Ej2^2
  V2<- Freq*V1
  V3<- sum(V2)
  #print(V2)
  Vstd<- (1/4)*(V3/(thetaupd[2]^4))

  #SS<- sum(Bj*Freq)/sum(Freq) 
  #S<- sqrt(SS)
  return(1/Vstd) 
}


VarstdNew(thetaupd=c(output1000m15EM[1],sqrt(output1000m15EM[2])),
       bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i])
#############################################################################################################