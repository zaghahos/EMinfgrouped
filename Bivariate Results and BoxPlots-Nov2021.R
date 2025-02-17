library(xtable)
##########################################################################################################

dataMLE_n50_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MLE Oct2021/MLE_n50_bin10.csv")

dataMLE_n100_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MLE Oct2021/MLE_n100_bin10.csv")

dataMLE_n300_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MLE Oct2021/MLE_n300_bin10.csv")

dataMLE_n600_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MLE Oct2021/MLE_n600_bin10.csv")

dataMLE_n1000_bin10<- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/MLE Oct2021/MLE_n1000_bin10.csv")

dataMLE_n50_bin10<- dataMLE_n50_bin10[,-1]

dataMLE_n100_bin10<- dataMLE_n100_bin10[,-1]

dataMLE_n300_bin10<- dataMLE_n300_bin10[,-1]

dataMLE_n600_bin10<- dataMLE_n600_bin10[,-1]

dataMLE_n1000_bin10<- dataMLE_n1000_bin10[,-1]

dataMLE_n50_bin10

dataMLE_n100_bin10

dataMLE_n300_bin10

dataMLE_n600_bin10


dataMLE_n1000_bin10
#################################################################################
muxMLE50b10<- sqrt(mean((dataMLE_n50_bin10[,1]-68)^2))
muxMLE100b10<- sqrt(mean((dataMLE_n100_bin10[,1]-68)^2))
muxMLE300b10<- sqrt(mean((dataMLE_n300_bin10[,1]-68)^2))
muxMLE600b10<- sqrt(mean((dataMLE_n600_bin10[,1]-68)^2))
muxMLE1000b10<- sqrt(mean((dataMLE_n1000_bin10[,1]-68)^2))

muxMLE50b10
muxMLE100b10
muxMLE300b10
muxMLE600b10
muxMLE1000b10
c(muxMLE50b10,muxMLE100b10,muxMLE300b10,muxMLE600b10,muxMLE1000b10)
###################################################################
muyMLE50b10<- sqrt(mean((dataMLE_n50_bin10[,2]-68)^2))
muyMLE100b10<- sqrt(mean((dataMLE_n100_bin10[,2]-68)^2))
muyMLE300b10<- sqrt(mean((dataMLE_n300_bin10[,2]-68)^2))
muyMLE600b10<- sqrt(mean((dataMLE_n600_bin10[,2]-68)^2))
muyMLE1000b10<- sqrt(mean((dataMLE_n1000_bin10[,2]-68)^2))

muyMLE50b10
muyMLE100b10
muyMLE300b10
muyMLE600b10
muyMLE1000b10
c(muyMLE50b10,muyMLE100b10,muyMLE300b10,muyMLE600b10,muyMLE1000b10)
#############################################################################
ssxMLE50b10<- sqrt(mean((dataMLE_n50_bin10[,3]-3)^2))
ssxMLE100b10<- sqrt(mean((dataMLE_n100_bin10[,3]-3)^2))
ssxMLE300b10<- sqrt(mean((dataMLE_n300_bin10[,3]-3)^2))
ssxMLE600b10<- sqrt(mean((dataMLE_n600_bin10[,3]-3)^2))
ssxMLE1000b10<- sqrt(mean((dataMLE_n1000_bin10[,3]-3)^2))

ssxMLE50b10
ssxMLE100b10
ssxMLE300b10
ssxMLE600b10
ssxMLE1000b10
c(ssxMLE50b10,ssxMLE100b10,ssxMLE300b10,ssxMLE600b10,ssxMLE1000b10)
###############################################################################
ssyMLE50b10<- sqrt(mean((dataMLE_n50_bin10[,4]-6)^2))
ssyMLE100b10<- sqrt(mean((dataMLE_n100_bin10[,4]-6)^2))
ssyMLE300b10<- sqrt(mean((dataMLE_n300_bin10[,4]-6)^2))
ssyMLE600b10<- sqrt(mean((dataMLE_n600_bin10[,4]-6)^2))
ssyMLE1000b10<- sqrt(mean((dataMLE_n1000_bin10[,4]-6)^2))

ssyMLE50b10
ssyMLE100b10
ssyMLE300b10
ssyMLE600b10
ssyMLE1000b10
c(ssyMLE50b10,ssyMLE100b10,ssyMLE300b10,ssyMLE600b10,ssyMLE1000b10)
################################################################################
rhoMLE50b10<- sqrt(mean((dataMLE_n50_bin10[,5]-0.4714045)^2))
rhoMLE100b10<- sqrt(mean((dataMLE_n100_bin10[,5]-0.4714045)^2))
rhoMLE300b10<- sqrt(mean((dataMLE_n300_bin10[,5]-0.4714045)^2))
rhoMLE600b10<- sqrt(mean((dataMLE_n600_bin10[,5]-0.4714045)^2))
rhoMLE1000b10<- sqrt(mean((dataMLE_n1000_bin10[,5]-0.4714045)^2))

rhoMLE50b10
rhoMLE100b10
rhoMLE300b10
rhoMLE600b10
rhoMLE1000b10
c(rhoMLE50b10,rhoMLE100b10,rhoMLE300b10,rhoMLE600b10,rhoMLE1000b10)
##################################################################################
method<- c("MLE EXAct","EM","MCEM")

n<- c(50,100,300,600,1000)

MXMLE<- c(muxMLE50b10,muxMLE100b10,muxMLE300b10,muxMLE600b10,muxMLE1000b10)

M1<- cbind(n,MXMLE)

MYMLE<- c(muyMLE50b10,muyMLE100b10,muyMLE300b10,muyMLE600b10,muyMLE1000b10)
M2<- cbind(n,MYMLE)

VXMLE<- c(ssxMLE50b10,ssxMLE100b10,ssxMLE300b10,ssxMLE600b10,ssxMLE1000b10)
M3<- cbind(n,VXMLE)


VYMLE<- c(ssyMLE50b10,ssyMLE100b10,ssyMLE300b10,ssyMLE600b10,ssyMLE1000b10)
M4<- cbind(n,VYMLE)


RHOMLE<- c(rhoMLE50b10,rhoMLE100b10,rhoMLE300b10,rhoMLE600b10,rhoMLE1000b10)
M5<- cbind(n,RHOMLE)

##########################################################################################
A1<- cbind(dataMLE_n50_bin10[,1],dataMLE_n100_bin10[,1],dataMLE_n300_bin10[,1],dataMLE_n600_bin10[,1],dataMLE_n1000_bin10[,1])
colnames(A1)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(A1)
B1<- cbind(dataMLE_n50_bin10[,2],dataMLE_n100_bin10[,2],dataMLE_n300_bin10[,2],dataMLE_n600_bin10[,2],dataMLE_n1000_bin10[,2])
colnames(B1)<- c("n=50","n=100","n=300","n=600","n=1000")

C1<- cbind(dataMLE_n50_bin10[,3],dataMLE_n100_bin10[,3],dataMLE_n300_bin10[,3],dataMLE_n600_bin10[,3],dataMLE_n1000_bin10[,3])
colnames(C1)<- c("n=50","n=100","n=300","n=600","n=1000")

D1<- cbind(dataMLE_n50_bin10[,4],dataMLE_n100_bin10[,4],dataMLE_n300_bin10[,4],dataMLE_n600_bin10[,4],dataMLE_n1000_bin10[,4])
colnames(D1)<- c("n=50","n=100","n=300","n=600","n=1000")

E1<- cbind(dataMLE_n50_bin10[,5],dataMLE_n100_bin10[,5],dataMLE_n300_bin10[,5],dataMLE_n600_bin10[,5],dataMLE_n1000_bin10[,5])
colnames(E1)<- c("n=50","n=100","n=300","n=600","n=1000")


############################################################################################################################
############################################################################################################################
############################################################################################################################
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
##########################################################################################################
muxEM50b10<- sqrt(mean((dataEM_n50_bin10[,1]-68)^2))
muxEM100b10<- sqrt(mean((dataEM_n100_bin10[,1]-68)^2))
muxEM300b10<- sqrt(mean((dataEM_n300_bin10[,1]-68)^2))
muxEM600b10<- sqrt(mean((dataEM_n600_bin10[,1]-68)^2))
muxEM1000b10<- sqrt(mean((dataEM_n1000_bin10[,1]-68)^2))

muxEM50b10
muxEM100b10
muxEM300b10
muxEM600b10
muxEM1000b10

c(muxEM50b10,muxEM100b10,muxEM300b10,muxEM600b10,muxEM1000b10)
##########################################################################################################
muyEM50b10<- sqrt(mean((dataEM_n50_bin10[,2]-68)^2))
muyEM100b10<- sqrt(mean((dataEM_n100_bin10[,2]-68)^2))
muyEM300b10<- sqrt(mean((dataEM_n300_bin10[,2]-68)^2))
muyEM600b10<- sqrt(mean((dataEM_n600_bin10[,2]-68)^2))
muyEM1000b10<- sqrt(mean((dataEM_n1000_bin10[,2]-68)^2))

muyEM50b10
muyEM100b10
muyEM300b10
muyEM600b10
muyEM1000b10

c(muyEM50b10,muyEM100b10,muyEM300b10,muyEM600b10,muyEM1000b10)
#######################################################################################################
vxEM50b10<- sqrt(mean((dataEM_n50_bin10[,3]-3)^2))
vxEM100b10<- sqrt(mean((dataEM_n100_bin10[,3]-3)^2))
vxEM300b10<- sqrt(mean((dataEM_n300_bin10[,3]-3)^2))
vxEM600b10<- sqrt(mean((dataEM_n600_bin10[,3]-3)^2))
vxEM1000b10<- sqrt(mean((dataEM_n1000_bin10[,3]-3)^2))

vxEM50b10
vxEM100b10
vxEM300b10
vxEM600b10
vxEM1000b10

c(vxEM50b10,vxEM100b10,vxEM300b10,vxEM600b10,vxEM1000b10)
#########################################################################
vyEM50b10<- sqrt(mean((dataEM_n50_bin10[,6]-6)^2))
vyEM100b10<- sqrt(mean((dataEM_n100_bin10[,6]-6)^2))
vyEM300b10<- sqrt(mean((dataEM_n300_bin10[,6]-6)^2))
vyEM600b10<- sqrt(mean((dataEM_n600_bin10[,6]-6)^2))
vyEM1000b10<- sqrt(mean((dataEM_n1000_bin10[,6]-6)^2))

vyEM50b10
vyEM100b10
vyEM300b10
vyEM600b10
vyEM1000b10

c(vyEM50b10,vyEM100b10,vyEM300b10,vyEM600b10,vyEM1000b10)
################################################################################
rhoEM50b10<- dataEM_n50_bin10[,4]/sqrt(dataEM_n50_bin10[,3]*dataEM_n50_bin10[,6])
rhoEM100b10<- dataEM_n100_bin10[,4]/sqrt(dataEM_n100_bin10[,3]*dataEM_n100_bin10[,6])
rhoEM300b10<- dataEM_n300_bin10[,4]/sqrt(dataEM_n300_bin10[,3]*dataEM_n300_bin10[,6])
rhoEM600b10<- dataEM_n600_bin10[,4]/sqrt(dataEM_n600_bin10[,3]*dataEM_n600_bin10[,6])
rhoEM1000b10<- dataEM_n1000_bin10[,4]/sqrt(dataEM_n1000_bin10[,3]*dataEM_n1000_bin10[,6])

cbind(rhoEM50b10,rhoEM100b10,rhoEM300b10,rhoEM600b10,rhoEM1000b10)

RMrhoEM50b10<- sqrt(mean((rhoEM50b10-0.4714045)^2)) 
RMrhoEM100b10<- sqrt(mean((rhoEM100b10-0.4714045)^2)) 
RMrhoEM300b10<- sqrt(mean((rhoEM300b10-0.4714045)^2)) 
RMrhoEM600b10<- sqrt(mean((rhoEM600b10-0.4714045)^2)) 
RMrhoEM1000b10<- sqrt(mean((rhoEM1000b10-0.4714045)^2)) 

c(RMrhoEM50b10,RMrhoEM100b10,RMrhoEM300b10,RMrhoEM600b10,RMrhoEM1000b10)
########################################################################################################
method<- c("MLE EXAct","EM","MCEM")

n<- c(50,100,300,600,1000)

MXEM<- c(muxEM50b10,muxEM100b10,muxEM300b10,muxEM600b10,muxEM1000b10)
MYEM<- c(muyEM50b10,muyEM100b10,muyEM300b10,muyEM600b10,muyEM1000b10)
VXEM<- c(vxEM50b10,vxEM100b10,vxEM300b10,vxEM600b10,vxEM1000b10)
VYEM<- c(vyEM50b10,vyEM100b10,vyEM300b10,vyEM600b10,vyEM1000b10)
RHOEM<- c(RMrhoEM50b10,RMrhoEM100b10,RMrhoEM300b10,RMrhoEM600b10,RMrhoEM1000b10)

N1<- cbind(n,MXEM)
N2<- cbind(n,MYEM)
N3<- cbind(n,VXEM)
N4<- cbind(n,VYEM)
N5<- cbind(n,RHOEM)
##########################################################################################################
A2<- cbind(dataEM_n50_bin10[,1],dataEM_n100_bin10[,1],dataEM_n300_bin10[,1],dataEM_n600_bin10[,1],dataEM_n1000_bin10[,1])
colnames(A2)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(A2)

B2<- cbind(dataEM_n50_bin10[,2],dataEM_n100_bin10[,2],dataEM_n300_bin10[,2],dataEM_n600_bin10[,2],dataEM_n1000_bin10[,2])
colnames(B2)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(A2)

C2<- cbind(dataEM_n50_bin10[,3],dataEM_n100_bin10[,3],dataEM_n300_bin10[,3],dataEM_n600_bin10[,3],dataEM_n1000_bin10[,3])
colnames(C2)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(A2)

D2<- cbind(dataEM_n50_bin10[,6],dataEM_n100_bin10[,6],dataEM_n300_bin10[,6],dataEM_n600_bin10[,6],dataEM_n1000_bin10[,6])
colnames(D2)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(A2)
E2<- cbind(rhoEM50b10,rhoEM100b10,rhoEM300b10,rhoEM600b10,rhoEM1000b10)
colnames(E2)<- c("n=50","n=100","n=300","n=600","n=1000")


#########################################################################################################
#pdf("BivariateMUX.pdf")
par(mfrow=c(3,2))
par(mar=c(4,4,4,4))

boxplot(A1,main="Estimation of Mean_X using MLE Exact",xlab="sample Size",
        ylab="Mean X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(67.5,68.5))
abline(h=68)

boxplot(B1,main="Estimation of Mean_Y using MLE Exact",xlab="sample Size",
        ylab="Mean Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(67,68.9))
abline(h=68)

boxplot(C1,main="Estimation of Variance_X using MLE Exact",xlab="sample Size",
        ylab="VAR X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1.5,4.5))
abline(h=3)

boxplot(D1,main="Estimation of Variance_Y using MLE Exact",xlab="sample Size",
        ylab="VAR Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(3.5,9.5))
abline(h=6)

boxplot(E1,main="Estimation of RHO using MLE Exact",xlab="sample Size",
        ylab="RHO",col=c("pink","pink","pink","pink","pink"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(0.1,0.65))
abline(h=0.4714045)

#dev.off()
#########################################################################################################################################
#pdf("BivariateMUX.pdf")
par(mfrow=c(3,2))
par(mar=c(4,4,4,4))

boxplot(A2,main="Estimation of Mean_X using EM Algorithm",xlab="sample Size",
        ylab="Mean X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(67.5,68.5))
abline(h=68)

boxplot(B2,main="Estimation of Mean_Y using EM Algorithm",xlab="sample Size",
        ylab="Mean Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(67,68.9))
abline(h=68)


boxplot(C2,main="Estimation of Variance_X using EM Algorithm",xlab="sample Size",
        ylab="VAR X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1.5,4.5))
abline(h=3)

boxplot(D2,main="Estimation of Variance_Y using EM Algorithm",xlab="sample Size",
        ylab="VAR Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(3.5,9.5))
abline(h=6)

boxplot(E2,main="Estimation of RHO using EM Algorithm",xlab="sample Size",
        ylab="RHO",col=c("pink","pink","pink","pink","pink"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(0.1,0.65))
abline(h=0.4714045)

#dev.off()
#######################################################################################################################################
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
#########################################################################################
##########################################################################################################
muxMCEM50b10<- sqrt(mean((dataMCEM_n50_bin10[,1]-68)^2))
muxMCEM100b10<- sqrt(mean((dataMCEM_n100_bin10[,1]-68)^2))
muxMCEM300b10<- sqrt(mean((dataMCEM_n300_bin10[,1]-68)^2))
muxMCEM600b10<- sqrt(mean((dataMCEM_n600_bin10[,1]-68)^2))
muxMCEM1000b10<- sqrt(mean((dataMCEM_n1000_bin10[,1]-68)^2))

muxMCEM50b10
muxMCEM100b10
muxMCEM300b10
muxMCEM600b10
muxMCEM1000b10

c(muxMCEM50b10,muxMCEM100b10,muxMCEM300b10,muxMCEM600b10,muxMCEM1000b10)
##########################################################################################################
muyMCEM50b10<- sqrt(mean((dataMCEM_n50_bin10[,2]-68)^2))
muyMCEM100b10<- sqrt(mean((dataMCEM_n100_bin10[,2]-68)^2))
muyMCEM300b10<- sqrt(mean((dataMCEM_n300_bin10[,2]-68)^2))
muyMCEM600b10<- sqrt(mean((dataMCEM_n600_bin10[,2]-68)^2))
muyMCEM1000b10<- sqrt(mean((dataMCEM_n1000_bin10[,2]-68)^2))

muyMCEM50b10
muyMCEM100b10
muyMCEM300b10
muyMCEM600b10
muyMCEM1000b10

c(muyMCEM50b10,muyMCEM100b10,muyMCEM300b10,muyMCEM600b10,muyMCEM1000b10)
##########################################################################################################
vxMCEM50b10<- sqrt(mean((dataMCEM_n50_bin10[,3]-3)^2))
vxMCEM100b10<- sqrt(mean((dataMCEM_n100_bin10[,3]-3)^2))
vxMCEM300b10<- sqrt(mean((dataMCEM_n300_bin10[,3]-3)^2))
vxMCEM600b10<- sqrt(mean((dataMCEM_n600_bin10[,3]-3)^2))
vxMCEM1000b10<- sqrt(mean((dataMCEM_n1000_bin10[,3]-3)^2))

vxMCEM50b10
vxMCEM100b10
vxMCEM300b10
vxMCEM600b10
vxMCEM1000b10

c(vxMCEM50b10,vxMCEM100b10,vxMCEM300b10,vxMCEM600b10,vxMCEM1000b10)
##########################################################################################################
vyMCEM50b10<- sqrt(mean((dataMCEM_n50_bin10[,6]-6)^2))
vyMCEM100b10<- sqrt(mean((dataMCEM_n100_bin10[,6]-6)^2))
vyMCEM300b10<- sqrt(mean((dataMCEM_n300_bin10[,6]-6)^2))
vyMCEM600b10<- sqrt(mean((dataMCEM_n600_bin10[,6]-6)^2))
vyMCEM1000b10<- sqrt(mean((dataMCEM_n1000_bin10[,6]-6)^2))

vyMCEM50b10
vyMCEM100b10
vyMCEM300b10
vyMCEM600b10
vyMCEM1000b10

c(vyMCEM50b10,vyMCEM100b10,vyMCEM300b10,vyMCEM600b10,vyMCEM1000b10)
##########################################################################################################
rhoMCEM50b10<- dataMCEM_n50_bin10[,4]/sqrt(dataMCEM_n50_bin10[,3]*dataMCEM_n50_bin10[,6])
rhoMCEM100b10<- dataMCEM_n100_bin10[,4]/sqrt(dataMCEM_n100_bin10[,3]*dataMCEM_n100_bin10[,6])
rhoMCEM300b10<- dataMCEM_n300_bin10[,4]/sqrt(dataMCEM_n300_bin10[,3]*dataMCEM_n300_bin10[,6])
rhoMCEM600b10<- dataMCEM_n600_bin10[,4]/sqrt(dataMCEM_n600_bin10[,3]*dataMCEM_n600_bin10[,6])
rhoMCEM1000b10<- dataMCEM_n1000_bin10[,4]/sqrt(dataMCEM_n1000_bin10[,3]*dataMCEM_n1000_bin10[,6])

cbind(rhoMCEM50b10,rhoMCEM100b10,rhoMCEM300b10,rhoMCEM600b10,rhoMCEM1000b10)

RMrhoMCEM50b10<- sqrt(mean((rhoMCEM50b10-0.4714045)^2)) 
RMrhoMCEM100b10<- sqrt(mean((rhoMCEM100b10-0.4714045)^2)) 
RMrhoMCEM300b10<- sqrt(mean((rhoMCEM300b10-0.4714045)^2)) 
RMrhoMCEM600b10<- sqrt(mean((rhoMCEM600b10-0.4714045)^2)) 
RMrhoMCEM1000b10<- sqrt(mean((rhoMCEM1000b10-0.4714045)^2)) 

c(RMrhoMCEM50b10,RMrhoMCEM100b10,RMrhoMCEM300b10,RMrhoMCEM600b10,RMrhoMCEM1000b10)
########################################################################################################
method<- c("MLE EXAct","EM","MCEM")

n<- c(50,100,300,600,1000)

MXMCEM<- c(muxMCEM50b10,muxMCEM100b10,muxMCEM300b10,muxMCEM600b10,muxMCEM1000b10)
MYMCEM<- c(muyMCEM50b10,muyMCEM100b10,muyMCEM300b10,muyMCEM600b10,muyMCEM1000b10)
VXMCEM<- c(vxMCEM50b10,vxMCEM100b10,vxMCEM300b10,vxMCEM600b10,vxMCEM1000b10)
VYMCEM<- c(vyMCEM50b10,vyMCEM100b10,vyMCEM300b10,vyMCEM600b10,vyMCEM1000b10)
RHOMCEM<- c(RMrhoMCEM50b10,RMrhoMCEM100b10,RMrhoMCEM300b10,RMrhoMCEM600b10,RMrhoMCEM1000b10)




#############################################################################       
###############################################################################################################
A3<- cbind(dataMCEM_n50_bin10[,1],dataMCEM_n100_bin10[,1],dataMCEM_n300_bin10[,1],dataMCEM_n600_bin10[,1],dataMCEM_n1000_bin10[,1])
colnames(A3)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(A3)

B3<- cbind(dataMCEM_n50_bin10[,2],dataMCEM_n100_bin10[,2],dataMCEM_n300_bin10[,2],dataMCEM_n600_bin10[,2],dataMCEM_n1000_bin10[,2])
colnames(B3)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(B3)

C3<- cbind(dataMCEM_n50_bin10[,3],dataMCEM_n100_bin10[,3],dataMCEM_n300_bin10[,3],dataMCEM_n600_bin10[,3],dataMCEM_n1000_bin10[,3])
colnames(C3)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(C3)

D3<- cbind(dataMCEM_n50_bin10[,6],dataMCEM_n100_bin10[,6],dataMCEM_n300_bin10[,6],dataMCEM_n600_bin10[,6],dataMCEM_n1000_bin10[,6])
colnames(D3)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(D3)
E3<- cbind(rhoMCEM50b10,rhoMCEM100b10,rhoMCEM300b10,rhoMCEM600b10,rhoMCEM1000b10)
colnames(E3)<- c("n=50","n=100","n=300","n=600","n=1000")
#boxplot(E3)
##############################################################################################################################

#pdf("BivariateMUX.pdf")
par(mfrow=c(3,2))
par(mar=c(4,4,4,4))

boxplot(A3,main="Estimation of Mean_X using MCEM Algorithm",xlab="sample Size",
        ylab="Mean X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(67.5,68.5))
abline(h=68)

boxplot(B3,main="Estimation of Mean_Y using MCEM Algorithm",xlab="sample Size",
        ylab="Mean Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7
        ,ylim=c(67,68.9))
abline(h=68)

boxplot(C3,main="Estimation of Variance_X using MCEM Algorithm",xlab="sample Size",
        ylab="VAR X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1.5,4.5))
abline(h=3)

boxplot(D3,main="Eystimation of Variance_Y using MCEM Algorithm",xlab="sample Size",
        ylab="VAR Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(3.5,9.5))
abline(h=6)

boxplot(E3,main="Estimation of RHO using MCEM Algorithm",xlab="sample Size",
        ylab="RHo",col=c("pink","pink","pink","pink","pink"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(0.1,0.65))
abline(h=0.4714045)

#dev.off()
#####################################################################################################



### RMSE TABLES#########
#######################################
method<- c("MLE EXAct","EM","MCEM")

n<- c(50,100,300,600,1000)


M1<- cbind(n,round(MXMLE,6))
M1
M2<- cbind(n,round(MYMLE,6))
M3<- cbind(n,round(VXMLE,6))
M4<- cbind(n,round(VYMLE,6))
M5<- cbind(n,round(RHOMLE,6))
M5
N1<- cbind(n,round(MXEM,6))
N2<- cbind(n,round(MYEM,6))
N3<- cbind(n,round(VXEM,6))
N4<- cbind(n,round(VYEM,6))
N5<- cbind(n,round(RHOEM,6))



P1<- cbind(n,round(MXMCEM,6))
P2<- cbind(n,round(MYMCEM,6))
P3<- cbind(n,round(VXMCEM,6))
P4<- cbind(n,round(VYMCEM,6))
P5<- cbind(n,round(RHOMCEM,6))


param1<- rep("mu_x1",5)
samplesize<- c(50,100,300,600,1000)
T1<- cbind(param1,samplesize, round(MXMLE,6),round(MXEM,6),round(MXMCEM,6))
T1
T1MX<- data.frame(T1)
colnames(T1MX)<- c("parameter","sample size","Exact MLE","EM","MCEM")
T1MX


#xtable(T1MX,digits = 6)
#########################################################
param2<- rep("mu_x2",5)

T2<- cbind(param2,samplesize, round(MYMLE,6),round(MYEM,6),round(MYMCEM,6))
T2
T2MY<- data.frame(T2)
colnames(T2MY)<- c("parameter","sample size","Exact MLE","EM","MCEM")

T2MY
#xtable(T2MY,digits=6)
#################################################################
param3<- rep("Var_x1",5)

T3<- cbind(param3,samplesize,round(VXMLE,6),round(VXEM,6),round(VXMCEM,6))
T3
T3VX<- data.frame(T3)
colnames(T3VX)<- c("parameter","sample size","Exact MLE","EM","MCEM")

T3VX
#xtable(T3VX,digits = 6)
#####################################################################       
param4<- rep("Var_x2",5)

T4<- cbind(param4,samplesize, round(VYMLE,6),round(VYEM,6),round(VYMCEM,6))
T4
T4VY<- data.frame(T4)
colnames(T4VY)<- c("parameter","sample size","Exact MLE","EM","MCEM")

T4VY
#xtable(T4VY,digits = 6)   
################################################
param5<- rep("RHO",5)

T5<- cbind(param5,samplesize, round(RHOMLE,6),round(RHOEM,6),round(RHOMCEM,6))
T5
T5RHO<- data.frame(T5)
colnames(T5RHO)<- c("parameter","sample size","Exact MLE","EM","MCEM")

T5RHO
#xtable(T5RHO,digits = 6)   

##############################################################################################
final_table<- rbind(T1MX,T2MY,T3VX,T4VY,T5RHO)
final_table
xtable(final_table)
##############################################################################################
#pdf("BivariateMUX.pdf")
par(mfrow=c(3,2))
par(mar=c(4,4,4,4))

boxplot(A1,main="Estimation of Mean_X using MLE Exact",xlab="sample Size",
        ylab="Mean X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(67.5,68.5))
abline(h=68)

boxplot(B1,main="Estimation of Mean_Y using MLE Exact",xlab="sample Size",
        ylab="Mean Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(67,68.9))
abline(h=68)

boxplot(C1,main="Estimation of Variance_X using MLE Exact",xlab="sample Size",
        ylab="VAR X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1.5,4.5))
abline(h=3)

boxplot(D1,main="Estimation of Variance_Y using MLE Exact",xlab="sample Size",
        ylab="VAR Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(3.5,9.5))
abline(h=6)

boxplot(E1,main="Estimation of RHO using MLE Exact",xlab="sample Size",
        ylab="RHO",col=c("pink","pink","pink","pink","pink"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(0.1,0.65))
abline(h=0.4714045)

#dev.off()
#########################################################################################################################################
#pdf("BivariateMUX.pdf")
par(mfrow=c(3,2))
par(mar=c(4,4,4,4))

boxplot(A2,main="Estimation of Mean_X using EM Algorithm",xlab="sample Size",
        ylab="Mean X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(67.5,68.5))
abline(h=68)

boxplot(B2,main="Estimation of Mean_Y using EM Algorithm",xlab="sample Size",
        ylab="Mean Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(67,68.9))
abline(h=68)


boxplot(C2,main="Estimation of Variance_X using EM Algorithm",xlab="sample Size",
        ylab="VAR X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1.5,4.5))
abline(h=3)

boxplot(D2,main="Estimation of Variance_Y using EM Algorithm",xlab="sample Size",
        ylab="VAR Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(3.5,9.5))
abline(h=6)

boxplot(E2,main="Estimation of RHO using EM Algorithm",xlab="sample Size",
        ylab="RHO",col=c("pink","pink","pink","pink","pink"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(0.1,0.65))
abline(h=0.4714045)

#dev.off()
########################################################################################################

#pdf("BivariateMUX.pdf")
par(mfrow=c(3,2))
par(mar=c(4,4,4,4))

boxplot(A3,main="Estimation of Mean_X using MCEM Algorithm",xlab="sample Size",
        ylab="Mean X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(67.5,68.5))
abline(h=68)

boxplot(B3,main="Estimation of Mean_Y using MCEM Algorithm",xlab="sample Size",
        ylab="Mean Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7
        ,ylim=c(67,68.9))
abline(h=68)

boxplot(C3,main="Estimation of Variance_X using MCEM Algorithm",xlab="sample Size",
        ylab="VAR X",col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1.5,4.5))
abline(h=3)

boxplot(D3,main="Eystimation of Variance_Y using MCEM Algorithm",xlab="sample Size",
        ylab="VAR Y",col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(3.5,9.5))
abline(h=6)

boxplot(E3,main="Estimation of RHO using MCEM Algorithm",xlab="sample Size",
        ylab="RHo",col=c("pink","pink","pink","pink","pink"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(0.1,0.65))
abline(h=0.4714045)

#dev.off()
#####################################################################################################
#pdf("BivariateMUX1.pdf")
pdf("BivariateMUX1-July2022.pdf")

par(mfrow=c(2,2))
par(mar=c(4,4,4,4))

boxplot(A1,main=expression(paste("Estimation of ", mu[x1], " using MLE Exact")),xlab="Sample size",
        ylab=expression(mu[x1]),col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(67,68.75))
abline(h=68)

boxplot(A2,main=expression(paste("Estimation of ", mu[x1], " using EM Algorithm")),xlab="Sample size",
        ylab=expression(mu[x1]),col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(67,68.75))
abline(h=68)

boxplot(A3,main=expression(paste("Estimation of ", mu[x1], " using MCEM Algorithm")),xlab="Sample size",
        ylab=expression(mu[x1]),col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(67,68.75))
abline(h=68)

dev.off()

#######################################################################################################################
#pdf("BivariateMUX2.pdf")
pdf("BivariateMUX2-July2022.pdf")

par(mfrow=c(2,2))
par(mar=c(4,4,4,4))

boxplot(B1,main=expression(paste("Estimation of ", mu[x2], " using MLE Exact")),xlab="Sample size",
        ylab=expression(mu[x2]),col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(66.75,69.25))
abline(h=68)

boxplot(B2,main=expression(paste("Estimation of ", mu[x2], " using EM Algorithm")),xlab="Sample size",
        ylab=expression(mu[x2]),col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(66.75,69.25))
abline(h=68)

boxplot(B3,main=expression(paste("Estimation of ", mu[x2], " using MCEM Algorithm")),xlab="Sample size",
        ylab=expression(mu[x2]),col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7
        ,ylim=c(66.75,69.25))
abline(h=68)


dev.off()
##########################################################################################################
#pdf("BivariatevarX1.pdf")
pdf("BivariatevarX1-July2022.pdf")

par(mfrow=c(2,2))
par(mar=c(4,4,4,4))

boxplot(C1,main=expression(paste("Estimation of ", sigma[x1]^2, " using MLE Exact")),xlab="Sample size",
        ylab=expression(sigma[x1]^2),col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,cex.lab=0.7
        ,ylim=c(1.5,5))
abline(h=3)


boxplot(C2,main=expression(paste("Estimation of ", sigma[x1]^2, " using EM Algorithm")),xlab="Sample size",
        ylab=expression(sigma[x1]^2),col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,
        cex.lab=0.7,
        ylim=c(1.5,5))
abline(h=3)


boxplot(C3,main=expression(paste("Estimation of ", sigma[x1]^2, " using MCEM Algorithm")),xlab="Sample size",
        ylab=expression(sigma[x1]^2),col=c("green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,
        cex.lab=0.7,
        ylim=c(1.5,5))
abline(h=3)

dev.off()
###########################################################################################################
pdf("BivariateVarX2.pdf")
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))

boxplot(D1,main=expression(paste("Estimation of ", sigma[x2]^2, " using MLE Exact")),xlab="Sample size",
        ylab=expression(sigma[x2]^2),col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7
        ,cex.lab=0.7
        ,ylim=c(3,10.5))
abline(h=6)


boxplot(D2,main=expression(paste("Estimation of ", sigma[x2]^2, " using EM Algorithm")),xlab="Sample size",
        ylab=expression(sigma[x2]^2),col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7
        ,cex.lab=0.7
        ,ylim=c(3,10.5))
abline(h=6)


boxplot(D3,main=expression(paste("Estimation of ", sigma[x2]^2, " using MCEM Algorithm")),xlab="Sample size",
        ylab=expression(sigma[x2]^2),col=c("light blue","light blue","light blue","light blue","light blue"),las=2,cex.axis=0.70,cex.main=0.7
        ,cex.lab=0.7
        ,ylim=c(3,10.5))
abline(h=6)

dev.off()
##########################################################################################################
#pdf("BivariateRHO.pdf")
pdf("BivariateRHO-July2022.pdf")

par(mfrow=c(2,2))
par(mar=c(4,4,4,4))


boxplot(E1,main=expression(paste("Estimation of ", rho, " using MLE Exact")),xlab="Sample size",
        ylab=expression(rho),col=c("pink","pink","pink","pink","pink"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(0.1,0.8))
abline(h=0.4714045)


boxplot(E2,main=expression(paste("Estimation of ", rho, " using EM Algorithm")),xlab="Sample size",
        ylab=expression(rho),col=c("pink","pink","pink","pink","pink"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(0.1,0.8))
abline(h=0.4714045)


boxplot(E3,main=expression(paste("Estimation of ", rho, " using MCEM Algorithm")),xlab="Sample size",
        ylab=expression(rho),col=c("pink","pink","pink","pink","pink"),las=2,cex.axis=0.70,cex.main=0.7,
        ylim=c(0.1,0.8))
abline(h=0.4714045)

dev.off()


