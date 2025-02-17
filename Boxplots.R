### BOXPLOTS###
### EM For MU###

bins8<- cbind(output50m8EM[,1],output100m8EM[,1],output300m8EM[,1],output600m8EM[,1],output1000m8EM[,1])
colnames(bins8)<- c("n=50","n=100","n=300","n=600","n=1000")
bins8

boxplot(bins8, )


bins15<- cbind(output50m15EM[,1],output100m15EM[,1],output300m15EM[,1],output600m15EM[,1],output1000m15EM[,1])
colnames(bins15)<- c("n=50","n=100","n=300","n=600","n=1000")
bins15


bins30<- cbind(output50m30EM[,1],output100m30EM[,1],output300m30EM[,1],output600m30EM[,1],output1000m30EM[,1])
colnames(bins30)<- c("n=50","n=100","n=300","n=600","n=1000")
bins30
boxplot(bins30)
A<- cbind(bins8,bins15,bins30)
boxplot(A)
abline(h=68)

#pdf("UniMUEM1.pdf")
boxplot(A,main="Estimation of MU using EM Algorithm",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70)
abline(h=68)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black")
#dev.off()


#pdf("MUEqual.pdf")
#par(mfrow=c(2,2))
#par(mar=c(4,4,4,4))
#dev.off

#outputMCEM1000m30
### MCEM For MU####
bins8MCEM<- cbind(outputMCEM50m8[,1],outputMCEM100m8[,1],outputMCEM300m8[,1],outputMCEM600m8[,1],outputMCEM1000m8[,1])
colnames(bins8MCEM)<- c("n=50","n=100","n=300","n=600","n=1000")
bins8MCEM

bins15MCEM<- cbind(outputMCEM50m15[,1],outputMCEM100m15[,1],outputMCEM300m15[,1],outputMCEM600m15[,1],outputMCEM1000m15[,1])
colnames(bins15MCEM)<- c("n=50","n=100","n=300","n=600","n=1000")
bins15MCEM

bins30MCEM<- cbind(outputMCEM50m30[,1],outputMCEM100m30[,1],outputMCEM300m30[,1],outputMCEM600m30[,1],outputMCEM1000m30[,1])
colnames(bins30MCEM)<- c("n=50","n=100","n=300","n=600","n=1000")
bins30MCEM

#pdf("UniMUMCEM1.pdf")
B<- cbind(bins8MCEM,bins15MCEM,bins30MCEM)
boxplot(B,main="Estimation of MU using MCEM Algorithm",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70)
abline(h=68)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black")


#dev.off()
### MLE EXACT for MU####

MLEbinned1000m8

bins8MLEbinned<- cbind(MLEbinned50m8[,1],MLEbinned100m8[,1],MLEbinned300m8[,1],MLEbinned600m8[,1],MLEbinned1000m8[,1])
colnames(bins8MLEbinned)<- c("n=50","n=100","n=300","n=600","n=1000")
bins8MLEbinned

bins15MLEbinned<- cbind(MLEbinned50m15[,1],MLEbinned100m15[,1],MLEbinned300m15[,1],MLEbinned600m15[,1],MLEbinned1000m15[,1])
colnames(bins15MLEbinned)<- c("n=50","n=100","n=300","n=600","n=1000")
bins15MLEbinned

bins30MLEbinned<- cbind(MLEbinned50m30[,1],MLEbinned100m30[,1],MLEbinned300m30[,1],MLEbinned600m30[,1],MLEbinned1000m30[,1])
colnames(bins30MLEbinned)<- c("n=50","n=100","n=300","n=600","n=1000")
bins30MLEbinned

C<- cbind(bins8MLEbinned,bins15MLEbinned,bins30MLEbinned)


#### MLE Unbinned For MU ####
#MLEunbinned1000m30

bins8MLEUnbinned<- cbind(MLEunbinned50m8[,1],MLEunbinned100m8[,1],MLEunbinned300m8[,1],MLEunbinned600m8[,1],MLEunbinned1000m8[,1])
colnames(bins8MLEUnbinned)<- c("n=50","n=100","n=300","n=600","n=1000")
bins8MLEUnbinned

bins15MLEUnbinned<- cbind(MLEunbinned50m15[,1],MLEunbinned100m15[,1],MLEunbinned300m15[,1],MLEunbinned600m15[,1],MLEunbinned1000m15[,1])
colnames(bins15MLEUnbinned)<- c("n=50","n=100","n=300","n=600","n=1000")
bins15MLEUnbinned

bins30MLEUnbinned<- cbind(MLEunbinned50m30[,1],MLEunbinned100m30[,1],MLEunbinned300m30[,1],MLEunbinned600m30[,1],MLEunbinned1000m30[,1])
colnames(bins30MLEUnbinned)<- c("n=50","n=100","n=300","n=600","n=1000")
bins30MLEUnbinned

D<- cbind(bins8MLEUnbinned,bins15MLEUnbinned,bins30MLEUnbinned)


#############################################################################################################
pdf("UnivariateMU.pdf")
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))

boxplot(D,main="Estimation of MU using MLE Ignoring grouping",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7)
abline(h=68)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)

boxplot(C,main="Estimation of MU using MLE EXACT",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.8)
abline(h=68)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)

boxplot(A,main="Estimation of MU using EM Algorithm",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.8)
abline(h=68)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)

boxplot(B,main="Estimation of MU using MCEM Algorithm",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.8)
abline(h=68)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)


dev.off()
#############################################################################################################
#############################################################################################################
### EM For Variance###

bins8V<- cbind(output50m8EM[,2],output100m8EM[,2],output300m8EM[,2],output600m8EM[,2],output1000m8EM[,2])
colnames(bins8V)<- c("n=50","n=100","n=300","n=600","n=1000")
bins8V


bins15V<- cbind(output50m15EM[,2],output100m15EM[,2],output300m15EM[,2],output600m15EM[,2],output1000m15EM[,2])
colnames(bins15V)<- c("n=50","n=100","n=300","n=600","n=1000")
bins15V


bins30V<- cbind(output50m30EM[,2],output100m30EM[,2],output300m30EM[,2],output600m30EM[,2],output1000m30EM[,2])
colnames(bins30V)<- c("n=50","n=100","n=300","n=600","n=1000")
bins30V
AA<- cbind(bins8V,bins15V,bins30V)
boxplot(AA)
abline(h=3.24)

### MCEM For Variance####
bins8MCEMV<- cbind(outputMCEM50m8[,2],outputMCEM100m8[,2],outputMCEM300m8[,2],outputMCEM600m8[,2],outputMCEM1000m8[,2])
colnames(bins8MCEMV)<- c("n=50","n=100","n=300","n=600","n=1000")
bins8MCEMV

bins15MCEMV<- cbind(outputMCEM50m15[,2],outputMCEM100m15[,2],outputMCEM300m15[,2],outputMCEM600m15[,2],outputMCEM1000m15[,2])
colnames(bins15MCEMV)<- c("n=50","n=100","n=300","n=600","n=1000")
bins15MCEMV

bins30MCEMV<- cbind(outputMCEM50m30[,2],outputMCEM100m30[,2],outputMCEM300m30[,2],outputMCEM600m30[,2],outputMCEM1000m30[,2])
colnames(bins30MCEMV)<- c("n=50","n=100","n=300","n=600","n=1000")
bins30MCEMV

BB<- cbind(bins8MCEMV,bins15MCEMV,bins30MCEMV)
boxplot(BB,main="Estimation of MU using MCEM Algorithm",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70,ylim=c(1,6))
abline(h=3.24)
#legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black")

### MLE EXACT for Variance####

#MLEbinned1000m8

bins8MLEbinnedV<- cbind(MLEbinned50m8[,2],MLEbinned100m8[,2],MLEbinned300m8[,2],MLEbinned600m8[,2],MLEbinned1000m8[,2])
colnames(bins8MLEbinnedV)<- c("n=50","n=100","n=300","n=600","n=1000")
bins8MLEbinnedV

bins15MLEbinnedV<- cbind(MLEbinned50m15[,2],MLEbinned100m15[,2],MLEbinned300m15[,2],MLEbinned600m15[,2],MLEbinned1000m15[,2])
colnames(bins15MLEbinnedV)<- c("n=50","n=100","n=300","n=600","n=1000")
bins15MLEbinnedV

bins30MLEbinnedV<- cbind(MLEbinned50m30[,2],MLEbinned100m30[,2],MLEbinned300m30[,2],MLEbinned600m30[,2],MLEbinned1000m30[,2])
colnames(bins30MLEbinnedV)<- c("n=50","n=100","n=300","n=600","n=1000")
bins30MLEbinnedV

CC<- cbind(bins8MLEbinnedV,bins15MLEbinnedV,bins30MLEbinnedV)

boxplot(CC,main="Estimation of MU using MCEM Algorithm",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70,ylim=c(1,6))
abline(h=3.24)

#### MLE Unbinned For Variance ####
#MLEunbinned1000m30

bins8MLEUnbinnedV<- cbind(MLEunbinned50m8[,2],MLEunbinned100m8[,2],MLEunbinned300m8[,2],MLEunbinned600m8[,2],MLEunbinned1000m8[,2])
colnames(bins8MLEUnbinnedV)<- c("n=50","n=100","n=300","n=600","n=1000")
bins8MLEUnbinnedV

bins15MLEUnbinnedV<- cbind(MLEunbinned50m15[,2],MLEunbinned100m15[,2],MLEunbinned300m15[,2],MLEunbinned600m15[,2],MLEunbinned1000m15[,2])
colnames(bins15MLEUnbinnedV)<- c("n=50","n=100","n=300","n=600","n=1000")
bins15MLEUnbinnedV

bins30MLEUnbinnedV<- cbind(MLEunbinned50m30[,2],MLEunbinned100m30[,2],MLEunbinned300m30[,2],MLEunbinned600m30[,2],MLEunbinned1000m30[,2])
colnames(bins30MLEUnbinnedV)<- c("n=50","n=100","n=300","n=600","n=1000")
bins30MLEUnbinnedV

DD<- cbind(bins8MLEUnbinnedV,bins15MLEUnbinnedV,bins30MLEUnbinnedV)

boxplot(DD,main="Estimation of MU using MCEM Algorithm",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70,ylim=c(1,6))
abline(h=3.24)
#############################################################################################################
#############################################################################
pdf("UnivariateVar.pdf")
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))

boxplot(DD,main="Estimation of Variance using MLE Ignoring Grouping",xlab="sample Size",
        ylab="Variance",col=c("blue","blue","blue","blue","blue",
                              "red","red","red","red","red",
                              "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1,6))
abline(h=3.24)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)


boxplot(CC,main="Estimation of Variance using MLE Exact",xlab="sample Size",
        ylab="Variance",col=c("blue","blue","blue","blue","blue",
                              "red","red","red","red","red",
                              "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1,6))
abline(h=3.24)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)




boxplot(AA,main="Estimation of Variance using EM Algorithm",xlab="sample Size",
        ylab="Variance",col=c("blue","blue","blue","blue","blue",
                              "red","red","red","red","red",
                              "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1,6))
abline(h=3.24)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)

boxplot(BB,main="Estimation of Variance using MCEM Algorithm",xlab="sample Size",
        ylab="Variance",col=c("blue","blue","blue","blue","blue",
                              "red","red","red","red","red",
                              "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.7,ylim=c(1,6))
abline(h=3.24)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)


dev.off()
########################################
pdf("UnivariateMU11.pdf")
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))

legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)

boxplot(A,main="Estimation of MU using EM Algorithm",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.8)
abline(h=68)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)

boxplot(B,main="Estimation of MU using MCEM Algorithm",xlab="sample Size",
        ylab="Mean",col=c("blue","blue","blue","blue","blue",
                          "red","red","red","red","red",
                          "green","green","green","green","green"),las=2,cex.axis=0.70,cex.main=0.8)
abline(h=68)
legend("topright",c("8 bins","15 bins","30 bins"),fill=c("blue","red","green"),border="black",cex=0.5)


dev.off()