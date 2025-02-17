#outputMCEM1000m30
##################### MEAN ##################################################################
####MCEM BIAS###########################################################################
MUBiasMCEM<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(MUBiasMCEM)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                        "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                        "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                        "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                        "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")
                        
MUBiasMCEM[,1]<- outputMCEM50m8[,1]-68
MUBiasMCEM[,2]<- outputMCEM50m15[,1]-68
MUBiasMCEM[,3]<- outputMCEM50m30[,1]-68

MUBiasMCEM[,4]<- outputMCEM100m8[,1]-68
MUBiasMCEM[,5]<- outputMCEM100m15[,1]-68
MUBiasMCEM[,6]<- outputMCEM100m30[,1]-68

MUBiasMCEM[,7]<- outputMCEM300m8[,1]-68
MUBiasMCEM[,8]<- outputMCEM300m15[,1]-68
MUBiasMCEM[,9]<- outputMCEM300m30[,1]-68

MUBiasMCEM[,10]<- outputMCEM600m8[,1]-68
MUBiasMCEM[,11]<- outputMCEM600m15[,1]-68
MUBiasMCEM[,12]<- outputMCEM600m30[,1]-68

MUBiasMCEM[,13]<- outputMCEM1000m8[,1]-68
MUBiasMCEM[,14]<- outputMCEM1000m15[,1]-68
MUBiasMCEM[,15]<- outputMCEM1000m30[,1]-68

MUaveBiasMCEM<- colMeans(MUBiasMCEM)
###MCEM RMSE ########################################################################
MUMSEMCEM<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(MUMSEMCEM)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                         "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                         "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                         "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                         "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

MUMSEMCEM[,1]<- (outputMCEM50m8[,1]-68)^2
MUMSEMCEM[,2]<- (outputMCEM50m15[,1]-68)^2
MUMSEMCEM[,3]<- (outputMCEM50m30[,1]-68)^2

MUMSEMCEM[,4]<- (outputMCEM100m8[,1]-68)^2
MUMSEMCEM[,5]<- (outputMCEM100m15[,1]-68)^2
MUMSEMCEM[,6]<- (outputMCEM100m30[,1]-68)^2

MUMSEMCEM[,7]<- (outputMCEM300m8[,1]-68)^2
MUMSEMCEM[,8]<- (outputMCEM300m15[,1]-68)^2
MUMSEMCEM[,9]<- (outputMCEM300m30[,1]-68)^2

MUMSEMCEM[,10]<- (outputMCEM600m8[,1]-68)^2
MUMSEMCEM[,11]<- (outputMCEM600m15[,1]-68)^2
MUMSEMCEM[,12]<- (outputMCEM600m30[,1]-68)^2

MUMSEMCEM[,13]<- (outputMCEM1000m8[,1]-68)^2
MUMSEMCEM[,14]<- (outputMCEM1000m15[,1]-68)^2
MUMSEMCEM[,15]<- (outputMCEM1000m30[,1]-68)^2

MUaveMSEMCEM<- colMeans(MUMSEMCEM)
MURMSEMCEM<- sqrt(colMeans(MUMSEMCEM))

###############################################################################################
########## VARIANCE ###########################################################################
####MCEM BIAS###########################################################################
VBiasMCEM<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(VBiasMCEM)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                         "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                         "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                         "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                         "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

VBiasMCEM[,1]<- outputMCEM50m8[,2]-6.25
VBiasMCEM[,2]<- outputMCEM50m15[,2]-6.25
VBiasMCEM[,3]<- outputMCEM50m30[,2]-6.25

VBiasMCEM[,4]<- outputMCEM100m8[,2]-6.25
VBiasMCEM[,5]<- outputMCEM100m15[,2]-6.25
VBiasMCEM[,6]<- outputMCEM100m30[,2]-6.25

VBiasMCEM[,7]<- outputMCEM300m8[,2]-6.25
VBiasMCEM[,8]<- outputMCEM300m15[,2]-6.25
VBiasMCEM[,9]<- outputMCEM300m30[,2]-6.25

VBiasMCEM[,10]<- outputMCEM600m8[,2]-6.25
VBiasMCEM[,11]<- outputMCEM600m15[,2]-6.25
VBiasMCEM[,12]<- outputMCEM600m30[,2]-6.25

VBiasMCEM[,13]<- outputMCEM1000m8[,2]-6.25
VBiasMCEM[,14]<- outputMCEM1000m15[,2]-6.25
VBiasMCEM[,15]<- outputMCEM1000m30[,2]-6.25

VaveBiasMCEM<- colMeans(VBiasMCEM)

###MCEM RMSE ########################################################################
VMSEMCEM<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(VMSEMCEM)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                        "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                        "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                        "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                        "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

VMSEMCEM[,1]<- (outputMCEM50m8[,2]-6.25)^2
VMSEMCEM[,2]<- (outputMCEM50m15[,2]-6.25)^2
VMSEMCEM[,3]<- (outputMCEM50m30[,2]-6.25)^2

VMSEMCEM[,4]<- (outputMCEM100m8[,2]-6.25)^2
VMSEMCEM[,5]<- (outputMCEM100m15[,2]-6.25)^2
VMSEMCEM[,6]<- (outputMCEM100m30[,2]-6.25)^2

VMSEMCEM[,7]<- (outputMCEM300m8[,2]-6.25)^2
VMSEMCEM[,8]<- (outputMCEM300m15[,2]-6.25)^2
VMSEMCEM[,9]<- (outputMCEM300m30[,2]-6.25)^2

VMSEMCEM[,10]<- (outputMCEM600m8[,2]-6.25)^2
VMSEMCEM[,11]<- (outputMCEM600m15[,2]-6.25)^2
VMSEMCEM[,12]<- (outputMCEM600m30[,2]-6.25)^2

VMSEMCEM[,13]<- (outputMCEM1000m8[,2]-6.25)^2
VMSEMCEM[,14]<- (outputMCEM1000m15[,2]-6.25)^2
VMSEMCEM[,15]<- (outputMCEM1000m30[,2]-6.25)^2

VaveMSEMCEM<- colMeans(VMSEMCEM)
VRMSEMCEM<- sqrt(colMeans(VMSEMCEM))


MUaveMSEMCEM<- colMeans(MUMSEMCEM)
MURMSEMCEM<- sqrt(colMeans(MUMSEMCEM))


#################################################################################################
#################################################################################################
##### EM ALgorithm ########################################################################################

###### MEAN ###############################################
#### BIAS #######################
MUBiasEM<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(MUBiasEM)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                         "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                         "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                         "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                         "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

MUBiasEM[,1]<- output50m8EM[,1]-68
MUBiasEM[,2]<- output50m15EM[,1]-68
MUBiasEM[,3]<- output50m30EM[,1]-68

MUBiasEM[,4]<- output100m8EM[,1]-68
MUBiasEM[,5]<- output100m15EM[,1]-68
MUBiasEM[,6]<- output100m30EM[,1]-68

MUBiasEM[,7]<- output300m8EM[,1]-68
MUBiasEM[,8]<- output300m15EM[,1]-68
MUBiasEM[,9]<- output300m30EM[,1]-68

MUBiasEM[,10]<- output600m8EM[,1]-68
MUBiasEM[,11]<- output600m15EM[,1]-68
MUBiasEM[,12]<- output600m30EM[,1]-68

MUBiasEM[,13]<- output1000m8EM[,1]-68
MUBiasEM[,14]<- output1000m15EM[,1]-68
MUBiasEM[,15]<- output1000m30EM[,1]-68

MUaveBiasEM<- colMeans(MUBiasEM)


###################################################################################
##### RMSE #########################################################################
MUMSEEM<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(MUMSEEM)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                       "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                       "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                       "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                       "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

MUMSEEM[,1]<- (output50m8EM[,1]-68)^2
MUMSEEM[,2]<- (output50m15EM[,1]-68)^2
MUMSEEM[,3]<- (output50m30EM[,1]-68)^2

MUMSEEM[,4]<- (output100m8EM[,1]-68)^2
MUMSEEM[,5]<- (output100m15EM[,1]-68)^2
MUMSEEM[,6]<- (output100m30EM[,1]-68)^2

MUMSEEM[,7]<- (output300m8EM[,1]-68)^2
MUMSEEM[,8]<- (output300m15EM[,1]-68)^2
MUMSEEM[,9]<- (output300m30EM[,1]-68)^2

MUMSEEM[,10]<- (output600m8EM[,1]-68)^2
MUMSEEM[,11]<- (output600m15EM[,1]-68)^2
MUMSEEM[,12]<- (output600m30EM[,1]-68)^2

MUMSEEM[,13]<- (output1000m8EM[,1]-68)^2
MUMSEEM[,14]<- (output1000m15EM[,1]-68)^2
MUMSEEM[,15]<- (output1000m30EM[,1]-68)^2

MUaveMSEEM<- colMeans(MUMSEEM)
MURMESEEM<- sqrt(colMeans(MUMSEEM))


###########################################################################################
##### VARIANCE #########################################################################
##### BIAS #################################################
VBiasEM<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(VBiasEM)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                       "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                       "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                       "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                       "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

VBiasEM[,1]<- output50m8EM[,2]-6.25
VBiasEM[,2]<- output50m15EM[,2]-6.25
VBiasEM[,3]<- output50m30EM[,2]-6.25

VBiasEM[,4]<- output100m8EM[,2]-6.25
VBiasEM[,5]<- output100m15EM[,2]-6.25
VBiasEM[,6]<- output100m30EM[,2]-6.25

VBiasEM[,7]<- output300m8EM[,2]-6.25
VBiasEM[,8]<- output300m15EM[,2]-6.25
VBiasEM[,9]<- output300m30EM[,2]-6.25

VBiasEM[,10]<- output600m8EM[,2]-6.25
VBiasEM[,11]<- output600m15EM[,2]-6.25
VBiasEM[,12]<- output600m30EM[,2]-6.25

VBiasEM[,13]<- output1000m8EM[,2]-6.25
VBiasEM[,14]<- output1000m15EM[,2]-6.25
VBiasEM[,15]<- output1000m30EM[,2]-6.25

VaveBiasEM<- colMeans(VBiasEM)


########################################################################################
##### RMSE ########################################################################
VMSEEM<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(VMSEEM)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                      "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                      "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                      "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                      "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

VMSEEM[,1]<- (output50m8EM[,2]-6.25)^2
VMSEEM[,2]<- (output50m15EM[,2]-6.25)^2
VMSEEM[,3]<- (output50m30EM[,2]-6.25)^2

VMSEEM[,4]<- (output100m8EM[,2]-6.25)^2
VMSEEM[,5]<- (output100m15EM[,2]-6.25)^2
VMSEEM[,6]<- (output100m30EM[,2]-6.25)^2

VMSEEM[,7]<- (output300m8EM[,2]-6.25)^2
VMSEEM[,8]<- (output300m15EM[,2]-6.25)^2
VMSEEM[,9]<- (output300m30EM[,2]-6.25)^2

VMSEEM[,10]<- (output600m8EM[,2]-6.25)^2
VMSEEM[,11]<- (output600m15EM[,2]-6.25)^2
VMSEEM[,12]<- (output600m30EM[,2]-6.25)^2

VMSEEM[,13]<- (output1000m8EM[,2]-6.25)^2
VMSEEM[,14]<- (output1000m15EM[,2]-6.25)^2
VMSEEM[,15]<- (output1000m30EM[,2]-6.25)^2

VaveMSEEM<- colMeans(VMSEEM)
VRMSEEM<- sqrt(colMeans(VMSEEM))


######################################################################################
################ MLE BINNED ##########################################################

##################### MEAN ##################################################################
#### BIAS###########################################################################
MUBiasMLEBin<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(MUBiasMLEBin)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                         "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                         "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                         "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                         "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

MUBiasMLEBin[,1]<- MLEbinned50m8[,1]-68
MUBiasMLEBin[,2]<- MLEbinned50m15[,1]-68
MUBiasMLEBin[,3]<- MLEbinned50m30[,1]-68

MUBiasMLEBin[,4]<- MLEbinned100m8[,1]-68
MUBiasMLEBin[,5]<- MLEbinned100m15[,1]-68
MUBiasMLEBin[,6]<- MLEbinned100m30[,1]-68

MUBiasMLEBin[,7]<- MLEbinned300m8[,1]-68
MUBiasMLEBin[,8]<- MLEbinned300m15[,1]-68
MUBiasMLEBin[,9]<- MLEbinned300m30[,1]-68

MUBiasMLEBin[,10]<- MLEbinned600m8[,1]-68
MUBiasMLEBin[,11]<- MLEbinned600m15[,1]-68
MUBiasMLEBin[,12]<- MLEbinned600m30[,1]-68

MUBiasMLEBin[,13]<- MLEbinned1000m8[,1]-68
MUBiasMLEBin[,14]<- MLEbinned1000m15[,1]-68
MUBiasMLEBin[,15]<- MLEbinned1000m30[,1]-68

MUaveBiasMLEBin<- colMeans(MUBiasMLEBin)


#######################################################################
#### RMSE ###########################################################################
MUMSEMLEBin<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(MUMSEMLEBin)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                           "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                           "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                           "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                           "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

MUMSEMLEBin[,1]<- (MLEbinned50m8[,1]-68)^2
MUMSEMLEBin[,2]<- (MLEbinned50m15[,1]-68)^2
MUMSEMLEBin[,3]<- (MLEbinned50m30[,1]-68)^2

MUMSEMLEBin[,4]<- (MLEbinned100m8[,1]-68)^2
MUMSEMLEBin[,5]<- (MLEbinned100m15[,1]-68)^2
MUMSEMLEBin[,6]<- (MLEbinned100m30[,1]-68)^2

MUMSEMLEBin[,7]<- (MLEbinned300m8[,1]-68)^2
MUMSEMLEBin[,8]<- (MLEbinned300m15[,1]-68)^2
MUMSEMLEBin[,9]<- (MLEbinned300m30[,1]-68)^2

MUMSEMLEBin[,10]<- (MLEbinned600m8[,1]-68)^2
MUMSEMLEBin[,11]<- (MLEbinned600m15[,1]-68)^2
MUMSEMLEBin[,12]<- (MLEbinned600m30[,1]-68)^2

MUMSEMLEBin[,13]<- (MLEbinned1000m8[,1]-68)^2
MUMSEMLEBin[,14]<- (MLEbinned1000m15[,1]-68)^2
MUMSEMLEBin[,15]<- (MLEbinned1000m30[,1]-68)^2

MUaveMSEMLEBin<- colMeans(MUMSEMLEBin)
MURMESMLEBin<- sqrt(colMeans(MUMSEMLEBin))


#######################################################################
#### VARIANCE##########################################################
##### BIAS #######
VBiasMLEBin<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(VBiasMLEBin)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                           "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                           "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                           "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                           "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

VBiasMLEBin[,1]<- MLEbinned50m8[,2]-6.25
VBiasMLEBin[,2]<- MLEbinned50m15[,2]-6.25
VBiasMLEBin[,3]<- MLEbinned50m30[,2]-6.25

VBiasMLEBin[,4]<- MLEbinned100m8[,2]-6.25
VBiasMLEBin[,5]<- MLEbinned100m15[,2]-6.25
VBiasMLEBin[,6]<- MLEbinned100m30[,2]-6.25

VBiasMLEBin[,7]<- MLEbinned300m8[,2]-6.25
VBiasMLEBin[,8]<- MLEbinned300m15[,2]-6.25
VBiasMLEBin[,9]<- MLEbinned300m30[,2]-6.25

VBiasMLEBin[,10]<- MLEbinned600m8[,2]-6.25
VBiasMLEBin[,11]<- MLEbinned600m15[,2]-6.25
VBiasMLEBin[,12]<- MLEbinned600m30[,2]-6.25

VBiasMLEBin[,13]<- MLEbinned1000m8[,2]-6.25
VBiasMLEBin[,14]<- MLEbinned1000m15[,2]-6.25
VBiasMLEBin[,15]<- MLEbinned1000m30[,2]-6.25

VaveBiasMLEBin<- colMeans(VBiasMLEBin)


###################################################################
##### RMSE #######
VMSEMLEBin<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(VMSEMLEBin)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                          "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                          "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                          "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                          "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

VMSEMLEBin[,1]<- (MLEbinned50m8[,2]-6.25)^2
VMSEMLEBin[,2]<- (MLEbinned50m15[,2]-6.25)^2
VMSEMLEBin[,3]<- (MLEbinned50m30[,2]-6.25)^2

VMSEMLEBin[,4]<- (MLEbinned100m8[,2]-6.25)^2
VMSEMLEBin[,5]<- (MLEbinned100m15[,2]-6.25)^2
VMSEMLEBin[,6]<- (MLEbinned100m30[,2]-6.25)^2

VMSEMLEBin[,7]<- (MLEbinned300m8[,2]-6.25)^2
VMSEMLEBin[,8]<- (MLEbinned300m15[,2]-6.25)^2
VMSEMLEBin[,9]<- (MLEbinned300m30[,2]-6.25)^2

VMSEMLEBin[,10]<- (MLEbinned600m8[,2]-6.25)^2
VMSEMLEBin[,11]<- (MLEbinned600m15[,2]-6.25)^2
VMSEMLEBin[,12]<- (MLEbinned600m30[,2]-6.25)^2

VMSEMLEBin[,13]<- (MLEbinned1000m8[,2]-6.25)^2
VMSEMLEBin[,14]<- (MLEbinned1000m15[,2]-6.25)^2
VMSEMLEBin[,15]<- (MLEbinned1000m30[,2]-6.25)^2

VaveMSEMLEBin<- colMeans(VMSEMLEBin)
VRMSEMLEBin<- sqrt(colMeans(VMSEMLEBin))


################################################################################################
########### MLE UNBINNED #######################################################################
BR1000m30
################# MEAN #########################################################################
#### BIAS ##########################
MUBiasMLEUNBin<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(MUBiasMLEUNBin)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                           "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                           "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                           "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                           "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

MUBiasMLEUNBin[,1]<- BR50m8[,1]
MUBiasMLEUNBin[,2]<- BR50m15[,1]
MUBiasMLEUNBin[,3]<- BR50m30[,1]

MUBiasMLEUNBin[,4]<- BR100m8[,1]
MUBiasMLEUNBin[,5]<- BR100m15[,1]
MUBiasMLEUNBin[,6]<- BR100m30[,1]

MUBiasMLEUNBin[,7]<- BR300m8[,1]
MUBiasMLEUNBin[,8]<- BR300m15[,1]
MUBiasMLEUNBin[,9]<- BR300m30[,1]

MUBiasMLEUNBin[,10]<- BR600m8[,1]
MUBiasMLEUNBin[,11]<- BR600m15[,1]
MUBiasMLEUNBin[,12]<- BR600m30[,1]

MUBiasMLEUNBin[,13]<- BR1000m8[,1]
MUBiasMLEUNBin[,14]<- BR1000m15[,1]
MUBiasMLEUNBin[,15]<- BR1000m30[,1]

MUaveBiasMLEUNBin<- colMeans(MUBiasMLEUNBin)
######################################################################
#### RMSE ##############################################################
MUMSEMLEUNBin<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(MUMSEMLEUNBin)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                             "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                             "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                             "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                             "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

MUMSEMLEUNBin[,1]<- BR50m8[,3]
MUMSEMLEUNBin[,2]<- BR50m15[,3]
MUMSEMLEUNBin[,3]<- BR50m30[,3]

MUMSEMLEUNBin[,4]<- BR100m8[,3]
MUMSEMLEUNBin[,5]<- BR100m15[,3]
MUMSEMLEUNBin[,6]<- BR100m30[,3]

MUMSEMLEUNBin[,7]<- BR300m8[,3]
MUMSEMLEUNBin[,8]<- BR300m15[,3]
MUMSEMLEUNBin[,9]<- BR300m30[,3]

MUMSEMLEUNBin[,10]<- BR600m8[,3]
MUMSEMLEUNBin[,11]<- BR600m15[,3]
MUMSEMLEUNBin[,12]<- BR600m30[,3]

MUMSEMLEUNBin[,13]<- BR1000m8[,3]
MUMSEMLEUNBin[,14]<- BR1000m15[,3]
MUMSEMLEUNBin[,15]<- BR1000m30[,3]

MUaveMSEMLEUNBin<- colMeans(MUMSEMLEUNBin)
MURMSEMLEUNBin<- sqrt(colMeans(MUMSEMLEUNBin))


#######################################################################################
########### VARIANCE #################################################################
#### BIAS ##########################
VBiasMLEUNBin<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(VBiasMLEUNBin)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                             "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                             "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                             "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                             "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

VBiasMLEUNBin[,1]<- BR50m8[,2]
VBiasMLEUNBin[,2]<- BR50m15[,2]
VBiasMLEUNBin[,3]<- BR50m30[,2]

VBiasMLEUNBin[,4]<- BR100m8[,2]
VBiasMLEUNBin[,5]<- BR100m15[,2]
VBiasMLEUNBin[,6]<- BR100m30[,2]

VBiasMLEUNBin[,7]<- BR300m8[,2]
VBiasMLEUNBin[,8]<- BR300m15[,2]
VBiasMLEUNBin[,9]<- BR300m30[,2]

VBiasMLEUNBin[,10]<- BR600m8[,2]
VBiasMLEUNBin[,11]<- BR600m15[,2]
VBiasMLEUNBin[,12]<- BR600m30[,2]

VBiasMLEUNBin[,13]<- BR1000m8[,2]
VBiasMLEUNBin[,14]<- BR1000m15[,2]
VBiasMLEUNBin[,15]<- BR1000m30[,2]

VaveBiasMLEUNBin<- colMeans(VBiasMLEUNBin)


######################################################################

#### RMSE ##############################################################
VMSEMLEUNBin<- matrix(rep(0,15*500),ncol=15,nrow=500)
colnames(VMSEMLEUNBin)<- c("n=50,bins=8","n=50,bins=15","n=50,bins=30",
                            "n=100,bins=8","n=100,bins=15","n=100,bins=30",
                            "n=300,bins=8","n=300,bins=15","n=300,bins=30",
                            "n=600,bins=8","n=600,bins=15","n=600,bins=30",
                            "n=1000,bins=8","n=1000,bins=15","n=1000,bins=30")

VMSEMLEUNBin[,1]<- BR50m8[,4]
VMSEMLEUNBin[,2]<- BR50m15[,4]
VMSEMLEUNBin[,3]<- BR50m30[,4]

VMSEMLEUNBin[,4]<- BR100m8[,4]
VMSEMLEUNBin[,5]<- BR100m15[,4]
VMSEMLEUNBin[,6]<- BR100m30[,4]

VMSEMLEUNBin[,7]<- BR300m8[,4]
VMSEMLEUNBin[,8]<- BR300m15[,4]
VMSEMLEUNBin[,9]<- BR300m30[,4]

VMSEMLEUNBin[,10]<- BR600m8[,4]
VMSEMLEUNBin[,11]<- BR600m15[,4]
VMSEMLEUNBin[,12]<- BR600m30[,4]

VMSEMLEUNBin[,13]<- BR1000m8[,4]
VMSEMLEUNBin[,14]<- BR1000m15[,4]
VMSEMLEUNBin[,15]<- BR1000m30[,4]

VaveMSEMLEUNBin<- colMeans(VMSEMLEUNBin)
VRMSEMLEUNBin<- sqrt(colMeans(VMSEMLEUNBin))


#################################################################################################
cbind(VRMSEMLEBin,VRMSEEM,VRMSEMCEM)
A<- rbind(VRMSEMLEBin,VRMSEEM,VRMSEMCEM)
A
library(xtable)
xtable(A,digits = 6)

VRMSEMLEUNBin[1:3]
cbind(VRMSEMLEUNBin[1:3],VRMSEMLEBin[1:3],VRMSEEM[1:3],VRMSEMCEM[1:3])

cbind(MUaveBiasMLEUNBin,MUaveBiasMLEBin,MUaveBiasEM,MUaveBiasMCEM)
cbind(MURMSEMLEUNBin,MURMESMLEBin,MURMESEEM,MURMSEMCEM)




a50bV<- c(VaveBiasMLEUNBin[1:3],VaveBiasMLEBin[1:3],VaveBiasEM[1:3],VaveBiasMCEM[1:3])
a50RV<- c(VRMSEMLEUNBin[1:3],VRMSEMLEBin[1:3],VRMSEEM[1:3],VRMSEMCEM[1:3])
a50V<- cbind(a50bV,a50RV)

a100bV<- c(VaveBiasMLEUNBin[4:6],VaveBiasMLEBin[4:6],VaveBiasEM[4:6],VaveBiasMCEM[4:6])
a100RV<- c(VRMSEMLEUNBin[4:6],VRMSEMLEBin[4:6],VRMSEEM[4:6],VRMSEMCEM[4:6])
a100V<- cbind(a100bV,a100RV)



a300bV<- c(VaveBiasMLEUNBin[7:9],VaveBiasMLEBin[7:9],VaveBiasEM[7:9],VaveBiasMCEM[7:9])
a300RV<- c(VRMSEMLEUNBin[7:9],VRMSEMLEBin[7:9],VRMSEEM[7:9],VRMSEMCEM[7:9])
a300V<- cbind(a300bV,a300RV)


a600bV<- c(VaveBiasMLEUNBin[10:12],VaveBiasMLEBin[10:12],VaveBiasEM[10:12],VaveBiasMCEM[10:12])
a600RV<- c(VRMSEMLEUNBin[10:12],VRMSEMLEBin[10:12],VRMSEEM[10:12],VRMSEMCEM[10:12])
a600V<- cbind(a600bV,a600RV)


a1000bV<- c(VaveBiasMLEUNBin[13:15],VaveBiasMLEBin[13:15],VaveBiasEM[13:15],VaveBiasMCEM[13:15])
a1000RV<- c(VRMSEMLEUNBin[13:5],VRMSEMLEBin[13:15],VRMSEEM[13:15],VRMSEMCEM[13:15])
a1000V<- cbind(a1000bV,a1000RV)




############################################################################
a50<- cbind(MUaveBiasMLEUNBin[1:3],MURMSEMLEUNBin[1:3],MUaveBiasMLEBin[1:3],MURMESMLEBin[1:3],
  MUaveBiasEM[1:3],MURMESEEM[1:3],MUaveBiasMCEM[1:3],MURMSEMCEM[1:3])
a50

a50bias<- c(MUaveBiasMLEUNBin[1:3],MUaveBiasMLEBin[1:3],MUaveBiasEM[1:3],MUaveBiasMCEM[1:3])
a50bias

a50RMSE<- c(MURMSEMLEUNBin[1:3],MURMESMLEBin[1:3],MURMESEEM[1:3],MURMSEMCEM[1:3])


cbind(VRMSEMLEUNBin[1:3],VRMSEMLEBin[1:3],VRMSEEM[1:3],VRMSEMCEM[1:3])

a50<- cbind(a50bias,a50RMSE)

#########################################
a100bias<- c(MUaveBiasMLEUNBin[4:6],MUaveBiasMLEBin[4:6],MUaveBiasEM[4:6],MUaveBiasMCEM[4:6])


a100RMSE<- c(MURMSEMLEUNBin[4:6],MURMESMLEBin[4:6],MURMESEEM[4:6],MURMSEMCEM[4:6])

a100<- cbind(a100bias,a100RMSE)
##################################################
a300bias<- c(MUaveBiasMLEUNBin[7:9],MUaveBiasMLEBin[7:9],MUaveBiasEM[7:9],MUaveBiasMCEM[7:9])


a300RMSE<- c(MURMSEMLEUNBin[7:9],MURMESMLEBin[7:9],MURMESEEM[7:9],MURMSEMCEM[7:9])

a300<- cbind(a300bias,a300RMSE)

#####################################################
a600bias<- c(MUaveBiasMLEUNBin[10:12],MUaveBiasMLEBin[10:12],MUaveBiasEM[10:12],MUaveBiasMCEM[10:12])


a600RMSE<- c(MURMSEMLEUNBin[10:12],MURMESMLEBin[10:12],MURMESEEM[10:12],MURMSEMCEM[10:12])

a600<- cbind(a600bias,a600RMSE)
######################################################################
a1000bias<- c(MUaveBiasMLEUNBin[13:15],MUaveBiasMLEBin[13:15],MUaveBiasEM[13:15],MUaveBiasMCEM[13:15])


a1000RMSE<- c(MURMSEMLEUNBin[13:15],MURMESMLEBin[13:15],MURMESEEM[13:15],MURMSEMCEM[13:15])

a1000<- cbind(a1000bias,a1000RMSE)
#########################################################################################

a50bV<- c(VaveBiasMLEUNBin[1:3],VaveBiasMLEBin[1:3],VaveBiasEM[1:3],VaveBiasMCEM[1:3])
a50RV<- c(VRMSEMLEUNBin[1:3],VRMSEMLEBin[1:3],VRMSEEM[1:3],VRMSEMCEM[1:3])
a50V<- cbind(a50bV,a50RV)

a100bV<- c(VaveBiasMLEUNBin[4:6],VaveBiasMLEBin[4:6],VaveBiasEM[4:6],VaveBiasMCEM[4:6])
a100RV<- c(VRMSEMLEUNBin[4:6],VRMSEMLEBin[4:6],VRMSEEM[4:6],VRMSEMCEM[4:6])
a100V<- cbind(a100bV,a100RV)



a300bV<- c(VaveBiasMLEUNBin[7:9],VaveBiasMLEBin[7:9],VaveBiasEM[7:9],VaveBiasMCEM[7:9])
a300RV<- c(VRMSEMLEUNBin[7:9],VRMSEMLEBin[7:9],VRMSEEM[7:9],VRMSEMCEM[7:9])
a300V<- cbind(a300bV,a300RV)


a600bV<- c(VaveBiasMLEUNBin[10:12],VaveBiasMLEBin[10:12],VaveBiasEM[10:12],VaveBiasMCEM[10:12])
a600RV<- c(VRMSEMLEUNBin[10:12],VRMSEMLEBin[10:12],VRMSEEM[10:12],VRMSEMCEM[10:12])
a600V<- cbind(a600bV,a600RV)


a1000bV<- c(VaveBiasMLEUNBin[13:15],VaveBiasMLEBin[13:15],VaveBiasEM[13:15],VaveBiasMCEM[13:15])
a1000RV<- c(VRMSEMLEUNBin[13:15],VRMSEMLEBin[13:15],VRMSEEM[13:15],VRMSEMCEM[13:15])
a1000V<- cbind(a1000bV,a1000RV)



################################################################################################

a50Res<- cbind(a50,a50V)

a100Res<- cbind(a100,a100V)

a300Res<- cbind(a300,a300V)
a600Res<- cbind(a600,a600V)
a1000Res<- cbind(a1000,a1000V)


a100Res
a300Res
a600Res
a1000Res


nc<- c("Method","Bins","Mean of Bias","RMSE","Mean of Bias","RMSE")
mm<- c(rep("MLE Ignore Grouping",3),rep("MLE Exact",3),rep("EM",3),rep("MCEM",3))
cl<- rep(c(8,15,30),4)

n50out<- cbind(mm,cl,round(a50Res,5))
colnames(n50out)<- nc
n50out

n100out<- cbind(mm,cl,round(a100Res,5))
colnames(n100out)<- nc

n300out<- cbind(mm,cl,round(a300Res,5))
colnames(n300out)<- nc

n600out<- cbind(mm,cl,round(a600Res,5))
colnames(n600out)<- nc

n1000out<- cbind(mm,cl,round(a1000Res,5))
colnames(n1000out)<- nc


xtable(n50out)

xtable(n100out)
xtable(n300out)
xtable(n600out)
xtable(n1000out)

#####################################################################################

Z<- c(MURMSEMLEUNBin[13:15],MURMESMLEBin[13:15],MURMESEEM[13:15],MURMSEMCEM[13:15])
Z
###############################################################################################
#### REQUIRED OUTPUT####
Z1<- c(MURMSEMLEUNBin[1],MURMSEMLEUNBin[4],MURMSEMLEUNBin[7],MURMSEMLEUNBin[10],MURMSEMLEUNBin[13])
Z1
Z2<- c(MURMSEMLEUNBin[2],MURMSEMLEUNBin[5],MURMSEMLEUNBin[8],MURMSEMLEUNBin[11],MURMSEMLEUNBin[14])
Z2
Z3<- c(MURMSEMLEUNBin[3],MURMSEMLEUNBin[6],MURMSEMLEUNBin[9],MURMSEMLEUNBin[12],MURMSEMLEUNBin[15])
Z3

Z4<- c(MURMESMLEBin[1],MURMESMLEBin[4],MURMESMLEBin[7],MURMESMLEBin[10],MURMESMLEBin[13])
Z4
Z5<- c(MURMESMLEBin[2],MURMESMLEBin[5],MURMESMLEBin[8],MURMESMLEBin[11],MURMESMLEBin[14])
Z5
Z6<- c(MURMESMLEBin[3],MURMESMLEBin[6],MURMESMLEBin[9],MURMESMLEBin[12],MURMESMLEBin[15])
Z6

Z7<- c(MURMESEEM[1],MURMESEEM[4],MURMESEEM[7],MURMESEEM[10],MURMESEEM[13])
Z7
Z8<- c(MURMESEEM[2],MURMESEEM[5],MURMESEEM[8],MURMESEEM[11],MURMESEEM[14])
Z8
Z9<- c(MURMESEEM[3],MURMESEEM[6],MURMESEEM[9],MURMESEEM[12],MURMESEEM[15])
Z9

Z10<- c(MURMSEMCEM[1],MURMSEMCEM[4],MURMSEMCEM[7],MURMSEMCEM[10],MURMSEMCEM[13])
Z10
Z11<- c(MURMSEMCEM[2],MURMSEMCEM[5],MURMSEMCEM[8],MURMSEMCEM[11],MURMSEMCEM[14])
Z11
Z12<- c(MURMSEMCEM[3],MURMSEMCEM[6],MURMSEMCEM[9],MURMSEMCEM[12],MURMSEMCEM[15])
Z12

XX1<- c(Z4,Z7,Z10)
XX2<- c(Z5,Z8,Z11)
XX3<- c(Z6,Z9,Z12)
methods<- c(rep("MLE Exact",5),rep("EM",5),rep("MCEM",5))
n<- rep(c(50,100,300,600,1000),3)
MeanRMSEinit<- cbind(XX1,XX2,XX3)
MeanRMSE<- cbind(methods,n,round(MeanRMSEinit,5))
colnames(MeanRMSE)<- c("Method","n","bins=8","bins=15","bins=30")
MeanRMSE
library(xtable)
xtable(MeanRMSE)
####################################################################################y
#Y1<- c(VRMSEMLEUNBin[4],VRMSEMLEUNBin[7],VRMSEMLEUNBin[10],VRMSEMLEUNBin[13])
#Y1
#Y2<- c(VRMSEMLEUNBin[5],VRMSEMLEUNBin[8],VRMSEMLEUNBin[11],VRMSEMLEUNBin[14])
#Y2
#Y3<- c(VRMSEMLEUNBin[6],VRMSEMLEUNBin[9],VRMSEMLEUNBin[12],VRMSEMLEUNBin[15])
#Y3

Y4<- c(VRMSEMLEBin[1],VRMSEMLEBin[4],VRMSEMLEBin[7],VRMSEMLEBin[10],VRMSEMLEBin[13])
Y4
Y5<- c(VRMSEMLEBin[2],VRMSEMLEBin[5],VRMSEMLEBin[8],VRMSEMLEBin[11],VRMSEMLEBin[14])
Y5
Y6<- c(VRMSEMLEBin[3],VRMSEMLEBin[6],VRMSEMLEBin[9],VRMSEMLEBin[12],VRMSEMLEBin[15])
Y6

Y7<- c(VRMSEEM[1],VRMSEEM[4],VRMSEEM[7],VRMSEEM[10],VRMSEEM[13])
Y7
Y8<- c(VRMSEEM[2],VRMSEEM[5],VRMSEEM[8],VRMSEEM[11],VRMSEEM[14])
Y8
Y9<- c(VRMSEEM[3],VRMSEEM[6],VRMSEEM[9],VRMSEEM[12],VRMSEEM[15])
Y9

Y10<- c(VRMSEMCEM[1],VRMSEMCEM[4],VRMSEMCEM[7],VRMSEMCEM[10],VRMSEMCEM[13])
Y10
Y11<- c(VRMSEMCEM[2],VRMSEMCEM[5],VRMSEMCEM[8],VRMSEMCEM[11],VRMSEMCEM[14])
Y11
Y12<- c(VRMSEMCEM[3],VRMSEMCEM[6],VRMSEMCEM[9],VRMSEMCEM[12],VRMSEMCEM[15])
Y12

WW1<- c(Y4,Y7,Y10)
WW2<- c(Y5,Y8,Y11)
WW3<- c(Y6,Y9,Y12)
methods<- c(rep("MLE Exact",5),rep("EM",5),rep("MCEM",5))
n<- rep(c(50,100,300,600,1000),3)
VRMSEinit<- cbind(WW1,WW2,WW3)
VarianceRMSE<- cbind(methods,n,round(VRMSEinit,5))
colnames(VarianceRMSE)<- c("Method","n","bins=8","bins=15","bins=30")
VarianceRMSE
library(xtable)
xtable(VarianceRMSE)
####################################################################################################
VMinit1<- cbind(XX1,WW1)
VM1<- cbind(methods,n,round(VMinit1,5))
colnames(VM1)<- c("Method","n","RMSE for Means","RMSE for Variances")

VM1

VMinit2<- cbind(XX2,WW2)
VM2<- cbind(methods,n,round(VMinit2,5))
colnames(VM2)<- c("Method","n","RMSE for Means","RMSE for Variances")

VM2

VMinit3<- cbind(XX3,WW3)
VM3<- cbind(methods,n,round(VMinit3,5))
colnames(VM3)<- c("Method","n","RMSE for Means","RMSE for Variances")

VM3

xtable(VM1)
xtable(VM2)
xtable(VM3)
