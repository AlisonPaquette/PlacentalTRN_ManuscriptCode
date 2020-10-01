#############################################################################
######## Simplified R Code for LASSO & Correlations On Server ###############3
##############################################################################

##Load Packages
library(randomForest)
library(foreach)
library(doParallel)
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataset_6192020.RData")
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TestingDataset_6192020.RData")

load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel1_Full.RData")





####### LASSO Regression ################
NCores<-20

cl <- makeCluster(NCores) #not to overload your computer #2o is too many=10 ok but check
registerDoParallel(cl)

##Write the Function


RF_Results<- function(Exprs_Train,Exprs_Test,TG,TFs){

  ##Set Up Training & Testing Data
  TG_Train<-as.numeric(Exprs_Train[TG,])
  TG_Test<-as.numeric(Exprs_Test[TG,])

  TFs_Train<-t(na.omit(as.matrix(Exprs_Train[TFs,]))) #can remove na.omit in real matrix
  TFs_Test<-t(na.omit(as.matrix(Exprs_Test[TFs,]))) #can remove na.omit in real matrix

  fit <- randomForest::randomForest( x = TFs_Train, y = TG_Train )
  TG_Predicted<- predict(fit, TFs_Test)
  x<-c((cor(TG_Predicted,TG_Test)^2),(cor(TG_Predicted,TG_Test)),(cor.test(TG_Predicted,TG_Test)$p.value),TG)
  names(x)<-c("RSquared","Cor","PVal","TG")
  x
}


# Parellelized Code
start<-c(NA,NA,NA,NA)

TRNFiltTest_RF<- foreach(j=1:length(names(TFsList)), .combine=rbind,.init=start) %dopar% {
#TRNFiltTest_RF<- foreach(j=1:20, .combine=rbind,.init=start) %dopar% {
  tempMatrix =RF_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
Model1_RF<-na.omit(TRNFiltTest_RF)
#stop cluster
stopCluster(cl)

save(Model1_RF,file="~/PlacentalTRN/IntermediateData/Model1RF.RData")



