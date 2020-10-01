############################################################################
######## Simplified R Code for LASSO & Correlations On Server ###############3
##############################################################################

##Load Packages
library(glmnet)
library(foreach)
library(doParallel)

load("/local/users/apaquett/PlacentalTRN/IntermediateData/TestingDataset_6192020.RData")
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataset_6192020.RData")

####### LASSO Regression ################


NCores<-20


cl <- makeCluster(NCores) #not to overload your computer #2o is too many=10 ok but check
registerDoParallel(cl)

##Write the Function


LASSO_Results<- function(Exprs_Train,Exprs_Test,TG,TFs){

  ##Set Up Training & Testing Data
  TG_Train<-as.numeric(Exprs_Train[TG,])
  TG_Test<-as.numeric(Exprs_Test[TG,])

  TFs_Train<-t(na.omit(as.matrix(Exprs_Train[TFs,]))) #can remove na.omit in real matrix
  TFs_Test<-t(na.omit(as.matrix(Exprs_Test[TFs,]))) #can remove na.omit in real matrix

  #The matrix of TFs is X, The Target gene is Y, the response variable
  cvfit = glmnet::cv.glmnet(TFs_Train, TG_Train,alpha=1)
  TG_Predicted<-predict(cvfit, newx = TFs_Test, s = "lambda.min")
  x<-c((cor(TG_Predicted,TG_Test)^2),(cor(TG_Predicted,TG_Test)),(cor.test(TG_Predicted,TG_Test)$p.value),TG)
  names(x)<-c("RSquared","Cor","PVal","TG")
  x
}



#### Model 1
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel1_Full.RData")

start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)

LASSOResults->Model1_LASSO

save(Model1_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/Model1Testing_LASSO.RData")

#Model 2
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel2_Full.RData")
start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)

LASSOResults->Model2_LASSO
save(Model2_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/Model2Testing_LASSO.RData")

#Model 3
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel3_Full.RData")
start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)

LASSOResults->Model3_LASSO

save(Model3_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/Model3Testing_LASSO.RData")


#Model 4
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel4_Full.RData")

start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)
LASSOResults->Model4_LASSO

save(Model4_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/Model4Testing_LASSO.RData")


#Model 5
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel5_Full.RData")
start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)

LASSOResults->Model5_LASSO
save(Model5_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/Model5Testing_LASSO.RData")

#Model 6
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel6_Full.RData")
start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)

LASSOResults->Model6_LASSO
save(Model6_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/Model6Testing_LASSO.RData")

#Model 7
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel7_REALBADCORS.RData")
start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)

LASSOResults->Model7_LASSO
save(Model7_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/Model7Testing_LASSO.RData")

#Model 8
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel8_FAKEGOODCORS.RData")
start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)

LASSOResults->Model8_LASSO
save(Model8_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/Model8Testing_LASSO.RData")

#Model 9
load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataModel9_FAKEBADCORS.RData")
start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)

LASSOResults->Model9_LASSO
save(Model9_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/Model9Testing_LASSO.RData")


################################################################################################
##  Run on Server: Rscript ServerCode_LASSO_TrainingIterations.R ##########################
################################################################################################

