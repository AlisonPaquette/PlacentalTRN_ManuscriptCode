############################################################################
######## Simplified R Code for LASSO & Correlations On Server ###############3
##############################################################################

##Load Packages
library(glmnet)
library(foreach)
library(doParallel)

load("/users/apaquett/PlacentalTRNTestingJuly2019/TestCode/TestingDataset363019.RData")
load("/users/apaquett/PlacentalTRNTestingJuly2019/TestCode/TrainingDataset363019.RData")


#load("/users/apaquett/PlacentalTRN_TestLASSO/Model1_Processed.RData")
#load("/users/apaquett/PlacentalTRN_TestLASSO/Model2_Processed.RData")
#load("/users/apaquett/PlacentalTRN_TestLASSO/Model3_Processed.RData")
#load("/users/apaquett/PlacentalTRN_TestLASSO/Model4_Processed.RData")
#load("/users/apaquett/PlacentalTRN_TestLASSO/Model5_Processed.RData")
#load("/users/apaquett/PlacentalTRN_TestLASSO/Model6_Processed.RData")
#load("/users/apaquett/PlacentalTRN_TestLASSO/Model7_Null_TFList.RData")
#load("/users/apaquett/PlacentalTRN_TestLASSO/Model8_Null_TFList.RData")
load("/users/apaquett/PlacentalTRN_TestLASSO/Model9_Null_TFList.RData")

#load("/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model2_Processed.RData")
#load("/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model3_Processed.RData")
#load("/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model4_Processed.RData")
#load("/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model5_Processed.RData")
#load("/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model6_Processed.RData")


####### LASSO Regression ################


NCores<-10

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


# Parellelized Code
start<-c(NA,NA,NA,NA)
#TRNFiltTest_Cor<- foreach(j=1:length(names(TRNFiltTest_TGsList)), .combine=rbind,.errorhandling='remove') %dopar% {
TRNFiltTest_LASSO<- foreach(j=1:length(names(TRNFiltTest_TFsList)), .combine=rbind,.init=start) %dopar% {
  tempMatrix =LASSO_Results(Exprs_Train,Exprs_Test,names(TRNFiltTest_TFsList)[j],TRNFiltTest_TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
TRNFiltTest_LASSO<-na.omit(TRNFiltTest_LASSO)
#stop cluster/users/apaquett/PlacentalTRNTestingJuly2019/Testing_August/
stopCluster(cl)

#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRN_TestLASSO/Model1_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRN_TestLASSO/Model2_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRN_TestLASSO/Model3_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRN_TestLASSO/Model4_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRN_TestLASSO/Model5_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRN_TestLASSO/Model6_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRN_TestLASSO/Model7Null_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRN_TestLASSO/Model8Null_LASSO.RData")
save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRN_TestLASSO/Model9Null_LASSO.RData")

#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model2_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model3_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model4_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model5_LASSO.RData")
#save(TRNFiltTest_LASSO,file="/users/apaquett/PlacentalTRNTestingJuly2019/Testing_UWTalk/Model6_LASSO.RData")
