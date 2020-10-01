############################################################################
######## Simplified R Code for LASSO & Correlations On Server ###############3
##############################################################################

##Load Packages
library(glmnet)
library(foreach)
library(doParallel)

load("/local/users/apaquett/PlacentalTRN/IntermediateData/TTDataset_6192020.RData")
load("/local/users/apaquett/PlacentalTRN/IntermediateData/HoldoutDataset_6192020.RData")

load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingTesting_ProcessedModelForHoldout.RData")

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
start<-c(NA,NA,NA,NA)
LASSOResults<- foreach(j=1:length(TFsList), .combine=rbind,.init=start,.errorhandling='remove') %dopar% {
  tempMatrix =LASSO_Results(Exprs_TT,Exprs_HO,names(TFsList)[j],TFsList[[j]])
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
LASSOResults<-na.omit(LASSOResults)

LASSOResults->Holdout_LASSO

save(Holdout_LASSO,file="/users/apaquett/PlacentalTRN/IntermediateData/HoldoutData_LASSO.RData")


####### Correlations ################
#Function: For each TF/TRN Network
Corr_Results<- function(Exprs_Train,Exprs_Test,TF,Genes){
  
  Corr<-data.frame(matrix(NA,nrow=length(Genes),ncol=6))
  colnames(Corr)<-c("Cor_Train","Cor_Test","P_Train","P_Test","Target_Gene","TF")
  for (i in 1:length(Genes)){
    
    if(length(intersect(rownames(Exprs_Train),Genes[i])) == 0){
      Corr[i,]<-NA
    } else {
      
      Training<-cor.test(as.numeric(Exprs_Train[TF,]),as.numeric(Exprs_Train[Genes[i],]))
      Testing<-cor.test(as.numeric(Exprs_Test[TF,]),as.numeric(Exprs_Test[Genes[i],]))
      Corr[i,1]<-Training$estimate
      Corr[i,2]<-Testing$estimate
      Corr[i,3]<-Training$p.value
      Corr[i,4]<-Testing$p.value
    }
    Corr$Target_Gene<-Genes
    Corr$TF<-rep(TF,length(Genes))
  }
  Corr
}


# Parellelized Code
Holdout_Cor<- foreach(j=1:length(names(TGsList)), .combine=rbind,.errorhandling='remove') %dopar% {
  tempMatrix = Corr_Results(Exprs_TT,Exprs_HO,names(TGsList)[j],TGsList[[j]]) #calling a function
  #do other things if you want
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}

#stop cluster
stopCluster(cl)


save(Holdout_Cor,file="/users/apaquett/PlacentalTRN/IntermediateData/HoldoutModel_CorrResults.RData")


