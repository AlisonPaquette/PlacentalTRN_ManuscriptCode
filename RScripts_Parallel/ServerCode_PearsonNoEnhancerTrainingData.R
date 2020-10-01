############################################################################
######## Simplified R Code for LASSO & Correlations On Server ###############
##############################################################################

##Load Packages
library(glmnet)
library(foreach)
library(doParallel)

load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataset_6192020.RData")
Exprs<-Exprs_Train

load("~/PlacentalTRN/PreprocessedData/NoEnhancerProcessedModel.RData")


####### LASSO Regression ################


NCores<-20

cl <- makeCluster(NCores) #not to overload your computer #2o is too many=10 ok but check
registerDoParallel(cl)

##Write the Function

####### Correlations ################
#Function: For each TF/TRN Network
Corr_Results<- function(Exprs,TG,TFs){
  
  Corr<-data.frame(matrix(NA,nrow=length(TFs),ncol=5))
  colnames(Corr)<-c("Cor","P","Target_Gene","TF","Rank")
  for (i in 1:length(TFs)){
      Training<-cor.test(as.numeric(Exprs[TG,]),as.numeric(Exprs[TFs[i],]))
      
      Corr[i,1]<-Training$estimate
      Corr[i,2]<-Training$p.value
    }
    Corr$Target_Gene<-rep(TG,length(TFs))
    Corr$TF<-TFs
    Corr<-Corr[order(abs(Corr$Cor),decreasing=T),]
    Corr$Rank<-1:length(TFs)
  Corr
}

NoEnhancer_Cor<- foreach(j=1:length(names(TFsList)), .combine=rbind,.errorhandling='remove') %dopar% {
  tempMatrix = Corr_Results(Exprs,names(TFsList)[j],TFsList[[j]]) #calling a function
  #do other things if you want
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}

#stop cluster
stopCluster(cl)


save(NoEnhancer_Cor,file="/users/apaquett/PlacentalTRN/IntermediateData/NoEnhancerCorrelationsTrainingData.RData")

