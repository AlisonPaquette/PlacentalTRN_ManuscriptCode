############################################################################
######## Simplified R Code for LASSO & Correlations On Server ###############
##############################################################################

##Load Packages
library(glmnet)
library(foreach)
library(doParallel)

load("/local/users/apaquett/PlacentalTRN/IntermediateData/TrainingDataset_6192020.RData")
Exprs<-Exprs_Train

load("~/PlacentalTRN/PreprocessedData/EnhancerProcessedModel.RData")


####### LASSO Regression ################


NCores<-20 #if nobody is on khaleesi

cl <- makeCluster(NCores) 
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

Enhancer_Cor<- foreach(j=1:length(names(TFsList)), .combine=rbind,.errorhandling='remove') %dopar% {
  tempMatrix = Corr_Results(Exprs,names(TFsList)[j],TFsList[[j]]) #calling a function
  tempMatrix 
}

#stop cluster
stopCluster(cl)


save(Enhancer_Cor,file="/users/apaquett/PlacentalTRN/IntermediateData/EnhancerCorrelationsTrainingData.RData")

