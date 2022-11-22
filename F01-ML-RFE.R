library(pROC)
library(caret)
require(glmnet)
require(tidyverse)
library(xgboost)
library(randomForest)
source("svm-ref.R")
RNGkind(sample.kind = "Rejection")
###calculate machine learning performance
ml.perf <- function(rfe_fit,X.test,Y.test)
{
  pred = predict(rfe_fit$fit,newdata = X.test,type="prob")
  
  rocm <- roc(Y.test,pred[,2])
  thre = coords(rocm, "best", ret = "threshold")
  
  if(is.data.frame(thre)){opt.thre = thre[1,]} else {opt.thre = as.numeric(thre)}
  
  #opt.thre = as.numeric(coords(rocm, "best", ret = "threshold",transpose = TRUE))
  pred.class <- ifelse(pred[,2] > opt.thre,1,0)
  dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
  T = confusionMatrix(dt)
  perf = c(T$overall[1],T$byClass[1:2])
  return(perf)
}

## run rfe model for each algorithm
runml <- function(j,X1)
{
  rfe_control = rfeControl(functions = caretFuncs, #caretFuncs here
                           method="cv",
                           number=10,saveDetails = TRUE)
  fit_control = trainControl(number = 10,method = "cv",
                             search="grid",classProbs = TRUE,summaryFunction = twoClassSummary)
  
  tunegrid1 <- expand.grid(.mtry=c(1:5,7,9)) ##rf
  tunegrid2 = expand.grid(alpha = seq(0.01,0.7,length.out = 20), lambda = seq(0.01,1,length.out = 25))
  tunegrid3 <- expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0.5,1))
  tunegrid4 <- expand.grid(k = c(5, 11, 21, 25))
  tunegrid5 =   expand.grid(
    nrounds =500,
    eta =  0.3,
    lambda = c(0.5,1),
    alpha = c(0,0.5,1)
  )
  
  
  done = FALSE
  while(done == FALSE)
  {
    while(j < 1000)
    {
      set.seed(252+7*j)
      train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
      X.train = data.frame(X1[train_ind,-1])
      Y.train = as.factor(X1[train_ind,1])
      X.test = data.frame(X1[-train_ind,-1])
      Y.test = as.factor(X1[-train_ind,1])
      
      set.seed(252+7*j)
      system.time(rfe_fit_rf <- rfe(x = X.train, y = as.factor(make.names(Y.train)),
                                    sizes = c(2:30,35,40), method = "rf",ntree = 1000, metric = 'ROC',
                                    rfeControl = rfe_control,trControl = fit_control,tuneGrid = tunegrid1))
      
      perf_rf = ml.perf(rfe_fit_rf,X.test,Y.test)
      final_rf = rfe_fit_rf$optVariables
      if(any(perf_rf==0))
      {
        j = j+100
        break;
      }
      
      set.seed(252+7*j)
      system.time(rfe_fit_nb <- rfe(x = X.train, y = as.factor(make.names(Y.train)),
                                    sizes = c(2:30,35,40), method = "nb",metric = 'ROC',
                                    rfeControl = rfe_control,trControl = fit_control,tuneGrid = tunegrid3))
      perf_nb = ml.perf(rfe_fit_nb,X.test,Y.test)
      final_nb = rfe_fit_nb$optVariables
      
      if(any(perf_nb==0))
      {
        j = j+100
        break;
      }
      set.seed(252+7*j)
      system.time(rfe_fit_knn <- rfe(x = X.train, y = as.factor(make.names(Y.train)),metric = 'ROC',
                                     sizes = c(2:30,35,40), method = "knn",
                                     rfeControl = rfe_control,trControl = fit_control,tuneGrid = tunegrid4))
      perf_knn = ml.perf(rfe_fit_knn,X.test,Y.test)
      final_knn = rfe_fit_knn$optVariables
      
      if(any(perf_knn==0))
      {
        j = j+100
        break;
      }
      set.seed(252+7*j)
      system.time(rfe_fit_xgb <- rfe(x = X.train, y = as.factor(make.names(Y.train)),metric = 'ROC',
                                     sizes = c(2:30,35,40), method = "xgbLinear",
                                     rfeControl = rfe_control,trControl = fit_control,tuneGrid = tunegrid5))
      final_xgb = rfe_fit_xgb$optVariables
      perf_xgb = ml.perf(rfe_fit_xgb,X.test[,final_xgb],Y.test)
      
      if(any(perf_xgb==0))
      {
        j = j+100
        break;
      }
      set.seed(252+7*j)
      system.time(rfe_fit_els <- rfe(x = X.train, y = as.factor(make.names(Y.train)),
                                     sizes = c(2:30,35,40), method = "glmnet",family = "binomial",metric = 'ROC',
                                     rfeControl = rfe_control,trControl = fit_control,tuneGrid = tunegrid2))
      
      perf_els = ml.perf(rfe_fit_els,X.test,Y.test)
      final_els = rfe_fit_els$optVariables
      
      if(any(perf_els==0))
      {
        j = j+100
        break;
      } 
      
      if(all(c(perf_rf,perf_els,perf_knn,perf_nb,perf_xgb)!= 0))
      {
        done = TRUE
        break;
      }
    }
    
  }
  
  dt <- rbind.data.frame(perf_rf,perf_els,perf_knn,perf_nb,perf_xgb)
  rownames(dt) = c("rf","Elastic net","knn","naive Bayes","xgboost")
  dt1 = list(rfe_fit_rf,rfe_fit_els,rfe_fit_knn,rfe_fit_nb,rfe_fit_xgb)
  return(list(performance = dt,allrfe = dt1))
}

#res = runml(j,X1=X.binary.f)

