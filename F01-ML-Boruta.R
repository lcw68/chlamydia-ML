library(Boruta)

ml.perf <- function(fit,X.test,final.select,Y.test)
{
  if(is.null(dim(X.test)))
  {
    X.test = data.frame(X.test)
    colnames(X.test) = final.select
    pred = predict(fit,newdata = X.test,type="prob")
    
  }else{
    pred = predict(fit,newdata = X.test,type="prob")
  }
  
  
  rocm <- roc(Y.test,pred[,2])
  aucm = as.numeric(rocm$auc)
  thre = coords(rocm, "best", ret = "threshold")
  
  if(is.data.frame(thre)){opt.thre = thre[1,]} else {opt.thre = as.numeric(thre)}
  
  #opt.thre = as.numeric(coords(rocm, "best", ret = "threshold",transpose = TRUE))
  pred.class <- ifelse(pred[,2] > opt.thre,1,0)
  dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
  T = confusionMatrix(dt)
  perf = c(T$overall[1],T$byClass[1:2],aucm)
  return(list(perf,rocm))
}

###runBoruta function
runburota <- function(j,X1)
{
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
  
  
  set.seed(252+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = data.frame(X1[train_ind,-1])
  Y.train = as.factor(X1[train_ind,1])
  X.test = data.frame(X1[-train_ind,-1])
  Y.test = as.factor(X1[-train_ind,1])
  set.seed(252+7*j)
  
  boruta.mod = Boruta(X.train,Y.train,doTrace = 2,pValue = 0.05,mcAdj = FALSE,maxRuns = 200)
  final.boruta <- TentativeRoughFix(boruta.mod)
  final.select = getSelectedAttributes(final.boruta, withTentative = F)
  if(length(final.select)==0)
  {
    perf_rf = perf_els = perf_knn = perf_nb = perf_xgb = rep(NA,4)
  }else{
    X.use = X1[train_ind,]
    X.use <- X.use %>% mutate(yinfect = factor(yinfect, 
                                               labels = make.names(levels(yinfect))))
    
    final.fit.rf <- train(formula(paste0("yinfect~",paste0(final.select,collapse = "+"))),
                          data = X.use,method = "rf",ntree = 1000,
                          trControl = fit_control,tuneGrid = tunegrid1,metric = "ROC")
    perf_rf = ml.perf(final.fit.rf,X.test[,final.select],final.select,Y.test)
    
    final.fit.nb <- train(formula(paste0("yinfect~",paste0(final.select,collapse = "+"))),
                          data = X.use,method = "nb",
                          trControl = fit_control,tuneGrid = tunegrid3,metric = "ROC")
    perf_nb = ml.perf(final.fit.nb,X.test[,final.select],final.select,Y.test)
    
    # final.fit.els <- train(formula(paste0("yinfect~",paste0(final.select,collapse = "+"))),
    #                       data = X.use,method = "glmnet",family = "binomial",
    #                       trControl = fit_control,tuneGrid = tunegrid2,metric = "ROC")
    # perf_els = ml.perf(final.fit.els,X.test[,final.select],final.select,Y.test)
    # 
    
    final.fit.knn <- train(formula(paste0("yinfect~",paste0(final.select,collapse = "+"))),
                           data = X.use,method = "knn",
                           trControl = fit_control,tuneGrid = tunegrid4,metric = "ROC")
    perf_knn = ml.perf(final.fit.knn,X.test[,final.select],final.select,Y.test)
    
    
    final.fit.xgb <- train(formula(paste0("yinfect~",paste0(final.select,collapse = "+"))),
                           data = X.use,method = "xgbLinear",
                           trControl = fit_control,tuneGrid = tunegrid5,metric = "ROC")
    perf_xgb = ml.perf(final.fit.xgb,X.test[,final.select],final.select,Y.test)
    
    
    #tp = list(selection = final.select,performance = rfperf)
  }
  dt <- rbind.data.frame(perf_rf[[1]],perf_knn[[1]],perf_nb[[1]],perf_xgb[[1]])
  rownames(dt) = c("rf","knn","naive Bayes","xgboost")
  colnames(dt) = c( "Accuracy", "Sensitivity", "Specificity","AUC")
  #dt1 = list(rfe_fit_rf,rfe_fit_els,rfe_fit_knn,rfe_fit_nb,rfe_fit_xgb)
  return(list(performance = dt,selection = final.select))
  #return(tp)
}