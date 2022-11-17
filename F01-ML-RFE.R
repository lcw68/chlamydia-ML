library(pROC)
library(caret)
require(glmnet)
require(tidyverse)
library(xgboost)
library(randomForest)
source("svm-ref.R")
RNGkind(sample.kind = "Rejection")

### interface function
ml.calculate <- function(j,X1,method)
{
  require("caret")
  require("glmnet")
  set.seed(232+7*j)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = data.frame(X1[train_ind,-1])
  Y.train = as.factor(X1[train_ind,1])
  X.use <- X1[train_ind,]
  if(method == "Elastic")
  {
    X.train = as.matrix(X1[train_ind,-1])
    Y.train = X1[train_ind,1]
    trctrl = trainControl(method = "cv", number = 10, search ="grid")
    alpha  <- seq(0.01,1,length.out = 40)
    lambda <- seq(0.01,1,length.out = 40)
    seed1= 232+7*j
    x = NA
    while(all(is.na(x)))
    {
      set.seed(seed1)
      dw_elnet <- train(X.train, Y.train, 
                        method = "glmnet", trControl = trctrl, family = "binomial", 
                        tuneGrid = expand.grid(alpha = alpha, lambda = lambda))
      
      imp.mat <- varImp(dw_elnet,lambda = dw_elnet$bestTune$lambda,alpha = dw_elnet$bestTune$alpha)$importance
      params = dw_elnet$bestTune
      x = imp.mat[,1]
      names(x) = rownames(imp.mat)
      seed1 = seed1 + 1
    }
    
    
    rank.score = sort(x,decreasing = TRUE)
    ge = 0
    for(k in 2:min(40,length(rank.score)))
    {
      sp =names(rank.score)[1:k]
      X.train <- X1[train_ind,sp]
      Y.train = as.factor(X1[train_ind,1])
      X.use <- X1[train_ind,c("yinfect",sp)]
      set.seed(1034+35*k)
      X.use <- X.use %>% mutate(yinfect = factor(yinfect, 
                                                 labels = make.names(levels(yinfect))))
      dw_elnet <- train(yinfect~.,data=X.use, 
                        method = "glmnet", trControl = trctrl, family = "binomial", 
                        tuneGrid = data.frame(params))
      ge[k] <- as.numeric(get_best_result(dw_elnet)[3])
      #ge[j] = unlist(calculate_num(j,sp,X,method = methods))
    }
    final.select <- names(rank.score)[1:(which.max(ge[-1])+1)]
    X.train = as.matrix(X1[train_ind,final.select])
    Y.train = X1[train_ind,1]
    X.test <- as.matrix(X1[-train_ind,final.select])
    Y.test = X1[-train_ind,1]
    alpha  <- seq(0.4,1,length.out = 50)
    lambda <- seq(0,1,length.out = 50)
    set.seed(232+7*j)
    dw_elnet1 <- train(X.train, Y.train, 
                       method = "glmnet", trControl = trctrl, family = "binomial", 
                       tuneGrid = expand.grid(alpha = alpha, lambda = lambda))
    
    model0 <- glmnet(X.train, Y.train, family = "binomial", 
                     alpha = dw_elnet1$bestTune$alpha, lambda = dw_elnet1$bestTune$lambda)
    
    pred <- predict(model0, newx = X.test, type = "response",
                    alpha = dw_elnet1$bestTune$alpha, lambda = dw_elnet1$bestTune$lambda)
    rocm <- roc(Y.test,pred[,1])
    thre = coords(rocm, "best", ret = "threshold")
    
    if(is.data.frame(thre)){opt.thre = thre[1,]} else {opt.thre = as.numeric(thre)}
    
    #opt.thre = as.numeric(coords(rocm, "best", ret = "threshold",transpose = TRUE))
    pred.class <- ifelse(pred[,1] > opt.thre,1,0)
    dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
    T = confusionMatrix(dt)
    perf = c(T$overall[1],T$byClass[1:2])
    return(list(score = rank.score, performance = perf, final.var = final.select))
  }
  else if (method == "rf")
  {
    #### rank var importance
    trctrl = trainControl(method = "cv", number = 10, search ="grid")
    tunegrid <- expand.grid(.mtry=c(1:7))
    set.seed(232+7*j)
    caret_res <- train(X.train, Y.train, method="rf", metric="Accuracy", 
                       tuneGrid=tunegrid, ntree = 1000, trControl=trctrl)
    #rfImp1 <- randomForest(yinfect ~ ., data = X1[train_ind,], ntree = 1000, importance = TRUE)
    rfs = varImp(caret_res)$importance
    mtry = caret_res$bestTune
    x = rfs[,1]
    names(x) = rownames(rfs)
    rank.score = sort(x,decreasing = TRUE)
    
    ### optimal number
    ge = 0
    for(k in 2:min(40,length(rank.score)))
    {
      sp =names(rank.score)[1:k]
      X.train <- X1[train_ind,sp]
      Y.train = as.factor(X1[train_ind,1])
      X.use <- X1[train_ind,c("yinfect",sp)]
      set.seed(1034+35*k)
      X.use <- X.use %>% mutate(yinfect = factor(yinfect, 
                                                 labels = make.names(levels(yinfect))))
      
      rftrain <- train(yinfect~.,data=X.use, 
                       method = "rf", trControl = trctrl, ntree = 100, 
                       tuneGrid = data.frame(mtry))
      ge[k] <- rftrain$results[,"Accuracy"]
      #ge[j] = unlist(calculate_num(j,sp,X,method = methods))
    }
    final.select <- names(rank.score)[1:(which.max(ge[-1])+1)]
    X.train = as.matrix(X1[train_ind,final.select])
    Y.train = X1[train_ind,1]
    X.test <- as.matrix(X1[-train_ind,final.select])
    Y.test = X1[-train_ind,1]
    
    tunegrid <- expand.grid(.mtry=seq(1,length(final.select),2))
    set.seed(232+7*j)
    
    trmod = train(X.train,Y.train,ntree=1000,method = "rf",tuneGrid = tunegrid,trControl = trctrl,metric='Accuracy')
    pred = predict(trmod,newdata = X.test,type = "prob")
    rocm <- roc(Y.test,pred[,2])
    thre = coords(rocm, "best", ret = "threshold")
    if(is.data.frame(thre)){opt.thre = thre[1,]} else {opt.thre = as.numeric(thre)}
    #opt.thre = as.numeric(coords(rocm, "best", ret = "threshold"))
    pred.class <- ifelse(pred[,2] > opt.thre,1,0)
    dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
    T = confusionMatrix(dt)
    perf = c(T$overall[1],T$byClass[1:2])
    return(list(score = rank.score, performance = perf, final.var = final.select))
  }
  else if(method == "xgboost")
  {
    # xgtrain <- xgb.DMatrix(data = as.matrix(X.train),label = as.numeric(Y.train)-1)
    # #trcon = trainControl(method="repeatedcv", number=5,repeats = 3 ,search="grid")
    # xgb1 = xgb.cv(data = xgtrain,eval_metric = "error",
    #               early_stopping_rounds= 50,nfold = 10, booster = "gblinear",nrounds=500, watchlist=list(data = xgtrain), objective = "binary:logistic",verbose = FALSE)
    # 
    # #xgb1 = xgb.train(data = xgtrain, booster = "gblinear",nrounds=500, watchlist=list(data = xgtrain), objective = "binary:logistic",verbose = FALSE)
    # xgi <- xgb.importance(model = xgb1)
    #xgi = xgi[match(colnames(X.train),xgi$Feature),]
    
    tunegrid <- expand.grid(
      nrounds =500,
      eta =  0.3,
      lambda = c(0.5,1),
      alpha = c(0,0.5,1)
    )
    trcon = trainControl(method="cv", number=10,search="grid")
    set.seed(232+7*j)
    trmod = train(X.train,Y.train,method="xgbLinear",trControl = trcon,tuneGrid = tunegrid)
    imps = varImp(trmod)$importance
    params = trmod$bestTune
    x = imps[,1]
    names(x) = rownames(imps)
    rank.score = sort(x,decreasing = TRUE)
    
    ### optimal number
    ge = 0
    for(k in 2:min(40,length(rank.score)))
    {
      sp =names(rank.score)[1:k]
      X.train <- X1[train_ind,sp]
      #Y.train = as.factor(X1[train_ind,1])
      set.seed(1034+35*k)
      xgtrain <- train(X.train,Y.train, 
                       method = "xgbLinear", trControl = trcon, tuneGrid = data.frame(params))
      ge[k] <- xgtrain$results[,"Accuracy"]
      #ge[j] = unlist(calculate_num(j,sp,X,method = methods))
    }
    final.select <- names(rank.score)[1:(which.max(ge[-1])+1)]
    X.train = as.matrix(X1[train_ind,final.select])
    Y.train = X1[train_ind,1]
    X.test <- as.matrix(X1[-train_ind,final.select])
    Y.test = X1[-train_ind,1]
    
    trmod = train(X.train,Y.train,method="xgbLinear",trControl = trcon,tuneGrid = tunegrid)
    pred = predict(trmod,newdata = X.test,type="prob")
    strc <- roc(Y.test,pred[,2])
    thre = coords(strc, "best", ret = "threshold")
    
    if(is.data.frame(thre)){opt.thre = thre[1,]} else {opt.thre = as.numeric(thre)}
    
    #opt.thre = as.numeric(coords(strc, "best", ret = "threshold",transpose = TRUE))
    pred.class <- ifelse(pred[,2] > opt.thre,1,0)
    dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
    T = confusionMatrix(dt)
    perf = c(T$overall[1],T$byClass[1:2])
    return(list(score = rank.score, performance = perf, final.var = final.select))
  }
  else if (method=="nb")
  {
    trcon = trainControl(method="cv", number=10,search="grid")
    tunegrid <- expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0,0.5,1.0))
    nbmod = train(X.train,Y.train,method="nb",trControl = trcon,tuneGrid = tunegrid)
    nbs = varImp(nbmod)$importance##importance
    params = nbmod$bestTune ## best parameter
    x = nbs[,1]
    names(x) = rownames(nbs)
    rank.score = sort(x,decreasing = TRUE)
    ### optimal number
    ge = 0
    for(k in 2:min(40,length(rank.score)))
    {
      sp =names(rank.score)[1:k]
      X.train <- X1[train_ind,sp]
      Y.train = as.factor(X1[train_ind,1])
      set.seed(1034+35*k)
      xgtrain <- train(X.train,Y.train, 
                       method = "nb", trControl = trcon, tuneGrid = data.frame(params))
      ge[k] <- xgtrain$results[,"Accuracy"]
      #ge[j] = unlist(calculate_num(j,sp,X,method = methods))
    }
    final.select <- names(rank.score)[1:(which.max(ge[-1])+1)]
    X.train = as.matrix(X1[train_ind,final.select])
    Y.train = X1[train_ind,1]
    X.test <- as.matrix(X1[-train_ind,final.select])
    Y.test = X1[-train_ind,1]
    trmod = train(X.train,Y.train,method="nb",trControl = trcon,tuneGrid = tunegrid)
    pred = predict(trmod,newdata = X.test,type="prob")
    rocm <- roc(Y.test,pred[,2])
    thre = coords(rocm, "best", ret = "threshold")
    
    if(is.data.frame(thre)){opt.thre = thre[1,]} else {opt.thre = as.numeric(thre)}
    
    #opt.thre = as.numeric(coords(rocm, "best", ret = "threshold",transpose = TRUE))
    pred.class <- ifelse(pred[,2] > opt.thre,1,0)
    dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
    T = confusionMatrix(dt)
    perf = c(T$overall[1],T$byClass[1:2])
    return(list(score = rank.score, performance = perf, final.var = final.select))
  }
  else if (method == "knn")
  {
    trcon = trainControl(method="cv", number=10,search="grid")
    tunegrid <- expand.grid(k = c(5, 11, 21, 25))
    knnmod = train(X.train, Y.train, method= "knn", trControl=trcon,tuneGrid = tunegrid)
    nbs = varImp(knnmod)$importance
    params = knnmod$bestTune
    x = nbs[,1]
    names(x) = rownames(nbs)
    rank.score = sort(x,decreasing = TRUE)
    ### optimal number
    ge = 0
    for(k in 2:min(40,length(rank.score)))
    {
      sp =names(rank.score)[1:k]
      X.train <- X1[train_ind,sp]
      Y.train = as.factor(X1[train_ind,1])
      set.seed(1034+35*k)
      knntrain <- train(X.train,Y.train, 
                        method = "knn", trControl = trcon, tuneGrid = data.frame(params))
      ge[k] <- knntrain$results[,"Accuracy"]
      #ge[j] = unlist(calculate_num(j,sp,X,method = methods))
    }
    final.select <- names(rank.score)[1:(which.max(ge[-1])+1)]
    X.train = as.matrix(X1[train_ind,final.select])
    Y.train = X1[train_ind,1]
    X.test <- as.matrix(X1[-train_ind,final.select])
    Y.test = X1[-train_ind,1]
    trmod = train(X.train,Y.train,method="knn",trControl = trcon,tuneGrid = tunegrid)
    pred = predict(trmod,newdata = X.test,type="prob")
    rocm <- roc(Y.test,pred[,2])
    thre = coords(rocm, "best", ret = "threshold")
    
    if(is.data.frame(thre)){opt.thre = thre[1,]} else {opt.thre = as.numeric(thre)}
    
    #opt.thre = as.numeric(coords(rocm, "best", ret = "threshold",transpose = TRUE))
    pred.class <- ifelse(pred[,2] > opt.thre,1,0)
    dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
    T = confusionMatrix(dt)
    perf = c(T$overall[1],T$byClass[1:2])
    return(list(score = rank.score, performance = perf, final.var = final.select))
  }
  else if (method == "svm-rfe")
  {
    nrows = length(train_ind)
    nfold = nrows
    folds = rep(1:nfold, len=nrows)[sample(nrows)]
    folds = lapply(1:nfold, function(x) which(folds == x)) ##each training data applied 5-fold svm-rfe
    set.seed(2011)
    results = lapply(folds, svmRFE.wrap, X1[train_ind,], k=1, halve.above=30)
    top.features = WriteFeatures(results, X1[train_ind,], save=FALSE)
    #x = top.features[order(top.features$FeatureID),c("AvgRank")]
    x = top.features$AvgRank
    #names(x) = colnames(X1[,-1])
    names(x) = top.features$FeatureName
    rank.score = sort(x,decreasing = FALSE)
    #featsweep = lapply(1:10, FeatSweep.wrap, results, X1[train_ind,])
    ge = 0
    for(k in 2:min(40,length(rank.score)))
    {
      sp =names(rank.score)[1:k]
      X.train <- X1[train_ind,sp]
      Y.train = as.factor(X1[train_ind,1])
      set.seed(1034+35*k)
      tuned = tune.svm(x=X.train,y = Y.train,gamma = 10^(-2:0),cost = 10^(0:3), tunecontrol=tune.control(cross=10))
      ge[k] <- tuned$best.performance
      #ge[j] = unlist(calculate_num(j,sp,X,method = methods))
    }
    final.select <- names(rank.score)[1:(which.min(ge[-1])+1)]
    X.train = as.matrix(X1[train_ind,final.select])
    Y.train = X1[train_ind,1]
    X.test <- as.matrix(X1[-train_ind,final.select])
    Y.test = X1[-train_ind,1]
    X.use <- X1[train_ind,c("yinfect",final.select)]
    
    model <- tune.svm(yinfect~.,data = X.use,gamma = 10^(-3:0),cost = 10^(0:3))
    model0<- svm(yinfect~.,data = X.use,gamma = as.numeric(model$best.parameters[1]),cost = as.numeric(model$best.parameters[2]),probability = TRUE)
    pr = predict(model0,newdata = X.test,probability=TRUE)
    pred = attr(pr,"probabilities")[,2]
    strc <- roc(Y.test,pred)
    thre = coords(strc, "best", ret = "threshold")
    
    if(is.data.frame(thre)){opt.thre = thre[1,]} else {opt.thre = as.numeric(thre)}
    
    #opt.thre = as.numeric(coords(strc, "best", ret = "threshold",transpose = TRUE))
    pred.class <- ifelse(pred > opt.thre,1,0)
    dt <- table(factor(pred.class,levels= c("1","0")),factor(Y.test,levels= c("1","0")))
    T = confusionMatrix(dt)
    perf = c(T$overall[1],T$byClass[1:2])
    return(list(score = rank.score, performance = perf, final.var = final.select))
    
  }
  
}
