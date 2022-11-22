args = commandArgs(trailingOnly=TRUE)  # passed from script
k <- as.numeric(args[1])
library(caret)

load("datafile.RData")



###run final model for each algorithm
var.imp <- function(seed0,X1,methods)
{
  set.seed(252+7*seed0)
  train_ind = createDataPartition(X1$yinfect,p=2/3,list = FALSE)[,1]
  X.train = data.frame(X1[train_ind,-1])
  Y.train = as.factor(X1[train_ind,1])
  trcon = trainControl(method="cv", number=10,search="grid")
  tunegrid2 = expand.grid(alpha = seq(0.01,0.7,length.out = 20), lambda = seq(0.01,1,length.out = 25))
  tunegrid3 <- expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0.5,1))
  tunegrid5 =   expand.grid(nrounds =500,eta =  0.3, lambda = c(0.5,1),
                            alpha = c(0,0.5,1))
  tunegrid1 <- expand.grid(.mtry=c(1:5,7,9)) ##rf
  tunegrid4 <- expand.grid(k = c(5, 11, 21, 25))
  if(methods == "Elastic")
  {
    set.seed(252+7*seed0)
    dw_elnet = train(X.train,Y.train,method="glmnet",trControl = trcon,tuneGrid = tunegrid2,family = "binomial")
    imp.mat <- varImp(dw_elnet,lambda = dw_elnet$bestTune$lambda,alpha = dw_elnet$bestTune$alpha)$importance
    x = imp.mat[,1]
    names(x) = rownames(imp.mat)
  }
  else if(methods == "rf")
  {
    set.seed(252+7*seed0)
    rfmod = train(X.train,Y.train,method="rf",trControl = trcon,tuneGrid = tunegrid1,ntree = 1000)
    rfs =  varImp(rfmod)$importance##importance
    x = rfs[,1]
    names(x) = rownames(rfs)
  }
  else if(methods == "nb")
  {
    set.seed(252+7*seed0)
    rfmod = train(X.train,Y.train,method="nb",trControl = trcon,tuneGrid = tunegrid3)
    rfs =  varImp(rfmod)$importance##importance
    x = rfs[,1]
    names(x) = rownames(rfs)
  }
  else if(methods == "xgboost")
  {
    set.seed(252+7*seed0)
    rfmod = train(X.train,Y.train,method="xgbLinear",trControl = trcon,tuneGrid = tunegrid5)
    rfs =  varImp(rfmod)$importance##importance
    x = rfs[,1]
    names(x) = rownames(rfs)
  }
  else if(methods == "knn")
  {
    set.seed(252+7*seed0)
    rfmod = train(X.train,Y.train,method="knn",trControl = trcon,tuneGrid = tunegrid4)
    rfs =  varImp(rfmod)$importance##importance
    x = rfs[,1]
    names(x) = rownames(rfs)
  }
  return(x)
}

calculate.num <- function(seed0,spp,X1,method)
{
  trcon = trainControl(method="cv", number=10,search="grid")
  tunegrid2 = expand.grid(alpha = seq(0.01,0.7,length.out = 20), lambda = seq(0.01,1,length.out = 25))
  tunegrid3 <- expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0.5,1))
  tunegrid5 =   expand.grid(nrounds =500,eta =  0.3, lambda = c(0.5,1),
                            alpha = c(0,0.5,1))
  tunegrid1 <- expand.grid(.mtry=c(1:5,7,9)) ##rf
  tunegrid4 <- expand.grid(k = c(5, 11, 21, 25))
  performance = numeric()
  for(j in 1:length(seed0))
  {
    set.seed(252+7*seed0[j])
    train_ind = createDataPartition(X1[,"yinfect"],p=2/3,list = FALSE)[,1]
    X.train <- X1[train_ind,spp]
    Y.train = as.factor(X1[train_ind,1])
    if(method == "Elastic")
    {
      set.seed(252+7*seed0[j])
      # X.use <- X.use %>% mutate(yinfect = factor(yinfect, 
      #                                            labels = make.names(levels(yinfect))))
      dw_elnet <- train(X.train,Y.train, 
                        method = "glmnet", trControl = trcon, family = "binomial", 
                        tuneGrid = tunegrid2)
      performance[j] <-  max(dw_elnet$results[,"Accuracy"],na.rm=TRUE)
    }
    else if(method == "rf")
    {
      set.seed(252+7*seed0[j])
      rftrain <- train(X.train,Y.train, 
                       method = "rf", trControl = trcon, tuneGrid = tunegrid1,ntree = 1000)
      performance[j] <- max(rftrain$results[,"Accuracy"],na.rm=TRUE)
    }
    else if(method == "xgboost")
    {
      set.seed(252+7*seed0[j])
      rftrain <- train(X.train,Y.train, 
                       method = "xgbLinear", trControl = trcon, tuneGrid = tunegrid5)
      performance[j] <- max(rftrain$results[,"Accuracy"],na.rm=TRUE)
    }
    else if(method == "knn")
    {
      set.seed(252+7*seed0[j])
      rftrain <- train(X.train,Y.train, 
                       method = "knn", trControl = trcon, tuneGrid = tunegrid4)
      performance[j] <- max(rftrain$results[,"Accuracy"],na.rm=TRUE)
    }
    else if(method == "nb")
    {
      set.seed(252+7*seed0[j])
      rftrain <- train(X.train,Y.train, 
                       method = "nb", trControl = trcon, tuneGrid = tunegrid3)
      performance[j] <- max(rftrain$results[,"Accuracy"],na.rm=TRUE)
    }
  }
  
  
  
  return(mean(performance,na.rm=TRUE))
}

optim.num <- function(seed0,rank.ave.score,X1,methods)
{
  ge = 0
  for(k in 2:min(35,length(rank.ave.score)))
  {
    sp =names(rank.ave.score)[1:k]
    ge[k] = calculate.num(seed0,sp,X1,method = methods)
  }
  final.select <- names(rank.ave.score)[1:(which.max(ge[-1])+1)]
  return(list(select = final.select,gaccu = ge))
}

rank.score = NULL

if(k==1)
{
  md ="rf"
}else if(k==2)
{
  md = "Elastic"
}else if(k == 3)
{
  md = "knn"
}else if(k==4)
{
  md = "nb"
}else if(k==5)
{   
  md = "xgboost"
}

seed0 = c()
for(j in 1:1005)
{
  #assign(paste0("res_",j),readRDS(paste0("asc_split",j,".RData")))
  if(file.exists(paste0("uninf_split",j,".RDS")))
  {
    res0 = readRDS(paste0("uninf_split",j,".RDS"))
    seed0 = c(seed0,res0$seed)
    print(res0$seed)
    res = var.imp(res0$seed,X1=XXX.binary.f,methods = md)
    rank.score = rbind(rank.score,res)
  }
  
}
dim(rank.score)
head(rank.score)
rank.ave.score =  sort(colMeans(rank.score,na.rm=TRUE),decreasing = TRUE)


#rank.ave.score.rf =  sort(colMeans(res.rf,na.rm=TRUE),decreasing = TRUE)
#rank.ave.score.xgb =  sort(colMeans(res.xgb,na.rm=TRUE),decreasing = TRUE)

### optimal number

opt = optim.num(seed0,rank.ave.score,X1=XXX.binary.f,methods = md)


save(opt,rank.ave.score,file = paste0(md,"-uninf-finalmodel.RData"))