args = commandArgs(trailingOnly=TRUE)  # passed from script
j <- as.numeric(args[1])

library(pROC)
library(caret)
require(glmnet)
require(tidyverse)
library(xgboost)
library(randomForest)
source("svm-ref.R")
RNGkind(sample.kind = "Rejection")
source("F01-ML-Boruta.R")

res = runburota(j,X1=XXX.binary.f)
saveRDS(res,file=paste0("uninf_split",j,".RDS"))