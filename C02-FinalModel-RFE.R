args = commandArgs(trailingOnly=TRUE)  # passed from script
k <- as.numeric(args[1])
library(caret)
###Take Uninfected as example

load("datafile.RData")
source("F02-FinalModel-RFE.R")
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

### optimal number

opt = optim.num(seed0,rank.ave.score,X1=XXX.binary.f,methods = md)
save(opt,rank.ave.score,file = paste0(md,"-uninf-finalmodel.RData"))