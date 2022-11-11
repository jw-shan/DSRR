# source
rm(list = ls())
library(hdm)
library(randomForest)
library(ggplot2)
library(grid)
library(glmnet)
library(Matrix)
library(MASS)
library(ncvreg)
library(splines)
library(glmnet)
library(crs)
library(np)
# library(causalDML)
source("estimation.R")
source("estimationLS.R")
source("estimate_nuisance.R")
source("DGP.R")
source("choose_bandwidth.R")
source("initialvalue.R")
source("nonparametricestimation.R")
source("CVlam.R")
source("CVr.R")


library(parallel)



# generate data
# SEED <- 2
n <- 3000
p <- 50
nsim <- 100
truevalue <- 1

cl <- makeCluster(50)
clusterExport(cl,ls())

print("Start!")

est_fun <- function(count){
  
  library(hdm)
  library(randomForest)
  library(ggplot2)
  library(grid)
  library(glmnet)
  library(Matrix)
  library(MASS)
  library(ncvreg)
  library(splines)
  library(glmnet)
  library(crs)
  library(np)
  
  SEED = count
  set.seed(SEED)
  
  data1 <- Datagen.nonln(SEED,n,p)
  # data1 <- Datagen.ln(SEED,n,p,3)
  
  
  Y  = data1$Y
  D  = data1$D
  X  = data1$X
  
  # data splitting
  N=length(D) 
  K=5
  inds.train.set = list()
  inds.eval.set = list()
  inds = 1:N
  for (i in 1:K) {
    inds.eval.set[[i]] <- sample(inds,floor(length(inds)/(K-i+1)))
    inds.train.set[[i]] <- setdiff(1:N,inds.eval.set[[i]])
    inds = setdiff(inds,inds.eval.set[[i]])
  }
  
  # estimate nuisance parameters
  ### initial
  mu1.hat <- rep(0,N)
  mu0.hat <- rep(0,N)
  ps.hat <- rep(0,N)
  
  ### eatimate
  for (i in 1:K) {
    
    ps.hat[inds.eval.set[[i]]]<-estimate_ps(D,X,inds.train.set[[i]],inds.eval.set[[i]])
    
    mu.hat<-estimate_mu(Y,D,X,inds.train.set[[i]],inds.eval.set[[i]]) 
    ### Extract regression function
    mu1.hat[inds.eval.set[[i]]]<-mu.hat[["mu1.hat"]]
    mu0.hat[inds.eval.set[[i]]]<-mu.hat[["mu0.hat"]]
  }
  
  DRSignal <- mu1.hat-mu0.hat + D * (Y - mu1.hat)/ ps.hat  - (1-D) * (Y - mu0.hat)/ (1-ps.hat) 
  
  return(DRSignal)
  
}


est  <- parSapply(cl,1:nsim,est_fun)
save.image("simu_para3.RData")

result <- matrix(nrow = 1, ncol = 4)
colnames(result)<-c("bias","stdev","RMSE","CR")
rownames(result)<-"naive"

Delta <- mean(est)
bias  <- Delta - truevalue
mse   <- 1/nsim*(sum((est-truevalue)^2))
stdev <- sqrt(1/nsim*(sum((est-Delta)^2)))
rmse<- sqrt(mse)

count<-0
for(j in 1:nsim){
  if(est[j]> truevalue-1.96*stdev & est[j]< truevalue+1.96*stdev){
    count<- count+1
    }
}
coverage_rate <- count/nsim
CR <- coverage_rate

result <- cbind(bias,stdev,rmse,CR)

result
save.image("simu_para3.RData")



