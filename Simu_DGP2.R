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


# generate data
SEED <- 2
n <- 3000
p <- 100

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


# # estimate nuisance parameters
# ### initial
# mu1.hat <- rep(0,N)
# mu0.hat <- rep(0,N)
# ps.hat <- rep(0,N)
# 
# ### eatimate
# for (i in 1:K) {
#   mu.hat<-estimate_mu(Y,D,X,inds.train.set[[i]],inds.eval.set[[i]],method = "rf") 
#   ### Extract regression function
#   mu1.hat[inds.eval.set[[i]]]<-mu.hat[["mu1.hat.eval"]]
#   mu0.hat[inds.eval.set[[i]]]<-mu.hat[["mu0.hat.eval"]]
#   
#   ps.hat[inds.eval.set[[i]]]<-estimate_ps(D,X,inds.train.set[[i]],inds.eval.set[[i]],method = "rf")
# }
# 
# muSignal <- mu1.hat-mu0.hat
# IPSignal <- D * Y / ps.hat  - (1-D) * Y / (1-ps.hat) 
# DRSignal <- mu1.hat-mu0.hat + D * (Y - mu1.hat)/ ps.hat  - (1-D) * (Y - mu0.hat)/ (1-ps.hat) 
# 
# 
# # CATE - project signal on X1
# list_kr_cates_inc = kr_cate(IPSignal,X[,1])
# plot(list_kr_cates_inc,z_label="X1",yrange=c(4,15))
# 
# 
# list_sr_cates_inc = spline_cate(DRSignal,X[,1])
# plot(list_sr_cates_inc,z_label="X1") +labs(title = "a") +geom_hline(aes(yintercept=10), colour="red", linetype="dashed")+geom_abline(intercept = 10,slope=1,color="blue")
# 



# for loop
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 3)))
# region <- function(row, col){
#   viewport(layout.pos.row = row, layout.pos.col = col)
# } 


method_ps = c("DSRR","rf") #j
method_mu = c("DSRR","rf") #k

for (j in 1:2) {
  for (k in 1:2) {
    
    
    # estimate nuisance parameters
    ### initial
    mu1.hat <- rep(0,N)
    mu0.hat <- rep(0,N)
    ps.hat <- rep(0,N)
    
    ### eatimate
    for (i in 1:K) {
      ps.hat[inds.eval.set[[i]]]<-estimate_ps(D,X,inds.train.set[[i]],inds.eval.set[[i]],method = method_ps[j])
      
      mu.hat<-estimate_mu(Y,D,X,inds.train.set[[i]],inds.eval.set[[i]],method = method_mu[k]) 
      ### Extract regression function
      mu1.hat[inds.eval.set[[i]]]<-mu.hat[["mu1.hat"]]
      mu0.hat[inds.eval.set[[i]]]<-mu.hat[["mu0.hat"]]
    }
    
    
    muSignal <- mu1.hat-mu0.hat
    IPSignal <- D * Y / ps.hat  - (1-D) * Y / (1-ps.hat) 
    DRSignal <- mu1.hat-mu0.hat + D * (Y - mu1.hat)/ ps.hat  - (1-D) * (Y - mu0.hat)/ (1-ps.hat) 
    
    mean(DRSignal)
    
    
    # plot 
    kr1  <- spline_cate(muSignal,X[,1])
    plt1 <- plot(kr1,z_label="X1",yrange=c(5,15)) + 
      geom_hline(aes(yintercept=10), colour="red", linetype="dashed")+geom_abline(intercept = 10,slope=1,color="blue", linetype="dashed")+
      labs(title = paste("ps-",as.character(method_ps[j]),"mu-",as.character(method_mu[k],sep = "") )   )
    
    kr2  <- spline_cate(IPSignal,X[,1])
    plt2 <- plot(kr2,z_label="X1")+ 
      geom_hline(aes(yintercept=10), colour="red", linetype="dashed")+geom_abline(intercept = 10,slope=1,color="blue", linetype="dashed")+
      labs(title = paste("ps-",as.character(method_ps[j]),"mu-",as.character(method_mu[k],as.character(p),sep = "") )   )
    
    kr3  <- spline_cate(DRSignal,X[,1])
    plt3 <- plot(kr3,z_label="X1")+ 
      geom_point(data = data.frame(x=seq(-4,4,0.05),y=seq(-4,4,0.05)^2),aes(x=x,y=y),color="blue",size=0.1)+
      geom_hline(aes(yintercept=1), colour="red", linetype="dashed")+
      labs(title = paste("DR,ps-",as.character(method_ps[j]),"mu-",as.character(method_mu[k],",p=",as.character(p),sep = "") )   )
    plt3
    
    print(plt1, vp = region(row = 2*(j-1)+k, col = 1))
    print(plt2, vp = region(row = 2*(j-1)+k, col = 2))
    print(plt3, vp = region(row = 2*(j-1)+k, col = 3))
    

  }
}

# # kernel
# for (j in 1:2) {
#   for (k in 1:2) {
#     
#     
#     # estimate nuisance parameters
#     ### initial
#     mu1.hat <- rep(0,N)
#     mu0.hat <- rep(0,N)
#     ps.hat <- rep(0,N)
#     
#     ### eatimate
#     for (i in 1:K) {
#       ps.hat[inds.eval.set[[i]]]<-estimate_ps(D,X,inds.train.set[[i]],inds.eval.set[[i]],method = method_ps[j])
#       
#       mu.hat<-estimate_mu(Y,D,X,inds.train.set[[i]],inds.eval.set[[i]],method = method_mu[k]) 
#       ### Extract regression function
#       mu1.hat[inds.eval.set[[i]]]<-mu.hat[["mu1.hat.eval"]]
#       mu0.hat[inds.eval.set[[i]]]<-mu.hat[["mu0.hat.eval"]]
#     }
#     
#     
#     muSignal <- mu1.hat-mu0.hat
#     IPSignal <- D * Y / ps.hat  - (1-D) * Y / (1-ps.hat) 
#     DRSignal <- mu1.hat-mu0.hat + D * (Y - mu1.hat)/ ps.hat  - (1-D) * (Y - mu0.hat)/ (1-ps.hat) 
#     
#     # plot 
#     kr1  <- kr_cate(muSignal,X[,1])
#     plt1 <- plot(kr1,z_label="X1",yrange=c(5,15)) + 
#       geom_hline(aes(yintercept=10), colour="red", linetype="dashed")+geom_abline(intercept = 10,slope=1,color="blue", linetype="dashed")+
#       labs(title = paste("ps-",as.character(method_ps[j]),"mu-",as.character(method_mu[k],sep = "") )   )
#     
#     kr2  <- kr_cate(IPSignal,X[,1])
#     plt2 <- plot(kr2,z_label="X1",yrange=c(0,20))+ 
#       geom_hline(aes(yintercept=10), colour="red", linetype="dashed")+geom_abline(intercept = 10,slope=1,color="blue", linetype="dashed")+
#       labs(title = paste("ps-",as.character(method_ps[j]),"mu-",as.character(method_mu[k],as.character(p),sep = "") )   )
#     
#     kr3  <- kr_cate(DRSignal,X[,1])
#     plt3 <- plot(kr3,z_label="X1",yrange=c(5,15))+ 
#       geom_hline(aes(yintercept=10), colour="red", linetype="dashed")+geom_abline(intercept = 10,slope=1,color="blue", linetype="dashed")+
#       labs(title = paste("DR,ps-",as.character(method_ps[j]),"mu-",as.character(method_mu[k],",p=",as.character(p),sep = "") )   )
#     
#     print(plt1, vp = region(row = 2*(j-1)+k, col = 1))
#     print(plt2, vp = region(row = 2*(j-1)+k, col = 2))
#     print(plt3, vp = region(row = 2*(j-1)+k, col = 3))
#     
#   }
# }
