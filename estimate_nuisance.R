



estimate_mu<-function(Y,D,X,inds.train,inds.eval,method="DSRR",...) {
  # This function estimates the conditional mean of Y given D and X
  
  
  if (method=="DSRR") {
    x = X[inds.train,]
    y = Y[inds.train]
    d = D[inds.train]
    
    id1=which(d==1)
    id0=which(d==0)
    
    y1=y[id1]
    y0=y[id0]
    x1=x[id1,]
    x0=x[id0,]
    
    # ==== estimate E(Y|D,X)
    
    K=5
    
    #result=CVr(y1,x1,K)
    #ropt=result[[1]]
    
    ropt=2
    r=2
    
    # result=CVlam(y1,x1,ropt,K)
    # lamopt=result[[1]]
    
    lamopt=1000
    
    result=estimationBini(y1,x1,ropt)
    Biniy1=result[[1]]
    
    result=estimation(y1,x1,ropt,lamopt,Biniy1)
    Bold=result[[1]]
    Bhat=Bold
    ind=which(apply(Bhat^2,1,sum)>0)
    xsels=x1[,ind]
    result=estimationLS(y1,xsels,ropt)
    Bhat=result[[1]]
    Bhat=as.matrix(Bhat)
    
    xsel=x1[,ind]
    uhat=xsel%*%Bhat
    uhat=as.matrix(uhat)
    
    xsel.eval=X[inds.eval,ind]
    uhat.eval=xsel.eval%*%Bhat
    uhat.eval=as.matrix(uhat.eval)
    
    result=choose_h_DIM(y1,uhat,ropt)
    band=result[[1]]
    
    result=nonparaestimate(y1,uhat,uhat.eval,band)
    mu1.hat=result[[1]]
    
    
    
    #result=CVr(y0,x0,K)
    #ropt=result[[1]]
    
    # result=CVlam(y0,x0,ropt,K)
    # lamopt=result[[1]]
    
    lamopt=500
    
    
    result=estimationBini(y0,x0,ropt)
    Biniy2=result[[1]]
    
    result=estimation(y0,x0,ropt,lamopt,Biniy2)
    Bold=result[[1]]
    Bhat=Bold
    ind=which(apply(Bhat^2,1,sum)>0)
    xsels=x0[,ind]
    xsels=as.matrix(xsels)
    result=estimationLS(y0,xsels,ropt)
    Bhat=result[[1]]
    Bhat=as.matrix(Bhat)
    
    xsel=x0[,ind]
    uhat=xsel%*%Bhat
    uhat=as.matrix(uhat)
    
    xsel.eval=X[inds.eval,ind]
    uhat.eval=xsel.eval%*%Bhat
    uhat.eval=as.matrix(uhat.eval)
    
    result=choose_h_DIM(y0,uhat,ropt)
    band=result[[1]]
    
    result=nonparaestimate(y0,uhat,uhat.eval,band)
    mu0.hat=result[[1]]
    
    
    } else if (method == "rf") {
      rf.data <- data.frame(Y=Y,D=D,X=X)[inds.train,]
      rf.reg <- randomForest(Y~.,rf.data)
      
      mu1.hat <- predict(rf.reg, newdata =data.frame(D=rep(1,length(inds.eval)),X=X[inds.eval,]))
      mu0.hat <- predict(rf.reg, newdata =data.frame(D=rep(0,length(inds.eval)),X=X[inds.eval,]))
    }
  
  
  return(list(mu1.hat=mu1.hat,mu0.hat=mu0.hat))
}


estimate_ps<-function(D,X,inds.train,inds.eval,method="DSRR",...){
  ## This function estimates the propensty score--the conditional mean of D given X
  ## output: the estimated valued in inds.eval set              
  
  if (method == "DSRR") {
    
    x = X[inds.train,]
    d = D[inds.train]
    
    # x = X[inds.train.set[[1]],]
    # d = D[inds.train.set[[1]]]
    
    K=5
    
    #result=CVr(d,x,K)
    #ropt=result[[1]]
    ropt=2
    r=2
    
    # result=CVlam(d,x,ropt,K)
    # lamopt=result[[1]]
    
    lamopt=220
    
    result=estimationBini(d,x,ropt)
    Bini=result[[1]] #初始值V(0)
    
    result=estimation(d,x,ropt,lamopt,Bini)
    Bold=result[[1]]
    Bhat=Bold
    
    # 去掉估计为0的x的部分再做一次回归。是否有必要？
    ind=which(apply(Bhat^2,1,sum)>0)
    xsel=x[,ind]
    result=estimationLS(d,xsel,ropt)
    Bhat=result[[1]]
    Bhat=as.matrix(Bhat)
    # 
    
    xsel=x[,ind]
    uhat=xsel%*%Bhat
    uhat=as.matrix(uhat) #用于估计条件期望的样本
    
    xsel.eval=X[inds.eval,ind]
    uhat.eval=xsel.eval%*%Bhat
    uhat.eval=as.matrix(uhat.eval) #估计完后做预测的点，这里和uhat是同样的点
    
    band=choose_h_DIM(d,uhat,ropt)[[1]]
    # band=result[[1]]
    
    ps.hat=nonparaestimate(d,uhat,uhat.eval,band)[[1]]
    # ps.hat=result[[1]]
    
    ps.hat[ps.hat==0] = 0.001
    ps.hat[ps.hat==1] = 1 - 0.001
    
    
    
  } else if (method == "rf") {
    
    rf.data <- data.frame(D=factor(D),X=X)[inds.train,]
    # rf.reg <- randomForest(D~.,rf.data,mtry=p/3,ntree=1000)
    rf.reg <- randomForest(D~.,rf.data)
    ps.hat <- predict(rf.reg, newdata =data.frame(X=X[inds.eval,]),type = "prob")[,2]
  }
  
# !!!!note that the order of ps.hat is corressponding to the order of inds.eval.
  
  return(ps.hat)
  
}
