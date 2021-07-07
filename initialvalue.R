# 返回V的初值V(0)
estimationBini<-function(D,x,r)
{ 
  n=nrow(x)
  dx=ncol(x)
  p=dx
  DD=D-mean(D)
  DD=DD%*%t(rep(1,1,length=(dx+1)))

  Dnew=cbind(1,x)*DD-rep(1,length=n)%*%t(cov(cbind(1,x),D)) #貌似是\tilde{W}

  #Bold=rbind(diag(r),matrix(0,nrow=(p-r),ncol=r))
  #Sig=ginv(t(x)%*%x) 
  #Lam=t(x)%*%(cbind(1,x)*DD)
  #betahat=Sig%*%Lam
  #UVD=svd(betahat)
  #Bold=(UVD$u)%*%diag(UVD$d)
  #Bold=Bold[,1:r]
  #Bold=as.matrix(Bold)
  #xB=x%*%Bold

  id=matrix(0,nrow=(p+1),ncol=p)

  jj=1
  fit=cv.glmnet (x, Dnew[,jj], family="gaussian")
  coef=coef(fit,s="lambda.min")
  coef=coef[2:length(coef)]
  coef=as.numeric(coef!=0) 
  lam=fit$lambda.min
  id[jj,]=coef

  for(jj in 2:(p+1)){
  fit=glmnet (x, Dnew[,jj], family="gaussian",lambda=lam)
  coef=predict(fit,type="coef")
  coef=coef[2:length(coef)]
  coef=as.numeric(coef!=0) 
  id[jj,]=coef
  }
  idd=apply(id,2,sum)
  idd=cbind(seq(1,p,length=p),idd)
  idd=idd[order(idd[,2],decreasing = TRUE),]
  idd=idd[which(idd[,2]!=0),]
  id1=idd[1:min(floor(sqrt(n)),nrow(idd)),1]

  
  xid1=x[,id1]
  Sig=ginv(t(xid1)%*%xid1) 
  Lam=t(xid1)%*%Dnew
  betahat=Sig%*%Lam
  UVD=svd(betahat)
  Bold1=(UVD$u)%*%diag(UVD$d)
  Bold1=Bold1[,1:r]
  Bold1=as.matrix(Bold1)
  Bold=matrix(0,nrow=p,ncol=r)
  Bold[id1,]=Bold1
  Bini=Bold
    return(list(Bini))
  }
      
