CVlam<-function(D,x,r,K)
{  
  n=nrow(x)
  dx=ncol(x)
  data=cbind(D,x)
  ind=1:n
  id=matrix(0,nrow=ceiling(n/K),ncol=K)
  id[,1]=sample(ind, ceiling(n/K), replace = FALSE, prob = NULL)
  for(i in 2:(K-1)){
  id[,i]=sample(ind[-as.vector(id[,1:(i-1)])], ceiling(n/K), replace = FALSE, prob = NULL)
  }
  
  idK=ind[-as.vector(id[,1:(K-1)])]
  id[1:length(idK),K]=idK

  lammin=1.5*sqrt(r*n*log(dx))
  lammax=3*sqrt(r*n*log(dx))
  nl=10
  lamseq=seq(lammin,lammax,length=nl)
  Ln=rep(0,length=nl)
  for(j in 1:nl){
  lam=lamseq[j]
  Lnlam=0
  for(i in 1:K){
  ki=as.vector(id[,-i])
  ki=sort(ki)
  ki=ki[which(ki!=0)]
  datrain=data[ki,]
  Dtrain=datrain[,1]
  xtrain=datrain[,-1]
  
  result=estimationBini(Dtrain,xtrain,r)
  Bini=result[[1]]
  result=estimation(Dtrain,xtrain,r,lam,Bini)
  Ahat=result[[2]]
  Bhat=result[[1]]
  #Dhat=result[[3]]
  
  ki=as.vector(id[,i])
  ki=sort(ki)
  ki=ki[which(ki!=0)]
  datest=data[ki,]
  Dtest=datest[,1]
  xtest=datest[,-1]
  
  ntest=nrow(xtest)
  dxtest=ncol(xtest)
  ptest=dxtest
  DDtest=Dtest-mean(Dtest)
  DDtest=DDtest%*%t(rep(1,1,length=(dxtest+1)))

  Dnewtest=cbind(1,xtest)*DDtest-rep(1,length=ntest)%*%t(cov(cbind(1,xtest),Dtest))
  Lnlam=Lnlam+sum(diag(t(Dnewtest-xtest%*%Bhat%*%t(Ahat))%*%(Dnewtest-xtest%*%Bhat%*%t(Ahat))))
  }
  Ln[j]=Lnlam
  }
  lamopt=lamseq[max(which(Ln==min(Ln)))]

  return(list(lamopt,Ln))
}