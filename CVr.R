
CVr<-function(D,x,K)
{  
  
  n=nrow(x)
  dx=ncol(x)
  data=cbind(D,x)
  ind=1:n
  ng=ceiling(n/K)
  id=matrix(0,nrow=ng,ncol=K)
  id[,1]=sample(ind, ng, replace = FALSE, prob = NULL)
  for(i in 2:(K-1)){
  id[,i]=sample(ind[-as.vector(id[,1:(i-1)])], ng, replace = FALSE, prob = NULL)
  }
  
  idK=ind[-as.vector(id[,1:(K-1)])]
  id[1:length(idK),K]=idK
 
  rseq=1:5
  Ln=rep(0,length=5)  
  for(j in 1:5){
  r=rseq[j]
  Lnlam=0
  
  result=CVlam(D,x,r,K) 
  lamopt=result[[1]]

  for(i in 1:K){
  ki=as.vector(id[,-i])
  ki=sort(ki)
  ki=ki[which(ki!=0)]
  datrain=data[ki,]
  Dtrain=datrain[,1]
  xtrain=datrain[,-1]
  
  result=estimationBini(Dtrain,xtrain,r)
  Bini=result[[1]]
  result=estimation(Dtrain,xtrain,r,lamopt,Bini)
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
  ropt=rseq[min(which(Ln==min(Ln)))]

  return(list(ropt,Ln))
}  
