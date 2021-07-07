estimationLS<-function(D,xsel,r)
{ 

  n=nrow(xsel)
  dx=ncol(xsel)
  p=dx
  DD=D-mean(D)
  DD=DD%*%t(rep(1,1,length=(dx+1)))

  Sig=solve(t(xsel)%*%xsel) 
  Lam=t(xsel)%*%(cbind(1,xsel)*DD)
  betahat=Sig%*%Lam
  UVD=svd(betahat)
  Bhat=UVD$u[,1:r]
     return(list(Bhat))
  }
      