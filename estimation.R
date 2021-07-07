
estimation<-function(D,x,r,lam,Bini)
{ 
  n=nrow(x)
  dx=ncol(x)
  p=dx
  DD=D-mean(D)
  DD=DD%*%t(rep(1,1,length=(dx+1)))

  Dnew=cbind(1,x)*DD-rep(1,length=n)%*%t(cov(cbind(1,x),D)) #\tilde{W}

  #Bold=rbind(diag(r),matrix(0,nrow=(p-r),ncol=r))
  #Sig=ginv(t(x)%*%x) 
  #Lam=t(x)%*%(cbind(1,x)*DD)
  #betahat=Sig%*%Lam
  #UVD=svd(betahat)
  #Bold=(UVD$u)%*%diag(UVD$d)
  #Bold=Bold[,1:r]
  #Bold=as.matrix(Bold)
  #xB=x%*%Bold
 
  Bold=Bini
  xB=x%*%Bold
  
  DV=t(Dnew)%*%xB
  DV=svd(DV)
  #pp=sqrt(DV$d)
  #pp=sqrt(n)
  pp=1
  Aold=(DV$u)%*%t(DV$v)/pp

  delta=2
  ep=10^(-2)
  step=0

  while((delta>ep)&&(step<=500)){
   step=step+1
     LB=0.5*sum(diag(t(Dnew-x%*%Bold%*%t(Aold))%*%(Dnew-x%*%Bold%*%t(Aold))))
     DLB=-t(x)%*%(Dnew-x%*%Bold%*%t(Aold))%*%Aold
     #LB=0.5*sum(t(Dnew%*%Aold-x%*%Bold)%*%(Dnew%*%Aold-x%*%Bold))
     #DLB=-t(x)%*%(Dnew%*%Aold-x%*%Bold)
     HBold=LB+lam*sum(sqrt(apply(Bold^2,1,sum)))

     eta=0.5
     tt=1/eta
     del=2
     while(del>0){
      tt=eta*tt
      Bup=Bold-tt*DLB
      Bupnorm=apply(Bup^2,1,sum)     
      Bupnorm=sqrt(Bupnorm)
      Bnorm=Bupnorm%*%t(rep(1,length=r))  
      lamt=lam*tt 
      Bnew=(Bup-lamt*Bup/Bnorm)*((Bupnorm>=lamt)%*%t(rep(1,length=r)))
      LBnew=0.5*sum(diag(t(Dnew-x%*%Bnew%*%t(Aold))%*%(Dnew-x%*%Bnew%*%t(Aold))))
      #LBnew=0.5*sum(t(Dnew%*%Aold-x%*%Bnew)%*%(Dnew%*%Aold-x%*%Bnew))     
      HBnew=LBnew+lam*sum(sqrt(apply(Bnew^2,1,sum)))
      del=HBnew-HBold+0.5*10^(-2)*sum(diag(t(Bnew-Bold)%*%(Bnew-Bold)))/tt
    }
      #Bnew=svd(Bnew)$u
      xB=x%*%Bnew
      DV=t(Dnew)%*%xB
      DV=svd(DV)
      Anew=(DV$u)%*%t(DV$v)/pp 
     LBnew=0.5*sum(diag(t(Dnew-x%*%Bnew%*%t(Anew))%*%(Dnew-x%*%Bnew%*%t(Anew))))
      #LBnew=0.5*sum(t(Dnew%*%Anew-x%*%Bnew)%*%(Dnew%*%Anew-x%*%Bnew))   
      HBnew=LBnew+lam*sum(sqrt(apply(Bnew^2,1,sum)))
     #delta=abs(HBnew-HBold)
    delta=sqrt(sum((Bold-Bnew)^2)/p)
     Aold=Anew
     Bold=Bnew

  }

#PBold=Bold%*%solve(t(Bold)%*%Bold)%*%t(Bold)
#B=rbind(cbind(beta1,beta2),matrix(0,nrow=(p-5),ncol=r))
#PB=B%*%solve(t(B)%*%B)%*%t(B)
#sum(diag(PB%*%PBold))/r

    return(list(Bold,Aold))
  }
      