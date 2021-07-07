nonparaestimate<-function(D,uhat,uhatf,band)
{ 
    muhat=npreg(tydat=D, txdat=uhat, exdat=uhatf,regtype="lc", bws=band)
    muhat=fitted(muhat)

  return(list(muhat))
}
