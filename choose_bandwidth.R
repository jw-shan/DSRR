#根据不同维度利用cv选取bandwidth
choose_h_DIM<-function(D,uhat,r)
{ 
  
  #result=crs(D~uhat[,1]+uhat[,2],cv="none",basis="tensor",degree=c(3,3),segments=c(1,1))
  
  if(r==1){
  result=npregbw(D~uhat[,1],regtype="lc") 
  }else if(r==2){
  result=npregbw(D~uhat[,1]+uhat[,2],regtype="lc")
  }else if(r==3){
  result=npregbw(D~uhat[,1]+uhat[,2]+uhat[,3],regtype="lc")
  }else if(r==4){
  result=npregbw(D~uhat[,1]+uhat[,2]+uhat[,3]+uhat[,4],regtype="lc")
  }else if(r==5){
  result=npregbw(D~uhat[,1]+uhat[,2]+uhat[,3]+uhat[,4]+uhat[,5],regtype="lc")
  }
  band=result[[1]]
  # muhat=fitted(npreg(exdat=uhatf, bws=result)) 
  return(list(band))
}

