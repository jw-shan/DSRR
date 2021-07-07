
library(MASS)


expit <- function(logodds){ 1/(1+exp(-logodds))}




# DGP1: linear with r=p1
Datagen.ln <- function(SEED,n,p,p1){
  set.seed(SEED)
  
  beta.true <- c(rep(1,p1),rep(0,p-p1))
  gamma.true <- c(rep(0.5,p1),rep(0,p-p1))
  
  X <- rnorm(n*p)
  X <- matrix(X,n,p)
  
  Y.0 <- rep(0,n)
  Y.1.true <- 10 + X %*% beta.true
  Y.1 <- Y.1.true +rnorm(n)
  Y.1 <- as.numeric(Y.1)
  
  
  U <- runif(n)
  D.true <- expit(X%*%gamma.true)
  D <- as.numeric(D.true>U)
  
  Y <- D*Y.1 + (1-D)*Y.0
  Y.true <- D*Y.1.true 
  
  
  return(list(Y=Y,D=D,X=X,D.true=D.true,Y.true=Y.true,beta.true=beta.true,gamma.true=gamma.true,ATE=10))
}

# DGP2: nonlinear with r=2
Datagen.nonln <- function(SEED,n,p){
  set.seed(SEED)
  
  ######AR(1) covariance matrix for x
  times <- 1:p
  rhoz <- 0.5
  sigmaz <- 1
  H <- abs(outer(times, times, "-"))
  sigX <- sigmaz*rhoz^H
  pp <- nrow(sigX)
  muX<-rep(0,p)
  
  X <- mvrnorm(n, muX, sigX)
  # X <- X-seq(1,1,length=n)%*%t(apply(X,2,mean)) # 每一列的x减去该列x的均值
  
  U <- runif(n)
  D.true <- expit( (X[,1]+X[,2])*(X[,3]+1)/2 ) # E(D|X)
  D <- as.numeric(D.true>U)
  
  # Y.true <- D + X[,1]^2 + X[,2]^2 #E(Y|X)
  
  Y.true <- (1+D) * X[,1]^2 + X[,2]^2 #E(Y|X)
  Y <- Y.true + + rnorm(n)
  
  return(list(Y=Y,D=D,X=X,D.true=D.true,Y.true=Y.true,ATE=1))
}


# data <- Datagen.asp(1,100,3)
