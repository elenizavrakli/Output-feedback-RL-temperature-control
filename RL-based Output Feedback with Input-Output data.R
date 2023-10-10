library(Matrix) 
library(e1071)
library(pracma)
library(devtools)
library(expm)
library(MASS)
library(tidyverse)
library(progress)
library(Rcpp)

#Noise function for data generation
noise_func <- function(m,t){
  noisevar=1.5
  noiseterm = sqrt(noisevar)*rnorm(1)+100*sin(rep(1.65,m)*t)+90*sin(runif(m,0,1.65)*t)+80*sin(runif(m,0,1.65)*t)+70*sin(runif(m,0,1.65)*t)+60*sin(runif(m,0,1.65)*t)+50*sin(runif(m,0,1.65)*t)+40*sin(runif(m,0,1.65)*t)+30*sin(runif(m,0,1.65)*t)+20*sin(runif(m,0,1.65)*t)+10*sin(runif(m,0,1.65)*t)
  return(noiseterm)
}

set.seed(1)

#System dimensions
n=6 #state length
m=7 #input length
p=5 #measurement length
N=6 #time horizon
l= N*m+N*p+p+m #length of augmented state z

#Weighting matrices for performance index
R <- diag(m) #input 
Q <- diag(p) #tracking
#Weighting matrices after Bayesian Optimisation
#R<-diag(c(0.14503672, 0.38862483, 0.02004028, 0.11571914, 0.09935187, 0.31611002, 0.0100000))
#Q<-diag(c(0.1740472, 0.05596686, 0.0100000, 0.0100000, 0.16019648))
C2 <- as.matrix(cbind(diag(p),-diag(p)))
Q2 <- t(C2)%*%Q%*%C2 

#Simulation parameter initialisation
n_dat <- 15000
n_iter <- 1000
pb_iter <- progress_bar$new(total =n_iter) #progress bar for training
gam <- 0.99 #discount factor
mu <- 0.0001 #Regularisation term

#Model matrices 
A <- matrix(c(0.992,0.0018,0,0,0,0,
              0.0023,0.9919,0.0043,0,0,0,
              0,-0.0042,1.0009,0.0024,0,0,
              0,0,0.0013,0.9979,0,0,
              0,0,0,0,0.9972,0,
              0,0,0,0,0,0.9953),
            nrow = n, byrow = TRUE) #state matrix

B <- matrix(c(1.0033,0,0,0,0,0,-0.2175,
              0,1.0460,0,0,0,0,-0.0788,
              0,0,1.0326,0,0,0,-0.0020,
              0,0,0,0.4798,0,0,-0.0669,
              0,0,0,0,0.8882,0,0.1273,
              0,0,0,0,0,1.1699,-0.1792),
            nrow = n, ncol = m, byrow = TRUE) #input matrix 

C <- matrix(c(0.992,0.00018,0,0,-0.0001,0,
              0.0023,1.3,0.0043,0,0,0,
              0,-0.0042,1.0109,0.0024,0,0.201,
              0,0,0.0013,0.989,0.00031,0.64,
              0,0,0,0,0.923,0.3),
            nrow = p, byrow = TRUE) #measurement matrix 

#augmented system matrices 
Ft <- diag(p)
Tmatr <-  as.matrix(bdiag(A,Ft))
B1 <- as.matrix(rbind(B,matrix(rep(0,n*m),n,m)))

#cost weighting matrices
C2 <- as.matrix(cbind(diag(p),-diag(p)))
Q2 <- t(C2)%*%Q%*%C2 

#Stabilising gain for the data generation. It gets applied to xaug
Kgain_data_sim <-  matrix(c(  0.7395,   -0.0076,   -0.0003,   -0.0264,    0.0194,   -0.0170,  
                              -0.0076,    0.7430 ,   0.0031,   -0.0093,    0.0068,   -0.0060, 
                              -0.0003,   -0.0033,    0.7599,    0.0021,    0.0002,   -0.0002,  
                              -0.0126,   -0.0042,    0.0016,    1.0971,    0.0092,   -0.0079, 
                              0.0171,    0.0058,    0.0002,    0.0170,    0.8179,    0.0108,  
                              -0.0198,   -0.0067,   -0.0002,   -0.0193,    0.0143,    0.6823, 
                              -0.1525,   -0.0519,   -0.0018,   -0.1412,    0.1091,   -0.0977),nrow=m,byrow=TRUE)

#Initialise all vectors used
x <- matrix(0,n,n_dat) #state
y <- matrix(0,p,n_dat) #measurement  
u <- matrix(0,m,n_dat) #input
r <- matrix(0,p,n_dat) #reference trajectory
yaug <- matrix(0,2*p,n_dat) #[measurement,reference]
z <- matrix(0,l,n_dat) #augmented state [\ubar_{k-1,k-N},\ybar_{k-1,k-N},r_{k_N},u_k]
ubar <- matrix(0,N*m,n_dat)
ybar <- matrix(0,N*p,n_dat)
cost <- matrix(0,n_dat,1)
lhsvi <- matrix(0,n_dat,l^2) #kronecker product term
Holdvi <- diag(l) #value function kernel matrix (what is being estimated in the training)
uest <- array(0,c(m,n_iter,n_dat)) #int estimate for training
zi <- matrix(0,l,n_dat) #augmented state [\ubar_{k-1,k-N},\ybar_{k-1,k-N},r_{k_N},uest_k]
costvi <- matrix(0,n_dat,1) #cost estimate inside the loop
errnorm <- matrix(0,n_iter,n_dat) #error norm for exiting the loop

#Initial values
x[,1]<- rep(50,n)
y[,1] <- rep(50,p)
#r[,1] <- c(190,180,170,165,160)
r[,1] <- c(150,160,170,175,180)
u[,1] <- -Kgain_data_sim%*%x[,1] + noise_func(m,1)
yaug[,1] <- c(y[,1],r[,1])
z[,1] <-  rnorm(l)
cost[1] <- t(yaug[,1])%*%Q2%*%yaug[,1]+t(u[,1])%*%R%*%u[,1]

for(k in 1:(n_dat-1)){
  x[,k+1] = A%*%x[,k] + B%*%u[,k]
  y[,k+1] = C%*%x[,k+1]
  r[,k+1] = Ft%*%r[,k]
  yaug[,k+1] = c(y[,k+1],r[,k+1])
  
  u[,k+1] = -Kgain_data_sim%*%x[,k+1] + noise_func(m,k+1)
  
  
  if(k<N){
    z[,k+1] <- rnorm(l)
  }else{
    for(bar in 0:(N-1)){
      ubar[(bar*m+1):((bar+1)*m),k+1] <- u[,k-bar]
      ybar[(bar*p+1):((bar+1)*p),k+1] <- y[,k-bar]
    }
    z[,k+1] <- c(ubar[,k+1],ybar[,k+1],r[,k-N+1],u[,k+1])
  }
  cost[k+1]=t(yaug[,k+1])%*%Q2%*%yaug[,k+1]+t(u[,k+1])%*%R%*%u[,k+1]
}

#Normalise data that will be used in training
z = z/1000
yaug = yaug/1000
ubar=ubar/1000
ybar=ybar/1000
r=r/1000
u=u/1000

#Kronecker product vectors
for(k in (N+1):n_dat){
  lhsvi[k,] = kronecker(t(z[,k]),t(z[,k]))
}

sourceCpp("summcrossprod.cpp")
summat=sumcrossprod(lhsvi,N,n_dat+l-1) #indices correct?????

leftmat = summat + mu*diag(l^2)
chol_decomp = chol(leftmat)
inverse_leftmat = chol2inv(chol_decomp)

for(i in 1:n_iter){
  pb_iter$tick()
  p0=Holdvi[(l-m+1):l,(l-m+1):l,drop=FALSE]
  pu=Holdvi[(l-m+1):l,1:(N*m),drop=FALSE]
  py=Holdvi[(l-m+1):l,(N*m+1):(N*m+N*p),drop=FALSE]
  pr=Holdvi[(l-m+1):l,(N*m+N*p+1):(N*m+N*p+p),drop=FALSE]
  
  for(k in N:(n_dat-1)){
    uest[,i,k+1] = -ginv(p0)%*%(cbind(pu,py,pr)%*%z[1:(l-m),k+1])
    zi[,k+1] = c(ubar[,k+1],ybar[,k+1],r[,k-N+1],uest[,i,k+1])
    costvi[k] = t(yaug[,k])%*%Q2%*%yaug[,k]+t(u[,k])%*%R%*%u[,k]+gam*t(zi[,k+1])%*%Holdvi%*%zi[,k+1]
  }
  
  sumvec = matrix(0,l^2,1)
  for(j in (N+1):(n_dat+l)){
    sumvec = sumvec + lhsvi[j,]*costvi[j]
  }
  
  rightvec = sumvec
  
  HLS_est = inverse_leftmat %*% rightvec  
  
  errnorm[i] = norm(matrix(HLS_est,nrow=l)-Holdvi)
  if(errnorm[i]<=0.001){
    break
  }
  
  Holdvi = matrix(HLS_est,nrow=l)
  p0=Holdvi[(l-m+1):l,(l-m+1):l,drop=FALSE]
  pu=Holdvi[(l-m+1):l,1:(N*m),drop=FALSE]
  py=Holdvi[(l-m+1):l,(N*m+1):(N*m+N*p),drop=FALSE]
  pr=Holdvi[(l-m+1):l,(N*m+N*p+1):(N*m+N*p+p),drop=FALSE]
  kgainfinal = ginv(p0)%*%cbind(pu,py,pr)
}
