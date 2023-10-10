#We perform Bayesian Optimisation on the weighting parameters Q and R for Observer-based output feedback
#The same code can be used for the data-driven output feedback case by replacing the code in the fitness function

library(rBayesianOptimization)
library(Matrix) 
library(tidyverse)
library(netcontrol)

set.seed(1)

fitness <- function(q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,r7){
  
  R <- diag(c(r1,r2,r3,r4,r5,r6,r7)) #input 
  Q <- diag(c(q1,q2,q3,q4,q5)) #tracking
  
  #System dimensions
  n=6 #state length
  m=7 #input length
  p=5 #measurement length
  
  layer_number =100
  #length of simulation
  
  tau =0.002 #small positive number for observer pole placement
  gamma = 0.99 #discount factor for perfomance index
  
  
  # Define system matrices 
  # State matrix
  A <- matrix(c(0.992,0.0018,0,0,0,0,
                0.0023,0.9919,0.0043,0,0,0,
                0,-0.0042,1.0009,0.0024,0,0,
                0,0,0.0013,0.9979,0,0,
                0,0,0,0,0.9972,0,
                0,0,0,0,0,0.9953),
              nrow = n, byrow = TRUE)
  # Input matrix
  B <- matrix(c(1.0033,0,0,0,0,0,-0.2175,
                0,1.0460,0,0,0,0,-0.0788,
                0,0,1.0326,0,0,0,-0.0020,
                0,0,0,0.4798,0,0,-0.0669,
                0,0,0,0,0.8882,0,0.1273,
                0,0,0,0,0,1.1699,-0.1792),
              nrow = n, ncol = m, byrow = TRUE)
  # Output matrix
  C <- matrix(c(0.992,0.00018,0,0,-0.0001,0,
                0.0023,1.3,0.0043,0,0,0,
                0,-0.0042,1.0109,0.0024,0,0.201,
                0,0,0.0013,0.989,0.00031,0.64,
                0,0,0,0,0.923,0.3),
              nrow = p, byrow = TRUE) #measurement matrix 
  
  #State matrix of augmented system
  Tmatr<- as.matrix(bdiag(A,diag(p))) #the reference signal is the same for all time points so F=I
  
  #Input matrix of augmented system
  B1<-as.matrix(rbind(B,matrix(rep(0,p*m),p,m)))
  
  #Weighting matrices for performance index
  
  C1<-as.matrix(cbind(diag(p),-diag(p)))
  
  Q1 <- cbind(rbind(t(C)%*%Q%*%C ,-Q%*%C),rbind(-t(C)%*%Q,Q))
  #Q1 <- t(C1)%*%Q%*%C1 
  
  
  # Initialise control policy randomly
  i=0
  Knew <- matrix(rnorm(m*(n+p),0,10),nrow=m)
  K <- matrix(rep(3,m*(n+p)),nrow=m)
  
  #Obtain optimal control policy through repeated solutions of the Lyapunov equation
  while(norm(K-Knew)>0.1 ){
    K = Knew
    P <- dlyap(gamma*(Tmatr-B1%*%K),Q1+t(K)%*%R%*%K)
    Knew <- gamma*ginv(R+gamma*t(B1)%*%P%*%B1)%*%t(B1)%*%P%*%Tmatr 
    i=i+1
  }
  
  #Initialize the state 
  x <- rep(20,n)
  xobs <-rep(50,n)
  xvec <- C%*%x #vector storing past states, used to plot the trajectories
  xvecobs <- C%*%x #vector storing past observer estimates 
  
  #Define optimal state to be tracked
  xopt <- c(150,160,170,175,180)
  
  #Define augmented state
  X<-c(xobs,xopt) 
  cost<- matrix(0,nrow=layer_number,ncol=1)
  u <- - Knew %*%X
  cost[1] = t(X)%*%Q1%*%X+t(u)%*%R%*%u
  for(i in 2:layer_number){
    #Calculate the optimal control
    u <- - Knew %*%X
    cost[i] = cost[i-1]+ gamma^(i-1)*(t(X)%*%Q1%*%X+t(u)%*%R%*%u) 
    #Simulation from the state space system, applying the optimal control 
    xnew<-A%*%x+B%*%u
    y=C%*%x
    
    #Equations for observer gain
    Psi <- C%*%t(C)+diag(tau,p)
    L <- A%*%t(C)%*%ginv(Psi)
    
    #State estimation from observer
    xobsnew <- A%*%xobs+ B%*%u + L%*%(y-C%*%xobs)
    
    #Update internal state and observer state
    x <- xnew
    xobs <- xobsnew
    X<-c(xobs,xopt) 
    
    #Save state and state estimate to corresponding vectors
    xvecobs<-cbind(xvecobs,C%*%xobsnew)
    xvec<-cbind(xvec,C%*%xnew)
  }
  QI=diag(p)
  QI2 <- cbind(rbind(t(C)%*%QI%*%C ,-QI%*%C),rbind(-t(C)%*%QI,QI))
  RI=diag(m)
  
  last_cost <-t(X)%*%QI2%*%X+t(u)%*%RI%*%u
  result <- list(Score=-last_cost/100000,Pred=0)
  
  return(result)
}

lower_lim=0.01
upper_lim_q=1
upper_lim_r=0.4

search_bounds <- list(q1=c(lower_lim,upper_lim_q),q2=c(lower_lim,upper_lim_q),q3=c(lower_lim,upper_lim_q),
                      q4=c(lower_lim,upper_lim_q),q5=c(lower_lim,upper_lim_q),
                      r1=c(lower_lim,upper_lim_r),r2=c(lower_lim,upper_lim_r),r3=c(lower_lim,upper_lim_r),
                      r4=c(lower_lim,upper_lim_r),r5=c(lower_lim,upper_lim_r),r6=c(lower_lim,upper_lim_r),
                      r7=c(lower_lim,upper_lim_r))

search_grid <- data.frame(q1=runif(50,lower_lim,upper_lim_q),q2=runif(50,lower_lim,upper_lim_q),
                          q3=runif(50,lower_lim,upper_lim_q),q4=runif(50,lower_lim,upper_lim_q),
                          q5=runif(50,lower_lim,upper_lim_q),
                          r1=runif(50,lower_lim,upper_lim_r),r2=runif(50,lower_lim,upper_lim_r),
                          r3=runif(50,lower_lim,upper_lim_r),r4=runif(50,lower_lim,upper_lim_r),
                          r5=runif(50,lower_lim,upper_lim_r),r6=runif(50,lower_lim,upper_lim_r),
                          r7=runif(50,lower_lim,upper_lim_r))

bayes_q_learning <- BayesianOptimization(FUN = fitness, bounds = search_bounds, 
                                         init_grid_dt = search_grid, init_points = 0, 
                                         n_iter = 100, acq = "ucb")

bayes_q_learning$Best_Par
