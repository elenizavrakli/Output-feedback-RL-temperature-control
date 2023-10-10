library(Matrix) 
library(tidyverse)
library(netcontrol)

set.seed(10)

#System dimensions
n=6 #state length
m=7 #input length
p=5 #measurement length

layer_number =100 #length of simulation

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
R = diag(m)
Q = diag(p)
#Weighting matrices obtained from Bayesian Optimisation
#Q= diag(c(0.9434381, 0.7622189, 0.5422118, 0.4198321,0.05144575))
#R= diag(c( 0.2997687, 0.2696552, 0.2805254, 0.09184189, 0.05358372, 0.2685642	, 0.3179601))
C1<-as.matrix(cbind(diag(p),-diag(p)))
Q1 <- cbind(rbind(t(C)%*%Q%*%C ,-Q%*%C),rbind(-t(C)%*%Q,Q))

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

#Initialize the system state and the observer states
x <- rep(20,n)
xobs <-rep(50,n)
xvec <- C%*%x #vector storing past states, used to plot the trajectories
xvecobs <- C%*%x #vector storing past observer estimates 

#Define optimal state to be tracked
xopt <- c(150,160,170,175,180)

#Define augmented state
X<-c(xobs,xopt) 

#Initialise vector to store cost terms
cost<- matrix(0,nrow=layer_number,ncol=1)

#Control for first time step and corresponding cost
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

#Plot state trajectories
xvec<-t(xvec)
xvec<-as.data.frame(cbind(1:(layer_number),xvec[1:layer_number,]))
colnames(xvec) <- c("time","measurement 1","measurement 2","measurement 3","measurement 4","measurement 5")
xvec[,c("time","measurement 1","measurement 2","measurement 3","measurement 4","measurement 5")] %>%
  #xvec[-1,c("layer","heater 1","heater 2","heater 3")] %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title = element_text(size = 20))+
  #ylim(75,200)+
  #theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())+
  labs(x='time',y='heater temperature evolution')
#ggsave("observer bo trajectories.pdf",width=10,height=4)

#Error at last step
norm(as.matrix(xvec[layer_number,-1]-xopt),type="2")

#Plot observer estimation error
obs_err<-as_data_frame(cbind(1:(layer_number),rowNorms(as.matrix(xvec[,-1]-t(xvecobs)),method="euclidean")))
colnames(obs_err)<-c("time","obs err")
obs_err[-1,c("time","obs err")] %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  theme(legend.position = "none", axis.title = element_text(size = 20))+
  labs(x='time',y='observer error')
#ggsave("observer error.pdf",width=10,height=4)

#Plot performance index
cost_plot<-as_data_frame(cbind(1:(layer_number),cost))
colnames(cost_plot)<-c("time","cost sum")
cost_plot[,c("time","cost sum")] %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  theme(legend.position = "none",axis.title = element_text(size = 20))+
  labs(x='time',y='sum of discounted costs')
#ggsave("observer, perf index.pdf",width=10,height=4)

#bocost=cost #to use when simulating with the weighting costs optimised through Bayesian Optimisation
#Plot comparison of performance indices before and after B.O.
#comp_cost_plot<-as_data_frame(cbind(1:101,cost[1:101],bocost[1:101]))
#colnames(comp_cost_plot)<-c("time","without BO","with BO")
#comp_cost_plot[,c("time","without BO","with BO")] %>%
#  gather(variable,value,-time) %>%
#  ggplot(aes(x=time,y=value,color=variable))+
#  geom_line(size=1)+
#  theme_classic()+
#  theme(axis.title = element_text(size = 20))+
#  #theme(legend.position = "none")+
#   labs(x='time',y='sum of discounted costs')
#ggsave("observer cost comparison.pdf",width=10,height=4)


