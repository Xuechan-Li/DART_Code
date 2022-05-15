### This is a simple sample for using the code, the results are not included in the paper

library("MASS")

setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/project_1/DART_Code/Simulation")


source("../DART_function.R")

source("../Find_grids.R")

########## Generate data: p-value and distance matrix ###########

# function to generate the hypotheses locations
nodes.2D.simu <- function(ntips){
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

# function to get the p-value
simdata <- function(n=90,nnodes,Dist0){
  ##arguments:
  #n: number of smaples
  #nnodes: number of hypothesis (denoted by m in the paper)
  #Dist0: Distantce matrix
  ##outputs:
  #stat: pvalues
  #beta1: true parameters of interest
  
 
  #setting parameters
  b11 <- 0
  b12 <- dnorm(Dist0[7,],0,0.1)*1
  b13 <- dnorm(Dist0[22,],0,1)*2-0.2
  b1 <- (b11+b13+b12)
  beta1 <- ifelse(b1>0.15,b1,0)
  beta1=beta1/2.5
  
  # generate observations
  X <- mvrnorm(n=n,mu=beta1,Sigma <- diag(rep(1,nnodes)))
  
  # get p-values
  T1 <- matrix(sqrt(n)*colMeans(X),1,nnodes)
  T1 <- qnorm(2*pnorm(abs(T1),lower.tail=FALSE),lower.tail=FALSE)
  return(list(stat=T1,beta1=beta1))
}

# Generate hypotheses locations and calculate the distance.
set.seed(123)
ntip=100;n=90                                                                                             
toy <- nodes.2D.simu(ntip) #hypotheses locations in a 2-D Eucleadian space.
Dist0 <- as.matrix(dist(toy)) #distance matrix
# Generate simulated data
data <- simdata(n=n,nnodes=100,Dist0=Dist0)

########## Use DART to conduct testing ###########

# Generate aggregation tree
Atree <- A.tree.mult2(Mgroup = 3, Dist0 = Dist0)
# Conduct multiple testing embedded in the aggregation tree
Test.result <- test.mult(alpha=0.05,Llist=Atree$Llist,Dist0=Dist0,T1=data$stat)
#Results: (Index of rejected nodes)
Test.result[[tpars$L]]


q(save="no")