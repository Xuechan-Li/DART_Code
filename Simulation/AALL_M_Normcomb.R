library(ape)
library(cubature)
library(foreach)
library(doParallel)
library("doRNG")
library(xtable)
library(tidyverse)
library("MASS")
library(survival)

source("../DART_function.R")
source("../FDRL_function.R")


# function to generate the hypotheses locations
nodes.2D.simu <- function(ntips){
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

# Generate datas (TEST STATSITICS)
simdata <- function(n=300,nnodes,Dist0,nullprop=0.01,SE=1,df=5,normalp=0.96){
  # SE1: simple
  # SE2: t
  # SE3: linear
  # SE4: cox
  
  #parameters set-up
  b11 <- pmax(dnorm(Dist0[156,],0,0.8)*3.4-0.8,0)
  b12 <- dnorm(Dist0[7,],0,0.05)*3
  b13 <- -dnorm(Dist0[90,],0,0.15)*2
  b1 <- (b11+b12)
  b1[c(300,200,400,100,500,600,700,800,900,1000)] <- 10
  
  beta1 <- ifelse(abs(b1)>0.15,b1,0)
  indc <- rbinom(nnodes,1,1-nullprop)
  selcb <- beta1[indc==0&beta1!=0]
  beta1[indc==0&beta1!=0]=0
  beta1[sample(which(beta1==0),length(selcb))]=selcb
  
  if(SE==1){
    beta1=beta1/5 # set the parameter scale to avoid too large/small power
    X <- mvrnorm(n=n,mu=beta1,Sigma <- diag(rep(1,nnodes)))
    T1 <- matrix(sqrt(n)*colMeans(X),1,nnodes)
    T1 <- qnorm(2*pnorm(abs(T1),lower.tail=FALSE),lower.tail=FALSE)
  }else if(SE==2){
    beta1=beta1/5
    X <- mvrnorm(n=n,mu=beta1,Sigma <- diag(rep(1,nnodes)))
    T1 <- matrix(sqrt(n)*colMeans(X),1,nnodes)
    norm <- rbinom(nnodes,1,normalp)
    T1 <- rt(nnodes,ncp=beta1*sqrt(n),df=df)*(1-norm)+T1*norm
    T1 <-  pnorm(-abs(T1))
    T1 <- qnorm(2*T1,lower.tail=FALSE)
    T1 <- matrix(T1,1,nnodes)
  }else if(SE==3){
    beta1=beta1/2.5
    beta0 <- rep(0.1,nnodes)
    beta2 <- rep(0.1,nnodes)
    
    trt <- rbinom(n,1,0.5)
    time <- runif(n,0.1,0.5)
    eij0 <- rnorm(n*nnodes)
    
    Y <- NULL
    for(j in 1:nnodes){
      Y <- cbind(Y,beta0[j]+beta1[j]*trt+beta2[j]*time+eij0[(j-1)*n+1:n])
    }
    
    X=cbind(1,"trt"=trt,"time"=time)
    sx <- solve(t(X)%*%X)
    beta.hat <- sx%*%t(X)%*%Y
    P <- X%*%sx%*%t(X)
    I <- diag(rep(1,n))
    sigma.hat <- t(rep(1,n))%*%((I-P)%*%Y)^2/n
    
    T1 <- beta.hat[2,]/(sqrt(sx[2,2]*sigma.hat))
    T1 <- qnorm(2*pnorm(abs(T1),lower.tail=FALSE),lower.tail=FALSE)
  }else{
    x1 <- rbinom(n,1,0.5)
    x2 <- runif(n,0.1,0.5) 
    #beta1=beta1/2.5
    beta1=beta1/2
    beta2 <- rep(0.1,nnodes)
    T1 <- NULL
    me <- NULL
    for(j in 1:nnodes){
      atime <- rexp(n,exp(x1*beta1[j]+beta2[j]*x2))
      ctime=runif(n,0,5)
      otime=pmin(atime,ctime)
      event=as.integer(atime<ctime)
      X1 <- x1
      X2 <- x2
      
      out <- coxph(Surv(otime,event)~X1+X2)
      me <- c(me,mean(event))
      T1 <- c(T1,qnorm(coef(summary(out))[1,5],lower.tail=FALSE))
    }
    T1 <- matrix(T1,1,nnodes)
  }
  
  return(list(stat=T1,beta1=beta1))
}


# Generate tree and corresponding aggregation tree
set.seed(123)
ntip=1000;n=300                                                                                             
toy <- nodes.2D.simu(ntip)

#[1]2
#[1]76  96 100

#[1]3
#[1]56 96

#[1]4
#[1]16

#[1]5
#[1]12

alphai=c(0.05,0.1,0.15,0.20)
sei=1:4
npi=c(0)
Mi=c(2:5)

grids <- list(c(76,96,110),c(56,96),c(16),c(12))
snm <- sqrt(n*log(ntip)*log(log(ntip)))
Atrees <- list()
for(M in Mi){
  Atrees[[M-1]] <- A.tree.mult(tree=toy,grids=grids[[M-1]]/snm,Mgroup=M)
}


t1 <- Sys.time()

nsim=200

cl=makeCluster(50)
registerDoParallel(cl)
rng <- RNGseq(nsim, 5555)

result_DART <- foreach(se=sei,.combine=rbind)%:%
  foreach(M=Mi,.combine=rbind)%:%
  foreach(i=c(1:nsim),r=rng[1:nsim],
          .combine=rbind,.packages=c("cubature","MASS","survival","tidyverse"))%dopar%{
            # set RNG seed
            rngtools::setRNG(r)
            Atree <- Atrees[[M-1]]
            out=mt.mult(n=300,alpha=alphai,
                    Llist=Atree$Llist,
                    Dist0=Atree$Dist0,SE=se,nullprop=0)
            out <- data.frame(out) %>% mutate(M=M)
            out}


save(result_DART,file="Results/Out_DART_DiffM.RData")


q(save='no')


