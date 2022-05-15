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
nsim=200

# function to generate the hypotheses locations
nodes.2D.simu <- function(ntips){
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

# Generate tree and corresponding aggregation tree
set.seed(123)
ntip=1000;n=300                                                                                             
toy <- nodes.2D.simu(ntip)
ceiling(log(ntip/80,3))
#grids <- c(8,22,56)
#grids <- c(38,58,72)
#grids <- c(22,32,60)
grids <- c(56,96)
snm <- sqrt(n*log(ntip)*log(log(ntip)))
t1 <- Sys.time()
Atree <- A.tree.mult(tree=toy,grids=grids/snm,Mgroup=3)
#Atree <- A.tree.mult2(tree=toy,Mgroup=3)
#save(Atree,file="Pars_1000.RData")
# t2 <- Sys.time()
# t2-t1
#load("Pars_1000.RData")
# set up replication times B.

##################### Abundance (SE3)########################
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
  m1 <- sum(beta1!=0)
  indc <- sample(which(beta1!=0),floor(nullprop*m1))
  indc0 <- sample(which(beta1==0),floor(nullprop*m1))
  beta1[indc0]=beta1[indc]
  beta1[indc]=0
  
  #indc <- rbinom(nnodes,1,1-nullprop)
  #selcb <- beta1[indc==0&beta1!=0]
  #beta1[indc==0&beta1!=0]=0
  #beta1[sample(which(beta1==0),length(selcb))]=selcb

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

cl=makeCluster(80)
registerDoParallel(cl)

alphai=c(0.05,0.1,0.15,0.20)
sei=1:4
npi=c(0,0.02,0.04,0.06,0.08,0.1)
rng <- RNGseq(nsim, 5555)

result_DART <- foreach(j=seq(alphai),.combine=rbind)%:%
  foreach(se=sei,.combine=rbind)%:%
  foreach(np=seq(npi),.combine=rbind)%:%
  foreach(i=c(1:nsim),r=rng[1:nsim],
          .combine=rbind,.packages=c("cubature","MASS","survival"))%dopar%{
    # set RNG seed
    rngtools::setRNG(r)
    mt.mult(n=300,alpha=alphai[j],
            Llist=Atree$Llist,
            Dist0=Atree$Dist0,SE=se,nullprop=npi[np])}

out_DART <- data.frame(result_DART) %>% group_by(alpha,Layers,SE,nullprop) %>% 
  summarize(FDP=mean(FDR),Sensitivity=mean(Power)) %>%
  mutate(Test= ifelse(Layers==1,"1 Layer",paste0(Layers," Layers")))

save(result_DART,file="Out_DART2.RData")

result_FDRL <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(se=sei,.combine=rbind)%:%
             foreach(np=seq(npi),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[1:nsim],
                     .combine=rbind,.packages=c("cubature","MASS","survival"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=300,alpha=alphai[j],SE=se,nullprop=npi[np],
                                Dist0=Atree$Dist0,k=5)})
stopCluster(cl)

out_FDRL <- data.frame(result_FDRL)%>%
  group_by(alpha,k,SE,nullprop) %>% 
  summarize(FDP=mean(fdr),Sensitivity=mean(power))%>%
  mutate(Test=paste0("FDRL,k=",k))

save(result_FDRL,file="Out_FDRL2.RData")



q(save='no')


