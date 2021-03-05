library(ape)
library(cubature)

library(foreach)
library(doParallel)
library("doRNG")
library(xtable)
library(tidyverse)
library("MASS")
library(LaplacesDemon)
library(survival)

source("../DART_function.R")
source("../FDRL_function.R")

# function to generate the tree
nodes.2D.simu <- function(ntips){
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

# Generate tree and corresponding aggregation tree
set.seed(123)
ntip=100;n=90                                                                                             
toy <- nodes.2D.simu(ntip)
grids <- 16
snm <- sqrt(n*log(ntip)*log(log(ntip)))
Atree <- A.tree.mult(toy,grids=grids/snm,Mgroup=Inf)

# set up replication times B.
nsim=20
##################### Abundance (SE4)########################
# Generate datas (TEST STATSITICS)
simdata <- function(n=500,nnodes,Dist0){
  # first assume \beta_{0j}=\beta_{2j}=\beta_{3j}
  
  trt <- rbinom(n,1,0.5)
  time <- runif(n,0.1,0.5) # only consider two time points
  
  scale1 <- scale2 <- 1
  
  b11 <- 0
  b12 <- dnorm(Dist0[7,],0,0.1)*1
  b13 <- dnorm(Dist0[22,],0,1)*2-0.2
  b1 <- (b11+b13+b12)
  beta1 <- ifelse(b1>0.15,b1,0)
  beta1=beta1*2
  
  beta0 <- rep(0.1,nnodes)
  beta2 <- rep(0.1,nnodes)
  beta3 <- rep(0,nnodes)
  
  eij0 <- rnorm(n*nnodes)
  
  Y <- NULL
  for(j in 1:nnodes){
    Y <- cbind(Y,beta0[j]+beta1[j]*trt+beta2[j]*time+
                 beta3[j]*time*trt+eij0[(j-1)*n+1:n])
  }
  
  X=cbind(1,"trt"=trt,"time"=time,"TT"=trt*time)
  
  sx <- solve(t(X)%*%X)
  beta.hat <- sx%*%t(X)%*%Y
  P <- X%*%sx%*%t(X)
  I <- diag(rep(1,n))
  sigma.hat <- t(rep(1,n))%*%((I-P)%*%Y)^2/(n-3)
  
  T1 <- beta.hat[2,]/(sqrt(sx[2,2]*sigma.hat))
  T1 <- qnorm(2*pnorm(abs(T1),lower.tail=FALSE),lower.tail=FALSE)
  
  return(list(stat=T1,beta1=beta1))
}



registerDoParallel(cores=15)
nsim=200
alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim*length(alphai), 5570)
result <- foreach(j=1:length(alphai),.combine=rbind)%:%foreach(i=c(1:nsim),r=rng[(j-1)*nsim+1:nsim],
                                                               .combine=rbind,.packages=c("cubature"))%dopar% {
                                                                 # set RNG seed
                                                                 rngtools::setRNG(r)
                                                                 mt.mult(n=90,alpha=alphai[j],
                                                                         Llist=Atree$Llist,
                                                                         Dist0=Atree$Dist0)}

write.csv(result,"abundance100_L3.csv")



result <- read.csv("abundance100_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

print(xtable(out,digits=4))


registerDoParallel(cores=15)
nsim=200
ki=c(2,3)
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=90,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})

result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]



outall <- rbind(out,out.fdrl)

write.csv(outall,"abundance_100_out.csv")




#####################nlap (SE2)############################

simdata <- function(n=90,nnodes,Dist0,normalp=0){
  scale1 <- scale2 <- 1
  
  b11 <- 0
  b12 <- dnorm(Dist0[7,],0,0.1)*1
  b13 <- dnorm(Dist0[22,],0,1)*2-0.2
  b1 <- (b11+b13+b12)
  beta1 <- ifelse(b1>0.15,b1,0)*sqrt(n)*0.6
  
  norm <- rbinom(nnodes,1,normalp)
  T1 <- rlaplace(nnodes,location=beta1)*(1-norm)+rnorm(nnodes,mean=beta1)*norm
  T1 <- plaplace(-abs(T1),location=0)*(1-normalp)+pnorm(abs(T1),mean=0,lower.tail=FALSE)*normalp
  T1 <- qnorm(2*T1,lower.tail=FALSE)
  T1 <- matrix(T1,1,nnodes)
  
  return(list(stat=T1,beta1=beta1))
}

registerDoParallel(cores=15)
alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim, 5555)

result <- foreach(i=c(1:nsim),r=rng,
                  .combine=rbind,.packages=c("cubature","LaplacesDemon"))%dopar% {
                    # set RNG seed
                    rngtools::setRNG(r)
                    mt.mult(n=90,alphai=alphai,
                            Llist=Atree$Llist,
                            Dist0=Atree$Dist0)}

write.csv(result,"nlap100_L3.csv")

result <- read.csv("nlap100_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

registerDoParallel(cores=15)
ki=c(2,3)
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature","LaplacesDemon"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=90,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})

result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]



outall <- rbind(out,out.fdrl)

write.csv(outall,"nlap_100_out.csv")

#####################simple (SE1)########################

simdata <- function(n=90,nnodes,Dist0){
  scale1 <- scale2 <- 1
  
  b11 <- 0
  b12 <- dnorm(Dist0[7,],0,0.1)*1
  b13 <- dnorm(Dist0[22,],0,1)*2-0.2
  b1 <- (b11+b13+b12)
  beta1 <- ifelse(b1>0.15,b1,0)
  beta1=beta1/2.5
  
  X <- mvrnorm(n=n,mu=beta1,Sigma <- diag(rep(1,nnodes)))
  
  T1 <- matrix(sqrt(n)*colMeans(X),1,nnodes)
  T1 <- qnorm(2*pnorm(abs(T1),lower.tail=FALSE),lower.tail=FALSE)
  return(list(stat=T1,beta1=beta1))
}

registerDoParallel(cores=15)
alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim, 5555)

result <- foreach(i=c(1:nsim),r=rng,
                  .combine=rbind,.packages=c("cubature","MASS"))%dopar% {
                    # set RNG seed
                    rngtools::setRNG(r)
                    mt.mult(n=90,alphai=alphai,
                            Llist=Atree$Llist,
                            Dist0=Atree$Dist0)}

write.csv(result,"simple100_L3.csv")



result <- read.csv("simple100_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

print(xtable(out,digits=4))


registerDoParallel(cores=15)
ki=c(2,3)
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature","MASS"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=90,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})

result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]



outall <- rbind(out,out.fdrl)

write.csv(outall,"simple_100_out.csv")


##################### 2T (SE3)##############################

simdata <- function(n=90,nnodes,Dist0,normalp=0,df=2){
  scale1 <- scale2 <- 1
  
  b11 <- 0
  b12 <- dnorm(Dist0[7,],0,0.1)*1
  b13 <- dnorm(Dist0[22,],0,1)*2-0.2
  b1 <- (b11+b13+b12)
  beta1 <- ifelse(b1>0.15,b1,0)
  
  beta1=beta1/6
  
  X <- mvrnorm(n=n,mu=beta1,Sigma <- diag(rep(1,nnodes)))
  T1 <- sqrt(n)*colMeans(X)
  norm <- rbinom(nnodes,1,normalp)
  T1 <- rt(nnodes,ncp=beta1*sqrt(n),df=df)*(1-norm)+T1*norm
  T1 <- pt(-abs(T1),ncp=beta1*sqrt(n),df=df)*(1-normalp)+pnorm(abs(T1),mean=0,lower.tail=FALSE)*normalp
  T1 <- qnorm(2*T1,lower.tail=FALSE)
  T1[is.na(T1)]=100
  T1 <- matrix(T1,1,nnodes)
  
  return(list(stat=T1,beta1=beta1))
}

registerDoParallel(cores=15)
alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim, 5555)

result <- foreach(i=c(1:nsim),r=rng,
                  .combine=rbind,.packages=c("cubature"))%dopar% {
                    # set RNG seed
                    rngtools::setRNG(r)
                    mt.mult(n=90,alphai=alphai,
                            Llist=Atree$Llist,
                            Dist0=Atree$Dist0)}

write.csv(result,"T2100_L3.csv")



result <- read.csv("T2100_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

head(out)

print(xtable(out,digits=4))


registerDoParallel(cores=15)
ki=c(2,3)
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=90,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})

result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]

outall <- rbind(out,out.fdrl)

write.csv(outall,"T2_100_out.csv")


########################## cox (SE5)##########################
simdata <- function(n=90,nnodes,Dist0){
  # first assume \beta_{0j}=\beta_{2j}=\beta_{3j}
  x1 <- rbinom(n*nnodes,1,0.5)
  x2 <- runif(n*nnodes,0.1,0.5) # only consider two time points
  
  scale1 <- scale2 <- 1
  
  b11 <- 0
  b12 <- dnorm(Dist0[7,],0,0.1)*1
  b13 <- dnorm(Dist0[22,],0,1)*2-0.2
  b1 <- (b11+b13+b12)
  beta1 <- ifelse(b1>0.15,b1,0)
  beta1=beta1*0.8
  
  beta2 <- rep(0.1,nnodes)
  T1 <- NULL
  me <- NULL
  for(j in 1:nnodes){
    atime <- rexp(n,exp(x1[(j-1)*n+1:n]*beta1[j]+
                          beta2[j]*x2[(j-1)*n+1:n]))
    ctime=runif(n,0,5)
    otime=pmin(atime,ctime)
    event=as.integer(atime<ctime)
    X1 <- x1[(j-1)*n+1:n]
    X2 <- x2[(j-1)*n+1:n]
    
    out <- coxph(Surv(otime,event)~X1+X2)
    me <- c(me,mean(event))
    T1 <- c(T1,qnorm(coef(summary(out))[1,5],lower.tail=FALSE))
  }
  T1 <- matrix(T1,1,nnodes)
  
  return(list(stat=T1,beta1=beta1))
}


registerDoParallel(cores=15)
alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim, 5555)

result <- foreach(i=c(1:nsim),r=rng,
                  .combine=rbind,.packages=c("cubature","survival"))%dopar% {
                    # set RNG seed
                    rngtools::setRNG(r)
                    mt.mult(n=90,alphai=alphai,
                            Llist=Atree$Llist,
                            Dist0=Atree$Dist0)}

write.csv(result,"cox100_L3.csv")



result <- read.csv("cox100_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

head(out)

print(xtable(out,digits=4))


registerDoParallel(cores=15)
ki=c(2,3)
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature","survival"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=90,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})

result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]



outall <- rbind(out,out.fdrl)

write.csv(outall,"cox_100_out.csv")


q(save='no')


