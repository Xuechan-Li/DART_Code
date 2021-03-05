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
grids <- c(8,22,56)
snm <- sqrt(n*log(ntip)*log(log(ntip)))
Atree <- A.tree.mult(toy,grids=grids/snm,Mgroup=3)

# set up replication times B.
nsim=200
##################### Abundance (SE4)########################
# Generate datas (TEST STATSITICS)
simdata <- function(n=300,nnodes,Dist0){
  # first assume \beta_{0j}=\beta_{2j}=\beta_{3j}
  
  trt <- rbinom(n,1,0.5)
  time <- runif(n,0.1,0.5) # only consider two time points
  
  scale1 <- scale2 <- 1
  
  b11 <- pmax(dnorm(Dist0[156,],0,0.8)*3.4-0.8,0)
  b12 <- dnorm(Dist0[7,],0,0.05)*3
  b13 <- -dnorm(Dist0[90,],0,0.15)*2
  b1 <- (b11+b12)
  #420,348
  b1[c(300,200,400,100,500,600,700,800,900,1000)] <- 10
  beta1 <- ifelse(abs(b1)>0.15,b1,0)
  
  beta1=beta1/1.2
  
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
  sigma.hat <- t(rep(1,n))%*%((I-P)%*%Y)^2/n
  
  T1 <- beta.hat[2,]/(sqrt(sx[2,2]*sigma.hat))
  T1 <- qnorm(2*pnorm(abs(T1),lower.tail=FALSE),lower.tail=FALSE)
  
  return(list(stat=T1,beta1=beta1))
}

cl=makeCluster(10)
registerDoParallel(cl)

alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim*length(alphai), 5570)
result <- foreach(j=1:length(alphai),.combine=rbind)%:%foreach(i=c(1:nsim),r=rng[(j-1)*nsim+1:nsim],
                                                               .combine=rbind,.packages=c("cubature"))%dopar% {
                                                                 # set RNG seed
                                                                 rngtools::setRNG(r)
                                                                 mt.mult(n=300,alpha=alphai[j],
                                                                         Llist=Atree$Llist,
                                                                         Dist0=Atree$Dist0)}

stopCluster(cl)
write.csv(result,"abundance_L3.csv")



result <- read.csv("abundance_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

print(xtable(out,digits=4))

cl=makeCluster(10)
registerDoParallel(cl)
ki=6:7
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=300,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})
stopCluster(cl)
result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]



outall <- rbind(out,out.fdrl)

write.csv(outall,"abundance_1000_out.csv")


#####################nlap (SE2)############################

simdata <- function(n=300,nnodes,Dist0,normalp=0){
  scale1 <- scale2 <- 1

  b11 <- pmax(dnorm(Dist0[156,],0,0.8)*3.4-0.8,0)
  b12 <- dnorm(Dist0[7,],0,0.05)*3
  b13 <- -dnorm(Dist0[90,],0,0.15)*2
  b1 <- (b11+b12)
  #420,348
  b1[c(300,200,400,100,500,600,700,800,900,1000)] <- 10
  beta1 <- ifelse(abs(b1)>0.15,b1,0)*sqrt(n)/6
  
  norm <- rbinom(nnodes,1,normalp)
  T1 <- rlaplace(nnodes,location=beta1)*(1-norm)+rnorm(nnodes,mean=beta1)*norm
  T1 <- plaplace(-abs(T1),location=0)*(1-normalp)+pnorm(abs(T1),mean=0,lower.tail=FALSE)*normalp
  T1 <- qnorm(2*T1,lower.tail=FALSE)
  T1 <- matrix(T1,1,nnodes)
  
  return(list(stat=T1,beta1=beta1))
}

# function to generate the tree
nodes.2D.simu <- function(ntips){
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

cl=makeCluster(10)
registerDoParallel(cl)
alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim, 5555)

result <- foreach(i=c(1:nsim),r=rng,
                  .combine=rbind,.packages=c("cubature","LaplacesDemon"))%dopar% {
                    # set RNG seed
                    rngtools::setRNG(r)
                    mt.mult(n=300,alphai=alphai,
                            Llist=Atree$Llist,
                            Dist0=Atree$Dist0)}
stopCluster(cl)
write.csv(result,"nlap_L3.csv")

result <- read.csv("nlap_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

cl=makeCluster(10)
registerDoParallel(cl)
ki=6:7
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature","LaplacesDemon"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=300,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})
stopCluster(cl)
result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]



outall <- rbind(out,out.fdrl)

write.csv(outall,"nlap_1000_out.csv")

#####################simple (SE1)########################

simdata <- function(n=300,nnodes,Dist0){
  scale1 <- scale2 <- 1
  
  b11 <- pmax(dnorm(Dist0[156,],0,0.8)*3.4-0.8,0)
  b12 <- dnorm(Dist0[7,],0,0.05)*3
  b13 <- -dnorm(Dist0[90,],0,0.15)*2
  b1 <- (b11+b12)
  b1[c(300,200,400,100,500,600,700,800,900,1000)] <- 10
  beta1 <- ifelse(abs(b1)>0.15,b1,0)
  
  beta1=beta1/7
  
  X <- mvrnorm(n=n,mu=beta1,Sigma <- diag(rep(1,nnodes)))
  
  T1 <- matrix(sqrt(n)*colMeans(X),1,nnodes)
  T1 <- qnorm(2*pnorm(abs(T1),lower.tail=FALSE),lower.tail=FALSE)
  return(list(stat=T1,beta1=beta1))
}

# function to generate the tree
nodes.2D.simu <- function(ntips){
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

cl=makeCluster(10)
registerDoParallel(cl)
alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim, 5555)

result <- foreach(i=c(1:nsim),r=rng,
                  .combine=rbind,.packages=c("cubature","MASS"))%dopar% {
                    # set RNG seed
                    rngtools::setRNG(r)
                    mt.mult(n=300,alphai=alphai,
                            Llist=Atree$Llist,
                            Dist0=Atree$Dist0)}
stopCluster(cl)
write.csv(result,"simple_L3.csv")



result <- read.csv("simple_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

print(xtable(out,digits=4))


cl=makeCluster(10)
registerDoParallel(cl)
ki=6:7
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature","MASS"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=300,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})
stopCluster(cl)
result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]



outall <- rbind(out,out.fdrl)

write.csv(outall,"simple_1000_out.csv")

##################### 2T (SE3) ##############################

simdata <- function(n=300,nnodes,Dist0,normalp=0,df=5){
  scale1 <- scale2 <- 1
  
  b11 <- pmax(dnorm(Dist0[156,],0,0.8)*3.4-0.8,0)
  b12 <- dnorm(Dist0[7,],0,0.05)*3
  b13 <- -dnorm(Dist0[90,],0,0.15)*2
  b1 <- (b11+b12)
  #420,348
  b1[c(300,200,400,100,500,600,700,800,900,1000)] <- 10
  beta1 <- ifelse(abs(b1)>0.15,b1,0)
  
  beta1=beta1/14
  
  X <- mvrnorm(n=n,mu=beta1,Sigma <- diag(rep(1,nnodes)))
  T1 <- matrix(sqrt(n)*colMeans(X),1,nnodes)
  norm <- rbinom(nnodes,1,normalp)
  T1 <- rt(nnodes,ncp=beta1*sqrt(n),df=df)*(1-norm)+T1*norm
  T1 <- pt(-abs(T1),ncp=beta1*sqrt(n),df=df)*(1-normalp)+pnorm(abs(T1),mean=0,lower.tail=FALSE)*normalp
  T1 <- qnorm(2*T1,lower.tail=FALSE)
  T1 <- matrix(T1,1,nnodes)
  
  return(list(stat=T1,beta1=beta1))
}

# function to generate the tree
nodes.2D.simu <- function(ntips){
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}


cl=makeCluster(10)
registerDoParallel(cl)
alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim, 5555)

result <- foreach(i=c(1:nsim),r=rng,
                  .combine=rbind,.packages=c("cubature","MASS"))%dopar% {
                    # set RNG seed
                    rngtools::setRNG(r)
                    mt.mult(n=300,alphai=alphai,
                            Llist=Atree$Llist,
                            Dist0=Atree$Dist0)}
stopCluster(cl)
write.csv(result,"T2_L3.csv")



result <- read.csv("T2_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

head(out)

print(xtable(out,digits=4))


cl=makeCluster(10)
registerDoParallel(cl)
ki=6:7
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature","MASS"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=300,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})
stopCluster(cl)
result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]

outall <- rbind(out,out.fdrl)

write.csv(outall,"T2_1000_out.csv")

########################## cox (SE5)##########################

simdata <- function(n=300,nnodes,Dist0){
  # first assume \beta_{0j}=\beta_{2j}=\beta_{3j}
  x1 <- rbinom(n,1,0.5)
  x2 <- runif(n,0.1,0.5) # only consider two time points
  
  scale1 <- scale2 <- 1
  
  b11 <- pmax(dnorm(Dist0[156,],0,0.8)*3.4-0.8,0)
  b12 <- dnorm(Dist0[7,],0,0.05)*3
  b13 <- -dnorm(Dist0[90,],0,0.15)*2
  b1 <- (b11+b12)
  #420,348
  b1[c(300,200,400,100,500,600,700,800,900,1000)] <- 10
  beta1 <- ifelse(abs(b1)>0.15,b1,0)
  beta1=beta1/3.5
  
  beta2 <- rep(0.1,nnodes)
  T1 <- NULL
  me <- NULL
  for(j in 1:nnodes){
    atime <- rexp(n,exp(x1*beta1[j]+
                          beta2[j]*x2))
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
  
  return(list(stat=T1,beta1=beta1))
}

# function to generate the tree
nodes.2D.simu <- function(ntips){
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

cl=makeCluster(10)
registerDoParallel(cl)
alphai=c(0.05,0.1,0.15,0.20)
rng <- RNGseq(nsim, 5555)

result <- foreach(i=c(1:nsim),r=rng,
                  .combine=rbind,.packages=c("cubature","survival"))%dopar% {
                    # set RNG seed
                    rngtools::setRNG(r)
                    mt.mult(n=300,alphai=alphai,
                            Llist=Atree$Llist,
                            Dist0=Atree$Dist0)}
stopCluster(cl)
write.csv(result,"cox_L3.csv")



result <- read.csv("cox_L3.csv")
result <- result[,-1]
result %>% group_by(alpha,Layers) %>% summarize(FDP=mean(FDR),Sensitivity=mean(Power))->out
out$Test <- ifelse(out$Layers==1,"1 Layer",paste0(out$Layers," Layers"))
out <- out[,-2]

head(out)

print(xtable(out,digits=4))


cl=makeCluster(10)
registerDoParallel(cl)
ki=6:7
rng <- RNGseq(nsim*length(alphai)*length(ki), 33344)
result <- (foreach(j=1:length(alphai),.combine=rbind)%:%
             foreach(i=1:length(ki),.combine=rbind)%:%
             foreach(ii=c(1:nsim),r=rng[((j-1)*length(alphai)+i-1)*nsim+1:nsim],
                     .combine=rbind,.packages=c("cubature","survival"))%dopar% {
                       # set RNG seed
                       rngtools::setRNG(r)
                       mt.fdrl2(n=300,alpha=alphai[j],
                                Dist0=Atree$Dist0,k=ki[i])})

stopCluster(cl)

result <- data.frame(result)
result %>% group_by(alpha,k) %>% summarize(FDP=mean(fdr),Sensitivity=mean(power))->out.fdrl


out.fdrl$Test <- paste0("FDRL,k=",out.fdrl$k)

ki <- table(out.fdrl$k)

out.fdrl <- out.fdrl[,-2]

outall <- rbind(out,out.fdrl)

write.csv(outall,"cox_1000_out.csv")

q(save='no')


