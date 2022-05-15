library(openxlsx)
library(dplyr)
library(vegan)
library(ape)
library(cubature)
library(xtable)
library(reshape2)
library(foreach)
library(doParallel)
library(data.table)
library(tidyverse)

source("../DART_function.R")
source("../FDRL_function.R")

load(file="clintest_new.RData")
load(file="Dist0_new.RData")


###### load Aggregation tree (Atree)
#load("Pars.RData")
Lmax=ceiling(log(nrow(Dist0)/50,3))
###### construct aggregation tree
set.seed(1234)
Dist <- Dist0
diag(Dist) <- 1000
md <- max(apply(Dist,1,min))
n= uniqueN(clintest$patid)
m <- length(nodes.index)
snm <- sqrt(n*log(m)*log(log(m)))
tpars <- find_pars(n=n,ntip=m,Dist=Dist0,Mgroup=3) 

t1 <- Sys.time()
Atree <- A.tree.mult(Dist0=Dist0,grids=tpars$grids,Mgroup=3)
t2 <- Sys.time()
t2-t1

###### Bootstrap sampling

boot <- function(clintest,Atree,alpha=0.1,seed,L=NULL){
  set.seed(seed)
  cin.index <- sample(1:nrow(clintest),size=nrow(clintest),replace=TRUE)
  
  clintest.train <- clintest[cin.index,]
  
  pvals <- NULL
  for(i in nodes.index){
    form2 <- reformulate(c("time","inout"),response=i)
    fm2 <- lm(form2, data=clintest.train)
    
    pvals <- c(pvals,summary(fm2)$coef[3,4])
  }
  
  #######################################################
  ######################## DART #########################
  #######################################################
  Rejs <- test.mult(Llist=Atree$Llist,Dist0=Dist0,
                    alpha=alpha,
                    T1=t(qnorm(pvals,lower.tail=FALSE)))
  
  detect <- rep(0,length(nodes.index))
  if(is.null(L)){
    Lmax=length(Atree$Llist)+1
  }else{
    Lmax=L
  }
  for(ell in rev(seq(Lmax))){
    detect[Rejs[[ell]]] <- ell
  }
  table(detect)
  
  
  #######################################################
  ####################### FDRL ##########################
  #######################################################
  fdrl.ii.3 <- fdrl.ii.2 <- fdrl.i.3 <- fdrl.i.2 <- rep(0,length(nodes.index))
  
  fdrl.i.2[test.fdrl(T1p=pvals,Dist0,k=3,alpha=alpha)] <- 1
  fdrl.i.3[test.fdrl(T1p=pvals,Dist0,k=5,alpha=alpha)] <- 1
  fdrl.ii.2[test.fdrl2(T1p=pvals,Dist0,k=3,alpha=alpha)] <- 1
  fdrl.ii.3[test.fdrl2(T1p=pvals,Dist0,k=5,alpha=alpha)] <- 1
  
  #write the testing result to a csv file
  outdf <- data.frame("ASV"=nodes.index,"pvals"=pvals,"Single.layer"=ifelse(detect==1,1,0),
                      "DART"=ifelse(detect>0,1,0),
                      "FDRL.I.2"=fdrl.i.2,"FDRL.I.3"=fdrl.i.3,
                      "FDRL.II.2"=fdrl.ii.2,"FDRL.II.3"=fdrl.ii.3)
  
  outdf.valid <- outdf
  rownames(outdf.valid) <- outdf.valid[,1]
  outdf.valid <- as.matrix(outdf.valid[,-1])
  
  return(outdf.valid)
}

cl=makeCluster(50)
registerDoParallel(cl)

nboot=200

Vresult=foreach(nboots=seq(nboot),.combine="+",
        .packages=c("cubature","MASS"))%dopar%boot(clintest=clintest,Atree=Atree,alpha=0.05,seed=nboots,L=2)
Vresult=Vresult/nboot

stopCluster(cl)
# validation output
save(Vresult,file="Vresult_new.RData")

q(save='no')