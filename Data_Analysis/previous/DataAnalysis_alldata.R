library(openxlsx)
library(dplyr)
library(vegan)
library(ape)
library(cubature)

source("../Find_grids.R")
source("../DART_function.R")
source("../FDRL_function.R")

load(file="Data.RData")#load the Data,contains 144 patients and 96 ASVs' relative abundances.
load(file="Dist0.RData")#load the distance matrix of the 96 ASVs
nodes.index <- colnames(clintest)[-(1:4)]
##################### Testing ###################
pvals <- NULL
for(i in nodes.index){
  form2 <- reformulate(c("time","inout"),response=i)
  fm2 <- lm(form2, data=clintest)
  
  pvals <- c(pvals,summary(fm2)$coef[3,4])
}

library(ggplot2)
Pvalues=pvals

#######################################################
######################## DART #########################
#######################################################
alpha=0.1 # set up desired FDR level


###### construct aggregation tree
set.seed(1234)
Dist <- Dist0
diag(Dist) <- 1000
md <- max(apply(Dist,1,min))
snm <- sqrt(144*log(96)*log(log(96)))
tpars <- find_pars(n=144,ntip=96,Dist=Dist0,Mgroup=3) # choose the tunning parameter
Atree <- A.tree.mult(grids=tpars$grids,Dist0=Dist0,Mgroup=tpars$M)

Rejs <- test.mult(Llist=Atree$Llist,Dist0=Dist0,
             T1=t(qnorm(pvals,lower.tail=FALSE)),alpha=alpha)

detect <- rep(0,length(nodes.index))
detect[Rejs[[1]]] <- 1
detect[Rejs[[2]]] <- 2
#detect <- factor(detect)


#write the testing result to a csv file
outdf <- data.frame("ASV"=nodes.index,"pvals"=pvals,"detect"=detect)
write.csv(outdf,"outdf_ANEW.csv")


#######################################################
####################### FDRL ##########################
#######################################################
ps <- fdrl.ii.3 <- fdrl.ii.2 <- fdrl.i.3 <- fdrl.i.2 <- rep(0,length(nodes.index))

fdrl.i.2[test.fdrl(T1p=pvals,Dist0,k=2,alpha=alpha)] <- 1
fdrl.i.3[test.fdrl(T1p=pvals,Dist0,k=3,alpha=alpha)] <- 1
fdrl.ii.2[test.fdrl2(T1p=pvals,Dist0,k=2,alpha=alpha)] <- 1
fdrl.ii.3[test.fdrl2(T1p=pvals,Dist0,k=3,alpha=alpha)] <- 1


#write the testing result to a csv file
outdf <- data.frame("ASV"=nodes.index,"pvals"=pvals,"detect1"=ifelse(detect==1,1,0),
                    "detect2"=ifelse(detect==2,1,0),
                    "FDRL.I.2"=fdrl.i.2,"FDRL.I.3"=fdrl.i.3,
                    "FDRL.II.2"=fdrl.ii.2,"FDRL.II.3"=fdrl.ii.3,
                    "ps"=0)
#outdf$detect1[54]=1
write.csv(outdf,"outdf_ANEW_ALL.csv")
