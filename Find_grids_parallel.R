library(ape)
library(cubature)
library(xtable)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyverse)
library(data.table)

source("../DART_function.R")

# # # # Generate location of hypotheses
nodes.2D.simu <- function(ntips){
  # function to generate the hypotheses locations
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

################################
#### Find_grids in Simulation
################################
Mgroup=3
set.seed(123)
ntip=1000;n=300;ncore=40;cm=100#56,96,M=2
locs <- nodes.2D.simu(ntip)
Dist <- NULL

t1 <- Sys.time()
snm <- sqrt(n*log(ntip)*log(log(ntip)))
Maxlayers <- ceiling(log(ntip/cm,Mgroup))
Maxlayers

if(is.null(Dist)){
  Dist <- as.matrix(dist(locs))
}
Dist0 <- Dist
diag(Dist) <- 10000*max(Dist)
md <- max(apply(Dist,1,min))

i.best <- 0
grids <- NULL

cl=makeCluster(ncore)
registerDoParallel(cl)
for(kk in 2:Maxlayers){
  if(i.best+4>md*snm){
    grids <- c(grids,i.best)
  }else{
    gridi <- seq(i.best+4,md*snm,4)
    
    outv=foreach(grid=gridi,.combine=rbind)%dopar%{
      Atree <- A.tree.mult(tree=locs,Dist0=Dist0,
                           grids=c(grids,grid)/snm,Mgroup=Mgroup)
      c("ct"=sum(rowSums(Atree$Llist[[kk-1]])>=2),"grid"=grid,"Md"=Atree$Md[[kk-1]])
    }
    
    outv1 <- data.frame(outv)%>%filter(!duplicated(ct))%>%top_n(1,ct)
    i.best <- outv1$grid
    grids <- c(grids,i.best)
    md <- outv1$Md
  }
}
stopCluster(cl)
t2 <- Sys.time()

Mgroup;grids


################################
#### Find_grids in Data analysis
################################

set.seed(123)

load(file="Dist0_new.RData")
load(file="clintest_new.RData")

Mgroup=3
Dist <- Dist0
ntip=nrow(Dist);n=456;ncore=40;cm=8#56,96,M=2

t1 <- Sys.time()
snm <- sqrt(n*log(ntip)*log(log(ntip)))
Maxlayers <- ceiling(log(ntip/cm,Mgroup))

Dist0 <- Dist
diag(Dist) <- 10000*max(Dist)
md <- max(apply(Dist,1,min))

i.best <- 0
grids <- NULL

cl=makeCluster(ncore)
registerDoParallel(cl)
for(kk in 2:Maxlayers){
  if(i.best+4>md*snm){
    grids <- c(grids,i.best)
  }else{
    gridi <- seq(i.best+4,md*snm,4)
    outv=foreach(grid=gridi,.combine=rbind)%dopar%{
      Atree <- A.tree.mult(tree=locs,Dist0=Dist0,
                           grids=c(grids,grid)/snm,Mgroup=Mgroup)
      c("ct"=sum(rowSums(Atree$Llist[[kk-1]])>=2),"grid"=grid,"Md"=Atree$Md[[kk-1]])
    }
    outv1 <- data.frame(outv)%>%filter(!duplicated(ct))%>%top_n(1,ct)
    i.best <- outv1$grid
    grids <- c(grids,i.best)
    md <- outv1$Md
  }
}
stopCluster(cl)
t2 <- Sys.time()
#281,311,315,317
Mgroup;grids

#283,312,316,318




#[1]2
#[1]76  96 100

#[1]3
#[1]56 96

#[1]4
#[1]16

#[1]5
#[1]12


