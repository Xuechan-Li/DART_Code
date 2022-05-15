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
nsim=30

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

grids <- c(56,96)
snm <- sqrt(n*log(ntip)*log(log(ntip)))
t1 <- Sys.time()
Atree <- A.tree.mult(tree=toy,grids=grids/snm,Mgroup=3)



q(save='no')


