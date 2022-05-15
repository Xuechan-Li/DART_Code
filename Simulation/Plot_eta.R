setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/project_1/DART_Code_last/Simulation")
rm(list=ls())

library(ape)
library(cubature)
library(foreach)
library(doParallel)
library("doRNG")
library(xtable)
library(tidyverse)
library("MASS")
library(survival)
library(ggplot2)

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
ceiling(log(ntip/50,3))
grids <- c(56,96)
snm <- sqrt(n*log(ntip)*log(log(ntip)))

Dist0 <- as.matrix(dist(toy))

b11 <- pmax(dnorm(Dist0[156,],0,0.8)*3.4-0.8,0)
b12 <- dnorm(Dist0[7,],0,0.05)*3
b13 <- -dnorm(Dist0[90,],0,0.15)*2
b1 <- (b11+b12)
b1[c(300,200,400,100,500,600,700,800,900,1000)] <- 10
eta1 <- ifelse(abs(b1)>0.15,b1,0)

data <- data.frame(toy,eta1=eta1)%>%
  mutate(logeta1=log(log(eta1+1))+0.005,
         etac=ifelse(eta1==0,"Null","Alternative"))
pp=ggplot(data,aes(x=d1,y=d2))+
  geom_point(aes(color=etac,size=log(eta1+1)+0.01),alpha=0.5)+
  theme_classic()+ labs(size = "L-eta",color="Feature")+
  scale_size(range=c(0.01,2))
pdf("Results/eta_plot.pdf",height=2.5,width=4)
print(pp)
dev.off()

q(save='no')


