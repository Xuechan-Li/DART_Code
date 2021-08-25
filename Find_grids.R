library(ape)
library(cubature)
library(xtable)

#source("DART_function.R")


find_pars <- function(locs=NULL,Dist=NULL,n,ntip,Mgroup=3,cm=30){
  #### function to find tunning parameters for aggregation tree construction.
  ## arguments:
  #locs: (matrix) hypotheses locations, not need to be defined is Distance matrix Dist is given.
  #n: (integer) number of samples per hypothesis.
  #ntip: (integer) number of hypothesis (also denoted by m in the paper).
  #Mgroup: (integer) Maximum cardinality.
  #cm: (integer) constant used for decide the maximum layer.
  ## outputs: (list) list of tunning parameters.
  
  snm <- sqrt(n*log(ntip)*log(log(ntip)))
  Maxlayers <- ceiling(log(ntip/cm,Mgroup))
  
  if(is.null(Dist)){
    Dist <- as.matrix(dist(locs))
  }
  Dist0 <- Dist
  diag(Dist) <- 10000*max(Dist)
  md <- max(apply(Dist,1,min))
  
  i.best <- 0
  grids <- md0 <- NULL
  
  for(kk in 2:Maxlayers){
    seq1 <- gi <- eqs <-  NULL
    i <- i.best+2
    eq <- 1
    md0 <- c(md0,md*(Mgroup^(kk-2)*2-1))
    
    while(i/snm<=md*(Mgroup^(kk-2)*2-1)&eq<=10){
      gridsi <- c(grids,i)
      Atree <- A.tree.mult(tree=locs,Dist0=Dist0,grids=gridsi/snm,Mgroup=Mgroup)
      seq1 <- c(seq1,i)
      gi <- c(gi,sum(rowSums(Atree$Llist[[length(gridsi)]])>=2))
      i=i+2
      eqs <- c(eqs,1)
      if(length(gi)>1){
        if(gi[length(gi)]<=gi[length(gi)-1]){
          eqs[length(gi)] <- eqs[length(gi)-1]+1
        }
      }
      eq <- eqs[length(eqs)]
    }
    i.best <- seq1[min(which(gi==max(gi)))]
    grids=c(grids,i.best)
    print(xtable(rbind(seq1,gi,eqs),digits=0))
  }
  
  return(list("grids"=grids/snm,"L"=Maxlayers,"M"=Mgroup))
}


# # # # Generate location of hypotheses
nodes.2D.simu <- function(ntips){
  # function to generate the hypotheses locations
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

# # # 
# # # # for M=3, m=1000, n=300
# set.seed(123)
# ntip=1000;n=300
# locs <- nodes.2D.simu(ntip)
# tpars <- find_pars(locs=locs,n=300,ntip=1000)
# tpars 
# # # # put the results: set tpars$grids as the grids argument in A.tree.mult function for aggregation tree construction
# 
# 
# # # # for M=3, m=100, n=90
# set.seed(123)
# ntip=100;n=90
# locs <- nodes.2D.simu(ntip)
# tpars <- find_pars(locs=locs,n=90,ntip=100)
# tpars 
# # # # put the results: set tpars$grids as the grids argument in A.tree.mult function for aggregation tree construction
# # # 
# # # # for M=Inf, m=1000, n=300
# set.seed(123)
# ntip=1000;n=300
# locs <- nodes.2D.simu(ntip)
# tpars <- find_pars(locs=locs,n=300,ntip=1000,M=Inf)
# tpars 
# # # # put the results: set tpars$grids as the grids argument in A.tree.mult function for aggregation tree construction
# # # 
# # # # for M=Inf, m=100, n=90
# set.seed(123)
# ntip=100;n=90
# locs <- nodes.2D.simu(ntip)
# tpars <- find_pars(locs=locs,n=90,ntip=100,M=Inf)
# tpars 
# # # # put the results: set tpars$grids as the grids argument in A.tree.mult function for aggregation tree construction
# # 
# # 
# #setwd("/Documents/Study/project_1/DART_Code_ALL/Data_Analysis/")
