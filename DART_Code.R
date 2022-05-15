#setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/project_1")

#setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/project_1/DART")

library(qlcMatrix)  
library(ape)
library(cubature)
library(xtable)

find.c <- function(TT,lTT,clast=10000,alpha=0.05,length.out=10000){
  # function for finding the rejection threshold on layer 1.
  m <- length(TT)
  mall <- sum(lTT)
  uplim <- 2.2*log(m)
  cl <- 0; cr <- min(uplim,2*clast)
  c0=seq(cl,cr,length.out=10000)
  TT <- matrix(TT,m,1)
  TTc0 <- apply(TT,1,'>',c0) #1000*199 matrix
  TTc0 <- TTc0%*%diag(lTT)
  
  estfdr1 <- mall*pchisq(c0,df=1,lower.tail=FALSE)
  estfdr2 <- (pmax(apply(TTc0,1,sum),1))
  
  indsmall <- which(estfdr1/estfdr2<=alpha)
  
  if(length(indsmall)==0){
    out=cr;indx=clast
  }else{
    indx=min(indsmall)
    out=c0[indx]
  }
  return(list(out=out,prevrej=c(estfdr1[indx],estfdr2[indx])))
}
#(pval2,J22[testindex],alpha=alpha,prevrej=prevrejs)

find.p <- function(Tp,lTT,clast=10000,alpha=0.05,length.out=10000,prevrej){
  #### function for finding the rejection threshold on layer \ell/
  m <- length(Tp)
  m0 <- sum(!is.na(Tp))
  mall <- sum(lTT[which(!is.na(Tp))])
  lowlim <- 1/m0^2
  cl <- lowlim; cr <- 1
  c0=seq(cl,cr,length.out=10000)
  Tp <- matrix(Tp,m,1)
  TTc0 <- apply(Tp,1,'<',c0) #1000*199 matrix
  TTc0[,which(is.na(TTc0[1,]))] <- 0
  
  if(is.null(dim(lTT))){
    TTc0 <- TTc0%*%lTT
  }else{
    TTc0 <- TTc0%*%diag(lTT)
  }
  
  estfdr1 <- prevrej[1]+mall*c0
  estfdr2 <- prevrej[2]+(pmax(apply(TTc0,1,sum,na.rm=TRUE),1))
  indsmall <- which(estfdr1/estfdr2<=alpha)
  if(length(indsmall)==0){
    out=cl;indx=1
  }else{
    indx=max(indsmall)
    out=c0[indx]
  }
  return(list(out=out,prevrej=c(estfdr1[indx],estfdr2[indx])))
}

make.group <- function(Dist,grids=0.5,miter=Inf,Mgroup=5){
  #### function for constructing the \ell layer in aggregation tree. Here, \ell= 2,...,L.
  ## arguments:
  # Dist: (matrix) distance matrix.
  # grids: (numeric) vector of maximum distance g^{(\ell)} on layers \ell\in\{2,...,L\}. 
  # miter: (integer) maximum number of iteration for finding the aggregated nodes. Should be greater than the total number of hypothesis m.
  # Mgroup: (integer) maximum cardinality of the child node set of the nodes in the aggregation tree.
  ## output: (list)
  # outmatrix: (matrix) The nodes on layer \ell presented in the form of matrix list.
  # Dist: (matrix) updated distance matrix on layer \ell.
  
  
  outmatrix0=diag(rep(1,nrow(Dist)))
  Dist1 <- Dist0 <- Dist
  # Dist0: original Distant matrix
  # Dist : Distant matrix for output
  # Dist1 : Distant matrix with large value
  infd <- max(Dist)*1000
  diag(Dist1) <- infd
  
  Grij <- NULL
  iter <- 0
  while(min(Dist1)<=grids&iter<=miter){
    maxdistn <- maxdist <- mM <- min(Dist1)
    argM <- which(Dist1 == mM, arr.ind = TRUE)
    gijk <- argM[1,]
    
    newgroupM <- sum(outmatrix0[gijk[1],]+outmatrix0[gijk[2],])
    
    if(newgroupM<=Mgroup){
      outmatrix0 <- rbind(outmatrix0,outmatrix0[gijk[1],]+outmatrix0[gijk[2],])
      outmatrix0 <- outmatrix0[-gijk,]
      
      newmax <- max(Dist[gijk,gijk])
      newrow <- qlcMatrix::rowMax(cbind(diag(Dist)[-gijk],Dist[-gijk,gijk],newmax))@x
      
      Dist <- rbind(cbind(Dist[-gijk,-gijk],newrow),c(newrow,newmax))
      Dist1 <- Dist
      diag(Dist1) <- infd
      
      if(newgroupM==Mgroup){
        fullgroup <- which(rowSums(outmatrix0)==Mgroup)
        Dist1[fullgroup,] <- infd
        Dist1[,fullgroup] <- infd
      }
    }else{
      Dist1[gijk[1],gijk[2]] <- Dist1[gijk[2],gijk[1]] <- infd
    }
    iter=iter+1
  }
  
  return(list(outmatrix=outmatrix0,Dist=Dist))
}


mdist.groups <- function(x1,x2,Dist){
  ####function to measure the distance between two groups.
  ####the distance is the maximum pair-wise distance among all nodes in the two groups.
  ## arguments:
  # x1: (vector) indexes of the nodes in the first group.
  # x2: (vector) indexes of the nodes in the second group.
  # Dist: (matrix) distance matrix.
  ## output: (numeric) distance betwee two groups.
  x <- na.omit(c(x1,x2))
  return(max(Dist[x,x],na.rm=TRUE))
}


A.tree.mult <- function(tree=NULL,grids,Dist0=NULL,Mgroup=3){
  #### function for aggregation tree construction
  ## arguments:
  # tree: (matrix) matrix containing the location information for the features linking to the hypotheses
  # Not needed if distance matrix Dist0 is ready.
  # grids: (vector) vector of maximum distance g^{(\ell)} on layers \ell=2,...,L. The maximum distance on the first layer is always 0. 
  # Dist0: (matrix) distance matrix
  # Mgroup: (integer) maximum cardinality of the child node set of the nodes in the aggregation tree.
  ## output: (list)
  # Llist: (list of matrixs) Aggregation tree presented in the form of matrix list.
  # Dist0: (matrix) distance matrix
  if(is.null(Dist0)){
    Dist0 <- Dist <- as.matrix(dist(tree))
  }else{
    Dist <- Dist0
  }
  m <- nrow(Dist0)
  # Layer 2
  Md <- Llist <- list()
  s=1
  
  # Layer \ell
  while(s<=length(grids)){
    outall <- make.group(Dist,grids=grids[s],Mgroup=Mgroup)
    Llist[[s]] <- outall$outmatrix
    Dist <- outall$Dist
    Dist1 <- Dist
    diag(Dist1) <- 10000*max(Dist)
    Md[[s]] <- max(apply(Dist1,1,min))
    s=s+1
  }
  return(list(Llist=Llist,Dist0=Dist0,Md=Md))
}


test.mult <- function(alpha=0.05,Llist,Dist0,T1){
  #### MAIN function of DART framwork.
  ## arguments:
  # alpha: (numeric) Vector of desired FDR
  # Llist: (list of matrixes) Aggregation tree presented in the form of matrix list.
  # Dist0: (matrix) Distance matrix
  # T1: (vector) Input Z-values transformed from p-values
  ## output: (list of vectors)
  # Index of rejected nodes on each layers
  Rej <- list()
  nnodes <- m <- length(c(T1))
  S=length(Llist)
  TT <- ifelse(T1==Inf,100,ifelse(T1==-Inf,-100,T1))
  
  pval1 <- pnorm(TT,lower.tail=FALSE)
  
  layer1.stats=find.p(pval1,rep(1,m),alpha=alpha,prevrej=c(0,0))
  
  cl1=layer1.stats$out
  prevrejs <- layer1.stats$prevrej
  
  Rej[[1]] <- rejs0 <- which(pval1<cl1)
  
  I2 <- diag(rep(1,m))
  I2[rejs0,rejs0] <- 0
  
  LI2 <- Llist[[1]]%*%I2
  
  # layer2
  J22 <- rowSums(LI2)
  testindex <- which(J22>=2) #index of elements that is in the random aggregation tree
  obsval2 <- LI2%*%t(TT)
  
  pval2 <- pnorm(obsval2[testindex]/sqrt(J22[testindex]),lower.tail=FALSE)
  
  layer2.stats=find.p(pval2,J22[testindex],alpha=alpha,prevrej=prevrejs)
  
  cl2=layer2.stats$out
  prevrejs <- layer2.stats$prevrej
  
  rejindex <- testindex[pval2<cl2]
  
  if(sum(pval2<cl2) %in% c(0,1)){
    rejs1 <- which(LI2[rejindex,]==1)
  }else{
    rejs1 <- which(colSums(LI2[rejindex,])==1)
  } 
  Rej[[2]] <- c(Rej[[1]],rejs1)
  
  rejs2 <- NULL
  s=2
  LI32 <- LI2
  while(s<=S){
    I3 <- rep(1,nrow(Llist[[s-1]]))
    I3[rejindex] <- 0
    I3[which(J22==0)] <- 0
    I3 <- diag(I3)
    
    LI3 <- Llist[[s]]%*%I3
    LI32 <- LI3%*%LI32
    
    J22 <- rowSums(LI3)
    
    testindex <- which(J22>=2) #index of elements that is in the random aggregation tree
    obsval3 <- LI32%*%t(TT)
    
    J32 <- rowSums(LI32)
    
    pval3 <- pnorm(obsval3[testindex]/sqrt(J32[testindex]),lower.tail=FALSE)
    
    layer3.stats=find.p(pval3,J32[testindex],alpha=alpha,prevrej=prevrejs)
    
    cl3 <- layer3.stats$out
    prevrejs <- layer3.stats$prevrej
    rejindex <- testindex[pval3<cl3]
    
    if(sum(pval3<cl3) %in% c(0,1)){
      rejs2 <- which(LI32[rejindex,]==1)
    }else{
      rejs2 <- which(colSums(LI32[rejindex,])==1)
    }
    Rej[[s+1]] <- c(Rej[[s]],rejs2)
    s=s+1
  }
  return(Rej)
}
