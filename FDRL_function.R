
####### Those functions are written based on 
##### Zhang, Chunming, Jianqing Fan, and Tao Yu. "Multiple testing via FDRl for large scale imaging data." Annals of statistics 39.1 (2011): 613.

p.med <- function(x,T1p,k){
  median(T1p[rank(x)<=k])
}

compare1 <- function(t,pmed){sum(pmed>=(1-t))}
compare2 <- function(t,pmed){sum(pmed>t)}

find.Q <- function(x,k,pn1,hatn1,j,pmed,T1p){
  nb <- x
  Aj=sample(1:hatn1,j)
  Dj=sample(1:k,k-j)
  
  Dj <- T1p[nb[Dj]]
  Aj <- T1p[pmed <= pn1][Aj]
  
  return(median(c(Dj,Aj)))
}

find.neighboor <- function(x,k){
  which(rank(x,ties.method="first")<=k)
  #sum(pmed[nei]<=pn1)
}

v0j <- function(x,pn1,pmed){
  sum(pmed[x]<=pn1)
}


hatG1 <- function(t,pmed){
  d=2*sum(pmed>0.5)+sum(pmed==0.5)
  n1 <- rep(0,length(t))
  if(sum(t<=0.5)>0){
    n1[t<=0.5] <- sapply(t[t<=0.5],compare1,pmed=pmed)/d
  }
  if(sum(t>0.5)>0){
    n1[t>0.5] <- 1-sapply(t[t>0.5],compare2,pmed=pmed)/d
  }
  return(n1)
}


find.tl1 <- function(pmed,alpha,lambda=0.1){
  t <- seq(0,1,0.0001)
  Wlambda <- sum(pmed>lambda)
  Rt <- length(pmed)-sapply(t,compare2,pmed=pmed)
  fdrl1 <- Wlambda*hatG1(t,pmed=pmed)/(pmax(Rt,1)*(1-hatG1(lambda,pmed=pmed)))
  
  return(list(tmin=t[max(which(fdrl1 <= alpha))],fdrl1=fdrl1))
}



hatG2 <- function(t,pmed,lambda=0.1,k=3,Dist0,T1p){
  hatn0 <- min(floor(sum(pmed>lambda)/(1-hatG1(lambda,pmed))),length(pmed)-1)
  hatn1 <- length(pmed)-hatn0
  
  pn1 <- pmed[rank(pmed,ties.method="first")==hatn1]
  
  neighboor <- apply(Dist0[pmed>pn1,],1,find.neighboor,k=k)
  n1v <- apply(neighboor,2,v0j,pn1=pn1,pmed=pmed)
  n1v2 <- c(n1v,0:(k-1))
  hattheta <- (table(n1v2)-1)/hatn0
  
  qjt <- NULL
  for(j in 1:(k-1)){
    if(j<=hatn1){
      pmed2 <- apply(neighboor[,n1v==0],2,find.Q,k=k,pn1=pn1,hatn1=hatn1,j=j,pmed=pmed,T1p=T1p)
      qjt <- cbind(qjt,
                   (length(pmed2)-sapply(t,compare2,pmed=pmed2))/table(n1v)[1])
    }else{
      qjt <- cbind(qjt,0)
    }
  }
  qjt <- cbind(hatG1(t,pmed),qjt)
  rowSums(qjt%*%diag(hattheta))
}


find.tl2 <- function(pmed,alpha,lambda=0.1,k=3,Dist0,T1p){
  t <- seq(0,1,0.0001)
  Wlambda <- sum(pmed>lambda)
  Rt <- length(pmed)-sapply(t,compare2,pmed=pmed)
  fdrl2 <- Wlambda*hatG2(t,pmed=pmed,k=k,lambda=lambda,Dist0=Dist0,T1p=T1p)/(pmax(Rt,1)*(1-hatG2(lambda,pmed=pmed,k=k,lambda=lambda,Dist0=Dist0,T1p=T1p)))
  index <- which(fdrl2 <= alpha)
  return(list(tmin=t[max(index)],fdrl2=fdrl2))
}


test.fdrl <- function(T1p,Dist0,k=5,lambda=0.1,alpha=0.05){
  pmed <- apply(Dist0,1,p.med,T1p=T1p,k=k)
  if(max(pmed)>=0.5){
    tl1 <- find.tl1(pmed,alpha=alpha)$tmin
    rej <- which(pmed<tl1)
  }else{
    rej <- NULL
  }
  return(rej)
}

mt.fdrl <- function(n=500,Dist0,k=5,lambda=0.1,alpha=0.05){
  nnodes=nrow(Dist0)
  data <- simdata(n=n,nnodes=nnodes,Dist0)
  true.result <- data$beta1
  T1 <- data$stat
  T1p <- pnorm(T1,lower.tail=FALSE)
  
  rej <- test.fdrl(T1p=T1p,Dist0=Dist0,k=k,lambda=lambda,alpha=alpha)
  
  testr3 <- truer <- rep(0,nnodes)
  testr3[rej] <- 1
  truer[true.result!=0] <- 1
  
  fdr <- sum((1-truer)*testr3)/max(sum(testr3),1)  
  power <- sum(truer*testr3)/sum(truer)
  
  return(c(alpha=alpha,k=k,fdr=fdr,power=power))
}


test.fdrl2 <- function(T1p,Dist0,k=5,lambda=0.1,alpha=0.05){
  pmed <- apply(Dist0,1,p.med,T1p=T1p,k=k)
  if(max(pmed)>=0.5){
    see <- find.tl2(pmed=pmed,alpha=alpha,Dist0=Dist0,k=k,T1p=T1p)
    tl2 <- see$tmin
    rej <- which(pmed<tl2)
  }else{
    rej=NULL
  }
  return(rej)
}

mt.fdrl2 <- function(n=500,Dist0,k=5,lambda=0.1,alpha=0.05,SE=1,nullprop=0){
  nnodes=nrow(Dist0)
  data <- simdata(n=n,nnodes=nnodes,Dist0,SE=SE,nullprop=nullprop)
  true.result <- data$beta1
  T1 <- data$stat
  T1p <- pnorm(T1,lower.tail=FALSE)
  
  rej <- test.fdrl(T1p=T1p,Dist0=Dist0,k=k,lambda=lambda,alpha=alpha)
  rej2 <- test.fdrl2(T1p=T1p,Dist0=Dist0,k=k,lambda=lambda,alpha=alpha)
  
  testr3 <- testr4 <- truer <- rep(0,nnodes)
  testr3[rej] <- 1
  testr4[rej2] <- 1
  truer[true.result!=0] <- 1
  
  fdr1 <- sum((1-truer)*testr3)/max(sum(testr3),1)  
  power1 <- sum(truer*testr3)/sum(truer)
  
  fdr2 <- sum((1-truer)*testr4)/max(sum(testr4),1)  
  power2 <- sum(truer*testr4)/sum(truer)
  
  out <-rbind(c(alpha,k,fdr1,power1,1,SE,nullprop),
        c(alpha,k+0.2,fdr2,power2,2,SE,nullprop))
  colnames(out) <- c("alpha","k","fdr","power","type","SE","nullprop")
  #out$SE=SE;out$nullprop=nullprop
  return(out)
}