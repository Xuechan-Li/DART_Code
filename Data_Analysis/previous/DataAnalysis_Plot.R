library("ggplot2")
library("ggdendro")
library(dendextend)
library(gridExtra)
library("reshape2")
library("grid")
library(reshape2)
library(scales)

source("../Find_grids.R")
source("../DART_function.R")
load(file="Dist0.RData")

outdf <- read.csv("outdf_ANEW.csv")
outdf2 <- read.csv("outdf_ANEW_ALL.csv")
colnames(outdf2)[1:9]<- c("X","ASV","pvals","BH","DART","FDRL.I.2",
                     "FDRL.I.3","FDRL.II.2","FDRL.II.3")
outdf2 <- outdf2[,1:9]
outdf2 <- outdf2[,c(9,8,7,6,4,5,3,2,1)]
outdf <- outdf[,-1]
outdf[,3] <- factor(outdf[,3])

rownames(Dist0) <- colnames(Dist0)

snm <- sqrt(144*log(96)*log(log(96)))
tpars <- find_pars(n=144,ntip=96,Dist=Dist0,Mgroup=3) # choose the tunning parameter
Atree <- A.tree.mult(grids=tpars$grids,Dist0=Dist0,Mgroup=tpars$M)

dsub=Dist0[which(Dist0==max(Dist0),arr.ind=TRUE)[1,],]
dsub1=which(dsub[1,]<dsub[2,])
dsub2=which(dsub[1,]>dsub[2,])

getseg <- function(llist){
##### function to draw the aggregation tree.
  Ind <- NULL
  Y <- X <- rep(0,ncol(llist))
  Y0 <- 1:ncol(llist)
  Yend <- Xend <- rep(1,ncol(llist))
  for(i in 1:nrow(llist)){
    ind=which(llist[i,]==1)
    if(ind[1] %in% dsub1){
      Y[ind] <- length(Ind)+1:length(ind)
      Ind=c(Ind,ind)
      Yend[ind] <- mean(Y[ind])
    }
  }
  for(i in 1:nrow(llist)){
    ind=which(llist[i,]==1)
    if(!(ind[1] %in% dsub1)){
      Y[ind] <- length(Ind)+1:length(ind)
      Ind=c(Ind,ind)
      Yend[ind] <- mean(Y[ind])
    }
  }
  data.frame(X,Y0,Y,Xend,Yend)
}

data <- cbind(getseg(Atree$Llist[[1]]),outdf)
colnames(data)[c(7,8)] <- c("-log(P-value)","Detection")
data[7] <- -log(data[7])

###### aggregation tree
p <- ggplot() +
  geom_segment(data         =  data,
               aes(x        =  X,
                   y        =  Y,
                   xend     =  Xend,
                   yend     =  Yend),
               lineend      =  "round",
               show.legend  =  FALSE)+labs(title="(a)")+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  geom_point(data=data,aes(x=X,y=Y,size=`-log(P-value)`))+
  scale_size(range=c(0.1,4))+coord_flip()

outdf2[,"X"] <- data$Y
outdf2 <- outdf2[order(outdf2[,"X"]),]
outdflong <- melt(outdf2,id=c("ASV","pvals","X"))
outdflong$Detection <- as.factor(outdflong$value)

###### Heatmap for detection
heatmap.detect <- ggplot(data = outdflong, aes(x = X, y = variable)) +
  geom_tile(aes(fill = Detection)) +
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0),labels=c("FDRL.I.2"=expression(list(FDR[L]~I,~k==2)),"FDRL.II.2"=expression(list(FDR[L]~II,~k==2)),"FDRL.I.3"=expression(list(FDR[L]~I,~k==3)),"FDRL.II.3"=expression(list(FDR[L]~II,~k==3))))+
  scale_fill_manual(name = "Test Result",
                    values = c("0" = "skyblue","1" = "darkorange"),
                    labels = c("Not Rejected","Rejected"))+
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        panel.border = element_blank()
        )

pdf("see.pdf")
print(heatmap.detect)
dev.off()
################# Histogram (from bootstrap) #################

Vresult <- read.csv("Vresult_ANEW.csv")
rownames(Vresult) <- Vresult[,1]
Vresult <- Vresult[,-c(1,9)]
Vresult[,-1] <- Vresult[,-1]
Vresult <- Vresult[,c(1,3,2,4,5,6,7)]
colnames(Vresult)<- c("pvals","DART","BH","FDR[L]~I*','*~k==2",
                      "FDR[L]~I*','*~k==3","FDR[L]~II*','*~k==2",
                      "FDR[L]~II*','*~k==3")
res2 <- melt(Vresult,id=c("pvals"))
res2$variable <- factor(res2$variable,levels=c("DART","BH","FDR[L]~I*','*~k==2",
                                               "FDR[L]~I*','*~k==3","FDR[L]~II*','*~k==2",
                                               "FDR[L]~II*','*~k==3"))
res3 <- res2[res2$value!=0,]
breaks=seq(0,0.8,0.1)

# histogram for rejection rates of ASVs.
h <- ggplot(res2,aes(x=value))+geom_histogram(color="black", fill="white",breaks=seq(0,0.9,0.1), aes(y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]))+
  facet_wrap(~variable,ncol=3,labeller=label_parsed)+ labs(title="(b)")+
  theme(strip.text.x.bottom=element_text(angle=45))+
  xlab("Rejection Rate (RR)")+ylab("Proportion of ASVs")



################Generate Figure 3 in DART manuscript#############
pdf("Dist_AT.pdf", width=8, height=8)
grid.newpage()
print(p,vp = viewport(x = 0.544, y = 0.87, width = 0.955, height = 0.2))
print(heatmap.detect, 
      vp = viewport(x = 0.49, y = 0.68, width = 0.99, height = 0.2))
print(h,vp = viewport(x = 0.49, y = 0.3, width = 0.8, height = 0.6))
dev.off()


################Generate Table 1 in DART manuscript#############
out <- NULL
for(i in 2:7){
  out <- rbind(out,c(sum(Vresult[,i]<=0.1),sum(Vresult[,i]>0.8))/sum(Vresult[,i]>0))
}
out <- round(out,2)
out <- cbind(colnames(Vresult)[-1],out)
colnames(out) <- c("Method","RR<=0.1","RR>0.8")


