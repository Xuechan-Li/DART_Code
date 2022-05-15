#setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/Micro_mediation/ALL_AT2D_100_1000m_Normcomb/Plot")

library("ggplot2")
library("ggpubr") 
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


Model <- c("simple","nlap","T2","abundance","cox")

#Names <- list("Normal","Linear Model","Cox Model","Laplace","Student T")
Names <- list("SE1","SE2")
File <- paste0(Model,"_100_out.csv")

GG.power <- GG.fdr <- list()

for(i in 1:length(File)){
  out <- read.csv(File[i])
  
  ot <- as.character(out$Test)
  out$Test <- as.character(out$Test)
  out$Test[which(substr(ot,nchar(ot)-1,nchar(ot))==".2")] <- 
    paste0(substr(ot[which(substr(ot,nchar(ot)-1,nchar(ot))==".2")],1,8)," II")
  out$Test[which(substr(ot,nchar(ot)-1,nchar(ot))!=".2" & substr(out$Test,1,1)=="F")] <-
    paste0(out$Test[which(substr(ot,nchar(ot)-1,nchar(ot))!=".2" & substr(out$Test,1,1)=="F")]," I")
  
  out <- out%>%filter(substr(Test,1,1)!="F")%>%mutate(Test=ifelse(Test=="1 Layer","No latent layer",
                                                                  "One latent layer"))

  solidls <- length(table(factor(out$Test)))
  
  out1 <- data.frame(out[,c(2,5,3)])
  out2 <- data.frame(out[,c(2,5,4)])
  colnames(out1) <- c("Desired_FDR","Test","FDP")
  colnames(out2) <- c("Desired_FDR","Test","Sensitivity")
  
  shape.value <- c(15,16)
  col.value <- c("black","red")
  
    gg.fdr <- (ggplot(out1, aes(x=Desired_FDR, y=FDP,group=Test))+
                 geom_point(aes(col=Test,shape=Test),size=15)+
                 scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                 scale_color_manual(values = col.value,labels=names(table(out$Test)))
               +geom_abline(intercept = 0, slope = 1,size=3)
               +ggtitle("(b)")
               #+labs(color = "Method (n=90,m=100)",shape="Method (n=90,m=100)")
               +xlim(0.05,0.22)+ylim(0,0.2)+theme_classic()+ theme(legend.position = "top",
                                                                   text = element_text(size=80)))
    
    gg.power <- (ggplot(out2, aes(x=Desired_FDR, y=Sensitivity,group=Test),size=3)+
                   geom_point(aes(col=Test,shape=Test),size=15)+
                   geom_line(aes(col=Test,linetype=Test),size=3)+
                   scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                   scale_color_manual(values = col.value,labels=names(table(out$Test)))+
                   scale_linetype_manual(values = c(rep(1,solidls)),
                                         labels = names(table(out$Test)))
                 +ggtitle("(c)")
                 +xlim(0.05,0.22)+ylim(0,0.5)+theme_classic()+ theme(legend.position = "top",
                                                                   text = element_text(size=80)))
 
  GG.fdr[[i]] <- gg.fdr+ xlab(expression(alpha))
  GG.power[[i]] <- gg.power+ xlab(expression(alpha))
}


# function to generate the tree
nodes.2D.simu <- function(ntips){
  d1 <- rnorm(ntips,0,2)
  d2 <- runif(ntips,0,4)
  return(cbind(d1,d2))
}

# Generate tree and corresponding aggregation tree
set.seed(123)
ntip=100;n=90                                                                                             
toy <- nodes.2D.simu(ntip)
grids <- 26
snm <- sqrt(n*log(ntip)*log(log(ntip)))
Atree <- A.tree.mult(toy,grids=grids/snm,Mgroup=3)
Dist0 <- Atree$Dist0

b11 <- 0
b12 <- dnorm(Dist0[7,],0,0.1)*1
b13 <- dnorm(Dist0[22,],0,1)*2-0.2
b1 <- (b11+b13+b12)
beta1 <- ifelse(b1>0.15,b1,0)
beta1=beta1*2


locplot <- data.frame(toy,beta1=beta1,Feature=ifelse(beta1==0,"Null","Alternative"))

outdot=ggplot(locplot,aes(x=d1,y=d2,color=Feature))+geom_point(size=10,alpha=0.8)+
  theme_classic()+
  xlab("Dimension 1")+ylab("Dimension 2")+
  ggtitle("(a)")+
  scale_color_manual(values =c("Orange","blue"),labels=c("Alternative","Null"))+
  theme(legend.position = "bottom",text = element_text(size=70))
  
# pdf.options(reset = TRUE, onefile = FALSE)
# pdf("ggplot_100_grant.pdf",width=5,height=5)
# print(outdot)
# dev.off()


pdf.options(reset = TRUE, onefile = FALSE)
pdf("ggplot_100_grant.pdf",width=17*3,height=17)
outgg <- ggarrange(GG.fdr[[1]],GG.power[[1]],nrow=1,ncol = 2,
                   common.legend=TRUE,
                  legend="bottom",align = "v")   
#print(outgg)
print(ggarrange(outdot,outgg,ncol=2,widths=c(1,2)))
dev.off()



