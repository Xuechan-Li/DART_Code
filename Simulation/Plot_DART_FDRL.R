#setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/Micro_mediation/ALL_AT2D_100_1000m_Normcomb/Plot")

library("ggplot2")
library("ggpubr") 

Model <- c("simple","nlap","T2","abundance","cox")

#Names <- list("Normal","Linear Model","Cox Model","Laplace","Student T")
Names <- list("SE1","SE2","SE3","SE4","SE5")
File <- paste0(Model,"_1000_out.csv")

GG.power <- GG.fdr <- list()

for(i in 1:length(File)){
  out <- read.csv(File[i])
  
  ot <- as.character(out$Test)
  out <- out[which(substr(ot,1,1)%in%c("1","2","3","4") |
                     substr(ot,6,8)%in%c("k=6","k=7")),]
  #out <- out[which(substr(ot,6,8)!="k=7"&substr(ot,1,1)!="5"),]
  out$Test <- ot <- as.character(out$Test)

  out$Test[which(substr(ot,nchar(ot)-1,nchar(ot))==".2")] <- 
    paste0(substr(ot[which(substr(ot,nchar(ot)-1,nchar(ot))==".2")],1,8)," II")
  out$Test[which(substr(ot,nchar(ot)-1,nchar(ot))!=".2" & substr(out$Test,1,1)=="F")] <-
    paste0(out$Test[which(substr(ot,nchar(ot)-1,nchar(ot))!=".2" & substr(out$Test,1,1)=="F")]," I")
  
  
  solidls <- length(table(factor(out$Test[substr(out$Test,1,1)!="F"])))
  dashls <- length(table(factor(out$Test[substr(out$Test,1,1)=="F"])))
  
  out1 <- data.frame(out[,c(2,5,3)])
  out2 <- data.frame(out[,c(2,5,4)])
  colnames(out1) <- c("Desired_FDR","Test","FDP")
  colnames(out2) <- c("Desired_FDR","Test","Sensitivity")
  
  shape.value <- c(15,16,17,18,25,3,7,4,10)
  col.value <- c("black","blue","green","red","purple","black","skyblue2",
                 "goldenrod2","brown4")
  #col.value <- 1:length(table(out$Test))+6
  #col.value[1] <- 6
  #col.value[length(col.value)] <- 5
  
  if(i==1){
    gg.fdr <- (ggplot(out1, aes(x=Desired_FDR, y=FDP,group=Test))+
                 geom_point(aes(col=Test,shape=Test),size=15)+
                 scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                 scale_color_manual(values = col.value,labels=names(table(out$Test)))
               +geom_abline(intercept = 0, slope = 1,size=3)
               +labs(color = "Method (n=300,m=1000)",shape="Method (n=300,m=1000)")
               +xlim(0.05,0.22)+ylim(0,0.4)+theme_classic()+ theme(legend.position = "top",
                                                                   text = element_text(size=100)))
    
    gg.power <- (ggplot(out2, aes(x=Desired_FDR, y=Sensitivity,group=Test),size=3)+
                   geom_point(aes(col=Test,shape=Test),size=15)+
                   geom_line(aes(col=Test,linetype=Test),size=3)+
                   scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                   scale_color_manual(values = col.value,labels=names(table(out$Test)))+
                   scale_linetype_manual(values = c(rep(1,solidls),rep(2,dashls)),
                                         labels = names(table(out$Test)))
                 +xlim(0.05,0.22)+ylim(0,1)+theme_classic()+ theme(legend.position = "none",
                                                                   text = element_text(size=100)))
  }else{
    gg.fdr <- (ggplot(out1, aes(x=Desired_FDR, y=FDP,group=Test))+
                 geom_point(aes(col=Test,shape=Test),size=15)+
                 scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                 scale_color_manual(values = col.value,labels=names(table(out$Test)))
               +geom_abline(intercept = 0, slope = 1,size=3)
               +xlim(0.05,0.22)+ylim(0,0.4)+theme_classic()+ theme(legend.position = "none",axis.title.y=element_blank(),
                                                                   text = element_text(size=100)))
    
    gg.power <- (ggplot(out2, aes(x=Desired_FDR, y=Sensitivity,group=Test),size=3)+
                   geom_point(aes(col=Test,shape=Test),size=15)+
                   geom_line(aes(col=Test,linetype=Test),size=3)+
                   scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                   scale_color_manual(values = col.value,labels=names(table(out$Test)))+
                   scale_linetype_manual(values = c(rep(1,solidls),rep(2,dashls)),
                                         labels = names(table(out$Test)))
                 +xlim(0.05,0.22)+ylim(0,1)+theme_classic()+ theme(legend.position = "none",axis.title.y=element_blank(),
                                                                   text = element_text(size=100)))
    
  }
  GG.fdr[[i]] <- gg.fdr+ xlab(expression(alpha))+ggtitle(Names[i])
  GG.power[[i]] <- gg.power+ xlab(expression(alpha))
}

pdf.options(reset = TRUE, onefile = FALSE)
pdf("ggplot_1000.pdf",width=15*5,height=15*2)
outgg <- ggarrange(GG.fdr[[1]],GG.fdr[[2]],GG.fdr[[3]],GG.fdr[[4]],GG.fdr[[5]],
                   GG.power[[1]],GG.power[[2]],GG.power[[3]],GG.power[[4]],
                   GG.power[[5]],nrow=2,ncol = 5,
                   common.legend = TRUE, legend="bottom",
                   align = "v")   
print(outgg)
dev.off()

outgg1000 <- outgg

GG.fdr1000 <- GG.fdr
GG.power1000 <- GG.power

################################
################################
################################


library("ggplot2")
library("ggpubr") 

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

  solidls <- length(table(factor(out$Test[substr(out$Test,1,1)!="F"])))
  dashls <- length(table(factor(out$Test[substr(out$Test,1,1)=="F"])))
  
  out1 <- data.frame(out[,c(2,5,3)])
  out2 <- data.frame(out[,c(2,5,4)])
  colnames(out1) <- c("Desired_FDR","Test","FDP")
  colnames(out2) <- c("Desired_FDR","Test","Sensitivity")
  
  shape.value <- c(15,16,25,3,4,10)
  col.value <- c("black","blue","darkorange4","red","purple","darkorange4")
  
  if(i==1){
    gg.fdr <- (ggplot(out1, aes(x=Desired_FDR, y=FDP,group=Test))+
                 geom_point(aes(col=Test,shape=Test),size=15)+
                 scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                 scale_color_manual(values = col.value,labels=names(table(out$Test)))
               +geom_abline(intercept = 0, slope = 1,size=3)
               +labs(color = "Method (n=90,m=100)",shape="Method (n=90,m=100)")
               +xlim(0.05,0.22)+ylim(0,0.4)+theme_classic()+ theme(legend.position = "top",
                                                                   text = element_text(size=100)))
    
    gg.power <- (ggplot(out2, aes(x=Desired_FDR, y=Sensitivity,group=Test),size=3)+
                   geom_point(aes(col=Test,shape=Test),size=15)+
                   geom_line(aes(col=Test,linetype=Test),size=3)+
                   scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                   scale_color_manual(values = col.value,labels=names(table(out$Test)))+
                   scale_linetype_manual(values = c(rep(1,solidls),rep(2,dashls)),
                                         labels = names(table(out$Test)))
                 +xlim(0.05,0.22)+ylim(0,1)+theme_classic()+ theme(legend.position = "none",
                                                                   text = element_text(size=100)))
  }else{
    gg.fdr <- (ggplot(out1, aes(x=Desired_FDR, y=FDP,group=Test))+
                 geom_point(aes(col=Test,shape=Test),size=15)+
                 scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                 scale_color_manual(values = col.value,labels=names(table(out$Test)))
               +geom_abline(intercept = 0, slope = 1,size=3)
               +xlim(0.05,0.22)+ylim(0,0.4)+theme_classic()+ theme(legend.position = "none",axis.title.y=element_blank(),
                                                                   text = element_text(size=100)))
    
    gg.power <- (ggplot(out2, aes(x=Desired_FDR, y=Sensitivity,group=Test),size=3)+
                   geom_point(aes(col=Test,shape=Test),size=15)+
                   geom_line(aes(col=Test,linetype=Test),size=3)+
                   scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                   scale_color_manual(values = col.value,labels=names(table(out$Test)))+
                   scale_linetype_manual(values = c(rep(1,solidls),rep(2,dashls)),
                                         labels = names(table(out$Test)))
                 +xlim(0.05,0.22)+ylim(0,1)+theme_classic()+ theme(legend.position = "none",axis.title.y=element_blank(),
                                                                   text = element_text(size=100)))
    
  }
  GG.fdr[[i]] <- gg.fdr+ xlab(expression(alpha))+ggtitle(Names[i])
  GG.power[[i]] <- gg.power+ xlab(expression(alpha))
}

pdf.options(reset = TRUE, onefile = FALSE)
pdf("ggplot_100.pdf",width=15*5,height=15*2)
outgg <- ggarrange(GG.fdr[[1]],GG.fdr[[2]],GG.fdr[[3]],GG.fdr[[4]],GG.fdr[[5]],
                   GG.power[[1]],GG.power[[2]],GG.power[[3]],GG.power[[4]],
                   GG.power[[5]],nrow=2,ncol = 5,
                   common.legend = TRUE, legend="bottom",
                   align = "v")   
print(outgg)
dev.off()




outgg100 <- outgg


pdf.options(reset = TRUE, onefile = FALSE)
pdf("ggplot_ALL.pdf",width=15*5,height=15*4)
outgg <- ggarrange(outgg100,outgg1000,nrow=2,ncol = 1,
                   common.legend = FALSE,
                   align = "v")   
print(outgg)
dev.off()