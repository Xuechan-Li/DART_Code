#setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/Micro_mediation/ALL_AT2D_100_1000m_Normcomb/Plot")

library("ggplot2")
library("ggpubr") 

Model <- c("simple","nlap","T2","abundance","cox")

#Names <- list("Normal","Linear Model","Cox Model","Laplace","Student T")
Names <- list("SE1","SE2","SE3","SE4","SE5")
File <- paste0(Model,"_1000_out.csv")
Fileinf <- paste0("INF_",File)
GG.power <- GG.fdr <- list()

for(i in 1:length(File)){
  out <- read.csv(File[i])
  outinf <- read.csv(Fileinf[i])
  
  ot <- as.character(out$Test)
  out <- out[which(substr(ot,1,1)=="4"),]
  otinf <- as.character(outinf$Test)
  outinf <- outinf[which(substr(otinf,1,1)=="4"),]
  
  out$Test <- paste0(as.character(out$Test),", M=5")
  outinf$Test <- paste0(as.character(outinf$Test),", M=Inf")
  out <- rbind(out,outinf)
  
  out1 <- data.frame(out[,c(2,5,3)])
  out2 <- data.frame(out[,c(2,5,4)])
  colnames(out1) <- c("Desired_FDR","Test","FDP")
  colnames(out2) <- c("Desired_FDR","Test","Sensitivity")
  
  shape.value <- c(15,16)
  col.value <- c("blue","goldenrod2")
  
  if(i==1){
    gg.fdr <- (ggplot(out1, aes(x=Desired_FDR, y=FDP,group=Test))+
                 geom_point(aes(col=Test,shape=Test),size=15)+
                 scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                 scale_color_manual(values = col.value,labels=names(table(out$Test)))
               +geom_abline(intercept = 0, slope = 1,size=3)
               +labs(color = "Method (n=500,m=1000)",shape="Method (n=500,m=1000)")
               +xlim(0.05,0.22)+ylim(0,0.4)+theme_classic()+ theme(legend.position = "top",
                                                                   text = element_text(size=100)))
    
    gg.power <- (ggplot(out2, aes(x=Desired_FDR, y=Sensitivity,group=Test),size=3)+
                   geom_point(aes(col=Test,shape=Test),size=15)+
                   geom_line(aes(col=Test,linetype=Test),size=3)+
                   scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                   scale_color_manual(values = col.value,labels=names(table(out$Test)))+
                   scale_linetype_manual(values = c(1,2),
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
                   scale_linetype_manual(values = c(1,2),
                                         labels = names(table(out$Test)))
                 +xlim(0.05,0.22)+ylim(0,1)+theme_classic()+ theme(legend.position = "none",axis.title.y=element_blank(),
                                                                   text = element_text(size=100)))
    
  }
  GG.fdr[[i]] <- gg.fdr+ xlab(expression(alpha))+ggtitle(Names[i])
  GG.power[[i]] <- gg.power+ xlab(expression(alpha))
}

pdf.options(reset = TRUE, onefile = FALSE)
pdf("ggplotM_1000.pdf",width=15*5,height=15*2)
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
Fileinf <- paste0("INF_",File)

GG.power <- GG.fdr <- list()

for(i in 1:length(File)){
  out <- read.csv(File[i])
  outinf <- read.csv(Fileinf[i])
  
  ot <- as.character(out$Test)
  out <- out[which(substr(ot,1,1)=="2"),]
  otinf <- as.character(outinf$Test)
  outinf <- outinf[which(substr(otinf,1,1)=="2"),]
  
  out$Test <- paste0(as.character(out$Test),", M=3")
  outinf$Test <- paste0(as.character(outinf$Test),", M=Inf")
  out <- rbind(out,outinf)
  
  out1 <- data.frame(out[,c(2,5,3)])
  out2 <- data.frame(out[,c(2,5,4)])
  colnames(out1) <- c("Desired_FDR","Test","FDP")
  colnames(out2) <- c("Desired_FDR","Test","Sensitivity")
  
  shape.value <- c(15,16)
  col.value <- c("blue","goldenrod2")
  
  if(i==1){
    gg.fdr <- (ggplot(out1, aes(x=Desired_FDR, y=FDP,group=Test))+
                 geom_point(aes(col=Test,shape=Test),size=15)+
                 scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                 scale_color_manual(values = col.value,labels=names(table(out$Test)))
               +geom_abline(intercept = 0, slope = 1,size=3)
               +labs(color = "Method (n=500,m=1000)",shape="Method (n=500,m=1000)")
               +xlim(0.05,0.22)+ylim(0,0.4)+theme_classic()+ theme(legend.position = "top",
                                                                   text = element_text(size=100)))
    
    gg.power <- (ggplot(out2, aes(x=Desired_FDR, y=Sensitivity,group=Test),size=3)+
                   geom_point(aes(col=Test,shape=Test),size=15)+
                   geom_line(aes(col=Test,linetype=Test),size=3)+
                   scale_shape_manual(values=shape.value,labels=names(table(out$Test)))+
                   scale_color_manual(values = col.value,labels=names(table(out$Test)))+
                   scale_linetype_manual(values = c(1,2),
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
                   scale_linetype_manual(values = c(1,2),
                                         labels = names(table(out$Test)))
                 +xlim(0.05,0.22)+ylim(0,1)+theme_classic()+ theme(legend.position = "none",axis.title.y=element_blank(),
                                                                   text = element_text(size=100)))
    
  }
  GG.fdr[[i]] <- gg.fdr+ xlab(expression(alpha))+ggtitle(Names[i])
  GG.power[[i]] <- gg.power+ xlab(expression(alpha))
}

pdf.options(reset = TRUE, onefile = FALSE)
pdf("ggplotM_100.pdf",width=15*5,height=15*2)
outgg <- ggarrange(GG.fdr[[1]],GG.fdr[[2]],GG.fdr[[3]],GG.fdr[[4]],GG.fdr[[5]],
                   GG.power[[1]],GG.power[[2]],GG.power[[3]],GG.power[[4]],
                   GG.power[[5]],nrow=2,ncol = 5,
                   common.legend = TRUE, legend="bottom",
                   align = "v")   
print(outgg)
dev.off()




outgg100 <- outgg


pdf.options(reset = TRUE, onefile = FALSE)
pdf("ggplotM_ALL.pdf",width=15*5,height=15*4)
outgg <- ggarrange(outgg100,outgg1000,nrow=2,ncol = 1,
                   common.legend = FALSE,
                   align = "v")   
print(outgg)
dev.off()
