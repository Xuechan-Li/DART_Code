#setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/project_1/DART_Code_last/Data_Analysis")
rm(list=ls())
library(openxlsx)
library(dplyr)
library(tidyverse)
library(data.table)
library(vegan)
library(ape)
library(cubature)
library(lme4)

#source("../Find_grids.R")
source("../DART_function.R")
source("../FDRL_function.R")
load("taxtab_rdp_rc.RData")

load(file="clintest_new.RData")
load(file="Dist0_new.RData")

Lmax=ceiling(log(nrow(Dist0)/80,3))
###### construct aggregation tree
set.seed(2424)
Dist <- Dist0
diag(Dist) <- 1000
md <- max(apply(Dist,1,min))
#n= uniqueN(clintest$patid)
n= uniqueN(clintest$sampid)
m <- length(nodes.index)
snm <- sqrt(n*log(m)*log(log(m)))

# t1 <- Sys.time()
# Atree <- A.tree.mult(Dist0=Dist0,grids=c(16,16)/snm,Mgroup=3)
# t2 <- Sys.time()
# t2-t1
# save(Atree,file="Pars.RData")

load("Pars.RData")

table(rowSums(Atree$Llist[[1]]))


# ASVs <- intersect(paste0("ASV",1:3000),colnames(clintest))
# clintest2 <- clintest %>% arrange(patid,time)%>%
#   filter(time<100)%>%
#   group_by(patid)%>% 
#   filter(row_number() %in% c(1,n()))%>%
#   mutate(time=diff(range(time)))%>%
#   #mutate(time=log(time))%>%
#   dplyr::select(c("patid","time","inout",ASVs))
# 
# clintest2a <- clintest2 %>% filter(row_number()==1) %>%
#   ungroup()%>%dplyr::select(-c("patid","time","inout"))
# clintest2b <- clintest2 %>% filter(row_number()==2) %>% 
#   ungroup()%>%dplyr::select(-c("patid","time","inout"))
# 
# ctasv <- clintest2b-clintest2a
# clintest2 <- cbind(clintest2%>%ungroup()%>%dplyr::select(c("patid","time","inout")),ctasv)

# Get p-value
load(file="clintest_new.RData")
set.seed(123)
Out <- NULL
clintest <- clintest%>%filter(time<=365)%>%
  mutate(inout=ifelse(inout=="Home Care",1,0))
  #mutate(time=log(time+1))
library(doParallel)
library(foreach)
cl=makeCluster(50)
registerDoParallel(cl)
t1 <- Sys.time()
Out <- foreach(i=nodes.index,.combine=rbind,.packages=c("lme4"))%dopar%
  {
    form0 <- reformulate(c("time","inout","(1|patid)"),response=i)
    form1 <- reformulate(c("time","inout","time*inout","(1|patid)"),response=i)
    fm0 <- lmer(form0, data=clintest)
    fm1 <- lmer(form1, data=clintest)
    cis <- confint(fm1)
    
    #sum of covariate confidence interval
    
    coefs <- summary(fm1)$coef
    cors <- vcov(fm1)[c(2,4),c(2,4)]
    
    coef23 <- coefs[2,1]+coefs[4,1]
    ci23 <- coef23+c(-1.96,1.96)*sqrt(sum(cors))

    c(coefs[2,1],cis[4,],
      coefs[3,1],cis[5,],
      coefs[4,1],cis[6,],
      coef23,ci23,
              anova(fm0,fm1)[2,8])
  }
stopCluster(cl)
t2 <- Sys.time()

# for(i in nodes.index){
#   form0 <- reformulate(c("time","(1|patid)"),response=i)
#   form1 <- reformulate(c("time*inout","(1|patid)"),response=i)
#   fm0 <- lmer(form0, data=clintest)
#   fm1 <- lmer(form1, data=clintest)
#   
#   outs <- c(summary(fm1)$coef[3,1],confint(fm1)[5,],
#            summary(fm1)$coef[4,1],confint(fm1)[6,],
#            anova(fm0,fm1)[2,8])
#   Out <- rbind(Out,outs)
#   }
# t2 <- Sys.time()

colnames(Out) <- c("Time","Time_lower","Time_upper","Group","Group_lower","Group_upper",
  "Interact","Interact_lower","Interact_upper",
  "Sum","Sum_lower","Sum_upper",
  "Pvals")
Out <- data.frame(Out)
#######################################################
######################## DART #########################
#######################################################
alpha=0.05 # set up desired FDR level

Rejs <- test.mult(Llist=Atree$Llist,Dist0=Dist0,
                  T1=t(qnorm(Out$Pvals,lower.tail=FALSE)),alpha=alpha)

detect <- rep(0,length(nodes.index))
#Lmax=rev(seq(length(Atree$Llist)+1))
Lmax=2
for(ell in rev(seq(Lmax))){
  detect[Rejs[[ell]]] <- ell
}
table(detect)

#write the testing result to a csv file
outdf <- data.frame("ASV"=nodes.index,"pvals"=Out$Pvals,"detect"=detect)
Out <- Out %>% mutate(ASV=nodes.index,"detect"=detect,pvals=Pvals)



save(Out,Atree,asvname,file="Results/outdf_ANEW.RData")


#1. 画图
#2. write the data analysis section




# ####################################################
# ############### Distance heatmap ###################
# ####################################################
# 
# # For distance heatmap
# Dist0df <- data.frame(Dist0)
# 
# Dist0df[data$Y,] <- Dist0df
# Dist0df[,data$Y] <- Dist0df
# rownames(Dist0df) <- colnames(Dist0df) <- paste0("ASV",1:96)
# 
# 
# Dist0df <- cbind("ASV"=rownames(Dist0df),"y"=seq(nrow(Dist0)),Dist0df)
# Dist0df <- melt(Dist0df,id=c("ASV","y"))
# Dist0df$variable <- as.character(Dist0df$variable)
# Dist0df$x <- as.integer(substr(Dist0df$variable,4,nchar(Dist0df$variable)))
# 
# # heatmap for distance matrix
# heatmap.dist <- ggplot(data = Dist0df, aes(x = x, y = y)) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_gradient2(name="Distance",low="green",mid="blue",high="white",midpoint=0.5) +
#   theme(axis.text.y = element_text(size = 20),
#         axis.text.x = element_text(size = 20,angle = 45),
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         legend.position = "left",
#         legend.title = element_text(size = 17),
#         legend.text = element_text(size = 15))+
#   scale_x_continuous(expand = c(0, 0))+
#   scale_y_continuous(expand = c(0, 0))
# 
# pdf("see2.pdf")
# heatmap.dist
# dev.off()
# 
