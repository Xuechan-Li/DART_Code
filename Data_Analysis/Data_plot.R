rm(list=ls())

setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/project_1/DART_Code_last/Data_Analysis/Results")
library("ggplot2")
library("ggdendro")
library(dendextend)
library(gridExtra)
library("reshape2")
library("grid")
library(reshape2)
library(scales)
library(tidyverse)
library(dplyr)
library(data.table)
library(xtable)
library(ggforce)
library(ggh4x)

#source("../Find_grids.R")
#source("../DART_function.R")
#load(file="Dist0_new.RData")
#load(file="clintest_new.RData")
load("outdf_ANEW.RData")
load("/Users/xuechan/Documents/Study/Rotation-Jichun/Micro_data/DART_Selected_seq/taxtab_rdp_rc.RData")

# genlecture <- data.table("Genus"=c("Lactobacillus","Streptococcus",
#   "Escherichia/Shigella","Pediococcus",
#   "Enterococcus","Staphylococcus"),"Literature"="Yes")

Out <- Out%>%right_join(asvname,by=c("ASV"="nodes.index"))%>%
  left_join(rownames_to_column(data.frame(rdptax_rc)),by=c("seq"="rowname"))%>%
  mutate(ASV=ifelse(is.na(detect),"ASVref",ASV),
         BH=ifelse(detect==1,1,0),
         DART=ifelse(detect>=1,1,0))

check <- Out%>%filter(DART+BH>0)%>%
  group_by(Genus,Family,Order,Class,Phylum)%>%
  summarize(BH=sum(BH),DART=sum(DART))%>%
  ungroup()%>%
  #inner_join(genlecture)%>%
  arrange(Order,Family, Genus)%>%
  dplyr::select(c("Phylum","Class","Order","Family","Genus","BH","DART"))

check2 <- Out%>%filter(DART+BH>0)%>%
  dplyr::select(c("Order","Family","Genus","BH","DART","pvals"))

print(xtable(check,digits=0),include.rownames=FALSE)
#Bacteroides https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7934839/
#Flavonifractor https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7218439/

outplot <- Out%>%filter(DART==1)%>%
  dplyr::select(-c("seq","BH","DART"))%>%
  pivot_longer(col=c("Time","Time_lower","Time_upper",
                     "Group","Group_lower","Group_upper",
                     "Interact","Interact_lower","Interact_upper",
                     "Sum","Sum_lower","Sum_upper"))%>%
  mutate(effect = colsplit(name,"_",c(1:2))[,1],
         stats = colsplit(name,"_",c(1:2))[,2],
         stats=ifelse(stats=="","value",stats))%>%
  dplyr::select(-c("name"))%>%
  pivot_wider(names_from=stats,values_from=c(value))%>%
  mutate(effect=ifelse(effect=="Interact","Group*Time",effect),
         Genus=ifelse(Genus=="Clostridium_XlVa","Clostridium",Genus))
#%>%
#  mutate(par=ifelse(effect=="Group","hat(theta)[1*','*i]",
#                    ifelse(effect=="Group*Time","hat(theta)[3*','*i]",
#                           "hat(theta)[2*','*i]")))



my_breaks <- function(y){ 
  if(max(y) >0.0004){
    c(-0.0004,0,0.0004)
  }else{
    c(-0.0002,0,0.0002)
  }
}

outplot2 <- outplot %>% 
  replace_na(list(Order="'-'",Family="'-'",Genus="'-'"))%>%
  filter(effect%in%c("Time","Sum"))%>%
  mutate(effect=ifelse(effect=="Time","SH","HC"))%>%
  mutate(Family=ifelse(Family=="Clostridiaceae_1","Clostridiaceae",Family),
         Genus=ifelse(Genus=="Clostridium_sensu_stricto","CSS",Genus))%>%
  mutate(Family=factor(Family,level=rev(sort(unique(Family)))))


pdf("Forest_plot.pdf",width=8, height=2.5)
ggp=ggplot(outplot2%>%filter(Genus!="'-'"), aes(x = ASV, y = value))+
  facet_nested(effect~Order+Family+Genus, scales = "free",labeller = label_parsed)+
  geom_hline(yintercept=0, linetype="dashed",color="grey")+
  geom_pointrange(aes(ymax = lower, ymin = upper),fatten=0.006,size=0.5)+
  force_panelsizes(cols = c(1.7,1,1,1.2,0.9,2,1.1,1.4)*2,respect = TRUE)+
  ylab("")+
  scale_y_continuous(breaks = my_breaks)+
  labs(tag = "Order:\n\nFamily:\n\nGenus:") +
  theme_classic()+
  theme(plot.tag.position = c(0.04, 0.82),
        plot.tag=element_text(size=7.5,face="bold"),
        strip.text.x = element_text(size=7),
        strip.text.y = element_text(size=7),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4), 
        strip.background = element_rect(color = "black", size = 0.4))
print(ggp)
dev.off()


pdf("Forest_plot_full.pdf",width=10, height=2.5)
ggp=ggplot(outplot2, aes(x = ASV, y = value))+
  facet_nested(effect~Order+Family+Genus, scales = "free",labeller = label_parsed)+
  geom_hline(yintercept=0, linetype="dashed",color="grey")+
  geom_pointrange(aes(ymax = lower, ymin = upper),fatten=0.006,size=0.5)+
  force_panelsizes(cols = c(0.6,1.7,1,1,1.2,0.9,1.3,2,1.1,0.2,1.4)*2,respect = TRUE)+
  ylab("")+
  scale_y_continuous(breaks = my_breaks)+
  labs(tag = "Order:\n\nFamily:\n\nGenus:") +
  #labs(tag = "Order:\n\nFamily:\n\nGenus:\n\n\n\n\n\n\n\n\n\n\n\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t #: CSS=Clostridium sensu stricto") +
  theme_classic()+
  theme(plot.tag.position = c(0.04, 0.82),
        plot.tag=element_text(size=7.5,face="bold"),
        strip.text.x = element_text(size=7),
        strip.text.y = element_text(size=7),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4), 
        strip.background = element_rect(color = "black", size = 0.4))
print(ggp)
dev.off()

# 
# pdf("Forest_plot.pdf",width=8, height=3)
# ggp=ggplot(outplot2, aes(x = ASV, y = value))+
# 
#   facet_nested(par~Order+Family+Genus, scales = "free",labeller = label_parsed)+
#   geom_hline(yintercept=0, linetype="dashed")+
#   geom_pointrange(aes(ymax = lower, ymin = upper))+
#   force_panelsizes(cols = c(0.6,1.2,1.2,1.7,1.5,0.3,1.4), 
#                    respect = TRUE)+
#   #scale_y_continuous(breaks = my_breaks)+
#   ylab("")+
#   labs(tag = "Order:\n\nFamily:\n\nGenus:") +
#   theme(plot.tag.position = c(0.04, 0.87),
#         plot.tag=element_text(size=8,face="bold"),
#         strip.text.x = element_text(size=8),
#         strip.text.y = element_text(size=8),
#         axis.text.x = element_text(size=8,angle=45,hjust=1),
#         axis.title.x=element_blank(),
#         axis.title.y=element_text(size=9),
#         panel.spacing=unit(0, "lines"),
#         panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
#         strip.background = element_rect(color = "black", size = 0.1))
# print(ggp)
# dev.off()
