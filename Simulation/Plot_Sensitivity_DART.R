setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/project_1/DART_Code_last/Simulation/Results")

rm(list=ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library('ggthemes')
library(ggh4x)
library(ggpubr)

#### Sensitivity for M

lower <- 0.25
upper <- 0.75

load("Out_DART_diffM.RData")

out_DART <- data.frame(result_DART)%>%
  rename(Sensitivity=Power,FDP=FDR)%>%
  pivot_longer(cols=c("Sensitivity","FDP"))%>%
  group_by(alpha,Layers,SE,M,name) %>% 
  summarize(
    Upper=quantile(value,upper,na.rm=TRUE),
    Lower=quantile(value,lower,na.rm=TRUE),
    value=mean(value,na.rm=TRUE)) %>%
  mutate(Test= ifelse(Layers==1,"BH",paste0(Layers," Layers")))%>%
  mutate(Test=ifelse(Layers==max(Layers),"DART",Test))


out <- out_DART%>%
  group_by(M)%>%
  filter(Layers==max(Layers))%>%
  ungroup()%>%
  mutate(M=as.character(M),aa=paste0("alpha*'='*",alpha))

dummy <- out%>%filter(name=="FDP")%>%mutate(intercept=alpha)
p1=ggplot(out,aes(x=SE, y=value,fill=M))+
  facet_grid(name~aa,scales="free_y",labeller = label_parsed)+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  geom_hline(data=dummy,aes(yintercept=intercept),linetype="dashed")+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        #panel.spacing=unit(0, "lines"),
        #panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 1))

pdf("Sensitivity_M.pdf",height=2.5,width=6)
print(p1)
dev.off()



####### Sensitivity smoothness

load("Out_DART2.RData")
load("Out_FDRL2.RData")

out_FDRL <- data.frame(result_FDRL)%>%
  rename(Sensitivity=power,FDP=fdr)%>%
  pivot_longer(cols=c("Sensitivity","FDP"))%>%
  group_by(alpha,k,SE,nullprop,name) %>% 
  summarize(Upper=quantile(value,upper,na.rm=TRUE),
            Lower=quantile(value,lower,na.rm=TRUE),
            value=mean(value,na.rm=TRUE))%>%
  mutate(k1=floor(k))%>%
  mutate(Test=ifelse(k1==k,"FDRL I","FDRL II"))%>%
  ungroup()

out_DART <- data.frame(result_DART) %>%
  filter(Layers!=2)%>%
  rename(Sensitivity=Power,FDP=FDR)%>%
  pivot_longer(cols=c("Sensitivity","FDP"))%>%
  group_by(alpha,Layers,SE,nullprop,name) %>% 
  summarize(Upper=quantile(value,upper,na.rm=TRUE),
            Lower=quantile(value,lower,na.rm=TRUE),
            value=mean(value,na.rm=TRUE))%>%
  mutate(Test=ifelse(Layers==1,"BH","DART"))%>%
  ungroup()

out <- out_DART%>%bind_rows(out_FDRL)%>%
  mutate(nullprop0=nullprop)%>%
  mutate(nullprop=paste0("tau*'='*",as.character(nullprop)),aa=paste0("alpha*'='*",alpha))


# out_FDRL <- out_FDRL%>%mutate(k1=floor(k))%>%
#   mutate(Test=ifelse(k1==k,"FDRL I","FDRL II"))
# #  mutate(Test=ifelse(k1==k,paste0("FDRL I, k=", k1),
# #                     paste0("FDRL II, k=", k1)))
# out_DART0 <- out_DART
# out <- out_DART0%>%filter(Layers!=2)%>%
#   mutate(Test=ifelse(Layers==1,"BH","DART"))%>%
#   bind_rows(out_FDRL)%>%
#   pivot_longer(cols=c("Sensitivity","FDP"))%>%
#   mutate(nullprop0=nullprop)%>%
#   mutate(nullprop=paste0("tau*'='*",as.character(nullprop)),aa=paste0("alpha*'='*",alpha))

dummy0 <- out%>%filter(name=="FDP",nullprop0==0)%>%mutate(intercept=alpha)
out0 <- out%>%filter(nullprop0==0)
compares <- ggplot(out0,aes(x=SE, y=value,fill=Test))+
  facet_grid(name~aa,labeller = label_parsed,scales="free")+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  geom_hline(data=dummy0,aes(yintercept=intercept),linetype="dashed")+
  theme_classic()+
  labs(title="A")+
  theme(axis.title.y=element_blank())


dummy <- out%>%filter(name=="FDP",nullprop0!=0)%>%mutate(intercept=alpha)
out_smooth <- out%>%filter(Test%in%c("DART"
                                     ,"BH"
                                     ,"FDRL I"
                                     ,"FDRL II"),nullprop0!=0)

s_fdp=ggplot(out_smooth%>%filter(name=="FDP"), 
          aes(x=SE, y=value,fill=Test))+
  facet_grid(nullprop~aa,labeller = label_parsed)+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  scale_y_continuous(breaks=c(0.0,0.2))+
  geom_hline(data=dummy,aes(yintercept=intercept),linetype="dashed")+
  theme_classic()+
  labs(title="B")+
  theme(axis.title.y=element_blank())

s_power=ggplot(out_smooth%>%filter(name=="Sensitivity"), 
          aes(x=SE, y=value,fill=Test))+
  facet_grid(nullprop~aa,labeller = label_parsed)+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  scale_y_continuous(breaks=c(0.0,0.4))+
  theme_classic()+
  theme(axis.title.y=element_blank())

pdf.options(reset = TRUE, onefile = FALSE)
pdf("Compare.pdf",height=6.3,width=5)
print(ggarrange(compares,s_fdp,common.legend=TRUE,ncol=1,heights = c(3, 4.5)))
dev.off()

pdf("nonnull_prop_power.pdf",height=3.5,width=6)
print(s_power)
dev.off()




