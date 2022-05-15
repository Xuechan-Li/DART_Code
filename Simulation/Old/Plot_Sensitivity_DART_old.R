setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/project_1/DART_Code_last/Simulation/Results")

rm(list=ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library('ggthemes')

load("Out_FDRL_findk.RData")
out_FDRL <- out_FDRL%>%mutate(k1=floor(k))%>%
  mutate(Test=ifelse(k1==k,paste0("FDRL I, k=", k1),
                     paste0("FDRL II, k=", k1)))%>%
  group_by(Test,SE)%>%
  summarize(FDP=mean(FDP),Sensitivity=mean(Sensitivity))

load("Out_DART.RData")
load("Out_FDRL.RData")
out_FDRL <- out_FDRL%>%mutate(k1=floor(k))%>%
  mutate(Test=ifelse(k1==k,paste0("FDRL I, k=", k1),
                     paste0("FDRL II, k=", k1)))
out_DART0 <- out_DART
out <- out_DART0%>%#filter(Layers%in%c())%>%
  bind_rows(out_FDRL)%>%
  pivot_longer(cols=c("Sensitivity","FDP"))%>%
  mutate(nullprop=as.character(nullprop))
dummy <- out%>%filter(name=="FDP")%>%mutate(intercept=0,slop=1)

p1=ggplot(out%>%filter(Test%in%c("3 Layers"
                                 ,"1 Layer"
                              ,"FDRL I, k=7"
                              ,"FDRL II, k=7"),
                       alpha==0.05,
                       name=="FDP"), 
          aes(x=Test, y=value,fill=nullprop))+
  facet_wrap(~SE,ncol=1)+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_hline(yintercept=0.05,linetype="dashed")+
  
  geom_point(aes(col=nullprop,shape=Test))+
  geom_abline(data=dummy,aes(intercept = intercept, slope = slop))+
  theme_tufte()+ theme(axis.line=element_line())

ggplot(out%>%filter(Test%in%c("3 Layers"
                              ,"1 Layer"
                              ),
                    nullprop==0), aes(x=alpha, y=value,group=Test))+
  facet_grid(name~SE,scales="free_y")+
  geom_point(aes(col=nullprop,shape=Test))+
  geom_abline(data=dummy,aes(intercept = intercept, slope = slop))+
  theme_tufte()+ theme(axis.line=element_line())

p1=ggplot(out%>%filter(Test%in%c("4 Layers"
                              #,"3 Layers"
                              #,"FDRL I, k=7"
                              ,"FDRL II, k=5"
                              ),nullprop==0), aes(x=alpha, y=value,group=Test))+
  facet_grid(name~SE,scales="free_y")+
  geom_point(aes(col=nullprop,shape=Test))+
  geom_abline(data=dummy,aes(intercept = intercept, slope = slop))+
  theme_tufte()+ theme(axis.line=element_line())

pdf("Sensitivity.pdf")
print(p1)
dev.off()
  #+scale_shape_manual(values=shape.value,labels=names(table(out$Test)))
  #+scale_color_manual(values = col.value,labels=names(table(out$Test)))

#+labs(color = "Method (n=500,m=1000)",shape="Method (n=500,m=1000)")
#+xlim(0.05,0.22)+ylim(0,0.4)+theme_classic()
#+theme(legend.position = "top",text = element_text(size=100)))
