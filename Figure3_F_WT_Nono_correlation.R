#!/usr/bin/Rscript --vanilla
setwd("~/github/Nono/")
rm(list=ls())
#=======WT Nono Set1 vs Set2=========
Nono=read.table("Set2_mESC_WT_Nono_ChIP.bin10")
Tet1=read.table("Set2_mESC_WT_Tet1_ChIP.bin10")

dd=merge(Nono,Tet1,by=c("V1","V2","V3"))
data=dd[,c(4,5)]
names(data)=c("Set1","Set2")
data=log10(data[data$Set1>0&data$Set2>0,])
nrow(data)
cor.test(data$Set1,data$Set2)
p=ggplot(data,aes(x=Set1,y=Set2))
p=p+geom_jitter(position="jitter",size=0.1,colour="blue",alpha=1/200)+geom_smooth(method="lm",colour="black",formula=y~x,size=0.5)
p=p+theme_classic(7)+xlab("log10(Set1_Nono)")+ylab("log10(Set2_Nono)")
p=p+annotate("text",color="black",label="R=0.80, P<2.2e-16", x=1.5, y=3.5,hjust = 0.5,size = 3)
p
#ggsave(p,file="Set2_WT_Nono_correlation_bin10.pdf",width=1.7,height=1.6)
