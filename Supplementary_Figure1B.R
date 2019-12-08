#!/usr/bin/Rscript --vanilla
setwd("~/github/Nono/")
rm(list=ls())
library(ggplot2)
library(reshape)
library(scales)
#==Set1===============================
Set1=read.table("set1_Expression_TPM.xls",sep="\t",header = T)
Set1[,2:13]=log2(Set1[,2:13]+1)
#D0====
cor.test(Set1$set1_RE_D0,Set1$set1_WT_D0)
set_d0=Set1[,c("set1_WT_D0","set1_RE_D0")]
p=ggplot(set_d0,aes(x=set1_WT_D0,y=set1_RE_D0))
p=p+geom_jitter(position="jitter",size=0.3,colour="blue",alpha=1/1)+geom_smooth(method="lm",colour="black",formula=y~x,size=0.8)
p=p+theme_classic(7)+xlab("log2(D0_WT TPM)")+ylab("log2(D0_NonoKO+WT TPM")
p
#D3====
cor.test(Set1$set1_RE_D3,Set1$set1_WT_D3)
set_D3=Set1[,c("set1_WT_D3","set1_RE_D3")]
p=ggplot(set_D3,aes(x=set1_WT_D3,y=set1_RE_D3))
p=p+geom_jitter(position="jitter",size=0.3,colour="blue",alpha=1/1)+geom_smooth(method="lm",colour="black",formula=y~x,size=0.8)
p=p+theme_classic(7)+xlab("log2(D3_WT TPM)")+ylab("log2(D3_NonoKO+WT TPM")
p
#D6====
cor.test(Set1$set1_RE_D6,Set1$set1_WT_D6)
set_D6=Set1[,c("set1_WT_D6","set1_RE_D6")]
p=ggplot(set_D6,aes(x=set1_WT_D6,y=set1_RE_D6))
p=p+geom_jitter(position="jitter",size=0.3,colour="blue",alpha=1/1)+geom_smooth(method="lm",colour="black",formula=y~x,size=0.8)
p=p+theme_classic(7)+xlab("log2(D6_WT TPM)")+ylab("log2(D6_NonoKO+WT TPM")
p
#D12====
cor.test(Set1$set1_RE_D12,Set1$set1_WT_D12)
set_D12=Set1[,c("set1_WT_D12","set1_RE_D12")]
p=ggplot(set_D12,aes(x=set1_WT_D12,y=set1_RE_D12))
p=p+geom_jitter(position="jitter",size=0.3,colour="blue",alpha=1/1)+geom_smooth(method="lm",colour="black",formula=y~x,size=0.8)
p=p+theme_classic(7)+xlab("log2(D12_WT TPM)")+ylab("log2(D12_NonoKO+WT TPM")
p

#==set2===============================
set2=read.table("set2_Expression_TPM.xls",sep="\t",header = T)
set2[,2:13]=log2(set2[,2:13]+1)
#D0====
cor.test(set2$set2_RE_D0,set2$set2_WT_D0)
set_d0=set2[,c("set2_WT_D0","set2_RE_D0")]
p=ggplot(set_d0,aes(x=set2_WT_D0,y=set2_RE_D0))
p=p+geom_jitter(position="jitter",size=0.3,colour="blue",alpha=1/1)+geom_smooth(method="lm",colour="black",formula=y~x,size=0.8)
p=p+theme_classic(7)+xlab("log2(D0_WT TPM)")+ylab("log2(D0_NonoKO+WT TPM")
p

#D3====
cor.test(set2$set2_RE_D3,set2$set2_WT_D3)
set_D3=set2[,c("set2_WT_D3","set2_RE_D3")]
p=ggplot(set_D3,aes(x=set2_WT_D3,y=set2_RE_D3))
p=p+geom_jitter(position="jitter",size=0.3,colour="blue",alpha=1/1)+geom_smooth(method="lm",colour="black",formula=y~x,size=0.8)
p=p+theme_classic(7)+xlab("log2(D3_WT TPM)")+ylab("log2(D3_NonoKO+WT TPM")
p

#D6====
cor.test(set2$set2_RE_D6,set2$set2_WT_D6)
set_D6=set2[,c("set2_WT_D6","set2_RE_D6")]
p=ggplot(set_D6,aes(x=set2_WT_D6,y=set2_RE_D6))
p=p+geom_jitter(position="jitter",size=0.3,colour="blue",alpha=1/1)+geom_smooth(method="lm",colour="black",formula=y~x,size=0.8)
p=p+theme_classic(7)+xlab("log2(D6_WT TPM)")+ylab("log2(D6_NonoKO+WT TPM")
p

#D12====
cor.test(set2$set2_RE_D12,set2$set2_WT_D12)
set_D12=set2[,c("set2_WT_D12","set2_RE_D12")]
p=ggplot(set_D12,aes(x=set2_WT_D12,y=set2_RE_D12))
p=p+geom_jitter(position="jitter",size=0.3,colour="blue",alpha=1/1)+geom_smooth(method="lm",colour="black",formula=y~x,size=0.8)
p=p+theme_classic(7)+xlab("log2(D12_WT TPM)")+ylab("log2(D12_NonoKO+WT TPM")
p

