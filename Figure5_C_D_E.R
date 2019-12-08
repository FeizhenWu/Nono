setwd("~/github/Nono")
rm(list=ls())
{
  #================================================================================
  # to retrieve promoter from UCSC
  #================================================================================
  promoterFromUCSC=function(UcscGenes,UpStream=-2000,DownStream=2000-1){
    if(missing(UcscGenes)){
      message("The UcscGenes(RangedData) is required. As default, it will retrieve hg19-refseq from ucsc.")
      UcscGenes=UcscTableQuery(myGenome='hg19',TableName='refGene')
    }
    message("Upstream and Downstream are ", UpStream, " and ", DownStream,", respectively. ")
    genes=as.data.frame(UcscGenes$data.frame)
    Tss=genes$txStart
    idx=genes$strand=="+"
    genes$start[idx]=Tss[idx]+UpStream
    genes$end[idx]=Tss[idx]+DownStream
    Tss=genes$txEnd
    idx=(genes$strand=="-")  
    genes$start[idx]=Tss[idx]-DownStream
    genes$end[idx]=Tss[idx]-UpStream
    
    message("Removing duplicated records ...")
    i=!duplicated(genes$name)
    genes=genes[i,]
    
    genes$size=genes$end-genes$start
    mypromoters=genes[,c("chrom","start","end","name","size","strand","name2")]
    names(mypromoters)=c("chr","start","end","RefseqID","size","strand","GeneSymbol")
    mypromoters=as.data.table(mypromoters)
    setkey(mypromoters, chr, start, end)
    return(mypromoters)
  }
  
  #================================================================================
  # to identify genes targeted by peaks on promoters
  #================================================================================
  Targeted_Genes_by_Peaks_OnPromoter=function(bedfile,promoters){
    if (missing(bedfile)) {
      stop("bedfile is required!")
    }
    if(missing(promoters)){
      message("Use the refseq genes of hg19 in the UCSC.")
      UcscRefseq=ucscTable(GenomeName="hg19",TableName="refGene")
      promoters=promoterFromUCSC(UcscGenes=UcscRefseq)
    }
    
    peaks<- fread(file=bedfile,header=FALSE)[,1:9]
    names(peaks)=c("chr","start","end","name","score","strand","foldchange","Neglog10Pvalue","Neglog10Qvalue")
    peaks=as.data.table(peaks)
    setkey(peaks, chr, start, end)
    
    idx=foverlaps(peaks, promoters, by.x=c("chr", "start", "end"),type="any", which=TRUE)
    idx=na.omit(idx)
    targeted_genes=promoters[idx$yid,]
    return(targeted_genes)
  }
  
  #================================================================================
  # Summarizes data.
  # Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
  #   data: a data frame.
  #   measurevar: the name of a column that contains the variable to be summariezed
  #   groupvars: a vector containing names of columns that contain grouping variables
  #   na.rm: a boolean that indicates whether to ignore NA's
  #   conf.interval: the percent range of the confidence interval (default is 95%)
  #================================================================================
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,conf.interval=.95, .drop=TRUE) 
  {
    require(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }
}
#=Load Data=======================
{
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  OrgDb="org.Mm.eg.db"
  #UcscRefseq=ucscTable(GenomeName="mm9",TrackName="RefSeq Genes",TableName="refGene")
  #save(UcscRefseq,file="mm9_UcscRefseq.RData")
  load("mm9_UcscRefseq.RData")
  promoters=promoterFromUCSC(UcscGenes=UcscRefseq)
  
  load("3NonoKO_D0_vs_3WT_D0.RData")
  
  Tet1_targeted_Genes=Targeted_Genes_by_Peaks_OnPromoter("Set1_Set2_Comm_Tet1_39041.narrowPeak",promoters)
  Tet1_targeted_Genes=as.character(unique(Tet1_targeted_Genes$GeneSymbol))
}

#=Tet1 change======================
{
  dd_Tet1=read.table("WT_NonoKO_Tet1_at_Promoter2K.txt",header = F)[,c(4,7,8)]
  nrow(dd_Tet1)
  names(dd_Tet1)=c("REFSEQ","WT","NonoKO")
  
  tmp=as.data.frame(bitr(dd_Tet1$REFSEQ, fromType="REFSEQ", toType="SYMBOL", OrgDb=OrgDb))
  dd_Tet1=merge(dd_Tet1,tmp,by="REFSEQ")
  dd_Tet1=dd_Tet1[!duplicated(dd_Tet1$SYMBOL),]
  dd_Tet1$group="No-Tet1-Bound"
  dd_Tet1$group[dd_Tet1$SYMBOL %in% Tet1_targeted_Genes]="Tet1-Bound"
  dd_Tet1$log2.fold_change=log2(dd_Tet1$NonoKO/dd_Tet1$WT)
  dd_Tet1=dd_Tet1[is.finite(dd_Tet1$log2.fold_change),]
  dd_Tet1=dd_Tet1[,c("group","log2.fold_change")]
  se=summarySE(dd_Tet1,measurevar = "log2.fold_change",groupvars = "group")
  se$group=factor(se$group,levels = c("Tet1-Bound","No-Tet1-Bound"))
  p=ggplot(se,aes(x=group,y=log2.fold_change,fill=group))+geom_bar(stat = "identity")+theme_classic(8)
  p=p+geom_errorbar(aes(ymin=log2.fold_change-se, ymax=log2.fold_change+se),  size=.3, width=.2, position=position_dodge(.9))+
    ylim(-0.2,0.2)+theme(legend.position="none")+xlab("mm9 Annotation Genes")+ylab("Average log2(NonoKO/WT)")+ggtitle("Tet1 Change")
  Figure5C=p+theme(panel.background = element_rect(fill='white', colour='black'),plot.title = element_text(hjust = 0.5,vjust = 0.5))
  Figure5C
}

#=5hmC change======================
{
  dd_5hmC=read.table("WT_NonoKO_hMeDIP_at_Promoter2K.txt",header = F)[,c(4,7,8)]
  names(dd_5hmC)=c("REFSEQ","WT","NonoKO")
  tmp=as.data.frame(bitr(dd_5hmC$REFSEQ, fromType="REFSEQ", toType="SYMBOL", OrgDb=OrgDb))
  dd_5hmC=merge(dd_5hmC,tmp,by="REFSEQ")
  dd_5hmC=dd_5hmC[!duplicated(dd_5hmC$SYMBOL),]
  
  dd_5hmC$group="No-Tet1-Bound"
  dd_5hmC$group[dd_5hmC$SYMBOL %in% Tet1_targeted_Genes]="Tet1-Bound"
  dd_5hmC$log2.fold_change=log2(dd_5hmC$NonoKO/dd_5hmC$WT)
  
  dd_5hmC=dd_5hmC[is.finite(dd_5hmC$log2.fold_change),]
  dd_5hmC=dd_5hmC[,c("group","log2.fold_change")]
  se=summarySE(dd_5hmC,measurevar = "log2.fold_change",groupvars = "group")
  se$group=factor(se$group,levels = c("Tet1-Bound","No-Tet1-Bound"))
  
  p=ggplot(se,aes(x=group,y=log2.fold_change,fill=group))+geom_bar(stat = "identity")+theme_classic(8)
  p=p+geom_errorbar(aes(ymin=log2.fold_change-se, ymax=log2.fold_change+se),  size=.3, width=.2, position=position_dodge(.9))+
    ylim(-0.5,0.3)+theme(legend.position="none")+xlab("mm9 Annotation Genes")+ylab("Average log2(NonoKO/WT)")+ggtitle("5hmC Change")
  Figure5D=p+theme(panel.background = element_rect(fill='white', colour='black'),plot.title = element_text(hjust = 0.5,vjust = 0.5))
  Figure5D
}

#==Figure 5E ======================
{
  dd=diffGeneEnrich$exp
  nrow(dd)
  dd=dd[is.finite(dd$log2.fold_change),]
  dd$group="No-Tet1-Bound"
  dd$group[dd$geneSymbol %in% Tet1_targeted_Genes]="Tet1-Bound"
  head(dd)
  dd=dd[,c("group","log2.fold_change")]
  se=summarySE(dd,measurevar = "log2.fold_change",groupvars = "group")
  se$group=factor(se$group,levels = c("Tet1-Bound","No-Tet1-Bound"))
  
  p=ggplot(se,aes(x=group,y=log2.fold_change,fill=group))+geom_bar(stat = "identity")+theme_classic(7)
  p=p+geom_errorbar(aes(ymin=log2.fold_change-se, ymax=log2.fold_change+se),  size=.3, width=.2, position=position_dodge(.9))+
    ylim(-0.3,.4)+theme(legend.position="none")+xlab("mm9 Annotation Genes")+ylab("Average log2(NonoKO/WT)")+ggtitle("Differential Expression")
  Figure5E=p+theme(panel.background = element_rect(fill='white', colour='black'),plot.title = element_text(hjust = 0.5,vjust = 0.5))
  Figure5E
}
