setwd("~/github/Nono/")
rm(list=ls())

library(ComplexHeatmap)
library(reshape)
library(ggplot2)
library(scales)
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
  
  DiffNum=function(diffGeneEnrich=diffGeneEnrich){
    diffGenes=diffGeneEnrich$diffgenes
    diffGenes$Reg="up"
    diffGenes$Reg[diffGenes$log2.fold_change<0]="down"
    diffnum=as.data.frame(table(diffGenes$Reg))
    return(diffnum)
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

#==Supplementary Figure 6 A===
load("3NonoKO_D0_vs_3WT_D0.RData")
WT_KO=diffGeneEnrich$diffgenes
WT_KO$Regulation="Up"
WT_KO$Regulation[WT_KO$log2.fold_change<0]="Down"
WT_KO$Tet1Bing="No"
WT_KO$Tet1Bing[WT_KO$gene %in% Tet1_targeted_Genes]="Yes"
SupFigure6A=table(WT_KO[,c("Regulation","Tet1Bing")])
SupFigure6A
chisq.test(SupFigure6A)

#===Supplementary Figure 6B===
BP=diffGeneEnrich$DwGeneEnrich$GOBPenrich@result

