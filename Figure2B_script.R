rm(list=ls())
#function
Nine_Square=function(diffGeneEnrich1,diffGeneEnrich2,myXlabel="",myYlabel="",org="hsa",
                     xMin=-4,xMax=4,yMin=-4,yMax=4,
                     cutoff=1.5,Fpkm_filter=F,miniFpkm1=0.5,miniFpkm2=0.5, 
                     labeledGenes=c(),label_xlim=c(NA,NA),label_ylim=c(NA,NA),label="text",
                     PerformEnrich=F,EnrichGroups="all"){
  #Data is:
  #Foldchange1 Foldchange2
  #                     Foldchange1 Foldchange2
  #0610007P14Rik              0.086        0.24
  #0610009B22Rik              0.026        0.19
  #0610009L18Rik             -0.077       -0.14
  #0610009O20Rik             -0.149        0.14
  #0610010B08Rik,Gm4724       0.000        0.00
  #0610010F05Rik             -0.339        0.16
  #myXlabel is corresponding to first column,myYlabel to second column
  #cutoff=1.5
  #labeledGenes=c("Pdcd1lg2","Ifng","H2-Ea-ps","Cxcr6")
  #xMin=-4
  #xMax=4
  #yMin=-4
  #yMax=4
  #myXlabel=""
  #myYlabel=""
  # for example
  #p=Nine_Square(py8119,hepa16,myXlabel="KO/WT with py8119",myYlabel="KO/WT with Hepa1-6",
  #              labeledGenes=labeledGenes,label_xlim = c(-6,-2),label_ylim=c(-6,-1),
  #              xMin = -6,xMax = 4,yMin = -6,yMax = 4,PerformEnrich=F,EnrichGroups="all",org="mmu")
  #=============================
  if(myXlabel==""){
    myXlabel=str_split_fixed(diffGeneEnrich1$title,"_Fpkm",2)[1]
  }
  if(myYlabel==""){
    myYlabel=str_split_fixed(diffGeneEnrich2$title,"_Fpkm",2)[1]
  }
  #=============================
  if(Fpkm_filter==TRUE){
    data=diffGeneEnrich1$exp
    idx=is.finite(data$log2.fold_change)&data$status=="OK"
    FoldChange1=data[data[,3]>miniFpkm1 & data[,4]>miniFpkm1 & idx,c(1,5)]
    names(FoldChange1)[2]="Foldchange1"
    
    data=diffGeneEnrich2$exp
    idx=is.finite(data$log2.fold_change)&data$status=="OK"
    FoldChange2=data[data[,3]>miniFpkm2 & data[,4]>miniFpkm2 & idx,c(1,5)]
    names(FoldChange2)[2]="Foldchange2"
  }else{
    data=diffGeneEnrich1$exp
    FoldChange1=data[,c(1,5)]
    names(FoldChange1)[2]="Foldchange1"
    
    data=diffGeneEnrich2$exp
    FoldChange2=data[,c(1,5)]
    names(FoldChange2)[2]="Foldchange2"
  }
  
  Data=merge(FoldChange1,FoldChange2,by="geneSymbol")
  row.names(Data)=Data$geneSymbol
  Data$geneSymbol=NULL
  idx=!(is.infinite(Data$Foldchange1)|is.infinite(Data$Foldchange2))
  Data=Data[idx,]
  #=============================
  upcutoff=log2(cutoff)
  dwcutoff=-log2(cutoff)
  idxA=Data$Foldchange1<dwcutoff & Data$Foldchange2>upcutoff
  idxB=abs(Data$Foldchange1)<upcutoff & Data$Foldchange2>upcutoff
  idxC=Data$Foldchange1>upcutoff & Data$Foldchange2>upcutoff
  
  idxD=Data$Foldchange1<dwcutoff & abs(Data$Foldchange2)<upcutoff
  idxE=abs(Data$Foldchange1)<upcutoff & abs(Data$Foldchange2)<upcutoff
  idxF=Data$Foldchange1>upcutoff & abs(Data$Foldchange2)<upcutoff
  
  idxG=Data$Foldchange1<dwcutoff & Data$Foldchange2<dwcutoff
  idxH=abs(Data$Foldchange1)<upcutoff & Data$Foldchange2<dwcutoff
  idxI=Data$Foldchange1>upcutoff & Data$Foldchange2<dwcutoff
  
  Data$group="no"
  Data$group[idxA]="A"
  Data$group[idxB]="B"
  Data$group[idxC]="C"
  Data$group[idxD]="D"
  Data$group[idxE]="E"
  Data$group[idxF]="F"
  Data$group[idxG]="G"
  Data$group[idxH]="H"
  Data$group[idxI]="I"
  Data2=Data[Data$group!="no",]
  #==================================
  Data2$Foldchange1[Data2$Foldchange1>xMax]=xMax
  Data2$Foldchange1[Data2$Foldchange1<xMin]=xMin
  Data2$Foldchange2[Data2$Foldchange2>yMax]=yMax
  Data2$Foldchange2[Data2$Foldchange2<yMin]=yMin
  Data2$GeneSymbol=row.names(Data2)
  
  #Convert gene symbol into ENTREZID 
  {
    tmp=Data2[,c("GeneSymbol","Foldchange1")]
    names(tmp)=c("SYMBOL","value")
    tmp1=geneIDconvertion(tmp,org)
    tmp1=tmp1$genes_values[,c(1,3)]
    names(tmp1)=c("GeneSymbol","ENTREZID")
    tmp1=tmp1[!duplicated(tmp1$GeneSymbol),]
    Data2=merge(Data2,tmp1,by="GeneSymbol",all.x = T)
    row.names(Data2)=Data2$GeneSymbol
    Data2=Data2[,c("Foldchange1","Foldchange2","group","GeneSymbol","ENTREZID")]
  }
  
  #===================================
  #plot 9 squares
  mycolour=c("others"="aliceblue",  "A"="magenta",  "C"="blue","B"="red", "D"="limegreen")
  p=ggplot(Data2,aes(x=Foldchange1,y=Foldchange2,colour=group,fill=group))+xlim(xMin-0.01,xMax+0.01)+ylim(yMin-0.01,yMax+0.01)
  p=p+geom_point(shape=".",alpha=1/1,size = 1)#+scale_color_manual(values=mycolour)
  p=p+geom_jitter(size = 1)
  p=p+theme_classic(8)
  p=p+geom_vline(xintercept = upcutoff,linetype = "dotted")
  p=p+geom_vline(xintercept = dwcutoff,linetype = "dotted")
  p=p+geom_hline(yintercept = upcutoff,linetype = "dotted") 
  p=p+geom_hline(yintercept = dwcutoff,linetype = "dotted") 
  #p=p+geom_vline(xintercept = ctrl_mean,linetype = "dotted")+geom_hline(yintercept = treat_mean,linetype = "dotted") 
  p=p+annotate("text",color="black",label=paste("A:",as.character(dim(Data2[Data2$group=="A",])[1]),sep=""), x=xMin/2, y=yMax/2,hjust = 0.5,size = 3)
  p=p+annotate("text",color="black",label=paste("B:",as.character(dim(Data2[Data2$group=="B",])[1]),sep=""), x=0,     y=yMax/2,hjust = 0.5,size = 3)
  p=p+annotate("text",color="black",label=paste("C:",as.character(dim(Data2[Data2$group=="C",])[1]),sep=""), x=xMax/2,  y=yMax/2,hjust = 0.5,size = 3)
  p=p+annotate("text",color="black",label=paste("D:",as.character(dim(Data2[Data2$group=="D",])[1]),sep=""), x=xMin/2, y=0,hjust = 0.5,size = 3)
  p=p+annotate("text",color="black",label=paste("E:",as.character(dim(Data2[Data2$group=="E",])[1]),sep=""), x=0,     y=0,hjust = 0.5,size = 3)
  p=p+annotate("text",color="black",label=paste("F:",as.character(dim(Data2[Data2$group=="F",])[1]),sep=""), x=xMax/2,  y=0,hjust = 0.5,size = 3)
  p=p+annotate("text",color="black",label=paste("G:",as.character(dim(Data2[Data2$group=="G",])[1]),sep=""), x=xMin/2,  y=yMin/2,hjust = 0.5,size = 3)
  p=p+annotate("text",color="black",label=paste("H:",as.character(dim(Data2[Data2$group=="H",])[1]),sep=""), x=0,     y=yMin/2,hjust = 0.5,size = 3)
  p=p+annotate("text",color="black",label=paste("I:",as.character(dim(Data2[Data2$group=="I",])[1]),sep=""), x=xMax/2, y=yMin/2,hjust = 0.5,size = 3)
  
  #p=p+geom_abline(slope=treat_down/ctrl_left, intercept=-0.05,linetype = "dotted")
  #p=p+geom_abline(slope=treat_down/ctrl_left, intercept=0.05,linetype = "dotted")
  p=p+xlab(myXlabel)+ylab(myYlabel)
  p=p+theme(legend.position="null")+theme(legend.title=element_blank())
  p=p+guides(col = guide_legend(ncol = 2,byrow =T))
  
  if(length(labeledGenes)>0){
    idx=Data2$GeneSymbol %in% labeledGenes
    Data3=Data2[idx,]
    if(label=="box"){
      p=p+geom_label_repel(data=Data3,aes(x=Foldchange1,y=Foldchange2,
                                          label = as.character(GeneSymbol)),xlim=label_xlim,ylim=label_ylim,
                           fontface = 'bold', color = 'white',size = 3,force=2,
                           box.padding = unit(0.5, "lines"),segment.color = 'grey50',
                           point.padding = unit(0.2, "lines")) 
    }else{
      p=p+geom_text_repel(data=Data3,aes(x=Foldchange1,y=Foldchange2,
                                         label = as.character(GeneSymbol)),xlim=label_xlim,ylim=label_ylim,
                          fontface = 'bold', color = 'black',size = 3,force=2,
                          box.padding = unit(0.5, "lines"),segment.color = 'grey50',
                          point.padding = unit(0.2, "lines"))
    }
  }
  
  results=list()
  results$nineSquare=p
  names(Data2)[grep("Foldchange1",names(Data2))]=myXlabel
  names(Data2)[grep("Foldchange2",names(Data2))]=myYlabel
  results$data=Data2
  #===================================
  #perform enrichment for groups
  if(PerformEnrich==T){
    if (org=="hsa"){
      suppressPackageStartupMessages(library(org.Hs.eg.db))
      OrgDb="org.Hs.eg.db"
    }else if(org=="mmu"){
      suppressPackageStartupMessages(library(org.Mm.eg.db))
      OrgDb="org.Mm.eg.db"
    }else if(org=="rno"){
      suppressPackageStartupMessages(library(org.Rn.eg.db))
      OrgDb="org.Rn.eg.db"
    }else{
      cat("Please specify species, according to http://www.genome.jp/kegg/catalog/org_list.html")
      exit()
    }
    
    Enrichment9square=list()
    data=p$data[,c("GeneSymbol","Foldchange1")]
    names(data)=c("SYMBOL","value")
    geneList=geneIDconvertion(data,org=org)
    geneList=geneList$geneList
    if(EnrichGroups=="all"){EnrichGroups = c("A","B","C","D","E","F","G","H","I")}
    for(k in 1:length(EnrichGroups)){
      group=EnrichGroups[k]
      dd=Data2[Data2$group==group,]
      genes=suppressWarnings(bitr(as.character(row.names(dd)), fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb, drop = TRUE))
      genes=genes$ENTREZID
      Enrichment9square[[k]]=enrichment_analysis(geneList,genes,plotTitle=paste(group,"_group",sep=""),org=org, output=T,type=c("GO","KEGG"),GOcatagory="BP")
      names(Enrichment9square)[k]=EnrichGroups[k]
    }
    results$Enrichment=Enrichment9square
  }  
  
  return(results)
}

geneIDconvertion=function(genes_values,org){
  # genes_values format
  #SYMBOL value
  #Fsbp    -9.3
  #Ntn5    -8.2
  #Btg3    -7.1
  #
  # or
  #        REFSEQ value
  #    NR_157155     1
  # NM_001001981     2
  #
  if (org=="hsa"){
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    OrgDb="org.Hs.eg.db"
  }else if(org=="mmu"){
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    OrgDb="org.Mm.eg.db"
  }else if(org=="rno"){
    suppressPackageStartupMessages(library(org.Rn.eg.db))
    OrgDb="org.Rn.eg.db"
  }else{
    cat("Please specify species, according to http://www.genome.jp/kegg/catalog/org_list.html")
    exit()
  }
  
  # to determine genename type: REFSEQ or SYMBOL
  idx=c(grep("NM_",genes_values[,1]),grep("NR_",genes_values[,1]),grep("NP_",genes_values[,1]),grep("XM_",genes_values[,1]),grep("XP_",genes_values[,1]),grep("XR_",genes_values[,1]))    
  if(length(idx)>0){
    names(genes_values)=c("REFSEQ","value")
  }else{
    names(genes_values)=c("SYMBOL","value")
  }
  orignType=names(genes_values)[1]
  genes=suppressWarnings(bitr(as.character(genes_values[,1]), fromType=orignType, toType="ENTREZID", OrgDb=OrgDb, drop = TRUE))
  genes_values_ENTREZID=merge(genes_values,genes,by=orignType)
  geneList=genes_values_ENTREZID$value
  names(geneList)=genes_values_ENTREZID$ENTREZID
  results=list()
  results$geneList=geneList
  results$genes_values=genes_values_ENTREZID
  return(results)
}
#============================================================
setwd("~/github/Nono/")
load("WT_D12_vs_D0.RData")
WT=diffGeneEnrich
load("NonoKO_D12_vs_D0.RData")
KO=diffGeneEnrich
rm(diffGeneEnrich)
load("Rescue_D12_vs_D0.RData")
RE=diffGeneEnrich
WT_KO=Nine_Square(WT,KO,myXlabel="log2(D12/D0) in WT",myYlabel="log2(D12/D0) in KO",org="mmu")
#============================================================
myGroup=WT_KO$data[WT_KO$data$group=="I"|WT_KO$data$group=="F",]
genes=row.names(myGroup)

RE1=RE
RE1$exp=RE1$exp[RE1$exp$geneSymbol %in%genes, ]
KO1=KO
KO1$exp=KO1$exp[KO1$exp$geneSymbol %in%genes, ]

RE1_KO1=Nine_Square(RE1,KO1,myXlabel="log2(D12/D0) in KO+WT",myYlabel="log2(D12/D0) in KO",org="mmu")
RE1_KO1$nineSquare
