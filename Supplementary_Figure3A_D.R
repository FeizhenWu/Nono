#!/usr/bin/Rscript --vanilla
setwd("~/github/Nono/")
rm(list=ls())
library(ggplot2)
library(reshape)
library(scales)
library(clusterProfiler)

#Functions===========================
{
  orgData=function(dd){
    dd[,2:5]=t(apply(dd[,2:5],1,rescale))
    dd=melt(dd,id="geneSymbol")
    dd1=str_split_fixed(dd$variable,"_",2)
    dd=cbind(dd,dd1)[,-2]
    names(dd)=c("geneSymbol","Exp","CellType","Day")
    return(dd)
  }
  DifferentialCourse=function(myAll,myTitle){
    Ave=aggregate(Exp ~ CellType+Day,data=myAll,FUN="median")
    Ave=Ave[order(Ave$CellType),]
    dodge <- position_dodge(width = 0.4)
    p=ggplot(myAll,aes(x=Day,y=Exp))
    p=p+geom_boxplot(position = dodge,na.rm=T,colour="gray80")
    p=p+facet_grid(CellType~.)
    p=p+stat_summary(fun.y=median, geom="point", size=2, color="red")
    p=p+theme_bw(8)+ggtitle(label = myTitle)+theme(plot.title = element_text(hjust=0.5))
    p=p+xlab("Differentiation Course")+ylab("Normalized FPKM")
    p
    return(p)
  }
  Nine_Square=function(diffGeneEnrich1,diffGeneEnrich2,myXlabel="",myYlabel="",org="hsa",
                       xMin=-4,xMax=4,yMin=-4,yMax=4,
                       cutoff=1.5,miniFpkm1=0.5,miniFpkm2=0.5, 
                       labeledGenes=c(),label_xlim=c(NA,NA),label_ylim=c(NA,NA),label="text",
                       PerformEnrich=F,EnrichGroups="all"){
    #=============================
    if(myXlabel==""){
      myXlabel=str_split_fixed(diffGeneEnrich1$title,"_Fpkm",2)[1]
    }
    if(myYlabel==""){
      myYlabel=str_split_fixed(diffGeneEnrich2$title,"_Fpkm",2)[1]
    }
    #=============================
    data=diffGeneEnrich1$exp
    FoldChange1=data[,c(1,5)]
    names(FoldChange1)[2]="Foldchange1"
    
    data=diffGeneEnrich2$exp
    FoldChange2=data[,c(1,5)]
    names(FoldChange2)[2]="Foldchange2"
    
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
    
    results=list()
    results$nineSquare=p
    names(Data2)[grep("Foldchange1",names(Data2))]=myXlabel
    names(Data2)[grep("Foldchange2",names(Data2))]=myYlabel
    results$data=Data2
    #===================================
    #perform enrichment for groups
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
  
  enrichment_analysis=function(geneList,genes,plotTitle="",org="hsa", output=F, 
                               type=c("GO"),GOcatagory=c("BP","CC","MF"), 
                               pvalueCutoff=1,qvalueCutoff=1,termNum=10,charLength=50,fontsize=0,figformat="pdf"){
    if (org=="hsa"){
      suppressPackageStartupMessages(library(org.Hs.eg.db))
      OrgDb=org.Hs.eg.db
      orgdb="org.Hs.eg.db"
    }else if(org=="mmu"){
      suppressPackageStartupMessages(library(org.Mm.eg.db))
      OrgDb=org.Mm.eg.db
      orgdb="org.Mm.eg.db"
    }else if(org=="rno"){
      suppressPackageStartupMessages(library(org.Rn.eg.db))
      OrgDb=org.Rn.eg.db  
      orgdb="org.Rn.eg.db"
    }else{
      cat("Please specify species, according to http://www.genome.jp/kegg/catalog/org_list.html")
      exit()
    }
    
    results=list()
    results$title=plotTitle
    data=data.frame(ENTREZID=names(geneList), value=geneList)
    tmp=as.data.frame(bitr(data$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb=orgdb))
    data=merge(data,tmp,by="ENTREZID")
    data$diff=data$ENTREZID %in% genes
    results$data=data
    
    type=sapply(type,toupper)
    geneList=geneList[order(geneList,decreasing = T)]
    if (any(type %in% "GO")){
      if (any(GOcatagory %in% "BP"))
        tryCatch({
          results$GOBPenrich <- enrichGO(gene     = genes,
                                         universe      = names(geneList),
                                         OrgDb         = OrgDb,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         minGSSize = 10, 
                                         maxGSSize = 500,
                                         pvalueCutoff  = pvalueCutoff,
                                         qvalueCutoff  = qvalueCutoff,
                                         readable      = TRUE)
          if (is.null(results$GOBPenrich)){
            results$GOBPenrich_plot=NULL
          }else{
            Rev=enrichment_plot(results$GOBPenrich@result,paste(plotTitle,"GOBP"),"blue",termNum,charLength,fontsize)
            results$GOBPenrich_plot = Rev$plot
            if (output){
              write.table(results$GOBPenrich@result,file = paste(plotTitle,"GO_BiologicalProcess_enrichment.txt",sep="_"), row.names = F,quote = F,sep="\t")
              ggsave(results$GOBPenrich_plot,filename = paste(plotTitle,"_GO_BiologicalProcess_enrichment.",figformat,sep=""),height = Rev$height/100,width = Rev$width/100)
            }
          }
        })
      if (any(GOcatagory %in% "CC"))
        tryCatch({
          results$GOCCenrich <- enrichGO(gene     = genes,
                                         universe      = names(geneList),
                                         OrgDb         = OrgDb,
                                         ont           = "CC",
                                         pAdjustMethod = "BH",
                                         minGSSize = 10, 
                                         maxGSSize = 500,
                                         pvalueCutoff  = pvalueCutoff,
                                         qvalueCutoff  = qvalueCutoff,
                                         readable      = TRUE)
          if (is.null(results$GOCCenrich)){
            results$GOCCenrich_plot=NULL
          }else{
            Rev=enrichment_plot(results$GOCCenrich@result,paste(plotTitle,"GOCC"),"blue",termNum,charLength,fontsize)
            results$GOCCenrich_plot=Rev$plot
            if (output){          
              write.table(results$GOCCenrich@result,file = paste(plotTitle,"GO_CellularComponent_enrichment.txt",sep="_"), row.names = F,quote = F,sep="\t")
              ggsave(results$GOCCenrich_plot,filename = paste(plotTitle,"_GO_CellularComponent_enrichment.",figformat,sep=""),height = Rev$height/100,width = Rev$width/100)
            }
          }
        })
      if (any(GOcatagory %in% "MF"))
        tryCatch({
          results$GOMFenrich <- enrichGO(gene     = genes,
                                         universe      = names(geneList),
                                         OrgDb         = OrgDb,
                                         ont           = "MF",
                                         pAdjustMethod = "BH",
                                         minGSSize = 10, 
                                         maxGSSize = 500,
                                         pvalueCutoff  = pvalueCutoff,
                                         qvalueCutoff  = qvalueCutoff,
                                         readable      = TRUE)
          if (is.null(results$GOMFenrich)){
            results$GOMFenrich_plot=NULL
          }else{
            Rev=enrichment_plot(results$GOMFenrich@result,paste(plotTitle,"GOMF"),"blue",termNum,charLength,fontsize)
            results$GOMFenrich_plot=Rev$plot
            if (output){
              write.table(results$GOMFenrich@result,file = paste(plotTitle,"GO_MolecularFunction_enrichment.txt",sep="_"), row.names = F,quote = F,sep="\t")
              ggsave(results$GOMFenrich_plot,filename = paste(plotTitle,"_GO_MolecularFunction_enrichment.",figformat,sep=""),height = Rev$height/100,width = Rev$width/100)
            }
          }
        })
    }
    return(results)
  }
  enrichment_plot=function(enrichment,plotTitle="",gridColour="blue",termNum=10,charLength=50,fontsize=8){
    #enrichment=results$KEGGenrich@result
    #plotTitle="";termNum=10;charLength=50;gridColour="blue"
    
    enrichment$pvalue=-log10(enrichment$pvalue)
    enrichment=enrichment[!is.na(enrichment$ID),]
    enrichment=enrichment[!duplicated(enrichment$Description),]
    if(nrow(enrichment)==0){return(NA)}
    
    if(nrow(enrichment)>=termNum){
      enrichment=enrichment[order(enrichment$pvalue,decreasing = T)[1:termNum],]
    }else{
      enrichment=enrichment[order(enrichment$pvalue,decreasing = T),]
    }
    enrichment=enrichment[order(enrichment$pvalue,decreasing = F),]
    {
      terms=as.character(enrichment$Description)
      terms=lapply(terms,function(x,k){
        x=as.character(x)
        x=paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="")
        if(nchar(x)>k){
          x=substr(x,start=1,stop=k)
          x=gsub("(\\w+)$", "...", x)
        }
        return(x)
      },charLength)
      enrichment$Description=do.call(rbind,terms)
    }
    
    longest=max(nchar(enrichment$Description))
    
    if(fontsize==0){fontsize=8}
    aa=enrichment$Description[!duplicated(enrichment$Description[,1]),1]
    enrichment$Name=factor(enrichment$Description,levels=aa)
    
    p1=ggplot(data=enrichment)
    p1=p1+geom_point(aes(x=pvalue,y=Name,size=Count,colour = factor(GeneRatio)))
    p1=p1+guides(color=F)
    p1=p1+scale_size_continuous(name="Gene Count")
    #p1=p1+guides(fill=guide_legend(title="Gene Count"))
    p1=p1+xlab(paste("-log10(Pvalue)",sep=""))+ylab("")+ggtitle(plotTitle)
    p1=p1+theme(panel.grid.major=element_line(colour=gridColour),
                panel.grid.minor=element_blank(), 
                panel.background = element_blank(),
                legend.title = element_text(size=fontsize,face="plain",colour ='black'),
                legend.text = element_text(hjust = 1,size=fontsize),
                legend.key.size = unit(1,"mm"),
                legend.background = element_rect(),
                plot.title = element_text(hjust = 0.5,size=fontsize+2),
                axis.line = element_blank(),
                axis.title = element_text(size=fontsize,face="plain",colour ='black'),
                axis.text.x=element_text(size=fontsize,face="plain",colour ='black'),
                axis.text.y=element_text(size=fontsize,face="plain",colour ='black'))
    Rev=list()
    Rev$plot = p1
    Rev$width = (longest+2)*5+350
    Rev$height = nrow(enrichment)*25+30
    return(Rev)
  }
}
#load Data===========================
{
  load("WT_D3_vs_D0.RData")
  dd=diffGeneEnrich$exp[,c(1,3,4)]
  load("WT_D6_vs_D3.RData")
  dd=merge(dd,diffGeneEnrich$exp[,c(1,4)],by="geneSymbol")
  load("WT_D12_vs_D0.RData")
  dd=merge(dd,diffGeneEnrich$exp[,c(1,4)],by="geneSymbol")
  WT=orgData(dd)
  load("NonoKO_D3_vs_D0.RData")
  dd=diffGeneEnrich$exp[,c(1,3,4)]
  load("NonoKO_D6_vs_D3.RData")
  dd=merge(dd,diffGeneEnrich$exp[,c(1,4)],by="geneSymbol")
  load("NonoKO_D12_vs_D0.RData")
  dd=merge(dd,diffGeneEnrich$exp[,c(1,4)],by="geneSymbol")
  NonoKO=orgData(dd)
  load("Rescue_D3_vs_D0.RData")
  dd=diffGeneEnrich$exp[,c(1,3,4)]
  load("Rescue_D6_vs_D3.RData")
  dd=merge(dd,diffGeneEnrich$exp[,c(1,4)],by="geneSymbol")
  load("Rescue_D12_vs_D0.RData")
  dd=merge(dd,diffGeneEnrich$exp[,c(1,4)],by="geneSymbol")
  rm(diffGeneEnrich)
  RE=orgData(dd)
  All=rbind(WT,NonoKO,RE)
  All$CellType=sub("RE","Rescure",All$CellType)
  All$Day=sub("D","Day",All$Day)
  All$CellType=factor(All$CellType,levels = c("WT","NonoKO","Rescure"))
  All$Day=factor(All$Day,levels = c("Day0","Day3","Day6","Day12"))
}
{
  load("WT_D3_vs_D0.RData")
  WT=diffGeneEnrich
  load("NonoKO_D3_vs_D0.RData")
  KO=diffGeneEnrich
  load("Rescue_D3_vs_D0.RData")
  RE=diffGeneEnrich
  rm(diffGeneEnrich)
}
#D3_vs_D0===========================
{

  #plot SupFigure 2 A left
  {
    WT_KO=Nine_Square(WT,KO,myXlabel="log2(D3/D0) in WT",myYlabel="log2(D3/D0) in KO",org="mmu")
    SupFigure2A_Left=WT_KO$nineSquare
    SupFigure2A_Left
  }
  
  #plot SupFigure 2 B left and SupFigure 3 A left-up
  {
    myGroup=WT_KO$data[WT_KO$data$group=="I"|WT_KO$data$group=="F",]
    genes=row.names(myGroup)
    RE1=RE
    RE1$exp=RE1$exp[RE1$exp$geneSymbol %in%genes, ]
    KO1=KO
    KO1$exp=KO1$exp[KO1$exp$geneSymbol %in%genes, ]
    Group_FI=Nine_Square(RE1,KO1,myXlabel="log2(D3/D0) in KO+WT",myYlabel="log2(D3/D0) in KO",org="mmu")
    SupFigure2B_Left=Group_FI$nineSquare
    #============================================================
    myGroup=Group_FI$data[Group_FI$data$group=="I"|Group_FI$data$group=="F",]
    genes=row.names(myGroup)
    myAll=All[All$geneSymbol %in% genes,]
    myTitle="D3 vs D0 Dynamic of Group F&I"
    
    SupFigure3A_LeftUp=DifferentialCourse(myAll,myTitle)
    SupFigure3A_LeftUp
  }
  
  #plot SupFigure 3 A left-down
  {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    OrgDb=org.Mm.eg.db
    genes=WT$exp[,c(1,5)]
    names(genes)=c("SYMBOL","value")
    geneList=geneIDconvertion(genes,org="mmu")
    
    myGroup=WT_KO$data[WT_KO$data$group=="I"|WT_KO$data$group=="F",]
    myTitle="Group F&I"
    genes=row.names(myGroup)
    genes=suppressWarnings(bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb, drop = TRUE))
    genes=genes$ENTREZID
    EN=enrichment_analysis(geneList$geneList,genes,plotTitle=myTitle,org="mmu",type=c("GO"),GOcatagory=c("BP"),output=T)
    SupFigure3A_LeftDown=EN$GOBPenrich_plot
    SupFigure3A_LeftDown
  }
  
  #plot SupFigure 3 A right
  {
    myGroup=WT_KO$data[WT_KO$data$group=="I"|WT_KO$data$group=="F",]
    genes=row.names(myGroup)
    WT1=WT
    WT1$exp=WT1$exp[WT1$exp$geneSymbol %in%genes, ]
    RE1=RE
    RE1$exp=RE1$exp[RE1$exp$geneSymbol %in%genes, ]
    KO1=KO
    KO1$exp=KO1$exp[KO1$exp$geneSymbol %in%genes, ]
    
    dd=merge(WT1$exp[c(1,5)],KO1$exp[c(1,5)],by="geneSymbol")
    names(dd)=c("geneSymbol","WT","NonoKO")
    dd=merge(dd,RE1$exp[c(1,5)],by="geneSymbol")
    names(dd)=c("geneSymbol","WT","NonoKO","Rescure")
    fc=log2(1.5)
    idx1=dd$WT>fc&dd$NonoKO<(-fc)&dd$Rescure>fc
    idx2=is.finite(dd$Rescure)
    dd1=dd[idx1&idx2,]
    row.names(dd1)=dd1$geneSymbol
    dd1$geneSymbol=NULL
    dd2=t(apply(dd1,1,scale))
    dd2=as.data.frame(dd2)
    names(dd2)=names(dd1)
    dd2=dd2[order(dd2$Rescure,decreasing = T),]
    dd2=dd2[1:50,]
    SupFigure3A_Right=pheatmap(dd2,cluster_cols = F,fontsize = 7,main="Top 50 Representative Rescured Genes")
    SupFigure3A_Right
  }
  
  #plot SupFigure 2 C left and SupFigure 3 D left-up
  {
    myGroup=WT_KO$data[WT_KO$data$group=="A"|WT_KO$data$group=="D",]
    genes=row.names(myGroup)
    RE1=RE
    RE1$exp=RE1$exp[RE1$exp$geneSymbol %in%genes, ]
    KO1=KO
    KO1$exp=KO1$exp[KO1$exp$geneSymbol %in%genes, ]
    Group_AD=Nine_Square(RE1,KO1,myXlabel="log2(D3/D0) in KO+WT",myYlabel="log2(D3/D0) in KO",org="mmu")
    SupFigure2C_Left=Group_AD$nineSquare
    
    myGroup=Group_AD$data[Group_AD$data$group=="A"|Group_AD$data$group=="D",]
    genes=row.names(myGroup)
    myAll=All[All$geneSymbol %in% genes,]
    myTitle="D3 vs D0 Dynamic of Group A&D"
    
    SupFigure3D_LeftUP=DifferentialCourse(myAll,myTitle)
    SupFigure3D_LeftUP 
  }

  #plot SupFigure 3 D left-down
  {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    OrgDb=org.Mm.eg.db
    genes=WT$exp[,c(1,5)]
    names(genes)=c("SYMBOL","value")
    geneList=geneIDconvertion(genes,org="mmu")
    
    myGroup=Group_AD$data[Group_AD$data$group=="A"|Group_AD$data$group=="D",]
    myTitle="Group A&D"
    genes=row.names(myGroup)
    genes=suppressWarnings(bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb, drop = TRUE))
    genes=genes$ENTREZID
    EN=enrichment_analysis(geneList$geneList,genes,plotTitle=myTitle,org="mmu",type=c("GO"),GOcatagory=c("BP"),output=T)
    SupFigure3D_LeftDown=EN$GOBPenrich_plot
    SupFigure3D_LeftDown
  }

  #plot SupFigure 3 D right
  {
    myGroup=WT_KO$data[WT_KO$data$group=="A"|WT_KO$data$group=="D",]
    genes=row.names(myGroup)
    WT1=WT
    WT1$exp=WT1$exp[WT1$exp$geneSymbol %in%genes, ]
    RE1=RE
    RE1$exp=RE1$exp[RE1$exp$geneSymbol %in%genes, ]
    KO1=KO
    KO1$exp=KO1$exp[KO1$exp$geneSymbol %in%genes, ]
    
    dd=merge(WT1$exp[c(1,5)],KO1$exp[c(1,5)],by="geneSymbol")
    names(dd)=c("geneSymbol","WT","NonoKO")
    dd=merge(dd,RE1$exp[c(1,5)],by="geneSymbol")
    names(dd)=c("geneSymbol","WT","NonoKO","Rescure")
    fc=log2(1.5)
    idx1=dd$WT<(-fc)&dd$NonoKO>fc&dd$Rescure<(-fc)
    idx2=is.finite(dd$Rescure)
    dd1=dd[idx1&idx2,]
    row.names(dd1)=dd1$geneSymbol
    dd1$geneSymbol=NULL
    dd2=t(apply(dd1,1,scale))
    dd2=as.data.frame(dd2)
    names(dd2)=names(dd1)
    dd2=dd2[order(dd2$WT,decreasing = T),]
    #dd2=dd2[1:50,]
    SupFigure3D_Right=pheatmap(dd2,cluster_cols = F,fontsize = 7,main="Top 50 Representative Rescured Genes")
    SupFigure3D_Right
  }
}
