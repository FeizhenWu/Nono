setwd("~/github/Nono/")
rm(list=ls())
#functions
{
  volcano_plot=function(diffGeneEnrich,boundary=3,labeledGenes=c(),UptopK=3,DwtopK=3,label="text",label_xlim=c(NA,NA),label_ylim=c(NA,NA)){
    #volcano_plot(diffGeneEnrich,topK=0)
    #boundary=3
    #labeledGenes=c()
    #topK=10
    #label="box"
    #label_xlim=c(NA,NA)
    #label_ylim=c(NA,NA)
    
    foldchange=log2(diffGeneEnrich$foldchange)
    pvalue=diffGeneEnrich$pvalue
    fpkmcutoff=diffGeneEnrich$fpkmcutoff
    tfpkm=diffGeneEnrich$fpkm_total
    
    dd=diffGeneEnrich$exp
    dd=dd[is.finite(dd$log2.fold_change)&dd$status=="OK",]
    dd=dd[(dd[,3]>=fpkmcutoff & dd[,4]>=fpkmcutoff),]
    dd=dd[(dd[,3]+dd[,4])>=tfpkm,]
    
    dd$pvalue[dd$pvalue<1e-10]=1e-10
    dd$Nlog10Pvalue=-log10(dd$pvalue)
    dd=dd[dd$Nlog10Pvalue!=0,]
    dd$geneSymbol=as.character(dd$geneSymbol)
    dd$cor="no"
    
    dd$log2.fold_change[dd$log2.fold_change>boundary]=boundary
    dd$log2.fold_change[dd$log2.fold_change<(-boundary)]=-boundary
    
    diffgenes=diffGeneEnrich$diffgenes
    UpGenes=as.character(diffgenes$gene_id[diffgenes$log2.fold_change>0])
    dd$cor[dd$geneSymbol %in% UpGenes]="up"
    DwGenes=as.character(diffgenes$gene_id[diffgenes$log2.fold_change<0])
    dd$cor[(dd$geneSymbol %in% DwGenes)]="down"
    
    mycolour=c("no"="gray",  "up"="red",  "down"="blue")
    p=ggplot(dd,aes(x=log2.fold_change,y=Nlog10Pvalue,colour=cor))
    p=p+geom_point(shape=".",alpha=1/1,size = 1)+scale_color_manual(values=mycolour)
    p=p+geom_jitter(size = 1)+theme_classic(8)
    p=p+xlab("log2(fold-change)")+ylab("-log10(Pvalue)")
    p=p+theme(legend.position="null")+theme(legend.title=element_blank())
    p=p+geom_vline(xintercept = foldchange,linetype = "dotted")
    p=p+geom_vline(xintercept = -foldchange,linetype = "dotted")
    p=p+geom_hline(yintercept = -log10(pvalue),linetype = "dotted") 
    p=p+annotate("text",color="black",label=paste("Down:",as.character(dim(dd[dd$cor=="down",])[1]),sep=""), x=-(boundary+foldchange)/2, y=max(dd$Nlog10Pvalue)*1.3,hjust = 0.5,size = 2.5)
    p=p+annotate("text",color="black",label=paste("Up:",as.character(dim(dd[dd$cor=="up",])[1]),sep=""), x=(boundary+foldchange)/2, y=max(dd$Nlog10Pvalue)*1.3,hjust = 0.5,size = 2.5)
    p=p+xlim(-boundary,boundary)+ylim(0,max(dd$Nlog10Pvalue)*1.4)
    
    if(length(labeledGenes)==0 && (UptopK+DwtopK)>0){
      dd1=dd[dd$log2.fold_change>0,]
      dd1=dd1[order(dd1$Nlog10Pvalue*abs(dd1$log2.fold_change),decreasing=T),]
      labeledGenes=dd1$geneSymbol[1:UptopK]
      
      dd1=dd[dd$log2.fold_change<0,]
      dd1=dd1[order(dd1$Nlog10Pvalue*abs(dd1$log2.fold_change),decreasing=T),]
      labeledGenes=c(labeledGenes,dd1$geneSymbol[1:DwtopK])
    }
    if(length(labeledGenes)>0){
      idx=dd$geneSymbol %in% labeledGenes
      dd1=dd[idx,]
      if(label=="box"){
        p=p+geom_label_repel(data=dd1,aes(x=log2.fold_change,y=Nlog10Pvalue,
                                          label = geneSymbol),xlim=label_xlim,ylim=label_ylim,
                             fontface = 'bold', color = 'black',size = 2,force=2,
                             box.padding = unit(0.5, "lines"),segment.color = 'grey50',
                             point.padding = unit(0.2, "lines")) 
      }else{
        p=p+geom_text_repel(data=dd1,aes(x=log2.fold_change,y=Nlog10Pvalue,
                                         label = geneSymbol),xlim=label_xlim,ylim=label_ylim,
                            fontface = 'bold', color = 'black',size = 2,force=2,
                            box.padding = unit(0.5, "lines"),segment.color = 'grey50',
                            point.padding = unit(0.2, "lines"))
      }
    }
    return(p)
  }
  TopN_diffgenes=function(diffGeneEnrich,geneNum=20,horizontal=T,downBoundary=-NA,upBoundary=NA,colourCutoff=NA){
    diffgenes=diffGeneEnrich$diffgenes
    diffgenes=diffgenes[order(abs(diffgenes$log2.fold_change),decreasing = T),]
    genes=as.character(diffgenes[1:geneNum,"gene_id"])
    p=GeneHeatmap(genes,diffGeneEnrich$Enrich)
    return(p)
  }
  GeneHeatmap=function(genes,Enrichment,horizontal=T,gene_Feature="",geneNum=NA,downBoundary=-NA,upBoundary=NA,colourCutoff=NA){
    #genes=c("Ccl1","Ccl6","Ccl9","Ccr2","Cxcl1","Cxcl9","Gbp4")
    #Enrichment=gMDSC_KO_H_vs_gMDSC_WT_H
    #gene_Feature="chemokine-related-genes
    #
    #geneheat=GeneHeatmap(genes, diffGeneEnrich$Enrich,gene_Feature ="chemokine-related-genes")
    #
    title = Enrichment$title
    if(gene_Feature==""){
      gene_Feature="geneList"
    }
    
    exp=Enrichment$data[Enrichment$data$SYMBOL %in% genes,c(3,2)]
    row.names(exp)=exp$SYMBOL
    
    if(!is.na(geneNum) && nrow(exp)>=geneNum){
      exp=exp[order(abs(exp$value),decreasing = T),]
      exp=exp[1:geneNum,]
    }
    
    exp=exp[order(exp$value,decreasing = T),]
    if(!is.na(downBoundary)){
      exp$value[exp$value<downBoundary]=downBoundary
    }
    if(!is.na(upBoundary)){
      exp$value[exp$value>upBoundary]=upBoundary
    }
    
    if(horizontal==T){
      exp1=t(exp$value) 
      colnames(exp1)=row.names(exp)
      rownames(exp1)=names(exp)[2]
      rownames(exp1)[1]="log2_FC"
      h=0.1+0.2+max(nchar(colnames(exp1)))*0.08
      pdf(file = paste(title,"@",gene_Feature,".pdf",sep=""), width = 0.8+ncol(exp1)*0.25,height = h,onefile=F)
    }else{
      exp1=data.frame(value=exp$value)
      row.names(exp1)=row.names(exp)
      names(exp1)=names(exp)[2]
      colnames(exp1)[1]="log2_FC"
      h=0.1+0.2+max(nchar(row.names(exp1)))*0.1
      pdf(file = paste(title,"@",gene_Feature,".pdf",sep=""), width = h,height = 1+nrow(exp1)*0.2,onefile=F)
    }
    
    if(is.na(colourCutoff)){
      bb=max(abs(exp1)+0.1)
      breaksList = seq(-bb,bb, by = 0.1)
    }else{
      bb=colourCutoff
      breaksList = seq(-bb,bb, by = 0.1)
    }
    p=pheatmap(exp1,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = T,
               show_colnames = T, show_rownames  = T, legend = F, border_color="white",
               number_format = "%.1f", number_color = "white", fontsize = 7,
               color = colorRampPalette(c("blue","navy", "firebrick3"))(length(breaksList)),
               breaks = breaksList, main = gene_Feature,
               fontsize_number =7)
    #p
    #dev.off()
    results=list()
    results$data=exp1
    results$pheatmap=p
    return(results)
  }
  Heatmap_KEGG_GO_Genes=function(pathways,Enrichment,type="KEGG",horizontal=T,geneNum=NA,excludedGenes=NA,downBoundary=-NA,upBoundary=NA,cutoff=NA)
  {
    #pathway="mmu03030"
    #title="gMDSC_KO_H_vs_WT_H"
    #Enrichment=gMDSC_KO_H_vs_gMDSC_WT_H
    #pathway="GO:0007186",Enrichment=diffGeneEnrich$Enrich,type="BP";horizontal=T;upBoundary=-2;downBoundary=2
    #GOBP5=Heatmap_KEGG_GO_Genes("GO:0007186",diffGeneEnrich$Enrich,type="BP",horizontal=T)
    #
    title = Enrichment$title
    if(type=="KEGG"){
      enrichedTerms=Enrichment$KEGGenrich@result
    }else if (type=="BP"){
      enrichedTerms=Enrichment$GOBPenrich@result
    }else if (type=="MF"){
      enrichedTerms=Enrichment$GOMFenrich@result
    }else if (type=="CC"){
      enrichedTerms=Enrichment$GOCCenrich@result
    }
    
    Exps=list()
    for(pathway in pathways){
      termName=pathway
      genes=enrichedTerms$geneID[enrichedTerms$ID==pathway]
      pathwayName = enrichedTerms$Description[enrichedTerms$ID==pathway]
      genes=str_split(genes,"/")[[1]]
      
      pathway=paste(gsub(":","",pathway),gsub(" ","_",pathwayName),sep="=")
      
      exp=Enrichment$data[Enrichment$data$SYMBOL %in% genes,c(3,2)]
      row.names(exp)=exp$SYMBOL
      
      if(length(na.omit(excludedGenes))>0){
        idx=exp$SYMBOL %in% excludedGenes
        exp=exp[!idx,]
      }
      
      if(!is.na(cutoff)){
        exp=exp[abs(exp$value)>cutoff,]
      }
      
      if(!is.na(geneNum) && nrow(exp)>=geneNum){
        exp=exp[order(abs(exp$value),decreasing = T),]
        exp=exp[1:geneNum,]
      }
      
      exp=exp[order(exp$value,decreasing = T),]
      if(!is.na(downBoundary)){
        exp$value[exp$value<downBoundary]=downBoundary
      }
      if(!is.na(upBoundary)){
        exp$value[exp$value>upBoundary]=upBoundary
      }
      
      names(exp)[2]="log2 (FC)"
      if(horizontal==T){
        exp1=t(exp$`log2 (FC)`) 
        colnames(exp1)=row.names(exp)
        rownames(exp1)=names(exp)[2]
        w=0.8+ncol(exp1)*0.25
        h=0.1+0.2+max(nchar(colnames(exp1)))*0.1
        if(nchar(pathwayName)*0.08>w){w=nchar(pathwayName)*0.08}
      }else{
        exp1=data.frame(log2_FC=exp$`log2 (FC)`)
        row.names(exp1)=row.names(exp)
        names(exp1)=names(exp)[2]
        w=0.1+0.2+max(nchar(row.names(exp1)))*0.1
        h=1+nrow(exp1)*0.2
        if(nchar(pathwayName)>h){h=nchar(pathwayName)}
      }
      
      bb=max(abs(exp1)+0.1)
      breaksList = seq(-bb,bb, by = 0.1)
      
      p=pheatmap(exp1,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = T,
                 show_colnames = T, show_rownames  = T, legend = F, border_color="white",
                 number_format = "%.1f", number_color = "white", fontsize = 7,
                 color = colorRampPalette(c("blue","navy", "firebrick3"))(length(breaksList)),
                 breaks = breaksList, main = pathwayName,
                 fontsize_number =7)
      p
      Exps[[termName]]=list()
      Exps[[termName]]$trimmed_value=exp1
      Exps[[termName]]$heatmap=p
    }
    return(Exps)
  }
}
load("3NonoKO_D0_vs_3WT_D0.RData")
Figure5A=volcano_plot(diffGeneEnrich)
Figure5A

pp1=TopN_diffgenes(diffGeneEnrich,30)
Figure5B_Top=pp1$pheatmap
Figure5B_Top

pp2=Heatmap_KEGG_GO_Genes("GO:0061564",diffGeneEnrich$Enrich,type="BP",geneNum=30)
Figure5B_Bottom=pp2$`GO:0061564`$heatmap
Figure5B_Bottom