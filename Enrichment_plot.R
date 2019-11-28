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
