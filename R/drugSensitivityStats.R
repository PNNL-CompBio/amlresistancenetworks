

#now we might wantt o examine some proteins
#' plotPRotsByMetric
#' plots basal protein expression by metric for a specific list of proteins
#' @param sens.data formatted data
#' @param genelist list of gene names to plot
#' @param metric either AUC or IC50
#' @param listname list
#' @return just plots
#' @export
plotProtsByMetric<-function(sens.data,genelist=c("BCL2","TRIM21","DUSP23"),
                            metric='AUC',listname='selected'){
  require(ggplot2)
  library(viridis)
  dat<-sens.data%>%
    subset(Gene%in%genelist)%>%
    subset(Metric==metric)
  
  ggplot2::ggplot(dat,ggplot2::aes(x=LogFoldChange,y=Value))+
    ggplot2::geom_point(ggplot2::aes(color=Condition))+
    ggplot2::geom_smooth(ggplot2::aes(color = Condition, fill = Condition), method = "lm") + 
    viridis::scale_fill_viridis(discrete=T)+
    viridis::scale_color_viridis(discrete=T,option='D')+
    ggplot2::ggtitle(paste("Protein expression vs.",metric))+
    ggplot2::facet_grid(Gene~.)
  
  ggplot2::ggsave(paste0('foldChangeVs',metric,listname,'prots.png'))
  
}
