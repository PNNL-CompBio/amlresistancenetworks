

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



#' plotall AUCS with clinical data
#' what are we measuring here
#' @param auc.data
#' @import pheatmap
#' @export
plotAllAUCs<-function(auc.data){
  auc.data%>%
    dplyr::select(`AML sample`,Condition,percAUC)%>%distinct()%>%
    ggplot(aes(x=percAUC,fill=Condition))+geom_histogram()+ theme(legend.position = "none")+
    ggtitle("Distribution of AUC values by drug")
  ggsave('percAUCdist.pdf')
  synapseStore('percAUCdist.pdf','syn22130776')
  auc.data%>%
    dplyr::select(`AML sample`,Condition,AUC)%>%distinct()%>%
    ggplot(aes(x=AUC,fill=Condition))+geom_histogram()+ theme(legend.position = "none")+
    ggtitle("Distribution of raw AUC values by drug")
  ggsave('AUCdist.pdf')
  synapseStore('AUCdist.pdf','syn22130776')
  
  library(pheatmap)
  pat.vars<-auc.data%>%
    dplyr::select(`AML sample`,gender,ageAtDiagnosis,vitalStatus,overallSurvival)%>%
    distinct()%>%
    tibble::column_to_rownames("AML sample")
  drug.vars<-auc.data%>%
    dplyr::select(Drug='Condition',family)%>%
    distinct()
  auc.mat<-auc.data%>%
    dplyr::select(`AML sample`,Drug='Condition',percAUC)%>%
    distinct()%>%
    tidyr::pivot_wider(names_from="AML sample",values_from="percAUC",values_fill=(list(percAUC=0)),
                       values_fn=list(percAUC=mean))%>%
    tibble::column_to_rownames('Drug')%>%
    as.matrix()
  
  pheatmap(auc.mat,annotation_col = pat.vars,filename = 'AUCheatmap.pdf',cellwidth = 10,cellheight = 10) 
  synapseStore('AUCheatmap.pdf','syn22130776')
}

