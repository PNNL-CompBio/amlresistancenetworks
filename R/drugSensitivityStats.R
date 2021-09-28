

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
#' @param pat.data Patient data frame
#' @param drug.metric 'AUC' is the primary drug metric
#' @param drug.column 'Condition' is default column, but we can also plot by 'family'
#' @import pheatmap
#' @export
plotAllAUCs<-function(auc.data,pat.data,drug.metric='AUC',drug.column='Condition'){

  library(wesanderson)
  
  pal<-wes_palette('Darjeeling1',length(unique(auc.data$Condition)),type='continuous')
    p1<-auc.data%>%
    dplyr::select('AML sample',drugCol=drug.column,drugMetric=drug.metric)%>%distinct()%>%
    ggplot(aes(x=drugMetric,fill=Condition))+geom_histogram()+ theme(legend.position = "none")+
    ggtitle(paste("Distribution of raw",drug.metric,"values by",drug.column))+scale_fill_manual(values=pal)
  print(p1)
  ggsave(paste0(drug.metric,drug.column,'dist.pdf'),p1)
 # synapseStore('AUCdist.pdf','syn22130776')
  
    library(pheatmap)
  dat.summ<-pat.data%>%group_by(`AML sample`)%>%
      summarize(RNA=if_else(all(mRNALevels==0),FALSE,TRUE),
               proteins=if_else(all(proteinLevels==0),FALSE,TRUE),
               mutations=if_else(all(geneMutations==0),FALSE,TRUE))%>%
    mutate(phosphoSite=`AML sample`%in%pat.phos$Sample)%>%
    left_join(select(auc.data,c(`AML sample`))%>%distinct())
  
  print(dat.summ)
  pat.vars<-auc.data%>%
    select(`AML sample`,overallSurvival)%>%
  #  dplyr::select(-c(Condition,percAUC,AUC,medAUC))%>%
#      `AML sample`,gender,ageAtDiagnosis,vitalStatus,overallSurvival)%>%
    distinct()%>%
    left_join(dat.summ)%>%
    tibble::column_to_rownames("AML sample")
  
  annote.colors<-lapply(names(select(pat.vars,-overallSurvival)),function(x) c(`FALSE`='darkgrey',`TRUE`='white'))
  names(annote.colors)<-setdiff(names(pat.vars),'overallSurvival')

  if('RNA'%in%names(pat.vars))
    pat.vars$RNA<-as.factor(pat.vars$RNA)
  if('mutations'%in%names(pat.vars))
    pat.vars$mutations<-as.factor(pat.vars$mutations)
  if('proteins'%in%names(pat.vars))
    pat.vars$proteins<-as.factor(pat.vars$proteins)
  if('phosphoSite'%in%names(pat.vars))
    pat.vars$phosphoSite<-as.factor(pat.vars$phosphoSite)
  
  drug.vars<-auc.data%>%
    dplyr::select(Drug='Condition')%>%
    distinct()
  
  pfn=list(0)
  names(pfn)=drug.metric
  pff=list(mean)
  names(pff)<-drug.metric
  
  auc.mat<-auc.data%>%
    dplyr::select(`AML sample`,Drug=drug.column,!!to.plot)%>%
    distinct()%>%
    tidyr::pivot_wider(names_from="AML sample",values_from=drug.metric,
                       values_fill=pfn,
                       values_fn=pff)%>%
    tibble::column_to_rownames('Drug')%>%
    as.matrix()
  pheatmap(auc.mat,annotation_col = pat.vars,clustering_distance_cols='correlation',
           clustering_method='ward',cellwidth = 10,cellheight = 10, color=pal,annotation_colors=annote.colors) 
  pheatmap(auc.mat,annotation_col = pat.vars,clustering_distance_cols='correlation',
           clustering_method='ward',filename = paste0(drug.metric,drug.column,'heatmap.pdf'),cellwidth = 10,cellheight = 10,
           color=pal,annotation_colors=annote.colors) 
 # synapseStore(paste0(drug.metric,'heatmap.pdf'),'syn22130776')
  return(pat.vars)
}

