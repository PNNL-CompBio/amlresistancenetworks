##lets compare the senstitivity of ex-vivo assays to drugs and drug combinations


library(amlresistancenetworks)
library(dplyr)

##first run gilteritinib data
sens.data<-readRDS(system.file('gilteritinibSensitivityData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(sens.data$Gene)


#' 
#' Get correlated proteins
#' Computes protein correlations with dose response metric, either AUC or IC50
#' returns list of drugs, proteins and correlation values
#' @param data
#' @param metric
#' @return data.frame
#' @require ggplot2
#' @require dplyr
#' @export
getCorrelatedProteins<-function(data,metric="AUC"){
  cor.stats<-subset(data,Metric==metric)%>%
    group_by(Condition,Gene)%>%
    dplyr::select(`AML sample`,LogFoldChange,Value)%>%distinct()%>%
    mutate(correlation=cor(LogFoldChange,Value))%>%
    dplyr::select(Condition,correlation)%>%
    distinct()
  ggplot2::ggplot(cor.stats)+
    ggplot2::geom_histogram(mapping=ggplot2::aes(x=correlation,fill=Condition),
                            position='dodge')+
    ggplot2::ggtitle(paste('Correlation between protein change and drug',metric))
  ggplot2::ggsave(paste0('drugProtein',metric,'correlation.png'))
  
  univ=unique(cor.stats$Gene)
  sapply(unique(cor.stats$Condition),function(cond){
    subset(cor.stats,Condition==cond)%>%
      dplyr::select(Gene,value='correlation')%>%
      amlresistancenetworks::computeGSEA(prot.univ = univ,prefix=paste(cond,'correlation with',metric))
  })
  
  cor.stats
  
}



toBin<-function(y){
  sapply(y,function(x){
  if(is.na(x)){
    return(0)
  }else if(!x){
    return(0)}
  else{
      return(1)
  }})
}

makeVenn<-function(cors,metric='AUC',thresh=0.5){
  library(ggplot2)
  library(UpSetR)
  posNeg=cors%>%
    mutate(IsPos=correlation>thresh,IsNeg=correlation<(-1*thresh))%>%
    subset(!is.na(Gene))

  pos<-posNeg%>%
    dplyr::select(-c(correlation,IsNeg))%>%
    tidyr::pivot_wider(names_from='Condition',values_from='IsPos',
                        names_prefix = 'Pos_',values_fn=list(IsPos=toBin))
    
  neg<-posNeg%>%  
    dplyr::select(-c(correlation,IsPos))%>%
    tidyr::pivot_wider(names_from='Condition',values_from='IsNeg',
                           names_prefix = 'Neg_',values_fn=list(IsNeg=toBin))
  tot<-pos%>%inner_join(neg)
  pdf(file=paste0('correlatedUpset_',metric,'thresh',thresh,'.pdf'),onefile=FALSE)
  print(upset(as.data.frame(tot),nsets=6,matrix.color='darkred'))
  dev.off()
}

#now we might wantt o examine some proteins
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

##first we compute all correlations
#auc.cors<-getCorrelatedProteins(sens.data,'AUC')
#ic50.cors<-getCorrelatedProteins(sens.data,'IC50')

##then we plot the overlap between the highly and uncorrelated
av<-makeVenn(auc.cors,'AUC',0.5)
iv<-makeVenn(ic50.cors,'IC50',0.5)

##lastly we look at correlated genes individually
sapply(unique(auc.cors$Condition),function(con){
    topAuc=auc.cors%>%
      group_by(Condition)%>%
      top_n(5,wt=correlation)%>%
      subset(Condition==con)%>%
      ungroup()%>%  
      dplyr::select(Gene)%>%unlist()%>%
      plotProtsByMetric(sens.data,genelist=.,metric='AUC',
                        listname=paste0(con,'posCors'))
    
    bottomAuc=auc.cors%>%
      group_by(Condition)%>%
      top_n(5,wt=rev(correlation))%>%
      subset(Condition==con)%>%
      ungroup()%>%  
      dplyr::select(Gene)%>%unlist()%>%
      plotProtsByMetric(sens.data,genelist=.,metric='AUC',
                        listname=paste0(con,'negCors'))
    
    
    topIC=ic50.cors%>%
      group_by(Condition)%>%
      top_n(5,wt=correlation)%>%
      subset(Condition==con)%>%
      ungroup()%>%  
      dplyr::select(Gene)%>%unlist()%>%
      plotProtsByMetric(sens.data,genelist=.,metric='IC50',listname=paste0(con,'posCors'))
    
    
    bottomIC=ic50.cors%>%
      group_by(Condition)%>%
      top_n(5,wt=rev(correlation))%>%
      subset(Condition==con)%>%
      ungroup()%>%  
      dplyr::select(Gene)%>%unlist()%>%
    plotProtsByMetric(sens.data,genelist=.,metric='IC50',listname=paste0(con,'negCors'))
    
})
