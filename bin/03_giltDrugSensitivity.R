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
  ggsave(paste0('drugProtein',metric,'correlation.png'))
  
  univ=unique(cor.stats$Gene)
  sapply(unique(cor.stats$Condition),function(cond){
    subset(cor.stats,Condition==cond)%>%
      dplyr::select(Gene,value='correlation')%>%
      amlresistancenetworks::computeGSEA(prot.univ = univ,prefix=paste(cond,'correlation with',metric))
  })
  
  cor.stats
  
}

auc.cors<-getCorrelatedProteins(sens.data,'AUC')
ic50.cors<-getCorrelatedProteins(sens.data,'IC50')


