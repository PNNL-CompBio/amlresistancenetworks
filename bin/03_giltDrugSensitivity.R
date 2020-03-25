##lets compare the senstitivity of ex-vivo assays to drugs and drug combinations


library(amlresistancenetworks)
library(dplyr)

##first run gilteritinib data
sens.data<-readRDS(system.file('gilteritinibSensitivityData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(sens.data$Gene)


#


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
