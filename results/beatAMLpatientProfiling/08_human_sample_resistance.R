##tests differential expression of proteins

library(amlresistancenetworks)
library(dplyr)
library(ggplot2)
# this is the sensitivity data
pat.data<-readRDS(system.file('patientMolecularData.Rds',package='amlresistancenetworks'))
pat.phos<-readRDS(system.file('patientPhosphoSampleData.Rds',package='amlresistancenetworks'))
pat.drugClin<-readRDS(system.file('patientDrugAndClinical.Rds', package='amlresistancenetworks'))



##' volcanoPlot
volcanoPlot<-function(result,title=""){
  library(ggplot2)
  ggplot(result,aes(x=condition_to_control,
                 y=-log10(p_val),shape=p_adj<0.05))+geom_point(aes(color=Drug))+ggtitle(title)
}


drugMolRegression<-function(){
  
  
}

drugMolCorrelation<-function(){
  
}

computeDiffExByDrug<-function(sens.data){
  #samps<-assignSensResSamps(sens.data,'AUC',sens.val,res.val)

  result<-sens.data%>%
    dplyr::select('Sample',Gene,cellLine='Drug',value="LogFoldChange",treatment='Status')%>%
    distinct()%>%
    subset(!is.na(cellLine))%>%
    group_by(cellLine)%>%
    group_modify(~ computeFoldChangePvals(.x,control=NA,conditions =c("Sensitive","Resistant")),keep=TRUE)%>%
    rename(Drug='cellLine')
    
  
  table(result%>%group_by(Drug)%>%subset(p_adj<0.1)%>%summarize(sigProts=n()))
  
  result
}


computeAUCCorVals<-function(clin.data,mol.data,mol.feature){
  tdat<-mol.data%>%select(Gene,`AML sample`,Mol=mol.feature)%>%
    subset(!is.na(Mol))%>%left_join(clin.data,by='AML sample')
  
  dcors<-tdat%>%select(Gene,Mol,Condition,AUC)%>%
      distinct()%>%
      group_by(Gene,Condition)%>%
      mutate(numSamps=n(),drugCor=cor(Mol,AUC))%>%
      select(Gene,Condition,numSamps,drugCor)%>%
      arrange(desc(drugCor))%>%subset(numSamps>10)
  
  ggplot(dcors,aes(x=drugCor))+geom_histogram(aes(fill=Condition))
  stats<-subset(dcors,abs(drugCor)>0.85)%>%
      ungroup()%>%
      group_by(Condition)%>%
      summarize(CorVals=n_distinct(Gene))
      
  return(dcors)
  
}
auc.dat<- subset(pat.drugClin,Metric%in%c('auc','AUC'))%>%
  select(`AML sample`,gender,ageAtDiagnosis,vitalStatus,overallSurvival,Condition,AUC='Value')

#auc.cor=computeCorVals(sens.data,metric='AUC')
#ic50.cor=computeCorVals(sens.data,metric='IC50')

