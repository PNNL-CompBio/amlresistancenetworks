##tests differential expression of proteins

library(amlresistancenetworks)
library(dplyr)
library(ggplot2)
# this is the sensitivity data
sens.data<-readRDS(system.file('patientProtSampleData.Rds',package='amlresistancenetworks'))


    #' assignSensResSamps
    #' @param sens.data sensitivity data
    #' @param metric such as AUc
    #' @param sens threshold
    #' @param res threshold
assignSensResSamps<-function(sens.data,metric='AUC',t.sens=0.25,t.res=0.75){

    res.samps<-sens.data%>%
        subset(Metric=='AUC')%>%
        dplyr::select(`AML sample`,Condition,Value)%>%
        distinct()%>%
        group_by(Condition)%>%
        mutate(Status=assignByQuant(Value,t.sens,t.res))%>%
       #dplyr::group_modify(~ assignByQuant(.x$Value,sens,res))
        select(`AML sample`,Condition,Status)
 #   print(res.samps)
    res.samps
}


##' volcanoPlot
volcanoPlot<-function(result,title=""){
  library(ggplot2)
  ggplot(result,aes(x=condition_to_control,
                 y=-log10(p_val),shape=p_adj<0.05))+geom_point(aes(color=Drug))+ggtitle(title)
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


computeCorVals<-function(sens.data,metric='AUC'){
  cors=sens.data%>%subset(Metric==metric)%>%
    group_by(Condition,Gene)%>%
    dplyr::select(LogFoldChange,Value)%>%
    summarize(corVal=cor(LogFoldChange,Value,method='spearman'),corSig=cor.test(LogFoldChange,Value,method='spearman')$p.value)%>%
    mutate(p_adj=p.adjust(corSig))
  
  return(cors)
  
}
sens.data%>%group_by(Drug,Status)%>%summarize(n_distinct(Sample))
res<-computeDiffExByDrug(sens.data)


#auc.cor=computeCorVals(sens.data,metric='AUC')
#ic50.cor=computeCorVals(sens.data,metric='IC50')

