##tests differential expression of proteins 

library(amlresistancenetworks)
library(dplyr)

sens.data<-readRDS(system.file('gilteritinibSensitivityData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(sens.data$Gene)

assign<-function(x,sens,res){
  sapply(x,function(y) ifelse(y<sens,'Sensitive',ifelse(y>res,'Resistant','Other')))
}

assignSensResSamps<-function(sens.data,metric='AUC',sens=100,res=100){
  
  res.samps<-sens.data%>%
    subset(Metric=='AUC')%>%
    dplyr::select(`AML sample`,Condition,Value)%>%
    distinct()%>%
    group_by(Condition)%>%rowwise()%>%
    mutate(Status=assign(Value,sens,res))%>%
    select(`AML sample`,Condition,Status)
  
  res.samps
}



volcanoPlot<-function(res){
  library(ggplot2)
  ggplot(res,aes(x=condition_to_control,
                 y=-log10(p_value),shape=p_adj<0.05))+geom_point(aes(color=Drug))
  
  
}

##first let's plot the AUC values and see what it looks lke
sens.data%>%
  dplyr::select(`AML sample`,Metric,Condition,Value)%>%
  subset(Metric=='AUC')%>%
  distinct()%>%
  ggplot(aes(x=`AML sample`,y=Value,fill=Condition))+
    geom_bar(stat='identity',position='dodge')+
    geom_hline(yintercept = 100)+
    geom_hline(yintercept = 150,linetype='dashed')+
    geom_hline(yintercept = 50,linetype='dashed')+
    ggtitle('AUC values by sample')

sens.data%>%
  dplyr::select(`AML sample`,Metric,Condition,Value)%>%
  subset(Metric=='IC50')%>%
  distinct()%>%
  ggplot(aes(x=`AML sample`,y=Value,fill=Condition))+
  geom_bar(stat='identity',position='dodge')+
 # geom_hline(yintercept = 100)+
#  geom_hline(yintercept = 150,linetype='dashed')+
#  geom_hline(yintercept = 50,linetype='dashed')+
  ggtitle('IC50 values by sample')


computeDiffExByThresh<-function(sens.data,sens,res){
  samps<-assignSensResSamps(sens.data,sens,res)

    new.df<-sens.data%>%
      subset(Metric=='AUC')%>%
      select('AML sample',Gene,value="LogFoldChange")%>%
      distinct()%>%
      left_join(samps)

  ##compute fol change and p-values for each drug
  result<-purrr::map_df(unique(new.df$Condition), function(cond){
    print(cond)
      new.df%>%subset(Condition==cond)%>%
      rename(cellLine='Condition',treatment='Status')%>%
      group_by(cellLine)%>%
      computeFoldChangePvals(control='Sensitive',conditions =c("Resistant"))%>%
      mutate(Drug=cond)
  })
  print(result%>%group_by(Drug)%>%subset(p_adj<0.05)%>%summarize(sigProts=n()))
  result
}

##now we can experiment with different thresholdls
reg2=computeDiffExByThresh(sens.data,sens=75,res=125)
reg=computeDiffExByThresh(sens.data,sens=100,res=100)
reg3=computeDiffExByThresh(sens.data,sens=50,res=150)
volcanoPlot(reg3)

plotProtsByMetric(sens.data,genelist=c("HLA-A","HLA-B","HLA-DRB1"),metric = 'AUC',listname='Immune Prots')
  
  