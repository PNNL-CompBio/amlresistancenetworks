##tests differential expression of proteins

library(amlresistancenetworks)
library(dplyr)
library(ggplot2)

syn <-synapseLogin()
sens.data<-syn$tableQuery("select * from  syn22156810")$asDataFrame()

# this is the sensitivity data
#sens.data<-readRDS(system.file('gilteritinibSensitivityData.Rds',package='amlresistancenetworks'))
#prot.univ<-unique(sens.data$Gene)

#'
#' helper function to assign values to sensitive or resistant depending on a
#' threshold provided
#' @param x is a vector of values
#' @param sens is the value below which things are deemed sensitive
#' @param res is the value above which cells are deemed resistant
assign<-function(x,sens,res){
    sapply(x,function(y) ifelse(y<sens,'Sensitive',
                         ifelse(y>res,'Resistant','Other')))
}

#'
#' helper function to assign values to sensitive or resistant depending on a
#' threshold provided
#' @param x is a vector of values
#' @param sens is the value below which things are deemed sensitive
#' @param res is the value above which cells are deemed resistant
assignByQuant<-function(x,t.sens,t.res){
    qs<-quantile(x,probs=c(t.sens,t.res))
   print(qs)
    result=sapply(x,function(y) ifelse(y<qs[[1]],'Sensitive',
                          ifelse(y>qs[[2]],'Resistant','Other')))
    print(paste("found",length(which(result=='Sensitive')),'sensitive and',length(which(result=='Resistant')),'resistant out of',length(x)))
 return(result)         
}


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



#' assignByVal
#' @param sens.data sensitivity data
#' @param metric such as AUc
#' @param sens threshold
#' @param res threshold
assignByVal<-function(sens.data,metric='AUC',r.val=100,s.val=100){
  
  res.samps<-sens.data%>%
    subset(Metric=='AUC')%>%
    dplyr::select(`AML sample`,Condition,Value)%>%
    distinct()%>%
    group_by(Condition)%>%
    mutate(Status=assign(Value,s.val,r.val))%>%
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

##first let's plot the AUC values and see what it looks lke
sens.data%>%
  dplyr::select(`AML sample`,Metric,Condition,Value)%>%
  subset(Metric=='AUC')%>%
  distinct()%>%
  ggplot(aes(x=reorder(`AML sample`,Value),y=Value,fill=Condition))+
    geom_bar(stat='identity',position='dodge')+
    geom_hline(yintercept = 100)+
  #  geom_hline(yintercept = 150,linetype='dashed')+
  #  geom_hline(yintercept = 50,linetype='dashed') +
    ggtitle('AUC values by sample')

sens.data%>%
  dplyr::select(`AML sample`,Metric,Condition,Value)%>%
  subset(Metric=='IC50')%>%
  distinct()%>%
  ggplot(aes(x=reorder(`AML sample`,Value),y=Value,fill=Condition))+
  geom_bar(stat='identity',position='dodge')+
#  geom_hline(yintercept = 100)+
#  geom_hline(yintercept = 150,linetype='dashed')+
#  geom_hline(yintercept = 50,linetype='dashed')+
  ggtitle('IC50 values by sample')


computeDiffExByVal<-function(sens.data,sens.val,res.val){
  samps<-assignByVal(sens.data,'AUC',sens.val,res.val)
  
  new.df<-sens.data%>%
    subset(Metric=='AUC')%>%
    select('AML sample',Gene,value="LogFoldChange")%>%
    distinct()%>%
    left_join(samps)
  
  ##compute fol change and p-values for each drug
  result<-purrr::map_df(unique(new.df$Condition), function(cond){
    # print(cond)
    new.df%>%subset(Condition==cond)%>%
      rename(cellLine='Condition',treatment='Status',Sample='AML sample')%>%
      group_by(cellLine)%>%
      computeFoldChangePvals(control='Sensitive',conditions =c("Resistant"))%>%
      mutate(Drug=cond)
  })
  
  table(result%>%group_by(Drug)%>%subset(p_adj<0.1)%>%summarize(sigProts=n()))
  
  result
}

computeDiffExByThresh<-function(sens.data,sens.val,res.val){
  samps<-assignSensResSamps(sens.data,'AUC',sens.val,res.val)

    new.df<-sens.data%>%
      subset(Metric=='AUC')%>%
      select('AML sample',Gene,value="LogFoldChange")%>%
      distinct()%>%
      left_join(samps)

  ##compute fol change and p-values for each drug
  result<-purrr::map_df(unique(new.df$Condition), function(cond){
   # print(cond)
      new.df%>%subset(Condition==cond)%>%
      rename(cellLine='Condition',treatment='Status',Sample='AML sample')%>%
      group_by(cellLine)%>%
      computeFoldChangePvals(control='Sensitive',conditions =c("Resistant"))%>%
      mutate(Drug=cond)
  })
  
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



computeMutInf<-function(sens.data,matric='AUC'){
  library(infotheo)
  gmat<-sens.data%>%
    dplyr::select(`AML sample`,LogFoldChange,names_from='AML sample',values_fn = list(LogFoldChange=mean))%>%
    subset(!is.na(Gene))%>%
    tibble::column_to_rownames('Gene')
   
}


auc.cor=computeCorVals(sens.data,metric='AUC')%>%mutate(Gene=unlist(Gene),metric='AUC')
#write.table(auc.cor,file='correlationWithAUC.csv',sep=',')

ic50.cor=computeCorVals(sens.data,metric='IC50')%>%mutate(Gene=unlist(Gene),metric='IC50')
#write.table(ic50.cor,file='correlationWithIC50.csv',sep=',')

total.cor<-rbind(auc.cor,ic50.cor)

reg1=computeDiffExByVal(sens.data,sens.val=100,res.val=100)
#write.table(reg1,file='results/singleComboResistance/diffExProteinsBySensRes.csv',sep=',')

##now we can experiment with different thresholdls
#reg2=computeDiffExByThresh(sens.data,sens=75,res=125)
#reg=computeDiffExByThresh(sens.data,sens=100,res=100)
reg5=computeDiffExByThresh(sens.data,sens.val=0.5,res.val=0.5)
reg4=computeDiffExByThresh(sens.data,sens.val=.4,res.val=.6)
reg3=computeDiffExByThresh(sens.data,sens.val=.30,res.val=.7)
reg2=computeDiffExByThresh(sens.data,sens.val=.2,res.val=.8)
#volcanoPlot(reg3)
library(cowplot)

cowplot::plot_grid(volcanoPlot(reg5,title='50/50 split'),volcanoPlot(reg4,title='40/60 split'),
        volcanoPlot(reg3,title='30/70 split'),volcanoPlot(reg2,title='20/80 split'))
ggsave('variousSplitsToDetermineDifferential.png')



#'plotByCor
#'
plotByCor<-function(result,auc.cor,title=""){
  library(ggplot2)
  result%>%left_join(rename(auc.cor,Drug='Condition'),by=c('Gene','Drug'))%>%
    ggplot(aes(x=condition_to_control,y=corVal,col=Drug))+geom_point()+ggtitle(title)
  
}

cowplot::plot_grid(plotByCor(reg5,auc.cor,title='50/50 split'),plotByCor(reg4,auc.cor,title='40/60 split'),
                   plotByCor(reg3,auc.cor,title='30/70 split'),plotByCor(reg2,auc.cor,title='20/80 split'))
ggsave('variousSplitsToDetermineCorrelation.png')

corstats=auc.cor%>%select(Gene,corVal,Condition)

doGSEA<-FALSE
if(doGSEA){
corstats%>%group_by(Gene)%>%
    summarize(value=mean(corVal),med=median(corVal))%>%
    arrange(value)%>%computeGSEA('meanExpressionAUCcor')

corstats%>%subset(Condition=='Venetoclax')%>%
  dplyr::rename(value=corVal)%>%arrange(value)%>%
  computeGSEA('VenetoclaxCorrelation')
  
corstats%>%subset(Condition=='Dora+Ven')%>%
  dplyr::rename(value=corVal)%>%arrange(value)%>%
  computeGSEA('DoraVenCorrelation')


corstats%>%subset(Condition=='Doramapimod')%>%
  dplyr::rename(value=corVal)%>%arrange(value)%>%
  computeGSEA('DoramapimodCorrelation')
}
plotProtsByMetric(sens.data,genelist=c("BCL2","BCL2A1","MCL1"),
                  metric = 'AUC',listname='Venetoclax Targets')

plotProtsByMetric(sens.data,genelist=c("IL1B","IL1RAP","MAPK14"),
                  metric = 'AUC',listname='Immune proteins')
plotProtsByMetric(sens.data,genelist=c("IL1B","IL1RAP","MAPK14"),
                  metric = 'IC50',listname='Immune proteins')