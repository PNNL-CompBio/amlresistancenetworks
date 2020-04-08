##phospho site analysis

library(amlresistancenetworks)
gilt.pdat<-readRDS(system.file('giltPhosphoSensData.Rds',package='amlresistancenetworks'))


library(KSEAapp)
library(dplyr)
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
   #   group_by(cellLine)%>%
      computeFoldChangePvals(control='Sensitive',conditions =c("Resistant"))%>%
      mutate(Drug=cond)
  })
  
  table(result%>%group_by(Drug)%>%subset(p_adj<0.1)%>%summarize(sigProts=n()))
  
  result
}

##' volcanoPlot
volcanoPlot<-function(result,title=""){
  library(ggplot2)
  ggplot(result,aes(x=condition_to_control,
                    y=-log10(p_val),shape=p_adj<0.05))+geom_point(aes(color=Drug))+ggtitle(title)
}

computeCorVals<-function(sens.data,metric='AUC'){
  cors=sens.data%>%subset(Metric==metric)%>%
    subset(!is.na(LogFoldChange))%>%
    subset(Gene!="")%>%
    group_by(Condition,Gene)%>%
    dplyr::select(LogFoldChange,Value)%>%
    summarize(corVal=cor(LogFoldChange,Value,method='spearman',use='pairwise.complete.obs'),
              corSig=cor.test(LogFoldChange,Value,method='spearman',use='pairwise.complete.obs')$p.value)%>%
    mutate(p_adj=p.adjust(corSig))
  
  return(cors)
  
}
gilt.pdat<-dplyr::rename(gilt.pdat,LogFoldChange='value')


gene.to.site<-dplyr::select(gilt.pdat,Gene,site,Peptide)%>%distinct()%>%
  dplyr::mutate(residue=stringr::str_replace(site,paste0(Gene,'-'),''))%>%
  dplyr::mutate(residue=stringr::str_replace_all(residue,"([STY])", ";\\1"))%>%
  dplyr::mutate(residue=stringr::str_replace(residue,"^;", ""))%>%
  dplyr::mutate(residue=stringr::str_replace_all(residue,"([sty])", ""))

#auc.cor<-gilt.pdat%>%
#  dplyr::filter(cellLine==!!cellLine)%>%
 # dplyr::select(Gene=site,treatment,Sample,LogFoldChange)%>%
#  amlresistancenetworks::computeFoldChangePvals(.,control,condition)%>%
#  dplyr::rename(site=Gene)%>%left_join(gene.to.site)

auc.cor=computeCorVals(gilt.pdat,metric='AUC')%>%left_join(gene.to.site)
ic50.cor=computeCorVals(gilt.pdat,metric='IC50')%>%left_join(gene.to.site)

reg5=computeDiffExByThresh(gilt.pdat,sens.val=0.5,res.val=0.5)
reg4=computeDiffExByThresh(gilt.pdat,sens.val=.4,res.val=.6)
reg3=computeDiffExByThresh(gilt.pdat,sens.val=.30,res.val=.7)
reg2=computeDiffExByThresh(gilt.pdat,sens.val=.2,res.val=.8)
#volcanoPlot(reg3)
library(cowplot)

cowplot::plot_grid(volcanoPlot(reg5,title='50/50 split'),volcanoPlot(reg4,title='40/60 split'),
                   volcanoPlot(reg3,title='30/70 split'),volcanoPlot(reg2,title='20/80 split'))
ggsave('variousSplitsToDetermineDifferentialPhos.png')



#'plotByCor
#'
plotByCor<-function(result,auc.cor,title=""){
  library(ggplot2)
  result%>%left_join(rename(auc.cor,Drug='Condition'),by=c('Gene','Drug'))%>%
    ggplot(aes(x=condition_to_control,y=corVal,col=Drug))+geom_point()+ggtitle(title)
  
}

cowplot::plot_grid(plotByCor(reg5,auc.cor,title='50/50 split'),plotByCor(reg4,auc.cor,title='40/60 split'),
                   plotByCor(reg3,auc.cor,title='30/70 split'),plotByCor(reg2,auc.cor,title='20/80 split'))
ggsave('variousSplitsToDetermineCorrelationPhos.png')

corstats=auc.cor
corstats%>%group_by(Gene,site,Peptide,residue)%>%
  summarize(value=mean(corVal),med=median(corVal))%>%
  arrange(value)%>%computeKSEA(prefix='meanExpressionAUCPhoscor')

corstats%>%subset(Condition=='Venetoclax')%>%
  dplyr::rename(value=corVal)%>%arrange(value)%>%
  computeKSEA(prefix='VenetoclaxPhosCorrelation')

corstats%>%subset(Condition=='Dora+Ven')%>%
  dplyr::rename(value=corVal)%>%arrange(value)%>%
  computeKSEA(prefix='DoraVenPhosCorrelation')


corstats%>%subset(Condition=='Doramapimod')%>%
  dplyr::rename(value=corVal)%>%arrange(value)%>%
  computeKSEA(prefix='DoramapimodPhosCorrelation')
