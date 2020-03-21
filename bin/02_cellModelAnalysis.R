##run through drug resistnace analysis


library(amlresistancenetworks)
library(dplyr)

##first run gilteritinib data
gilt.data<-readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(gilt.data$Gene)


early.data<-gilt.data%>%
  subset(treatment%in%(c('None','Early Gilteritinib')))%>%
  dplyr::select(Gene,Sample,cellLine,treatment,value)

late.data<-gilt.data%>%
  subset(treatment%in%(c('Early Gilteritinib','Late Gilteritinib')))%>%
  dplyr::select(Gene,Sample,cellLine,treatment,value)


plotDataByCondition<-function(data,control='None',condition='Early Gilteritinib',
                              cellLine='MOLM14',doNetworks=TRUE){
  
  total.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(data,cellLine==cellLine),
                                                                  control,condition)
                                                                    
  
  ##first plot FGF2
  genes.with.values=total.mean.diffs%>%
    ungroup()%>%
    subset(Condition==condition)%>%
    dplyr::select(Gene,value=condition_to_control,p_adj)
  
  total.lig=computeGSEA(genes.with.values,prot.univ,prefix=paste(cellLine,gsub(' ','',control),'vs',gsub(' ','',condition),sep='_'))
  
  if(doNetworks){

    
    lig.network<-computeProteinNetwork(sig.vals=subset(genes.with.values,p_adj<0.05),
                                       all.vals=genes.with.values,nrand=1000)
    RCy3::createNetworkFromIgraph(lig.network,title=paste(cellLine,gsub(' ','',control),'vs',gsub(' ','',condition),sep='_'))
    lig.network
  }else{
    return(NULL)
  }
}

dn=FALSE #no networks for one run

#molm14
m.early.net<-plotDataByCondition(gilt.data,control='None',condition='Early Gilteritinib',cellLine='MOLM14',doNetwork=dn)

m.late.early<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MOLM14',doNetwork=dn)

m.late.net<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='None',cellLine='MOLM14',doNetwork=dn)


#mv411
v.early.net<-plotDataByCondition(gilt.data,control='None',condition='Early Gilteritinib',cellLine='MV411',doNetwork=dn)

v.late.early<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MV411',doNetwork=dn)

v.late.net<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='None',cellLine='MV411',doNetwork=dn)

##now do networks

#molm.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(early.data, cellLine=='MOLM14'),
#                                                               control='None',
#                                                               conditions=c("FL","FGF2"))
#mv411.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(early.data, cellLine=='MV411'),
#                                                                control='None',
#                                                                conditions=c("FL","FGF2"))


#quiz.data<-readRDS()