##run through drug resistnace analysis


library(amlresistancenetworks)
library(dplyr)

##first run gilteritinib data
gilt.data<-readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(gilt.data$Gene)


early.data<-gilt.data%>%
  subset(treatment%in%(c('None','Early Gilteritinib')))%>%
  dplyr::select(Gene,Sample,cellLine,treatment,value,ligand)

late.data<-gilt.data%>%
  subset(treatment%in%(c('Early Gilteritinib','Late Gilteritinib')))%>%
  dplyr::select(Gene,Sample,cellLine,treatment,value,ligand)



#'this is our main analysis function  
#'it computes differential expression betwee control and various conditions
#'runs GSEA, plots, and creates a network
plotDataByCondition<-function(data,control='None',condition=c('Early Gilteritinib'),
                              cellLine='MOLM14',doNetworks=TRUE,prefix=''){
  
  total.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(data,cellLine==cellLine),
                                                                  control,condition)
  
  ##iterate over all conditions
  genes.with.values=purrr::map_df(condition,function(cond){
      diff.res<-total.mean.diffs%>%
          ungroup()%>%
          subset(Condition==cond)%>%
          dplyr::select(Gene,value=condition_to_control,p_adj,Condition)
  
      total.lig=amlresistancenetworks::computeGSEA(diff.res,
                                                     prefix=paste(prefix,cellLine,gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
  
    if(doNetworks){
    
      lig.network<-computeProteinNetwork(sig.vals=subset(diff.res,p_adj<0.05),
                                       all.vals=diff.res,nrand=1000)
      RCy3::createNetworkFromIgraph(lig.network,title=paste(prefix,cellLine,
                                                            gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
      }
    return(total.lig%>%mutate(Condition=cond))      
      })
  
  return(genes.with.values%>%mutate(CellLine=cellLine,Control=control))
}

dn=FALSE #no networks for one run


#####here we run the results
#quiz.data<-readRDS()

#molm14
#m.vs.parental<-plotDataByCondition(gilt.data,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MOLM14',doNetwork=dn)

#m.late.early<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MOLM14',doNetwork=dn)


#mv411
#v.vs.parental<-plotDataByCondition(gilt.data,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MV411',doNetwork=dn)

#v.late.early<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MV411',doNetwork=dn)


sapply(c("FGF2","FLT3"),function(lig){
  m.early.net<-plotDataByCondition(subset(gilt.data,ligand%in%c('None',lig)),
                                   control='None',
                                   condition='Early Gilteritinib',
                                   cellLine='MOLM14',doNetwork=dn,prefix=lig)
  v.early.net<-plotDataByCondition(subset(gilt.data,ligand%in%c('None',lig)),
                                   control='None',condition='Early Gilteritinib',
                                   cellLine='MV411',doNetwork=dn,prefix=lig)
})

