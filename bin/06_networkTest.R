
library(amlresistancenetworks)
library(dplyr)
library(RCy3)
gilt.data<-readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
gilt.pdat<-readRDS(system.file('giltPhosphoData.Rds',package='amlresistancenetworks'))



##first, test out combining bulk and regular proteomics
createCombinedGraph<-function(bulk.data,
                              phospho.data,control='None',
                              condition=c('Early Gilteritinib','Late Gilteritinib'),
                              cellLine='',
                              pvalThresh=0.01,prefix=''){
  total.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(bulk.data,cellLine==cellLine),
                                                                  control,condition)
  phospho.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(phospho.data,cellLine==cellLine),
                                                                    control,condition)
 
  
  all.nets=purrr::map(condition,function(cond){
    diff.res<-total.mean.diffs%>%
      ungroup()%>%
      subset(Condition==cond)%>%
      dplyr::select(Gene,value=condition_to_control,p_adj,Condition)
    
    diff.phos<-phospho.mean.diffs%>%
      ungroup()%>%
      subset(Condition==cond)%>%
      dplyr::select(Gene,value=condition_to_control,p_adj,Condition)
    
    overlapping<-intersect(unlist(subset(diff.phos,p_adj<pvalThresh)%>%dplyr::select(Gene)),unlist(subset(diff.res,p_adj<pvalThresh)%>%dplyr::select(Gene)))
  print(overlapping)
    print(paste('removing',length(overlapping),'genes in both phospho and bulk data'))
    
    sig.res<-diff.res%>%
        subset(p_adj<pvalThresh)%>%
        subset(!Gene%in%overlapping)
    sig.phos<-diff.phos%>%
      subset(p_adj<pvalThresh)%>%
      subset(!Gene%in%overlapping)
    
    lig.network<-computeProteinNetwork(sig.vals=rbind(sig.res,sig.phos),
                                       all.vals=diff.res,phos.vals=diff.phos,nrand=100)
    #RCy3::createNetworkFromIgraph(lig.network,title=paste(prefix,cellLine,
    #                                                      gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
    
    return(lig.network)
  }) 
  
  return(all.nets)
}

o.res<-createCombinedGraph(bulk.data=gilt.data,phospho.data=gilt.pdat,prefix='combined',pvalThresh = 0.01,cellLine='MOLM14')

v.res<-createCombinedGraph(bulk.data=gilt.data,phospho.data=gilt.pdat,prefix='combined',pvalThresh = 0.01,cellLine='MV411')

##test out multinetwork approach

createMultiGraph<-function(data,control='None',
                           condition=c('Early Gilteritinib','Late Gilteritinib'),
                           cellLine=''){
  
  
  total.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(data,cellLine==cellLine),
                                                                control,condition)

##iterate over all conditions
all.nets=purrr::map(condition,function(cond){
  diff.res<-total.mean.diffs%>%
    ungroup()%>%
    subset(Condition==cond)%>%
    dplyr::select(Gene,value=condition_to_control,p_adj,Condition)
  
  
    lig.network<-computeProteinNetwork(sig.vals=subset(diff.res,p_adj<0.05),
                                       all.vals=diff.res,nrand=1000)
  
  return(lig.network)
    })


    }