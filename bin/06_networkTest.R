##test out multinetwork approach

createMultiGraph<-function(data,control='None',
                           condition=c('Early Gilteritinib','Late Gilteritinib'),
                           cellLine=c('MV411','MOLM14')){
  
  
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