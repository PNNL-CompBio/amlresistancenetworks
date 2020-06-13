##phospho site analysis

library(amlresistancenetworks)
gilt.pdat<-querySynapseTable('syn22156809')#readRDS(system.file('giltPhosphoData.Rds',package='amlresistancenetworks'))


library(KSEAapp)
library(dplyr)

#'this is our main analysis function  
#'it computes differential expression betwee control and various conditions
#'runs GSEA, plots, and creates a network
plotPhosphoByCondition<-function(data,control='None',condition=c('Early Gilteritinib'),
                              cellLine='MOLM14',doNetworks=TRUE,prefix='phos'){
  library(KSEAapp)
  
  gene.to.site<-dplyr::select(data,Gene,site,Peptide)%>%distinct()%>%
      dplyr::mutate(residue=stringr::str_replace(site,paste0(Gene,'-'),''))%>%
      dplyr::mutate(residue=stringr::str_replace_all(residue,"([STY])", ";\\1"))%>%
    dplyr::mutate(residue=stringr::str_replace(residue,"^;", ""))%>%
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([sty])", ""))
    
  total.mean.diffs<-data%>%
    dplyr::filter(cellLine==!!cellLine)%>%
    dplyr::select(Gene=site,treatment,Sample,value,cellLine)%>%
    amlresistancenetworks::computeFoldChangePvals(.,control,condition)%>%
    dplyr::rename(site=Gene)%>%left_join(gene.to.site)
  
  print(head(total.mean.diffs))  
  
  ##first plot FGF2
  genes.with.values=purrr::map_df(condition,function(cond){
    diff.res<-total.mean.diffs%>%
      ungroup()%>%
      subset(Condition==cond)%>%
      dplyr::select(Gene,value=condition_to_control,p_adj,Condition,residue,Peptide)
    total.lig=amlresistancenetworks::computeKSEA(diff.res,
                                                 prefix=paste(prefix,cellLine,gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
    
    if(doNetworks){
      lig.network<-computeProteinNetwork(sig.vals=subset(diff.res,p_adj<0.005)%>%subset(abs(value)>1),
                                         all.vals=diff.res,nrand=1000,beta=.5)
      RCy3::createNetworkFromIgraph(lig.network,title=paste(prefix,cellLine,
                                                            gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
    }
    diff.res%>%mutate(Condition=cond)      
  })
  
  return(genes.with.values%>%mutate(Control=control,cellLine=cellLine))
}


dn=FALSE
#molm14
m.vs.parental<-plotPhosphoByCondition(gilt.pdat,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MOLM14',doNetwork=dn)

m.late.early<-plotPhosphoByCondition(gilt.pdat,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MOLM14',doNetwork=dn)


#mv411
v.vs.parental<-plotPhosphoByCondition(gilt.pdat,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MV411',doNetwork=dn)

v.late.early<-plotPhosphoByCondition(gilt.pdat,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MV411',doNetwork=dn)

dn=FALSE


sapply(c("FGF2","FLT3"),function(lig){
#  m.early.net<-plotPhosphoByCondition(subset(gilt.pdat,ligand%in%c('None',lig)),
 #                                  control='None',
 #                                  condition=c('Early Gilteritinib','Late Gilteritinib'),
#                                   cellLine='MOLM14',doNetwork=dn,prefix=lig)
#  v.early.net<-plotPhosphoByCondition(subset(gilt.pdat,ligand%in%c('None',lig)),
 #                                  control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),
#                                   cellLine='MV411',doNetwork=dn,prefix=lig)
  m.early.net<-plotPhosphoByCondition(subset(gilt.pdat,ligand%in%c('None',lig)),
                                      control='Early Gilteritinib',
                                      condition=c('Late Gilteritinib'),
                                      cellLine='MOLM14',doNetwork=dn,prefix=lig)
  v.early.net<-plotPhosphoByCondition(subset(gilt.pdat,ligand%in%c('None',lig)),
                                      control='Early Gilteritinib',
                                      condition=c('Late Gilteritinib'),
                                      cellLine='MV411',doNetwork=dn,prefix=lig)
})