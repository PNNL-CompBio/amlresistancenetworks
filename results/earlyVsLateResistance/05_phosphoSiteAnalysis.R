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
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([sty])", ""))%>%
    mutate(site=unlist(site))
  
  gene.to.site$Gene=unlist(gene.to.site$Gene)
    
  total.mean.diffs<-data%>%
    dplyr::filter(cellLine==!!cellLine)%>%
    dplyr::select(Gene=site,treatment,Sample,value,cellLine)%>%
    amlresistancenetworks::computeFoldChangePvals(.,control,condition)%>%
    mutate(site=unlist(as.character(Gene)))%>%select(-Gene)%>%left_join(gene.to.site,by='site')
  
  print(head(total.mean.diffs))  
  
  ##first plot FGF2
  genes.with.values=purrr::map_df(condition,function(cond){
    diff.res<-total.mean.diffs%>%
      ungroup()%>%
      subset(Condition==cond)%>%
      dplyr::select(Gene,value=condition_to_control,p_adj,Condition,residue,Peptide)
    
    total.lig=amlresistancenetworks::computeKSEA(diff.res,
                                                 prefix=paste(prefix,cellLine,gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
    
  print(total.lig)
        if(doNetworks){
      lig.network<-computeProteinNetwork(sig.vals=subset(diff.res,p_adj<0.005)%>%subset(abs(value)>1),
                                         all.vals=diff.res,nrand=1000,beta=.5)
      RCy3::createNetworkFromIgraph(lig.network,title=paste(prefix,cellLine,
                                                            gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
    }
    total.lig%>%mutate(Condition=paste0(cond,'_',cellLine))      
  })
  
  return(genes.with.values%>%mutate(Control=control,cellLine=cellLine))
}


#' build networks from data frame
#' @param data.res
#' @param gene.col
#' @param weight.col
#' @param condition.col
#' @return network list?
runNetworksFromDF<-function(data,gene.col='Kinase.Gene',
                            weight.col='aveSubstrateLog2FC',
                            condition.col='Condition',extra.col=c('Substrate.Gene','Source','log2FC'),
                            signif=0.05){
  res = data%>%
    # dplyr::select(cond=condition.col,value=weight.col,Gene=gene.col,p.value)%>%
    mutate(signif=p.value<signif)%>%
    dplyr::select(c(condition.col,weight.col,gene.col,'signif',extra.col))%>%distinct()%>%
    dplyr::rename(cond=condition.col,value=weight.col,Gene=gene.col)%>%
    group_by(cond)%>%
    dplyr::select(c('cond','Gene','value',extra.col,'signif'))%>%
    group_map(~ amlresistancenetworks::computeProteinNetwork(.x),.keep=TRUE)
  return(res)
}


dn=FALSE
#molm14
m.vs.parental<-plotPhosphoByCondition(gilt.pdat,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MOLM14',doNetwork=dn)

m.late.early<-plotPhosphoByCondition(gilt.pdat,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MOLM14',doNetwork=dn)
#runNetworksFromDF(m.vs.parental)
runNetworksFromDF(m.late.early)

if(FALSE){

#mv411
v.vs.parental<-plotPhosphoByCondition(gilt.pdat,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MV411',doNetwork=dn)

v.late.early<-plotPhosphoByCondition(gilt.pdat,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MV411',doNetwork=dn)
runNetworksFromDF(v.vs.parental)
runNetworksFromDF(v.late.early)
dn=FALSE


sapply(c("FGF2","FLT3"),function(lig){
  m.early.net<-plotPhosphoByCondition(subset(gilt.pdat,ligand%in%c('None',lig)),
                                   control='None',
                                   condition=c('Early Gilteritinib','Late Gilteritinib'),
                                   cellLine='MOLM14',doNetwork=dn,prefix=lig)
  v.early.net<-plotPhosphoByCondition(subset(gilt.pdat,ligand%in%c('None',lig)),
                                   control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),
                                   cellLine='MV411',doNetwork=dn,prefix=lig)
  #m.early.net<-plotPhosphoByCondition(subset(gilt.pdat,ligand%in%c('None',lig)),
  #                                    control='Early Gilteritinib',
  #                                    condition=c('Late Gilteritinib'),
  #                                    cellLine='MOLM14',doNetwork=dn,prefix=lig)
  #v.early.net<-plotPhosphoByCondition(subset(gilt.pdat,ligand%in%c('None',lig)),
  #                                    control='Early Gilteritinib',
  #                                    condition=c('Late Gilteritinib'),
  #                                    cellLine='MV411',doNetwork=dn,prefix=lig)
})

}