##run through drug resistnace analysis


library(amlresistancenetworks)
library(dplyr)

##first run gilteritinib data
gilt.data<-readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(gilt.data$Gene)


early.data<-gilt.data%>%
  subset(treatment%in%(c('None','Early Gilteritinib')))%>%
  dplyr::select(Gene,Sample,CellLine,ligand,value)%>%
  rename(ligand='treatment')

plotDataByLigand<-function(data,ligand='FL'){
  
  total.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(data,control='None',conditions=c("FL","FGF2"))
  
  ##first plot FGF2
  genes.with.values=total.mean.diffs%>%
    ungroup()%>%
    subset(Condition==ligand)%>%
    dplyr::select(Gene,value=condition_to_control,p_adj)
  
  total.lig=computeGSEA(genes.with.values,prot.univ,prefix=paste('Combined',ligand,sep='_'))
  
  
  lig.network<-computeProteinNetwork(sig.vals=subset(genes.with.values,p_adj<0.05),
                                     all.vals=genes.with.values,nrand=1000)
  RCy3::createNetworkFromIgraph(lig.network,title=paste(ligand,'Network'))
  lig.network
}


plotDataByLigand(early.data,"FGF2")
plotDataByLigand(early.data,"FL")

##now do networks

#molm.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(early.data, CellLine=='MOLM14'),
#                                                               control='None',
#                                                               conditions=c("FL","FGF2"))
#mv411.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(early.data, CellLine=='MV411'),
#                                                                control='None',
#                                                                conditions=c("FL","FGF2"))


#quiz.data<-readRDS()