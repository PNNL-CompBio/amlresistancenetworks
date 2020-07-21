##run through drug resistnace analysis


library(amlresistancenetworks)
library(dplyr)

##first run gilteritinib data
gilt.data<-querySynapseTable('syn22156807')%>%
  filter(!is.na(Gene))
#  readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(gilt.data$Gene)


early.data<-gilt.data%>%
  subset(treatment%in%(c('None','Early Gilteritinib')))%>%
  dplyr::select(Gene,Sample,cellLine,ligand,value)%>%
  dplyr::rename(treatment='ligand')

late.data<-gilt.data%>%
  subset(treatment%in%(c('None','Late Gilteritinib')))%>%
  dplyr::select(Gene,Sample,cellLine,ligand,value)%>%
  dplyr::rename(treatment='ligand')


plotDataByLigand<-function(fulldata,ligand='FLT3',treatment='Early',doNetworks=FALSE){
  
  
  for(cl in c("MOLM14",'MV411')){
    data <-subset(fulldata,cellLine==cl)  
    total.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(data,control='None',conditions=c("FLT3","FGF2"))
  
  ##first plot FGF2
  genes.with.values=total.mean.diffs%>%
    ungroup()%>%
    subset(Condition==ligand)%>%
    dplyr::select(Gene,value=condition_to_control,p_adj)
  
  total.lig=amlresistancenetworks::computeGSEA(genes.with.values,prefix=paste(treatment,'Treatment',cl,ligand,sep='_'))
  
  if(doNetworks){
    lig.network<-computeProteinNetwork(sig.vals=subset(genes.with.values,p_adj<0.05),
                                     all.vals=genes.with.values,nrand=1000)
    RCy3::createNetworkFromIgraph(lig.network,title=paste(treatment,'Treatment',cl,ligand,'Network'))
    lig.network
  }else{
    return(NULL)
  }
  }
}

fgf2.net=plotDataByLigand(early.data,"FGF2",treatment='Early')
fl.net=plotDataByLigand(early.data,"FLT3",treatment='Early')
late.fgf2.net=plotDataByLigand(late.data,"FGF2",treatment='Late')
late.fl.net=plotDataByLigand(late.data,"FLT3",treatment='Late')


##now do networks

#molm.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(early.data, cellLine=='MOLM14'),
#                                                               control='None',
#                                                               conditions=c("FL","FGF2"))
#mv411.mean.diffs<-amlresistancenetworks::computeFoldChangePvals(subset(early.data, cellLine=='MV411'),
#                                                                control='None',
#                                                                conditions=c("FL","FGF2"))

#quiz.data<-readRDS()