##time course analysis

library(amlresistancenetworks)
library(ggplot2)
library(dplyr)


updateControls<-function(prot.data){
  treats=subset(prot.data,!treatment%in%c('no treatment','No treatment'))
  conts=subset(prot.data,treatment%in%c('no treatment','No treatment'))%>%
    mutate(timePoint='0 hr')
  ret=rbind(treats,conts)%>%mutate(LogFoldChange=as.numeric(as.character(LogFoldChange)))
  
  return(ret)}

protSample=readRDS(system.file('timeCourseData.Rds',package='amlresistancenetworks'))%>%
  updateControls()
phosSample=readRDS(system.file('timeCoursePhosphoData.Rds',package='amlresistancenetworks'))%>%
  updateControls()


doDiffEx<-function(prot.data,cellLine,control,treatment){
  
}

clusterProtsAcrossTreatments<-function(prot.data){
  ##group data by treatment, cell line
  #cluster genes together, how many clusters do we get? 
  
  #return prots and clusters
  clusts<-prot.data%>%dplyr::select(-c(Sample))%>%
    subset(!is.na(treatment))%>%
    subset(!treatment%in%c('no treatment','No treatment'))%>%
    group_by(treatment,cellLine)%>%
    group_modify(~ getCoClusteredProteins(.x),keep=TRUE)
  return(clusts)
}



doEnrichmentOfClusters<-function(allClusts){

  res<-allClusts%>%
    dplyr::select(cellLine,treatment,cluster,Gene)%>%
    group_by(cellLine,treatment,cluster)%>%
    dplyr::group_modify(~ doRegularGo(.x$Gene))
  
}


doClust=FALSE
if(doClust){
  allClusts<-clusterProtsAcrossTreatments(protSample)
  clustEnrich<-doEnrichmentOfClusters(allClusts)
  clustSum=clustEnrich%>%group_by(Description,cellLine,treatment)%>%summarize(nClusts=n_distinct(cluster))
  ggsave('DistinctTermsPerCluster.png')
  
  ##which terms are unique to a specific cluster
  clustSum%>%ggplot()+geom_histogram(aes(x=nClusts,fill=treatment),position='dodge')
library(ggaluvial)
  ##not really impressed by these, let's try something else. 
}


doDiffEx=TRUE
if(doDiffEx){
  #lets try this!
  cellLines=c('CMK','HL60','K562','MOLM')
  dres<-protSample%>%    
    subset(!is.na(treatment))%>%
    rename(treatment='drug')%>%
    rename(timePoint='treatment')%>%
    rename(LogFoldChange='value')%>%
    mutate(Gene=as.character(Gene))
  
  cl.res=cellLines%>%
    purrr::map_df(~computeFoldChangePvals(subset(dres,cellLine==.x),control='0 hr',
                                                  conditions=c("30 min","3 hr","16 hr")))
                                     
}