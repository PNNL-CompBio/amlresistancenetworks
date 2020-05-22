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

createUpDownTable<-function(diffexRes,lfcThresh=1.0,prefix=''){
  
  
  summ<-diffexRes%>%
    dplyr::select(-c(Control,p_val,p_adj))%>%
    #subset(abs(condition_to_control)>1)
    mutate(direction=if_else(condition_to_control>lfcThresh,'Up',if_else(condition_to_control< (-1*lfcThresh),'Down','No Change')))
  
  ##fix the factors
  summ$direction<-factor(summ$direction,levels=c("Down","No Change","Up"))
  summ$Condition<-factor(summ$Condition,levels=c('30 min','3 hr','16 hr'))
  
  lps<-lapply(unique(diffexRes$cellLine),function(cl){
    rs=subset(summ,cellLine==cl)
    changing<-subset(rs,direction%in%c("Up","Down"))%>%dplyr::select(Gene)

    goterms <- doRegularGo(unique(changing$Gene),prefix)
    write.csv(goterms,file=paste0('results/cellLineTimeCourse/',prefix,cl,'changingGenes.csv'))
      p<-rs%>%subset(Gene%in%changing$Gene)%>%
        ggplot(aes(x=Condition,stratum=direction,alluvium=Gene,fill=direction,label=Condition))+
        geom_flow(stat='alluvium',lode.guidance='frontback')+
        geom_stratum()+
        ggtitle(paste("Changing",prefix,"proteins in",cl))+
        theme_minimal()+
        viridis::scale_fill_viridis(3,discrete=T)
    
    })
   
   cowplot::plot_grid(plotlist=lps,nrow=2)
  
  }
doClust=FALSE
if(doClust){
  allClusts<-clusterProtsAcrossTreatments(protSample)
  clustEnrich<-doEnrichmentOfClusters(allClusts)
  clustSum=clustEnrich%>%group_by(Description,cellLine,treatment)%>%summarize(nClusts=n_distinct(cluster))
  ggsave('DistinctTermsPerCluster.png')
  
  ##which terms are unique to a specific cluster
  clustSum%>%ggplot()+geom_histogram(aes(x=nClusts,fill=treatment),position='dodge')

  ##not really impressed by these, let's try something else. 
}

library(ggalluvial)
doDiffEx=TRUE
if(doDiffEx){
  #lets try this!
  cellLines=c('CMK','HL60','K562','MOLM')
  dres<-protSample%>%    
    subset(!is.na(treatment))%>%
    subset(Gene!="")%>%
    dplyr::rename(drug='treatment')%>%
    dplyr::rename(treatment='timePoint')%>%
    dplyr::rename(value='LogFoldChange')%>%
    mutate(Gene=as.character(Gene))
  
  cl.res=cellLines%>%
    purrr::map_df(~ amlresistancenetworks::computeFoldChangePvals(subset(dres,cellLine==.x),control='0 hr',
                                                  conditions=c("30 min","3 hr","16 hr"),
                                                  doManual=TRUE)%>%
                    mutate(cellLine=.x))
  
  plots<-createUpDownTable(cl.res)
  ggsave('results/cellLineTimeCourse/changingProtsInCellLines.png')
  
  pres<-phosSample%>%
    subset(!is.na(treatment))%>%
    subset(Gene!="")%>%
    dplyr::rename(drug='treatment')%>%
    dplyr::rename(treatment='timePoint')%>%
    dplyr::rename(value='LogFoldChange')%>%
    mutate(Gene=as.character(Gene))
  
  pl.res=cellLines%>%
    purrr::map_df(~ amlresistancenetworks::computeFoldChangePvals(subset(pres,cellLine==.x),control='0 hr',
                                                                  conditions=c("30 min","3 hr","16 hr"),
                                                                  doManual=TRUE)%>%
                    mutate(cellLine=.x))
  
  plots<-createUpDownTable(pl.res,prefix='phoshpo')
  ggsave('results/cellLineTimeCourse/changingPhosphoProtsInCellLines.png')
  
                                     
}