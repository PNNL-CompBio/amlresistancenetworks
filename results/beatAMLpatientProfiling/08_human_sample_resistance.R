##tests differential expression of proteins

library(amlresistancenetworks)
library(dplyr)
library(ggplot2)
# this is the sensitivity data
pat.data<-querySynapseTable('syn22156868')#readRDS(system.file('patientMolecularData.Rds',package='amlresistancenetworks'))
pat.phos<-querySynapseTable("syn22156830")#readRDS(system.file('patientPhosphoSampleData.Rds',package='amlresistancenetworks'))
pat.drugClin<-querySynapseTable("syn22156866")#readRDS(system.file('patientDrugAndClinical.Rds', package='amlresistancenetworks'))
drug.class<-querySynapseTable("syn22156956")%>%rename(Condition='inhibitor')
pat.drugClin<-querySynapseTable("syn22156866")%>%
  mutate(Condition=unlist(Condition))%>%
  left_join(drug.class,by='Condition')


#' working on this function still
drugMolRegression<-function(clin.data,mol.data,mol.feature){
  tdat<-mol.data%>%select(Gene,`AML sample`,Mol=mol.feature)%>%
    subset(!is.na(Mol))%>%left_join(clin.data,by='AML sample')
  
  dcors<-tdat%>%select(Gene,Mol,Condition,AUC)%>%
    distinct()%>%
    group_by(Gene,Condition)%>%
    purrr::pmap_df(lm(x=Mol,y=AUC))
  
  
}

#' how do we visualize the correlations of each drug/gene pair?
#' @param cor.res
#' @param cor.thresh 
plotCorrelationsByDrug<-function(cor.res,cor.thresh){
  ##for each drug class - what is the distribution of correlations broken down by data type

  do.p<-function(dat,cor.thresh){
    print(head(dat))
    fam=dat$family[1]
    #fam=dat%>%dplyr::select(family)%>%unlist()
    #fam=fam[1]
    fname=paste0(fam,'_correlations.png')
    p1<-ggplot(dat,aes(y=Condition,x=drugCor))+
      geom_density_ridges_gradient(aes(fill=feature,alpha=0.5))+
      scale_fill_viridis_d()+ggtitle(paste('Correlation with',fam))
    ##for each drug, how many genes have a corelation over threshold
    p2<-subset(dat,abs(drugCor)>cor.thresh)%>%
      ungroup()%>%
      group_by(Condition,feature,family)%>%
      summarize(CorVals=n_distinct(Gene))%>%ggplot(aes(x=Condition,y=CorVals,fill=feature))+
        geom_bar(stat='identity',position='dodge')+
      scale_fill_viridis_d()+ggtitle(paste("Correlation >",cor.thresh))
    
    cowplot::plot_grid(p1,p2,nrow=2) 
    
    ggsave(fname)
    return(fname)
  }
  
  famplots<-all.cors%>%split(all.cors$family)%>%purrr::map(do.p,cor.thresh)                                                  
  
  lapply(famplots,synapseStore,'syn22130776')
  
}


computeDiffExByDrug<-function(sens.data){
  #samps<-assignSensResSamps(sens.data,'AUC',sens.val,res.val)

  result<-sens.data%>%
    dplyr::select('Sample',Gene,cellLine='Drug',value="LogFoldChange",treatment='Status')%>%
    distinct()%>%
    subset(!is.na(cellLine))%>%
    group_by(cellLine)%>%
    group_modify(~ computeFoldChangePvals(.x,control=NA,conditions =c("Sensitive","Resistant")),keep=TRUE)%>%
    rename(Drug='cellLine')
    
  
  table(result%>%group_by(Drug)%>%subset(p_adj<0.1)%>%summarize(sigProts=n()))
  
  result
}


computeAUCCorVals<-function(clin.data,mol.data,mol.feature){
  tdat<-mol.data%>%select(Gene,`AML sample`,Mol=mol.feature)%>%
    subset(!is.na(Mol))%>%left_join(clin.data,by='AML sample')
  
  dcors<-tdat%>%select(Gene,Mol,Condition,AUC,family)%>%
      distinct()%>%
      group_by(Gene,Condition)%>%
      mutate(numSamps=n(),drugCor=cor(Mol,AUC))%>%
      select(Gene,Condition,numSamps,drugCor)%>%distinct()%>%
      arrange(desc(drugCor))%>%subset(numSamps>10)%>%
    mutate(feature=mol.feature)
  
  ggplot(dcors,aes(x=drugCor))+geom_histogram(aes(fill=Condition))

      
  return(dcors)
  
}
auc.dat<- subset(pat.drugClin,Metric%in%c('auc','AUC'))%>%
  select(`AML sample`,gender,ageAtDiagnosis,vitalStatus,overallSurvival,family,Condition,AUC='Value')%>%
  mutate(Condition=unlist(Condition))

all.cors<-purrr::map_df(list(mRNA='transcriptCounts',
                           protein='LogFoldChange',
                           gene='Tumor VAF'),~ computeAUCCorVals(auc.dat,pat.data,.x))

plotCorrelationsByDrug(all.cors,cor.thresh=0.8)
