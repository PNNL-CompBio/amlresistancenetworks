##just looad all the data for this project in a single file

library(amlresistancenetworks)
library(dplyr)
library(tidyr)



#'summarize dataset
#'@param auc.data - AUC data and clinical
#'@param mol.data -molelcular data
#'@import pheatmap
#'@import dplyr
plotAllPatients<-function(auc.data,pat.data){
  library(gridExtra)
  numDrugs=auc.data%>%
    group_by(`AML sample`)%>%
    summarize(numDrugs=n_distinct(Condition))
  pat.df<-pat.data%>%
    group_by(`AML sample`)%>%
    summarize(RNA=any(mRNALevels!=0),mutations=any(geneMutations!=0),proteins=any(proteinLevels!=0))%>%
    right_join(numDrugs)
  pdf('patientSummaryTab.pdf',height=11)
  grid.table(pat.df)
  dev.off()
  #ggsave("PatientSummaryTab.pdf")
  return(pat.df)
  
}



print("loading molecular data")
pat.data<-querySynapseTable('syn22172602')

pat.data<-pat.data%>%rename(proteinLevels='LogFoldChange')%>%
  rename(mRNALevels='transcriptCounts')%>%
  rename(geneMutations='Tumor VAF')%>%
  mutate(Gene=unlist(Gene))

pat.phos<-querySynapseTable("syn22156830")#readRDS(system.file('patientPhosphoSampleData.Rds',package='amlresistancenetworks'))

print("Loading patient variables and drug response")
#readRDS(system.file('patientMolecularData.Rds',package='amlresistancenetworks'))

drug.class<<-querySynapseTable("syn22156956")%>%
  dplyr::rename(Condition='inhibitor')%>%
  mutate(Condition=unlist(Condition))%>%
  mutate(family=unlist(family))

#drug response
pat.drugClin<-querySynapseTable("syn22170540")%>%
  mutate(Condition=unlist(Condition))%>%
  left_join(drug.class,by='Condition')

print("Fixing Vargatef mispelling")
pat.drugClin$Condition<-sapply(pat.drugClin$Condition,function(x){
  ifelse(x=='Vargetef','Vargatef',x)})
pat.drugClin$Condition<-sapply(pat.drugClin$Condition,function(x){
  ifelse(x=='Doramapimod','Doramapimod (BIRB 796) ',x)})


print("Reformating AUC data")
clin.dat<-pat.drugClin%>%
  # mutate(family=unlist(family))%>%
  dplyr::select(`AML sample`,gender,ageAtDiagnosis,vitalStatus,overallSurvival,Condition)%>%
  distinct()

auc.dat<- subset(pat.drugClin,Metric%in%c('auc','AUC'))%>%
  dplyr::select(`AML sample`,Condition,AUC='Value')%>%
  distinct()%>%
  group_by(`AML sample`,Condition)%>%
  mutate(meanAUC=mean(AUC))%>%
  ungroup()%>%
  mutate(Condition=unlist(Condition))%>%
  group_by(Condition)%>%
  mutate(medAUC=median(AUC))%>%
  mutate(percAUC=100*AUC/medAUC)%>%
  ungroup()%>%
  left_join(clin.dat)%>%
  select(-AUC)%>%
  rename(AUC='meanAUC')

drug.combos<-unique(auc.dat$Condition[grep(" - |\\+",auc.dat$Condition)])
print("Removing drug combinations")
auc.dat<-subset(auc.dat,!Condition%in%drug.combos)



##reduce dims
print("Getting Latent Variables")
lv.df<-querySynapseTable('syn22274890')


##getting kinase
print('Getting kinase estiamtes')
pat.kin <-mapPhosphoToKinase(pat.phos)



