##just looad all the data for this project in a single file

library(amlresistancenetworks)
library(dplyr)
library(tidyr)



#'summarize dataset
#'@param auc.data - AUC data and clinical
#'@param mol.data -molelcular data
#'@import pheatmap
#'@import dplyr
plotAllPatients<-function(auc.data,pat.data,pphos){
  library(gridExtra)
  numDrugs=auc.data%>%
    group_by(`AML sample`)%>%
    summarize(numDrugs=n_distinct(Condition))
  pat.df<-pat.data%>%
    group_by(`AML sample`)%>%
    summarize(RNA=any(mRNALevels!=0),mutations=any(geneMutations!=0),proteins=any(proteinLevels!=0))%>%
    right_join(numDrugs)
  
  pat.df<-pphos%>%group_by(`Sample`)%>%summarize(phosphoSites=any(LogFoldChange!=0))%>%
    rename(`AML sample`='Sample')%>%right_join(pat.df)
    
  pdf('patientSummaryTab.pdf',height=11)
  grid.table(pat.df)
  dev.off()
  return(pat.df)
  
}


print("loading molecular data")
pat.data<-querySynapseTable('syn22172602')

pat.data<-pat.data%>%rename(proteinLevels='LogFoldChange')%>%
  rename(mRNALevels='transcriptCounts')%>%
  rename(geneMutations='Tumor VAF')%>%
  mutate(Gene=unlist(Gene))

soraf.prot<-querySynapseTable("syn22314121")%>%
  subset(Treatment=='Vehicle')%>%
  subset(`Cell number`>=10000000)%>%
  dplyr::select(Gene,LogFoldChange,`AML sample`)%>%distinct()%>%
  full_join(pat.data,by=c('AML sample','Gene'))%>%
  rowwise()%>%
  mutate(proteinLevels=max(proteinLevels,LogFoldChange,na.rm=T))%>%
  select(-LogFoldChange)%>%
  mutate(mRNALevels=tidyr::replace_na(mRNALevels,0))%>%
  mutate(geneMutations=tidyr::replace_na(geneMutations,0))%>%
 # mutate(countMetric=unlist(countMetric))%>%
  distinct()
pat.data<-soraf.prot

non.soraf.prot<-subset(pat.data,!`AML sample`%in%soraf.prot$`AML sample`)

###replace values in original table....
print("Getting phosphosite data")
pat.phos<-querySynapseTable("syn22156830")
pat.phos$site<-unlist(pat.phos$site)
pat.phos$Gene<-unlist(pat.phos$Gene)

extra.phos<-querySynapseTable("syn22156814")%>%
  dplyr::select(Gene,site,Peptide,LogFoldChange='value',Sample="AML sample")

soraf.phos<-querySynapseTable("syn22314122")%>%
  subset(Treatment=='Vehicle')%>%
  subset(`Cell number`>=10000000)%>%
  dplyr::select(Gene,site,Peptide,LogFoldChange,Sample="AML sample")%>%distinct()

soraf.phos$site<-unlist(soraf.phos$site)
soraf.phos$Gene<-unlist(soraf.phos$Gene)

extra.phos$site <-unlist(extra.phos$site)
extra.phos$Gene <-unlist(extra.phos$Gene)
pat.phos<-rbind(pat.phos,unique(extra.phos),soraf.phos)

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
  ifelse(x=='Doramapimod','Doramapimod (BIRB 796)',x)})
pat.drugClin$Condition<-sapply(pat.drugClin$Condition,function(x){
  ifelse(x=='Gilteritinib','Gilteritinib (ASP-2215)',x)})

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

res<-plotAllPatients(auc.dat,pat.data,pat.phos)

numSens<-auc.dat%>%
  group_by(Condition)%>%
  subset(AUC<100)%>%summarize(numSens=n())
fracSens<-auc.dat%>%group_by(Condition)%>%
  summarize(nSamps=n())%>%
  left_join(numSens)%>%mutate(fracSens=numSens/nSamps)

threshold<-0.10
withSens=subset(fracSens,fracSens>threshold)%>%
  subset(numSens>1)
include<-union(withSens$Condition,'Gilteritinib (ASP-2215)')
auc.dat<-subset(auc.dat,Condition%in%include)

drug.combos<-unique(auc.dat$Condition[grep(" - |\\+",auc.dat$Condition)])
print("Removing drug combinations")
auc.dat<-subset(auc.dat,!Condition%in%drug.combos)

##reduce dims
print("Getting Latent Variables")
lv.df<-querySynapseTable('syn22274890')

##getting kinase
print('Getting kinase estimates')
pat.kin <-mapPhosphoToKinase(pat.phos)

##now get network dsitances
print('Getting network distances')
pat.net<-querySynapseTable("syn22343177")
filtered=pat.net%>%mutate(same=(`Hypha 1`==`Hypha 2`))%>%
  filter(same)%>%
  subset(net2_type=='community')
mut.nets<-filtered%>%subset(`Hypha 1`=='mutations')%>%
  select(Community='Network 2', distance,`AML sample`='Network 1')%>%distinct()
prot.nets<-filtered%>%subset(`Hypha 1`=='proteomics')%>%
  select(Community='Network 2', distance,`AML sample`='Network 1')%>%distinct()
  
