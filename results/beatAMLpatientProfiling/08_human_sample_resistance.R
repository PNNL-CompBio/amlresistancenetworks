##tests differential expression of proteins

library(amlresistancenetworks)
library(dplyr)
library(ggplot2)
# this is the sensitivity data
pat.data<-querySynapseTable('syn22168781')#readRDS(system.file('patientMolecularData.Rds',package='amlresistancenetworks'))
pat.phos<-querySynapseTable("syn22156830")#readRDS(system.file('patientPhosphoSampleData.Rds',package='amlresistancenetworks'))
#   pat.drugClin<-querySynapseTable("syn22156866")
    #readRDS(system.file('patientDrugAndClinical.Rds', package='amlresistancenetworks'))
drug.class<-querySynapseTable("syn22156956")%>%
  rename(Condition='inhibitor')%>%
  mutate(Condition=unlist(Condition))%>%
  mutate(family=unlist(family))

pat.drugClin<-querySynapseTable("syn22168720")%>%
  mutate(Condition=unlist(Condition))%>%
  left_join(drug.class,by='Condition')



#'plotAll patients
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
      left_join(numDrugs)
  pdf('patientSummaryTab.pdf')
  grid.table(pat.df)
  dev.off()
  #ggsave("PatientSummaryTab.pdf")

}


print("Fixing mispelling")
pat.drugClin$Condition<-sapply(pat.drugClin$Condition,function(x){
  ifelse(x=='Vargetef','Vargatef',x)})


print("Reformating AUC data")
clin.dat<-pat.drugClin%>%
  mutate(family=unlist(family))%>%
  select(`AML sample`,gender,ageAtDiagnosis,vitalStatus,overallSurvival,family,Condition)%>%
  distinct()

auc.dat<- subset(pat.drugClin,Metric%in%c('auc','AUC'))%>%
  dplyr::select(`AML sample`,Condition,AUC='Value')%>%
  distinct()%>%
  mutate(Condition=unlist(Condition))%>%
  group_by(Condition)%>%
  mutate(medAUC=median(AUC))%>%
  mutate(percAUC=100*AUC/medAUC)%>%
  ungroup()%>%
  left_join(clin.dat)


pat.data<-pat.data%>%rename(proteinLevels='LogFoldChange')%>%
  rename(mRNALevels='transcriptCounts')%>%
  rename(geneMutations='Tumor VAF')

#summarizing AUC data
plotAllAUCs(auc.dat)
plotAllPatients(auc.dat,pat.data)

print("Getting predictors")
all.preds<-purrr::map_df(list(mRNA='mRNALevels',
                             protein='proteinLevels',
                             gene='geneMutations'),~ drugMolRegression(auc.dat,
                                                                  pat.data,
                                                                  .x,category='Condition'))

full.tab<-all.preds%>%
  rename(Condition='var')%>%
  ungroup()%>%
  left_join(drug.class)%>%
  mutate(noGenes=(numGenes==0))

class.size<-drug.class%>%group_by(family)%>%
  summarize(size=n())
bigger.fams<-subset(class.size,size>1)

print("Plotting predictors")

p<- full.tab%>%
    subset(family%in%bigger.fams$family)%>%
    ggplot(aes(x=family,y=MSE,fill=noGenes))+
    geom_boxplot()+
    facet_grid(~Molecular)+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle('Mean Squared Error of predictors by features')
ggsave('MSEpreds.png')


if(FALSE){
  print("Getting correlations")
  all.cors<-purrr::map_df(list(mRNA='mRNALevels',
                           protein='proteinLevels',
                           gene='geneMutations'),~ computeAUCCorVals(auc.dat,pat.data,.x))

  print('Plotting')
  plotCorrelationsByDrug(all.cors,cor.thresh=0.8)
}