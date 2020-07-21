##tests differential expression of proteins

library(amlresistancenetworks)
library(dplyr)
library(ggplot2)
# this is the sensitivity data
#molecular data
pat.data<-querySynapseTable('syn22172602')#readRDS(system.file('patientMolecularData.Rds',package='amlresistancenetworks'))
#phoshpo data
pat.phos<-querySynapseTable("syn22156830")#readRDS(system.file('patientPhosphoSampleData.Rds',package='amlresistancenetworks'))
#   pat.drugClin<-querySynapseTable("syn22156866")
    #readRDS(system.file('patientDrugAndClinical.Rds', package='amlresistancenetworks'))


drug.class<-querySynapseTable("syn22156956")%>%
  dplyr::rename(Condition='inhibitor')%>%
  mutate(Condition=unlist(Condition))%>%
  mutate(family=unlist(family))

#drug response
pat.drugClin<-querySynapseTable("syn22170540")%>%
  mutate(Condition=unlist(Condition))%>%
  left_join(drug.class,by='Condition')


print("Fixing mispelling")
pat.drugClin$Condition<-sapply(pat.drugClin$Condition,function(x){
  ifelse(x=='Vargetef','Vargatef',x)})

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
      right_join(numDrugs)
  pdf('patientSummaryTab.pdf',height=11)
  grid.table(pat.df)
  dev.off()
  #ggsave("PatientSummaryTab.pdf")
  return(pat.df)

}




print("Reformating AUC data")
clin.dat<-pat.drugClin%>%
  mutate(family=unlist(family))%>%
  dplyr::select(`AML sample`,gender,ageAtDiagnosis,vitalStatus,overallSurvival,family,Condition)%>%
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

drug.combos<-unique(auc.dat$Condition[grep(" - |\\+",auc.dat$Condition)])
print("Removing drug combinations")
auc.dat<-subset(auc.dat,!Condition%in%drug.combos)

pat.data<-pat.data%>%rename(proteinLevels='LogFoldChange')%>%
  rename(mRNALevels='transcriptCounts')%>%
  rename(geneMutations='Tumor VAF')%>%
  mutate(Gene=unlist(Gene))

#summarizing AUC data
pat.summary<-plotAllPatients(auc.dat,pat.data)
auc.dat<-auc.dat%>%left_join(pat.summary)
plotAllAUCs(auc.dat,'percAUC')
plotAllAUCs(auc.dat,'AUC')


class.size<-drug.class%>%group_by(family)%>%
  summarize(size=n())
bigger.fams<-subset(class.size,size>1)

if(FALSE){
  print("Getting correlations")
  auc.only<-auc.dat%>%select(`AML sample`,Condition,AUC)%>%distinct()
  
  all.cors<-purrr::map_df(list(mRNA='mRNALevels',
                               protein='proteinLevels',
                               gene='geneMutations'),~ computeAUCCorVals(auc.only,pat.data,.x))
  
  print('Plotting')
  all.cors<-all.cors%>%
    left_join(drug.class)
  red.cors<-all.cors%>%
    subset(family%in%bigger.fams$family)
  plotCorrelationsByDrug(red.cors,cor.thresh=0.8)
}

if(FALSE){
print("Getting Regression predictors")
all.preds<-purrr::map_df(list(mRNA='mRNALevels',
                             protein='proteinLevels',
                             gene='geneMutations'),~ drugMolRegression(auc.dat,
                                                                  pat.data,
                                                                  .x,category='Condition'))




full.preds<-drugMolRegression(auc.dat,pat.data,c('mRNALevels','proteinLevels','geneMutations'))
##now how do we visualize this? 

}

print("Getting Random Forest predictors")
all.preds<-purrr::map_df(list(mRNA='mRNALevels',
                              protein='proteinLevels',
                              gene='geneMutations'),~ drugMolRandomForest(auc.dat,
                                                                        pat.data,
                                                                        .x,category='Condition'))




full.preds<-drugMolRandomForest(auc.dat,pat.data,c('mRNALevels','proteinLevels','geneMutations'))
##now how do we visualize this? 


comb.preds<-rbind(all.preds,full.preds%>%mutate(Molecular='allThree'))

full.tab<-comb.preds%>%
  rename(Condition='var')%>%
  ungroup()%>%
  left_join(drug.class)%>%
  mutate(hasPredictiveGenes=(numGenes!=0))

min.Preds<-full.tab%>%
  group_by(Condition)%>%
  summarize(minType=Molecular[which.min(MSE)])

print("Plotting predictors")

p<- full.tab%>%
    subset(family%in%bigger.fams$family)%>%
    ggplot(aes(x=family,fill=hasPredictiveGenes))+
    geom_bar(position='dodge')+
    facet_grid(Molecular~.)+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle('Number of drugs predicted')
ggsave('allMSEpreds.png',width=8,units='in')

p<- full.tab%>%
  subset(hasPredictiveGenes==TRUE)%>%
  subset(family%in%bigger.fams$family)%>%
  ggplot(aes(x=family,y=MSE,fill=Molecular))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle('Mean Squared Error of predictors with features')
ggsave('allMSEGenePreds.png',width=10,units='in')

red.tab<-full.tab%>%
  subset(hasPredictiveGenes==TRUE)%>%
  subset(family%in%bigger.fams$family)

for(fam in unique(red.tab$family)){
  drtab<-subset(red.tab,family==fam)
  ggplot(drtab,aes(x=Condition,y=MSE,fill=Molecular))+geom_bar(position='dodge',stat='identity')+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      ggtitle(paste0('Mean Squared Error of ',fam,'response'))
    ggsave(paste0(fam,'_MSEGenePreds.png'),width=10,units='in')
  }



