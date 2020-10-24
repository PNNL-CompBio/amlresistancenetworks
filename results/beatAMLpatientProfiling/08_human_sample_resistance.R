##tests differential expression of proteins

library(amlresistancenetworks)
library(dplyr)
library(ggplot2)
# this is the sensitivity data
if(!exists('dataLoaded')){
#  source('beatAMLdata.R')
  amlresistancenetworks::loadBeatAMLData()
  dataLoaded=TRUE
}
class.size<-drug.class%>%group_by(family)%>%
  summarize(size=n())
bigger.fams<-subset(class.size,size>1)


#summarizing AUC data
pat.summary<-plotAllPatients(auc.dat,pat.data,pat.phos)
auc.dat<-auc.dat%>%left_join(pat.summary)

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
  reg.preds<-purrr::map_df(list(mRNA='mRNALevels',
                               protein='proteinLevels',
                               gene='geneMutations'),~ drugMolRegression(auc.dat,
                                                                    pat.data,
                                                                    .x,category='Condition'))
  



#  full.preds<-drugMolRegression(auc.dat,pat.data,c('mRNALevels','proteinLevels','geneMutations'))
##now how do we visualize this? 

}

if(FALSE){
  print("Getting Random Forest predictors")
  rf.preds<-purrr::map_df(list(mRNA='mRNALevels',
                                protein='proteinLevels',
                                gene='geneMutations'),~ drugMolRandomForest(auc.dat,
                                                                          pat.data,
                                                                        .x,category='Condition'))



 # full.preds<-drugMolRandomForest(auc.dat,pat.data,
#                                  c('mRNALevels','proteinLevels','geneMutations'))
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
  
  plotPredictions(full.tab)

}

plotPredictions<-function(full.tab){
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

p2<-subset(all.preds,method!='RandomForest')%>%subset(numFeatures>0)%>%subset(Molecular%in%c('geneMutations','Phosphosite','proteinLevels','mRNALevels'))%>%ggplot(aes(x=numFeatures,MSE,col=var,shape=Molecular))+geom_point()+facet_grid(~method)
ggsave(p2,file='numFeaturesByDrug.png')
}


