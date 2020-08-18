##predict drug response from kinase averages

library(amlresistancenetworks)

if(!exists('dataLoaded')){
  source('beatAMLdata.R')
  dataLoaded=TRUE
}
pat.kin <-mapPhosphoToKinase(pat.phos)

cor.results<-computeAUCCorVals(auc.dat,rename(pat.kin,Gene='Kinase',`AML sample`='Sample'),'meanLFC')

reg.results=pat.kin%>%rename(Gene='Kinase',`AML sample`='Sample')%>%
  drugMolRegression(auc.dat,.,'meanLFC')
rf.results=pat.kin%>%rename(Gene='Kinase',`AML sample`='Sample')%>%
  drugMolRandomForest(auc.dat,.,'meanLFC')


