##run PLIER analysis on data

library(amlresistancenetworks)
library(PLIER)
library(dplyr)
library(tidyr)
#rc2path = "https://ndownloader.figshare.com/files/10881866"

##download these behemoth separately
#https://figshare.com/articles/recount_rpkm_RData/5716033/4

#then we want to looad this file
plier.results <- readRDS("recount_PLIER_model.RDS")

#'
#'Collado-Torres L, Nellore A, Kammers K, Ellis SE, Taub MA, Hansen KD, Jaffe AE, Langmead B and Leek JT (2017). "Reproducible RNA-seq analysis using recount2." Nature Biotechnology. doi: 10.1038/nbt.3838
#'And the PLIER preprint, if you use the PLIER model:
#'  Mao W, Harmann B, Sealfon SC, Zaslavsky E, and Chikina M (2017). "Pathway-Level Information ExtractoR (PLIER) for gene expression data." bioRxiv. doi: 10.1101/116061
#'  
#

pat.data<-querySynapseTable('syn22172602')#readRDS(system.file('patientMolecularData.Rds',package='amlresistancenetworks'))

#drug response
pat.drugClin<-querySynapseTable("syn22170540")%>%
  mutate(Condition=unlist(Condition))%>%
  left_join(drug.class,by='Condition')


print("Fixing mispelling")
pat.drugClin$Condition<-sapply(pat.drugClin$Condition,function(x){
  ifelse(x=='Vargetef','Vargatef',x)})


print("Reformating AUC data")
clin.dat<-pat.drugClin%>%
 # mutate(family=unlist(family))%>%
  dplyr::select(`AML sample`,gender,ageAtDiagnosis,vitalStatus,overallSurvival,Condition)%>%
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

exprs.mat<-pat.data%>%
  select(`AML sample`,'Gene','transcriptCounts')%>%
  pivot_wider(values_from=transcriptCounts, names_from=`AML sample`,
              values_fn=list(transcriptCounts=mean),
              values_fill=list(transcriptCounts=0.01))%>% #no zeroes allowed!
  tibble::column_to_rownames('Gene')%>%as.matrix()

source("../../../multi-plier/util/plier_util.R")
pat.recount.b <- GetNewDataB(exprs.mat = exprs.mat,
                             plier.model = plier.results)



###now that all the data are loaded we can compute the correlation, regression, and random forest values

lv.df<-pat.recount.b%>%
    as.data.frame()%>%
    tibble::rownames_to_column("Gene")%>%
    pivot_longer(-Gene,names_to='AML sample',values_to="LV")


reg.results<-drugMolRegression(auc.dat,lv.df,'LV')
rf.results<-drugMolRandomForest(auc.dat,lv.df,'LV')
cor.results<-computeAUCCorVals(auc.dat,lv.df,'LV')


