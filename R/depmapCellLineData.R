##new script to pull and process depmap data



loadDepMapMolecularData<-function(){
  require(dplyr)
  require(tidyr)
  #proteomic data
  syn<-synapseLogin()
  
  
#  prot.samps<-read_xlsx(syn$get('syn22975117')$path,sheet=2)%>%
#    subset(`Tissue of Origin`%in%c('Acute Myeloid Leukemia',
#                                   'Haematopoietic and Lymphoid Tissue'))%>%select(`Cell Line`)%>%distinct()
  prot<-read.csv(syn$get('syn22975106')$path)
  prot<-prot[,grep('HAEMAT',colnames(prot))]
  
  #transcript data
  exp<-read.csv(syn$get('syn22975095')$path)

  #genomic dat
  
  #phospho data
  

}

loadDepMapClinicalDrugData<-function(){
  require(dplyr)
  require(tidyr)
  require(readxl)
  syn<-synapseLogin()
  samps<-read.csv(syn$get('syn22978057')$path)%>%
    subset(primary_disease=='Leukemia')
  
  drug.data<-read.csv(syn$get('syn22975062')$path,check.names=F)
  dn<-colnames(drug.data)
  dn[1]<-'Sample'
  colnames(drug.data)<-dn
  
  drug.dat<-drug.data%>%subset(Sample%in%samps$DepMap_ID)%>%
    as.data.frame()%>%
    pivot_longer(cols=setdiff(dn,'Sample'),
                 names_to='DrugDose',values_to='foldChange')
  
}

loadDepMapData<-function(){
  loadDepMapClinicalDrugData()
  loadDepMapMolecularData()
}

