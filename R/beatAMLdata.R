#'
#'Specialized set of commands that specifically load beatAML
#'data - it's pretty big
##just looad all the data for this project in a single file


#'summarize dataset
#'@param auc.data - AUC data and clinical
#'@param mol.data -molelcular data
#'@import pheatmap
#'@import dplyr
#'@import gridExtra
plotAllPatients<-function(auc.data,pat.data,pphos){
  library(gridExtra)
  numDrugs=auc.data%>%
    group_by(`AML sample`)%>%
    summarize(numDrugs=n_distinct(Condition))
  pat.df<-pat.data%>%
    group_by(`AML sample`)%>%
    summarize(RNA=any(mRNALevels!=0),mutations=any(geneMutations!=0),
              proteins=any(proteinLevels!=0))%>%
    right_join(numDrugs)
  
  pat.df<-pphos%>%group_by(`Sample`)%>%summarize(phosphoSites=any(LogFoldChange!=0))%>%
    rename(`AML sample`='Sample')%>%right_join(pat.df)
    
  pdf('patientSummaryTab.pdf',height=11)
  grid.table(pat.df)
  dev.off()
  return(pat.df)
  
}



#' loadBeatAMLMolecularData
#' @import dplyr
#' @import tidyr
loadBeatAMLMolecularData<-function(){
  require(tidyr)
  require(dplyr)
  print("loading molecular data")
  orig.data<-querySynapseTable('syn22172602')

  orig.data<-orig.data%>%rename(proteinLevels='LogFoldChange')%>%
    rename(mRNALevels='transcriptCounts')%>%
    rename(geneMutations='Tumor VAF')%>%
    mutate(Gene=unlist(Gene))%>%rowwise()%>%
    mutate(binaryMutations=ifelse(geneMutations==0,0,1))
  
  pat.data<<-querySynapseTable("syn22314121")%>%
    subset(Treatment=='Vehicle')%>%
    subset(`Cell number`>=10000000)%>%
    dplyr::select(Gene,LogFoldChange,`AML sample`)%>%distinct()%>%
    full_join(orig.data,by=c('AML sample','Gene'))%>%
    rowwise()%>%
    mutate(proteinLevels=max(proteinLevels,LogFoldChange,na.rm=T))%>%
    select(-LogFoldChange)%>%
    mutate(mRNALevels=tidyr::replace_na(mRNALevels,0))%>%
    mutate(geneMutations=tidyr::replace_na(geneMutations,0))%>%
        mutate(binaryMutations=tidyr::replace_na(binaryMutations,0))%>%
 # mutate(countMetric=unlist(countMetric))%>%
   distinct()
  
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
  pat.phos<<-rbind(pat.phos,unique(extra.phos),soraf.phos)
  
  ##getting kinase
  print('Getting kinase estimates')
  pat.kin <<-mapPhosphoToKinase(pat.phos)

  
}

#' looadBeatAMLClinicalDrugData
#' This function looads the clinical parameters and the
#' drug dosage data for each patient
#' @param threshold is the fraction of patients with an AUC under
#' 100, required to be included
#' @import dplyr
#' @import tidyr
loadBeatAMLClinicalDrugData<-function(threshold=0.10){
  print("Loading patient variables and drug response")

  require(dplyr)
  require(tidyr)
  
  mut.status<<-querySynapseTable("syn23538858")%>%
    mutate(status=tolower(status))%>%
    pivot_wider(values_from='status',names_from='variant')
  
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
    rename(AUC='meanAUC')%>%
    left_join(mut.status)
  
    numSens<-auc.dat%>%
      group_by(Condition)%>%
      subset(AUC<100)%>%summarize(numSens=n())
    
    fracSens<-auc.dat%>%group_by(Condition)%>%
      summarize(nSamps=n())%>%
      left_join(numSens)%>%mutate(fracSens=numSens/nSamps)
    
    withSens=subset(fracSens,fracSens>threshold)%>%
      subset(numSens>1)
    include<-union(withSens$Condition,'Gilteritinib (ASP-2215)')
    auc.dat<-subset(auc.dat,Condition%in%include)
  
    drug.combos<-unique(auc.dat$Condition[grep(" - |\\+",auc.dat$Condition)])
    print("Removing drug combinations")
    auc.dat<<-subset(auc.dat,!Condition%in%drug.combos)
  
  
}
  
  
getNetworksAndLVs<-function(){

  ##reduce dims
  print("Getting Latent Variables")
  lv.df<-querySynapseTable('syn22274890')
  
  ##now get network dsitances
  print('Getting network distances')
  #pat.net<-querySynapseTable("syn22343177")
  pat.net <-querySynapseTable("syn23448907")
  filtered = pat.net%>% 
    filter(net2_type=='community')%>%
    filter(hyp1=='patients')%>%
    filter(hyp2=='panCan')
  
#  filtered=pat.net%>%mutate(same=(`Hypha 1`==`Hypha 2`))%>%
#    filter(same)%>%
#    subset(net2_type=='community')
#  mut.nets<<-filtered%>%subset(`Hypha 1`=='mutations')%>%
#    select(Community='Network 2', distance,`AML sample`='Network 1')%>%distinct()
  prot.nets<<-filtered%>%select(Community='net2',distance,`AML sample`='net1')
    #filtered%>%subset(`Hypha 1`=='proteomics')%>%
    #select(Community='Network 2', distance,`AML sample`='Network 1')%>%distinct()
}


#'loadBeatAMLData
#'General function that calls all beatAML data into memory, required for
#'analysis code to work
#'@export
loadBeatAMLData<-function(){
  getNetworksAndLVs()
  loadBeatAMLMolecularData()
  loadBeatAMLClinicalDrugData()
  res<-plotAllPatients(auc.dat,pat.data,pat.phos)
  
  full.pats<<-res%>%
    rowwise()%>%
    mutate(fullData=(as.character(RNA)=="TRUE" && as.character(proteins)=="TRUE"
                     && as.character(mutations)=="TRUE" && as.character(phosphoSites)=="TRUE"))%>%
    subset(fullData==TRUE)%>%
    select(`AML sample`)
  
}