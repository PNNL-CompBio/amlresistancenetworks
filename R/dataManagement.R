##load data from files into rdata files

#' reads in metadata from excel spreadsheet and tidies it
#' @require readxl
#' @require dplyr
#' 

readAndTidyMetadata<-function(){
  metadata<-readxl::read_xlsx('../../Projects/CPTAC/ex12_data/LabelingMetadata_exp12.xlsx')
 
  #update sample
  metadata$Sample<-sapply(metadata$`Sample ID`,function(x) paste('Sample',x,sep='_'))
  
  metadata<-metadata%>%
      tidyr::separate(col='Sample Info',into=c('treatment','cellLine','flask','ligand'),sep=', ',fill="right")%>%
      dplyr::select(Sample,cellLine,treatment,flask,ligand)
  
  #update cell line
  cl<-subset(metadata,cellLine=='10M')%>%dplyr::select(Sample,treatment)%>%mutate(cellLine=stringr::str_replace(treatment,'parental ',''))%>%
      dplyr::select(Sample,cellLine)%>%mutate(treatment='None',ligand='None')
 
  lig.dat<-metadata%>%subset(!is.na(ligand))%>%
    dplyr::select(c(flask,cellLine,ligand))
  
  fin=metadata%>%
    subset(cellLine!='10M')%>%
    dplyr::select(cellLine,treatment,flask,Sample)%>%
    left_join(lig.dat,by=c('flask','cellLine'))%>%
    dplyr::select(-flask)%>%
    rbind(cl)
  
  
  return(fin)
  
}


##load data from files into rdata files

#' reads in metadata from excel spreadsheet and tidies it
#' @require readxl
#' @require dplyr
#' @require tidyr
#' 

readAndTidyQuizMetadata<-function(){
  metadata<-readxl::read_xlsx('../../Projects/CPTAC/ex12_data/quizartLabelingMetadata.xlsx',skip=3)
  
  res<-metadata%>%
      dplyr::rowwise()%>%
    dplyr::mutate(SampleTemp=paste(Plex,`Sample #`,sep='.'))%>%
    dplyr::mutate(Sample=stringr::str_replace(SampleTemp,'Sample ',''))%>%
      dplyr::select(Sample,`Sample Name`,Ligand)
  
  fin<-res%>%
    tidyr::separate(col=`Sample Name`,into=c('tmpCellLine','Flask'),sep=', ')%>%
    dplyr::mutate(treatment='None')
  
  
  pfin<-subset(fin,is.na(Flask))%>%
    dplyr::select(Sample,treatment)%>%
    dplyr::mutate(cellLine='MOLM14',Ligand='None')
  
  lig<-fin%>%
    dplyr::select(Flask,tmpLigand='Ligand')%>%distinct()%>%
    subset(!is.na(tmpLigand))%>%
    rowwise()%>%
    mutate(Ligand=stringr::str_replace(tmpLigand,'FLT3 ligand','FLT3'))%>%
    dplyr::select(Flask,Ligand)


  tnfin<-subset(fin,!is.na(Flask))%>%
      dplyr::select(Sample,Flask)%>%
      mutate(cellLine='MOLM14')%>%
      left_join(lig,by='Flask')%>%
    dplyr::select(Sample,cellLine,treatment='Ligand')%>%mutate(treatment='Quizartinib')
  
  fin<-rbind(tnfin,pfin)
  
  return(fin)
  
}

#' reads and tidies protein measurements
#' joins with metadata
#' saves tidied data and returns
#' @export
#' @require dplyr
#
readAndTidyProtMeasures<-function(){
  metadata<-readAndTidyMetadata()
  dat<-read.table('../../projects/CPTAC/ex12_data/ptrc_ex12_global_median_centered_with_genes.txt',sep='\t',header=T)
  library(tidyr)
  gilteritinib.data<-tidyr::pivot_longer(dat,-c(Protein,Gene),"Sample")%>%
      dplyr::left_join(metadata,by='Sample')
saveRDS(gilteritinib.data,file='inst/gilteritinibData.Rds')
  return(gilteritinib.data)

    }



#' reads and tidies protein measurements from quizartinib treatment
#' joins with metadata
#' saves tidied data and returns
#' @export
#' @require dplyr
#' @require tidyr
#
readAndTidyQuizProtMeasures<-function(){
  metadata<-readAndTidyQuizMetadata()
  dat<-read.table('../../Projects/CPTAC/ex12_data/ptrc_exp12_global_round2_with_genes.txt',sep='\t',header=T)
  library(tidyr)
  quiz.data<-tidyr::pivot_longer(dat,-c(Entry_name,Gene),"Sample")%>%
      dplyr::left_join(metadata,by='Sample')
  saveRDS(quiz.data,file='inst/quizartinibData.Rds')
  return(quiz.data)
  
}



#' reads and tidies experiment 11 metadata
#' @require dplyr
#' @require readxl
#'
readAndTidySensMetadata<-function(){
  samp.names<-readxl::read_xlsx('../../Projects/CPTAC/ex11_data/Agarwal_PNNL_CPs_Final list 012919.xlsx')%>%
    dplyr::select(c(`Specimen ID`,Barcode))%>%
    dplyr::rename(`AML sample`='Specimen ID')%>%
    distinct()
  
  metadata<-readxl::read_xlsx('../../Projects/CPTAC/ex11_data/Agarwal_PNNL_CPs_Final list 012919_mod032020.xlsx',sheet=1,skip=1)%>%
    subset(!is.na(`AML sample`))%>%
    subset(`AML sample`!='AML sample')
  
 # clin.data<-metadata%>%
#    tidyr::pivot_longer(cols=c(Sex,`ELN Risk`,`Disease status`,`FLT3-ITD`,NPM1,comment),names_to='Clinical Variable',values_to='Clinical Value')%>%
#    dplyr::select(Age,`Clinical Variable`,`Clinical Value`,`AML sample`)%>%
#    distinct()
  ic50<-metadata%>%
    tidyr::pivot_longer(cols=c(4,5,6),names_to='IC50 Condition',values_to='Value')%>%
    dplyr::select(`AML sample`,`IC50 Condition`,Value)%>%
    tidyr::separate(`IC50 Condition`,into=c("Metric","Condition"),sep=' ')%>%
    distinct()
  auc<-metadata%>%
    tidyr::pivot_longer(cols=c(7,8,9),names_to='AUC Condition',values_to='Value')%>%
    dplyr::select(`AML sample`,`AUC Condition`,Value)%>%
    tidyr::separate(`AUC Condition`,into=c("Metric","Condition"),sep=' ')%>%
    distinct()
  # 
  # pd.ic50<-metadata%>%
  #   tidyr::pivot_longer(cols=c(8,9,10),names_to='IC50 Condition',values_to='Value')%>%
  #   dplyr::select(`AML sample`,`IC50 Condition`,Value)%>%
  #   tidyr::separate(`IC50 Condition`,into=c("Condition","Metric"),sep=' ')%>%
  #   distinct()
  # pd.auc<-metadata%>%
  #   tidyr::pivot_longer(cols=c(11,12,13),names_to='AUC Condition',values_to='Value')%>%
  #   dplyr::select(`AML sample`,`AUC Condition`,Value)%>%
  #   tidyr::separate(`AUC Condition`,into=c("Condition","Metric"),sep=' ')%>%
  #   distinct()
  # 
  
  full.metadata<-rbind(auc,ic50)%>%
    #left_join(clin.data,by='AML sample')%>%
    left_join(samp.names)
    
  return(full.metadata)
}

#' reads and tidies bulk proteomic data for experiment 11
#' @require dplyr
#' @require tidyr
#' @export
#' 
readAndTidySensProtMeasure<-function(){
  metadata<-readAndTidySensMetadata()
  dat<-read.table('../../Projects/CPTAC/ex11_data/ptrc_ex11_kurtz_2plex_global_d2_with_genes.txt',sep='\t',header=T)%>%
    tidyr::pivot_longer(cols=c(3:ncol(dat)),names_to='Sample', values_to='LogFoldChange')%>%
    mutate(Barcode=as.numeric(stringr::str_replace(Sample,"X","")))%>%
    dplyr::select(Barcode,Gene, LogFoldChange)%>%
    dplyr::left_join(metadata,by='Barcode')
  
  drugSensData<-dat
  saveRDS(drugSensData,file='inst/gilteritinibSensitivityData.Rds')
  
  dat
}