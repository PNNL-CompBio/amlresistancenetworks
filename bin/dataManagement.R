##load data from files into rdata files
#this file should only be modified rarely when new data are to be added

library(amlresistancenetworks)
library(dplyr)

#' reads in metadata from excel spreadsheet and tidies it
#' @require readxl
#' @require dplyr
#' @require reticulate
#'

readAndTidyMetadata<-function(){
  syn=synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22136312')$path)
  #%>%
  #  stringr::str_replace()

  #update sample
  metadata$Sample<-sapply(metadata$`Sample ID`,function(x) paste('Sample',x,sep='_'))

  metadata<-metadata%>%
      tidyr::separate(col='Sample Info',into=c('treatment','cellLine','flask','ligand'),
                      sep=', ',fill="right")%>%
      dplyr::select(Sample,cellLine,treatment,flask,ligand)

  #update cell line
  cl<-subset(metadata,cellLine=='10M')%>%
    dplyr::select(Sample,treatment)%>%
    mutate(cellLine=stringr::str_replace(treatment,'parental ',''))%>%
      dplyr::select(Sample,cellLine)%>%
    mutate(treatment='None',ligand='None')

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

readAndTidyQuizMetadata<-function(){
  syn=synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22136311')$path,skip=3)

  res<-metadata%>%
      dplyr::rowwise()%>%
    dplyr::mutate(SampleTemp=paste(Plex,`Sample #`,sep='.'))%>%
    dplyr::mutate(Sample=stringr::str_replace(SampleTemp,'Sample ',''))%>%
      dplyr::select(Sample,`Sample Name`,Ligand)

  fin<-res%>%
    tidyr::separate(col=`Sample Name`,into=c('tmpCellLine','Flask'),sep=', ')%>%
    dplyr::mutate(treatment='None')
  
  
  ligs<-fin%>%select(Flask,Ligand)%>%
    subset(!is.na(Ligand))%>%distinct()
  
    
  pfin <- fin%>%
    select(Sample,tmpCellLine,Flask)%>%
    left_join(ligs)%>%
    replace_na(list(Ligand='None'))%>%
    select(-Flask)%>%
    rowwise()%>%
    mutate(cellLine=stringr::str_replace(tmpCellLine,' Parental .+',''))%>%
    select(-tmpCellLine)
  
    return(pfin)

}

#' reads and tidies protein measurements
#' joins with metadata
#' saves tidied data and returns
#' @require dplyr
readAndTidyQuizProtMeasures<-function(){
  metadata<-readAndTidyQuizMetadata()
  syn=synapseLogin()
  dat<-read.table(syn$get('syn22136310')$path,sep='\t',header=T)
  library(tidyr)
  quiz.data<-tidyr::pivot_longer(dat,-c(Entry_name,Gene),"Sample")%>%
    subset(Gene!="")%>%
    dplyr::left_join(metadata,by='Sample')%>%
    mutate(value=tidyr::replace_na(value,0))%>%
    mutate(Gene=as.character(Gene))

  synTableStore(quiz.data,'Quizartinib Resistance Proteomics Data')
  
  return(quiz.data)
}



#' reads and tidies phosphoprotein measurements
#' joins with metadata
#' saves tidied data and 
#' 
#' @require dplyr
#
readAndTidyQuizPhosphoProtMeasures<-function(){
  metadata<-readAndTidyQuizMetadata()
  syn=synapseLogin()
  dat<-read.table(syn$get('syn23407326')$path,sep='\t',header=T)
  library(tidyr)
  quiz.phospho.data<-dat%>%
    tidyr::pivot_longer(-c(Entry_name,Gene,site,Peptide),"Sample")%>%
    dplyr::left_join(metadata,by='Sample')%>%
    mutate(value=tidyr::replace_na(value,0))%>%
    mutate(site=as.character(site))%>%
    mutate(Gene=as.character(Gene))
  # subset(!is.na(value))
  
  synTableStore(quiz.phospho.data,'Quizartinib Resistance Phosphoproteomics Data')
  return(quiz.phospho.data)
  
}

#' reads and tidies protein measurements
#' joins with metadata
#' saves tidied data and returns
#' @require dplyr
#
readAndTidyProtMeasures<-function(){
  metadata<-readAndTidyMetadata()
  syn=synapseLogin()
  dat<-read.table(syn$get('syn22136301')$path,sep='\t',header=T)
  library(tidyr)
  gilteritinib.data<-tidyr::pivot_longer(dat,-c(Protein,Gene),"Sample")%>%
      dplyr::left_join(metadata,by='Sample')%>%
    mutate(value=tidyr::replace_na(value,0))
  synTableStore(gilteritinib.data,'Gilteritinib Resistance Proteomics Data')

  return(gilteritinib.data)
}



#' reads and tidies phosphoprotein measurements
#' joins with metadata
#' saves tidied data and returns
#' @require dplyr
#
readAndTidyPhosphoProtMeasures<-function(){
  metadata<-readAndTidyMetadata()
  syn=synapseLogin()
  dat<-read.table(syn$get('syn22136309')$path,sep='\t',header=T)
  library(tidyr)
  gilt.phospho.data<-dat%>%
    tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),"Sample")%>%
    dplyr::left_join(metadata,by='Sample')%>%
    mutate(value=tidyr::replace_na(value,0))
   # subset(!is.na(value))

  synTableStore(gilt.phospho.data,'Gilteritinib Resistance Phosphoproteomics Data')
  return(gilt.phospho.data)

}


#' reads and tidies uncorrected phosphoprotein measurements
#' joins with metadata
#' saves tidied data and returns
#' @require dplyr
#
readAndTidyUncorrectedQuizPhosphoMeasures<-function(){
 metadata<-readAndTidyQuizMetadata()
  syn=synapseLogin()
  dat<-read.table(syn$get('syn24305006')$path,sep='\t',header=T)
  library(tidyr)
  gilt.phospho.data<-dat%>%
    tidyr::pivot_longer(-c(Entry_name,Gene,site,Peptide),"Sample")%>%
    dplyr::left_join(metadata,by='Sample')%>%
    mutate(value=tidyr::replace_na(value,0))
   # subset(!is.na(value))

  synTableStore(gilt.phospho.data,'Quizartinib Resistance Phosphoproteomics Unnormalized')
  return(gilt.phospho.data)

}




#' reads and tidies phosphototeomic data for experimnet 11
#' @require dplyr
#' @require tidyr
readAndTidySensPhosMeasures<-function(){
  metadata<-readAndTidySensMetadata()
  syn<-synapseLogin()
  dat<-read.table(syn$get('syn22130837')$path,sep='\t',header=T)
  gilt.sens.pdata<-dat%>%
    tidyr::pivot_longer(-c(Entry,Gene,site,Peptide,ids),"Sample")%>%
    dplyr::mutate(Barcode=as.numeric(stringr::str_replace(Sample,"X","")))%>%
    dplyr::left_join(metadata,by='Barcode')%>%
    mutate(value=tidyr::replace_na(value,0))

  synTableStore(gilt.sens.pdata,'Drug Combination Phosphoproteomic Data')
  return(gilt.sens.pdata)

}

#' reads and tidies bulk proteomic data for experiment 11
#' @require dplyr
#' @require tidyr
#'
readAndTidySensProtMeasure<-function(){
  metadata<-readAndTidySensMetadata()
  syn<-synapseLogin()
  dat<-read.table(syn$get('syn22130842')$path,sep='\t',header=T)

  drugSensData<-dat%>%tidyr::pivot_longer(cols=c(3:ncol(dat)),names_to='Sample', values_to='LogFoldChange')%>%
    dplyr::mutate(Barcode=as.numeric(stringr::str_replace(Sample,"X","")))%>%
    dplyr::select(Barcode,Gene, LogFoldChange)%>%
    mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))%>%
    dplyr::left_join(metadata,by='Barcode')
  synTableStore(drugSensData,'Drug Combination Proteomic Data')

  return(drugSensData)
}


#' reads and tidies experiment 11 metadata
#' @require dplyr
#' @require readxl
readAndTidySensMetadata<-function(){
  syn=synapseLogin()
  samp.names<-readxl::read_xlsx(syn$get('syn22130839')$path)%>%
    dplyr::select(c(`Specimen ID`,Barcode))%>%
    dplyr::rename(`AML sample`='Specimen ID')%>%
    distinct()

  metadata<-readxl::read_xlsx(syn$get('syn22130840')$path,skip=1)%>%
    subset(!is.na(`AML sample`))%>%
    subset(`AML sample`!='AML sample')

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

  full.metadata<-rbind(auc,ic50)%>%
    #left_join(clin.data,by='AML sample')%>%
    left_join(samp.names)
  return(full.metadata)
}



#' @import readxl
#'
getTimeCourseMetadata<-function(){
  syn<-synapseLogin()                                
    samp.data<-readxl::read_xlsx(syn$get('syn22130821')$path)%>%
        dplyr::select(Sample='Tube name',cellLine='Cell line',treatment='Treated with',
                      timePoint='Time point')
    return(samp.data)
}

#' getTimecourseData
#' @require dplyr
getTimeCourseData<-function(){
    library(dplyr)
    syn=synapseLogin()
    metadata<-getTimeCourseMetadata()
    dat<-read.csv2(syn$get("syn22130819")$path,
                   sep='\t',header=T,stringsAsFactors = FALSE)
    timeCourseData<-dat%>%
      tidyr::pivot_longer(cols=c(4:ncol(dat)),names_to='sample', 
                          values_to='LogFoldChange')%>%
    #dplyr::mutate(specId=stringr::str_replace(sample,".","-"))%>%
    dplyr::mutate(Sample=stringr::str_replace(sample,stringr::fixed("."),"-"))%>%
    dplyr::select(Sample,Gene, LogFoldChange)%>%
    left_join(metadata)%>%
      mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
      subset(!is.na(LogFoldChange))
  synTableStore(timeCourseData,'Cell Line Time Course Proteomics')

}

#' @import dplyr
getTimeCoursePhosphoData<-function(){
  library(dplyr)
  syn<-synapseLogin()
  metadata<-getTimeCourseMetadata()
  timeCoursePhospho<-read.csv2(syn$get('syn22130820')$path,
                 sep='\t',header=T,stringsAsFactors = F)%>%
    tidyr::pivot_longer(-c(Entry,Gene,site,Peptide,ids,Entry.name,Protein),"Sample",
                        values_to='LogFoldChange')%>%
    dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
    dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
    dplyr::select(Sample,Gene, site,Peptide,LogFoldChange)%>%
    mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))%>%
    left_join(metadata)

  synTableStore(timeCoursePhospho,'Cell Line Time Course Phosphoproteomics')

}


#' getCytokine data
#' gets cytokine sensitivity data and metadata in one call
#' @import dplyr
#' @import readxl
getCytokineSensData<-function(){
  library(dplyr)
  syn<-synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22175391')$path)%>%
    mutate(sample=stringr::str_replace(`Sample name (DMS)`,'PTRC_Ex15_',''))%>%
    dplyr::select(sample,CellType,TimePoint,Treatment)%>%
    subset(!is.na(sample))



  metadata$CellType<-sapply(metadata$CellType,function(x) {
    switch(x,Late_M='Late MOLM-13',Late_MR='Late MOLM-13 Tr Resistant',
           M='MOLM-13',MR='MOLM-13 Tr Resistant')})

  #  ifelse(x=='Late_M','Late MOLM-13',ifelse(x=='Late_MR','Late MOLM-13 Tr Resistant',x))})
  metadata$Treatment<-unlist(sapply(metadata$Treatment,function(x){
    switch(x,M='MCP-1',T='Trametinib',TW='Trametinib Withdrawn',none="none",`T+M`='Trametinib+MCP-1')}))

  pdat<-read.csv2(syn$get('syn22862628')$path,sep='\t')%>%
    tidyr::pivot_longer(-c(Protein,Gene),names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)%>%
    subset(sample!='Peptide')

  phdat<-read.csv2(syn$get('syn22862617')$path,sep='\t')%>%
    tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),
                        names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)

  pdat2<-read.csv2(syn$get('syn22173207')$path,sep='\t')%>%
    tidyr::pivot_longer(-c(Protein,Gene),names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)

    phdat2<-read.csv2(syn$get('syn22173201')$path,sep='\t')%>%
      tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),
                          names_to='tsamp',values_to='LogRatio')%>%
      mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
      dplyr::select(-tsamp)%>%
      left_join(metadata)

    fullp<-rbind(pdat,pdat2)
    fullph<-rbind(phdat,phdat2)
    synTableStore(rbind(pdat2,pdat),'Cytokine-induced Drug Sensitivity Proteomics')
    synTableStore(rbind(phdat2,phdat),'Cytokine-induced Drug Sensitivity Phospho-proteomics')

}


updateCytokinePhospho<-function(){
    library(dplyr)
  syn<-synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22175391')$path)%>%
    mutate(sample=stringr::str_replace(`Sample name (DMS)`,'PTRC_Ex15_',''))%>%
    dplyr::select(sample,CellType,TimePoint,Treatment)%>%
    subset(!is.na(sample))



  metadata$CellType<-sapply(metadata$CellType,function(x) {
    switch(x,Late_M='Late MOLM-13',Late_MR='Late MOLM-13 Tr Resistant',
           M='MOLM-13',MR='MOLM-13 Tr Resistant')})

  #  ifelse(x=='Late_M','Late MOLM-13',ifelse(x=='Late_MR','Late MOLM-13 Tr Resistant',x))})
  metadata$Treatment<-unlist(sapply(metadata$Treatment,function(x){
    switch(x,M='MCP-1',T='Trametinib',TW='Trametinib Withdrawn',none="none",`T+M`='Trametinib+MCP-1')}))

    phdat<-read.csv2(syn$get('syn24305135')$path,sep='\t')%>%
    tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),
                        names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)

      phdat2<-read.csv2(syn$get('')$path,sep='\t')%>% ##TBD
      tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),
                          names_to='tsamp',values_to='LogRatio')%>%
      mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
      dplyr::select(-tsamp)%>%
      left_join(metadata)
    synTableStore(rbind(phdat2,phdat),'Cytokine-induced Drug Sensitivity Phosphoproteomics Unnormalized')

 
}

