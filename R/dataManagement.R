##load data from files into rdata files
#this file should only be modified rarely when new data are to be added



#' reads in metadata from excel spreadsheet and tidies it
#' @require readxl
#' @require dplyr
#' @require reticulate
#'

readAndTidyMetadata<-function(){
  syn=synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22136312')$path)%>%
    stringr::string_replace()

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
readAndTidyQuizProtMeasures<-function(){
  metadata<-readAndTidyQuizMetadata()
  syn=synapseLogin()
  dat<-read.table(syn$get('syn22136310')$path,sep='\t',header=T)
  library(tidyr)
  quiz.data<-tidyr::pivot_longer(dat,-c(Protein,Gene),"Sample")%>%
    dplyr::left_join(metadata,by='Sample')%>%
    mutate(value=tidyr::replace_na(value,0))

  synTableStore(quiz.data,'Quizartinib Resistance Proteomics Data')
  
  return(quiz.data)
}



#' reads and tidies phosphoprotein measurements
#' joins with metadata
#' saves tidied data and returns
#' @export
#' @require dplyr
#
readAndTidyQuizPhosphoProtMeasures<-function(){
  metadata<-readAndTidyQuizMetadata()
  syn=synapseLogin()
  dat<-read.table(syn$get('syn23407326')$path,sep='\t',header=T)
  library(tidyr)
  quiz.phospho.data<-dat%>%
    tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),"Sample")%>%
    dplyr::left_join(metadata,by='Sample')%>%
    mutate(value=tidyr::replace_na(value,0))
  # subset(!is.na(value))
  
  synTableStore(quiz.phospho.data,'Quizartinib Resistance Phosphoproteomics Data')
  return(quiz.phospho.data)
  
}

#' reads and tidies protein measurements
#' joins with metadata
#' saves tidied data and returns
#' @export
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
#' @export
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

#' reads and tidies protein measurements from quizartinib treatment
#' joins with metadata
#' saves tidied data and returns
#' @export
#' @require dplyr
#' @require tidyr
#
readAndTidyQuizProtMeasures<-function(){
  metadata<-readAndTidyQuizMetadata()
  library(tidyr)
  quiz.data<-tidyr::pivot_longer(dat,-c(Entry_name,Gene),"Sample")%>%
      dplyr::left_join(metadata,by='Sample')
  return(quiz.data)

}


#' reads and tidies experiment 11 metadata
#' @require dplyr
#' @require readxl
#'
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

#' reads and tidies phosphototeomic data for experimnet 11
#' @require dplyr
#' @require tidyr
#' @export
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
#' @export
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

##################################BEATAML PATIENT samPles##################

getPatientTranscript<-function(patientlist){
  library(dplyr)
  library(readxl)
  syn=synapseLogin()
  #we dont need the RPKM because we have the CPM
  exp.3.megafile=syn$get('syn22130786')$path
  patient.cpm<-readxl::read_xlsx(exp.3.megafile,sheet='Table S9-Gene Counts CPM')%>%
    tidyr::pivot_longer(-c(Gene,Symbol),names_to='patient',values_to='transcriptCounts')%>%
    mutate(countMetric='CPM')%>%
    select(-Gene)%>%
    rename(Gene='Symbol',`AML sample`='patient')

  subset(patient.cpm, `AML sample`%in%patientlist)

}

getPatientVariants<-function(patientlist){
  library(dplyr)
  library(readxl)
  syn=synapseLogin()
  exp.3.megafile=syn$get('syn22130786')$path
  gene.var<-readxl::read_xlsx(exp.3.megafile,sheet='Table S7-Variants for Analysis')%>%
    subset(labId%in%patientlist)%>%
    select(labId,t_vaf,symbol)%>%distinct()%>%
    #tidyr::pivot_longer(c(t_vaf,n_vaf),names_to='Metric',values_to='Value')%>%
    rename(`AML sample`='labId',`Tumor VAF`='t_vaf',Gene='symbol')
  return(gene.var)
}

#'add in wave3/4 genetics
#'@import readxl
#'@import dplyr
addWave34GeneticsToTable<-function(){
  syn=synapseLogin()
  newtab<-readxl::read_xlsx(syn$get('syn22130786')$path,sheet=9)%>%
    dplyr::select(Gene=symbol,`AML sample`=original_id,`Tumor VAF`=t_vaf)%>%distinct()
  newtab

}


#' this grabs drug responses from three different sources, all on synapse
#' @param patientlist it will only return drugs for these patients
#' @return data frame of drug response values
#' @import dplyr
#' @import readxl
getPatientDrugResponses<-function(patientlist){
  library(dplyr)
  library(readxl)
  syn=synapseLogin()
    exp.3.megafile=syn$get('syn22130786')$path

    #this file has the primary drug responses
  dose.response<-readxl::read_xlsx(exp.3.megafile,sheet='Table S10-Drug Responses')%>%
    tidyr::pivot_longer(c(ic50,auc),names_to='Metric',values_to='Value')%>%
    dplyr::rename(`AML sample`='lab_id',Condition='inhibitor')

  #here is an additional set of data
  other.data<-readAndTidySensMetadata()%>%
    dplyr::select(-Barcode)

  #even more data
  extra.files<-c('syn22170222','syn22170223','syn22170225','syn22170226','syn22170227')
  more.data<-purrr::map_df(extra.files,function(x){
    fname=syn$get(x)$path
    if(length(grep('xlsx',fname))>0)
      ex<-readxl::read_xlsx(fname)
    else
      ex<-readxl::read_xls(fname)

    ex%>%dplyr::select(`AML sample`='Specimen: Lab ID',
                    Condition='Inhibitor Panel Definition: Drug',
                    AUC='Probit Interpretation: Area Under Curve',
                    IC50='Probit Interpretation: IC50')%>%
      tidyr::pivot_longer(c(IC50,AUC),names_to='Metric',values_to='Value')
  })
  comb.response<-rbind(dose.response,other.data,more.data)
  return(subset(comb.response,`AML sample`%in%patientlist))

}



#' getCuratedMutstatus
#' @export
#' @import dplyr
#' @import readxl
#' @import tidyr
getCuratedMutStatus<-function(){
  library(dplyr)
  syn<-synapseLogin()
  new.file<-syn$get('syn23538805')$path
#  exp.3.megafile=syn$get('syn22130786')$path
  
  patData<-readxl::read_xlsx(new.file,sheet='wv1to4')%>%
    dplyr::select(`AML sample`='labId','FLT3-ITD',NPM1,`FLT3-MUT`='FLT3')%>%distinct()%>%
      tidyr::pivot_longer(-`AML sample`,names_to='variant',values_to='status')%>%
      tidyr::replace_na(list(status='not_tested'))

  
# patData<-readxl::read_xlsx(exp.3.megafile,sheet='Clinical Summary')%>%
#  subset(labId%in%unlist(patients))%>%
#  select(`AML sample`='labId',contains("FLT3"))%>%
#  distinct()%>%
#   pivot_longer(-`AML sample`,values_to='FLT3_Status',names_to='FLT3_Assay')%>%
#   separate(FLT3_Status,into=c("FLT3_Status","FLT3_Details"),sep=' ',extra='merge')
  # mutate(Value=tidyr::replace_na(Value,0))
 synTableStore(patData,'BeatAML Pilot Mutation Status')

}
 
#' getPatientMetadata
#' @export
#' @import dplyr
#' @import readxl
getPatientMetadata<-function(synid='syn22170540'){
    require(dplyr)
#    synapser::synLogin()
  syn<-synapseLogin()
    exp.3.megafile=syn$get('syn22130786')$path

  patients<-readxl::read_xlsx(exp.3.megafile,sheet='Sample Summary')%>%
    dplyr::select('Specimen ID')%>%distinct()
  patients=c(unlist(patients),'16-01254')

  drugs<-getPatientDrugResponses(unlist(patients))

 patData<-readxl::read_xlsx(exp.3.megafile,sheet='Clinical Summary')%>%
   subset(labId%in%unlist(patients))%>%
   select(`AML sample`='labId',gender,ageAtDiagnosis, priorMalignancyType,
          vitalStatus,overallSurvival,causeOfDeath)%>%
          distinct()%>%
          left_join(drugs)%>%
   mutate(Value=tidyr::replace_na(Value,0))
#   subset(!is.na(Value))

 if(!is.null(synid)){##we need to delete rows and re-uplooad
   synTableUpdate(patData,synid)
 }else{
  synTableStore(patData,'BeatAML Pilot Drug and Clinical Data')
 }

 return(patData)

}

#' store drug class for now..
#' @import dplyr
storeDrugClassInfo<-function(){
  syn=synapseLogin()
  beat.samps<-syn$get('syn22130788')$path
  drug.class<-readxl::read_xlsx(beat.samps,sheet='Table S11-Drug Families')%>%
    distinct()
  synTableStore(drug.class,'Drug Classes')
}

# get all molecular data for patient baselines
#'@export
getPatientMolecularData<-function(synid='syn22172602'){
  syn=synapseLogin()
    exp.3.megafile=syn$get('syn22130786')$path

  patients<-readxl::read_xlsx(exp.3.megafile,sheet='Sample Summary')%>%
    dplyr::select('Specimen ID')%>%distinct()
  patients=c(unlist(patients),'16-01254')
  rna<-getPatientTranscript(unlist(patients))
  variants<-getPatientVariants(unlist(patients))
  moreVars<-addWave34GeneticsToTable()
  variants<-rbind(variants,moreVars)
  prots<-getPatientBaselines()
  patientMolecularData<-prots%>%
    full_join(rna,by=c('AML sample','Gene'),na_matches="never")%>%
    full_join(variants,by=c('AML sample','Gene'),na_matches="never")%>%
    mutate(transcriptCounts=tidyr::replace_na(transcriptCounts,0))%>%
    mutate(`Tumor VAF`=tidyr::replace_na(`Tumor VAF`,0))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))

  if(!is.null(synid)){##we need to delete rows and re-uplooad
    synTableUpdate(patientMolecularData,synid)
    }else{
    synTableStore(patientMolecularData,'BeatAML Pilot Molecular Data')
  }

  return(patientMolecularData)
}

#' get proteomic data for beataml samples
#' @import dplyr
getPatientBaselines<-function(){
  library(dplyr)
  syn=synapseLogin()
  dat<-read.table(syn$get('syn22130778')$path,sep='\t',header=T)#'../../Projects/CPTAC/exp_3'/PTRC_baseline_global_std_ref_with_genes.txt,sep='\t',header=T)
  patientProtSamples<-dat%>%tidyr::pivot_longer(cols=c(5:ncol(dat)),names_to='Sample', 
                                                values_to='LogFoldChange')%>%
    dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
    dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
    dplyr::select(`AML sample`='Sample',Gene, LogFoldChange)%>%
    mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))

  #this data is stored elsewhere
  metadata<-readAndTidySensMetadata()%>%
    dplyr::select(`AML sample`,Barcode)%>%
    distinct()%>%
    mutate(Barcode=as.character(Barcode))

  dat2<-read.table(syn$get('syn22130842')$path,sep='\t',header=T)
  secondData<-dat2%>%
    tidyr::pivot_longer(cols=c(3:ncol(dat2)),names_to='Sample', 
                        values_to='LogFoldChange')%>%
    dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
    dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
    dplyr::select(Barcode='Sample',Gene, LogFoldChange)%>%
    mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))%>%
    left_join(metadata)%>%
    select(-Barcode)

  return(rbind(patientProtSamples,secondData))
}


#' @export
#' @import dplyr
getPatientPhosphoBaselines<-function(){
    library(dplyr)
    syn<-synapseLogin()                          
    dat<-read.table(syn$get('syn22130779')$path,sep='\t',header=T)

    patientPhosphoSamples<-dat%>%
        tidyr::pivot_longer(-c(Entry,Gene,site,Peptide,ids,Entry_name),"Sample",
                            values_to='LogFoldChange')%>%
        dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
        dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
        dplyr::select(Sample,Gene, site,Peptide,LogFoldChange)%>%
      mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
      mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))

    synTableStore(patientPhosphoSamples,'BeatAML Pilot Phosphoproteomics')

    return(patientPhosphoSamples)
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
#' @export
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

#' @export
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
#' @export
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



#' getSorafenibSamples
#' Gets sorafenib trreated samples and stores to synaspe tables
#' @import dplyr
#' @import readxl
#' @import stringr
#' @export
getSorafenibSamples<-function(){
  library(dplyr)
  syn<-synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22314059')$path)%>%
    dplyr::select(specimen='Specimen ID',`FLT ITD`, `Sorafenib sensitivity`,`Sorafenib IC50`, `Cell number`,`Treatment`)%>%
    rowwise()%>%mutate(`AML sample`=stringr::str_split_fixed(specimen,stringr::fixed('.'),2)[1])%>%
    subset(!is.na(Treatment))

  pdat<-read.csv2(syn$get('syn22313435')$path,sep='\t')%>%
    dplyr::select(-ids)%>%
    tidyr::pivot_longer(-Gene,values_to='LogFoldChange',names_to='specIds')%>%
    mutate(LogFoldChange=as.numeric(as.character(LogFoldChange)))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))%>%
    mutate(specimen=unlist(stringr::str_replace(specIds,'X','')))%>%
    mutate(specimen=stringr::str_replace(specimen,'\\.','-'))%>%
    left_join(metadata)%>%
    select(-c(specIds,specimen))%>%distinct()%>%
    subset(!is.na(`AML sample`))

  phdat<-read.csv2(syn$get('syn22313433')$path,sep='\t')%>%
  dplyr::select(-Entry,ids)%>%
    tidyr::pivot_longer(-c(Gene,site,Peptide),values_to='LogFoldChange',names_to='specIds')%>%
    mutate(LogFoldChange=as.numeric(as.character(LogFoldChange)))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))%>%
    mutate(specimen=stringr::str_replace(specIds,'X',''))%>%
    mutate(specimen=stringr::str_replace(specimen,'\\.','-'))%>%
    left_join(metadata)%>%select(-c(specIds,specimen))%>%distinct()%>%subset(!is.na(`AML sample`))

  synTableStore(pdat,'Sorafenib treated proteomics')
  synTableStore(phdat,'Sorafenib treated Phospho-proteomics')
}
