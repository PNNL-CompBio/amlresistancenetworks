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



#' reads and tidies phosphoprotein measurements
#' joins with metadata
#' saves tidied data and returns
#' @export
#' @require dplyr
#
readAndTidyPhosphoProtMeasures<-function(){
  metadata<-readAndTidyMetadata()
  dat<-read.table('../../projects/CPTAC/ex12_data/ptrc_ex12_phospho_median_centered_stoich_with_sites.txt',sep='\t',header=T)
  library(tidyr)
  gilt.phospho.data<-dat%>%
    tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),"Sample")%>%
    dplyr::left_join(metadata,by='Sample')

    saveRDS(gilt.phospho.data,file='inst/giltPhosphoData.Rds')
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
  dat<-read.table('../../Projects/CPTAC/ex11_data/ptrc_ex11_kurtz_2plex_phospho_d2_stoich_with_sites.txt',sep='\t',header=T)
  gilt.sens.pdata<-dat%>%
    tidyr::pivot_longer(-c(Entry,Gene,site,Peptide,ids),"Sample")%>%
    dplyr::mutate(Barcode=as.numeric(stringr::str_replace(Sample,"X","")))%>%
    dplyr::left_join(metadata,by='Barcode')

  saveRDS(gilt.sens.pdata,file='inst/giltPhosphoSensData.Rds')
  return(gilt.sens.pdata)

}

#' reads and tidies bulk proteomic data for experiment 11
#' @require dplyr
#' @require tidyr
#' @export
#'
readAndTidySensProtMeasure<-function(){
  metadata<-readAndTidySensMetadata()
  dat<-read.table('../../Projects/CPTAC/ex11_data/ptrc_ex11_kurtz_2plex_global_d2_with_genes.txt',sep='\t',header=T)

  drugSensData<-dat%>%tidyr::pivot_longer(cols=c(3:ncol(dat)),names_to='Sample', values_to='LogFoldChange')%>%
    dplyr::mutate(Barcode=as.numeric(stringr::str_replace(Sample,"X","")))%>%
    dplyr::select(Barcode,Gene, LogFoldChange)%>%
    dplyr::left_join(metadata,by='Barcode')

  saveRDS(drugSensData,file='inst/gilteritinibSensitivityData.Rds')
  return(drugSensData)
}

##################################BEATAML PATIENT samPles
#


getPatientTranscript<-function(patientlist){
  library(dplyr)
  library(readxl)
  synapser::synLogin()
  #we dont need the RPKM because we have the CPM
  #patient.rpkm<-readxl::read_xlsx(exp.3.megafile,sheet='BeatAML S8 Gene Counts RPKM')%>%
  #  tidyr::pivot_longer(-c(Gene,Symbol),names_to='patient',values_to='counts')%>%
                                        #  mutate(countMetric='RPKM')
  exp.3.megafile=synapser::synGet('syn22130786')$path
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
  synapser::synLogin()
  exp.3.megafile=synapser::synGet('syn22130786')$path
  gene.var<-readxl::read_xlsx(exp.3.megafile,sheet='Table S7-Variants for Analysis')%>%
    subset(labId%in%patientlist)%>%
    select(labId,t_vaf,symbol)%>%
    #tidyr::pivot_longer(c(t_vaf,n_vaf),names_to='Metric',values_to='Value')%>%
    rename(`AML sample`='labId',`Tumor VAF`='t_vaf',Gene='symbol')
  return(gene.var)
}

getPatientDrugResponses<-function(patientlist){
  library(dplyr)
  library(readxl)
  synapser::synLogin()
  exp.3.megafile=synapser::synGet('syn22130786')$path

  dose.response<-readxl::read_xlsx(exp.3.megafile,sheet='BeatAML S10 Drug Responses')%>%
    tidyr::pivot_longer(c(ic50,auc),names_to='Metric',values_to='Value')%>%
    dplyr::rename(`AML sample`='lab_id',Condition='inhibitor')

  other.data<-readAndTidySensMetadata()%>%
    dplyr::select(-Barcode)

  comb.response<-rbind(dose.response,other.data)
  return(subset(comb.response,`AML sample`%in%patientlist))

}
#' getPatientMetadata
#' @export
#' @require dplyr
#' @require readxl
getPatientMetadata<-function(){
    require(dplyr)
    synapser::synLogin()
    exp.3.megafile=synapser::synGet('syn22130786')$path

  patients<-readxl::read_xlsx(exp.3.megafile,sheet='Sample Summary')%>%
    dplyr::select('Specimen ID')%>%distinct()

  drugs<-getPatientDrugResponses(unlist(patients))
 patData<-readxl::read_xlsx(exp.3.megafile,sheet='Clinical Summary')%>%
   subset(labId%in%unlist(patients))%>%
   select(`AML sample`='labId',gender,ageAtDiagnosis, priorMalignancyType,
          vitalStatus,overallSurvival,causeOfDeath)%>%distinct()%>%left_join(drugs)

    saveRDS(patData,file='inst/patientDrugAndClinical.Rds')
    synStore(File('inst/patientDrugAndClinical.Rds',parentId='syn22130776'))
 return(patData)

}

# get all molecular data for patient baselines
#'@export
getPatientMolecularData<-function(){
  patients<-readxl::read_xlsx(exp.3.megafile,sheet='Sample Summary')%>%
    dplyr::select('Specimen ID')%>%distinct()

  rna<-getPatientTranscript(unlist(patients))
  variants<-getPatientVariants(unlist(patients))
  prots<-getPatientBaselines()
  patientMolecularData<-rna%>%left_join(variants,by=c('AML sample','Gene'))%>%
    left_join(prots,by=c('AML sample','Gene'))

  saveRDS(patientMolecularData,file='inst/patientMolecularData.Rds')
  synStore(File('inst/patientMolecularData.Rds',parentId='syn22130776'))
  return(patientMolecularData)
}

#' get proteomic data for beataml samples
getPatientBaselines<-function(){
  library(dplyr)
#  metadata<-getPatientMetadata()
  dat<-read.table('../../Projects/CPTAC/exp_3/PTRC_baseline_global_std_ref_with_genes.txt',sep='\t',header=T)
  patientProtSamples<-dat%>%tidyr::pivot_longer(cols=c(5:ncol(dat)),names_to='Sample', values_to='LogFoldChange')%>%
    dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
    dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
    dplyr::select(`AML sample`='Sample',Gene, LogFoldChange)

  #saveRDS(patientProtSamples,file='inst/patientProtSampleData.Rds')
  return(patientProtSamples)
}


#' @export
#' @require synapser
#' @require dplyr
getPatientPhosphoBaselines<-function(){
  library(dplyr)
 # metadata<-getPatientMetadata()
  dat<-read.table('../../Projects/CPTAC/exp_3/PTRC_baseline_phospho_std_ref_with_sites_stoich.txt',sep='\t',header=T)
  patientPhosphoSamples<-dat%>%
    tidyr::pivot_longer(-c(Entry,Gene,site,Peptide,ids,Entry_name),"Sample",values_to='LogFoldChange')%>%
    dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
    dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
    dplyr::select(Sample,Gene, site,Peptide,LogFoldChange)

  saveRDS(patientPhosphoSamples,file='inst/patientPhosphoSampleData.Rds')
  return(patientPhosphoSamples)
}

#' @require synapser
#'
getTimeCourseMetadata<-function(){
  samp.data<-readxl::read_xlsx('../../Projects/CPTAC/exp_2/Pilot2_treated cell_TMT labels.xlsx')%>%
    dplyr::select(Sample=`Tube name`,cellLine=`Cell line`,treatment=`Treated with`,timePoint=`Time point`)
  return(samp.data)
}
#'@export
getTimeCourseData<-function(){
  library(dplyr)
  metadata<-getTimeCourseMetadata()
  dat<-read.csv2('../../Projects/CPTAC/exp_2/PTRC_pilot2_global_with_genes.txt',
                  sep='\t',header=T)
  timeCourseData<-dat%>%tidyr::pivot_longer(cols=c(4:ncol(dat)),names_to='sample', values_to='LogFoldChange')%>%
    #dplyr::mutate(specId=stringr::str_replace(sample,".","-"))%>%
    dplyr::mutate(Sample=stringr::str_replace(sample,stringr::fixed("."),"-"))%>%
    dplyr::select(Sample,Gene, LogFoldChange)%>%
    left_join(metadata)

  saveRDS(timeCourseData,file='inst/timeCourseData.Rds')
  return(timeCourseData)

}

#' @export
getTimeCoursePhosphoData<-function(){
  library(dplyr)
  metadata<-getTimeCourseMetadata()
  timeCoursePhospho<-read.csv2('../../Projects/CPTAC/exp_2/PTRC_pilot2_phospho_stoich_with_sites.txt',
                 sep='\t',header=T)%>%
    tidyr::pivot_longer(-c(Entry,Gene,site,Peptide,ids,Entry.name,Protein),"Sample",values_to='LogFoldChange')%>%
    dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
    dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
    dplyr::select(Sample,Gene, site,Peptide,LogFoldChange)%>%
    left_join(metadata)
  saveRDS(timeCoursePhospho,file='inst/timeCoursePhosphoData.Rds')
}
