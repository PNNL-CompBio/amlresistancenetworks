##load data from files into rdata files

#' reads in metadata from excel spreadsheet and tidies it
#' @require readxl
#' 

readAndTidyMetadata<-function(){
  metadata<-readxl::read_xlsx('../../Projects/CPTAC/LabelingMetadata_exp12.xlsx')
 
  #update sample
  metadata$Sample<-sapply(metadata$`Sample ID`,function(x) paste('Sample',x,sep='_'))
  
  metadata<-metadata%>%
      tidyr::separate(col='Sample Info',into=c('treatment','cellLine','flask','ligand'),sep=', ',fill="right")%>%
      dplyr::select(Sample,cellLine,treatment,flask,ligand)
  
  #update cell line
  cl<-subset(metadata,cellLine=='10M')%>%select(Sample,treatment)%>%mutate(cellLine=stringr::str_replace(treatment,'parental ',''))%>%select(Sample,cellLine)%>%mutate(treatment='None',ligand='None')
 
  lig.dat<-metadata%>%subset(!is.na(ligand))%>%select(c(flask,cellLine,ligand))
  
  fin=metadata%>%
    subset(cellLine!='10M')%>%
    dplyr::select(cellLine,treatment,flask,Sample)%>%
    left_join(lig.dat,by=c('flask','cellLine'))%>%
    select(-flask)%>%
    rbind(cl)
  
  
  return(fin)
  
}


##load data from files into rdata files

#' reads in metadata from excel spreadsheet and tidies it
#' @require readxl
#' 

readAndTidyQuizMetadata<-function(){
  metadata<-readxl::read_xlsx('../../Projects/CPTAC/quizartLabelingMetadata.xlsx',skip=3)
  
  res<-metadata%>%
      rowwise()%>%
      mutate(SampleTemp=paste(Plex,`Sample #`,sep='.'))%>%
      mutate(Sample=stringr::str_replace(SampleTemp,'Sample ',''))%>%
      dplyr::select(Sample,`Sample Name`,Ligand)
  
  fin<-res%>%separate(`Sample Name`,into=c(CellLine,Flask),by=',')
  return(fin)
  
}

#' reads and tidies protein measurements
#' joins with metadata
#' saves tidied data and returns
#' @export
#
readAndTidyProtMeasures<-function(){
  metadata<-readAndTidyMetadata()
  dat<-read.table('../../ex12_data/ptrc_ex12_global_median_centered_with_genes.txt',sep='\t',header=T)
  library(tidyr)
  gilteritinib.data<-tidyr::pivot_longer(dat,-c(Protein,Gene),"Sample")%>%dplyr::left_join(metadata,by='Sample')
saveRDS(gilteritinib.data,file='inst/gilteritinibData.Rds')
  return(gilteritinib.data)

    }



#' reads and tidies protein measurements from quizartinib treatment
#' joins with metadata
#' saves tidied data and returns
#' @export
#
readAndTidyQuizProtMeasures<-function(){
  metadata<-readAndTidyQuizMetadata()
  dat<-read.table('../../Projects/CPTAC/ptrc_exp12_global_round2_with_genes.txt',sep='\t',header=T)
  library(tidyr)
  quiz.data<-tidyr::pivot_longer(dat,-c(Entry_name,Gene),"Sample")%>%dplyr::left_join(metadata,by='Sample')
  saveRDS(gilteritinib.data,file='inst/quizartinibData.Rds')
  return(gilteritinib.data)
  
}

