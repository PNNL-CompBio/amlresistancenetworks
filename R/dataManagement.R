##load data from files into rdata files

# Reads in metadata from excel spreadsheet and tidies it
# @requires(readxl)
#
readAndTidyMetadata<-function(){
  metadata<-readxl::read_xlsx('../../LabelingMetadata_exp12.xlsx')
 
  #update sample
  metadata$Sample<-sapply(metadata$`Sample ID`,function(x) paste('Sample',x,sep='_'))
  
  #update cell line
  metadata$CellLine<-rep("Reference",nrow(metadata))
  metadata$CellLine[grep('MOLM14',metadata$`Sample Info`)]<-'MOLM14'
  metadata$CellLine[grep('MV411',metadata$`Sample Info`)]<-'MV411'
  
  #update treatment
  metadata$treatment<-rep("None",nrow(metadata))
  metadata$treatment[grep('Early',metadata$`Sample Info`)]<-'Early Gilteritinib'
  metadata$treatment[grep('Late',metadata$`Sample Info`)]<-'Late Gilteritinib'
  
  
  #update ligand
  metadata$ligand<-rep('None',nrow(metadata))
  metadata$ligand[grep('FGF2',metadata$`Sample Info`)]<-'FGF2'
  metadata$ligand[grep('FLT3',metadata$`Sample Info`)]<-'FL'
  
  return(metadata)
  
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
