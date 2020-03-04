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



#' Reads in data frame and computes fold change
#' @param data.frame
#' @return tidied data frame with p-values and fold change
#' @import dplyr
#' @export
computeFoldChangePvals<-function(gilt.data,
                                  treatmentCond='Early Gilteritnib',
                                 ligands=c("FL","FGF2")){
  
  data<-gilt.data%>%
    subset(treatment%in%(c('None',treatmentCond)))%>%
    dplyr::select(CellLine,Gene,value,ligand)%>%
    subset(!is.na(value))
  

  #remove those data that have only 1 replicate  
  reps<-data%>%
    group_by(CellLine,Gene,ligand)%>%
    mutate(reps=length(value))%>%
    subset(reps>1)%>%
    ungroup()
  
  res<-lapply(ligands,function(lig){
    #compute the significance for FL-modulaetd changes in both cell lines
    vals=reps%>%
      group_by(CellLine,Gene)%>%
      subset(ligand%in%c('None',lig))%>%
      mutate(ligand=as.factor(ligand))%>%
      mutate(conds=n_distinct(ligand))%>%
      subset(conds>1)
    
    signif=vals%>%
      group_modify(~infer::t_test(.,value~ligand,order=c(lig,"None")),keep=TRUE)%>%
      dplyr::select(c(CellLine,Gene,p_value))%>%
      ungroup()%>%
      group_by(CellLine)%>%
      mutate(p_adj=p.adjust(p_value))%>%
      ungroup()%>%
      rename(p_value=paste(lig,'pVal',sep='_'),
             p_adj=paste(lig,'pAdj',sep='_'))
    
    
  })
  ##now compute the significance for FGF2 treated cells
  

  mean.diffs<-data%>%group_by(CellLine,Gene,ligand)%>%
    summarize(meanVal=mean(value))%>%
    tidyr::pivot_wider(names_from=ligand,values_from = meanVal)%>%
    mutate(FGF2_to_parent=FGF2-None,FL_to_parent=FL-None)%>%
    full_join(fl.signif,by=c('CellLine','Gene'))%>%
    full_join(fgf.signif,by=c('CellLine','Gene'))
}