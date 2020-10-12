

##eval sigs on depmap
#dm<-query(eh,'depmap')

##first get eh <- ExperimentHub()
#all the samples, filtering for AML 


library(amlresistancenetworks)
formatCTRPData<-function(){
  syn=synapseLogin()
  cells<-read.csv(syn$get('syn22990030')$path,sep='\t')
  drugs<-read.csv(syn$get('syn22990031')$path,sep='\t')
 
  cell_line<-read.csv(syn$get('syn22978057')$path)%>%
    dplyr::select(Sample='stripped_cell_line_name',sex,age)
  
   mapping<-read.csv(syn$get('syn22990033')$path,sep='\t')%>%
    left_join(cells,by='master_ccl_id')
  
  auc<-read.csv(syn$get('syn22989927')$path,sep='\t')%>%
    left_join(drugs,by='master_cpd_id')%>%
    left_join(mapping,by='experiment_id')
  auc%>%subset(ccle_primary_hist=='haematopoietic_neoplasm')%>%
    dplyr::select(Sample='ccl_name',Condition='cpd_name',auc='area_under_curve')%>%
    left_join(cell_line)%>%
    tidyr::pivot_longer(cols='auc',names_to='Metric',values_to='Value')
  
}


formatSangerData<-function(){
  syn=synapseLogin()
  cell_line<-read.csv(syn$get('syn22978057')$path)%>%
    dplyr::select(stripped_cell_line_name,sex,COSMIC_ID='COSMICID',age)
  sanger<-read.csv(syn$get('syn22990888')$path)%>%
    left_join(cell_line,by='COSMIC_ID')%>%
    dplyr::select(Condition='DRUG_NAME',Sample='stripped_cell_line_name',auc,log2.ic50,ec50,sex,age)%>%
    tidyr::pivot_longer(cols=c('auc','log2.ic50','ec50'),
                        names_to='Metric',values_to='Value')
  return(sanger)
  
}

getStoreCellLineData<-function(){
  #' PRimary function that pulls all the data together. 
  #' Keeping for documentation but should only be run onece!
  
  drugDat<-rbind(formatSangerData()%>%dplyr::mutate(source='Sanger'),
                 formatCTRPData()%>%dplyr::mutate(source='CTRP'))%>%
    tidyr::replace_na(list(age=0,Value=0))
  library("depmap")
  library("ExperimentHub")
  eh <- ExperimentHub()
  
  #then get sdrug, gene, protein, x
  sampdata<-eh[['EH3086']]%>%
    subset(lineage=='leukemia')
  #moresamps<-eh[['EH2266']]
#Taking the latest one for now
#expr<-eh[['EH2264']]%>%
#  subset(depmap_id%in%sampdata$depmap_id)
#moreexpr<-eh[['EH2554']]%>%
#  subset(depmap_id%in%sampdata$depmap_id)
  
evenmoreexp<-eh[['EH3084']]%>%
  subset(depmap_id%in%sampdata$depmap_id)%>%
  left_join(sampdata)%>%
  dplyr::select(Gene='gene_name',transcriptCounts='expression',
                fullsamp='stripped_cell_line_name')%>%
  tidyr::separate(fullsamp,'_',into=c("Sample","rest"))%>%
  dplyr::select(-c(rest))

#muts<-eh[['EH2265']]%>%
#  subset(depmap_id%in%sampdata$depmap_id)
#moremuts<-eh[['EH2555']]%>%
#  subset(depmap_id%in%sampdata$depmap_id)
evenmoremuts<-eh[['EH3085']]%>%
  subset(depmap_id%in%sampdata$depmap_id)%>%left_join(sampdata)%>%
  dplyr::select(Gene='gene_name',Sample='stripped_cell_line_name',var_annotation)%>%
  group_by(Gene,Sample)%>%
  dplyr::mutate(numMuts=n())%>%dplyr::select(-var_annotation)


#now the proteomics
prots<-read.csv(syn$get('syn22975106')$path)
hemcols<-colnames(prots)[grep('HAEMA',colnames(prots))]
othercols<-c('Gene_Symbol')
prot.data<-prots%>%dplyr::select(c(othercols,hemcols))%>%
  tidyr::pivot_longer(cols=hemcols,values_to='LogFoldChange',names_to='fullsamp')%>%
  dplyr::rename(Gene='Gene_Symbol')%>%
  tidyr::separate(fullsamp,'_',into=c("Sample","rest"))%>%
  dplyr::select(-c(rest))


all.mol<-full_join(evenmoremuts,evenmoreexp,by=c('Gene','Sample'))%>%
  full_join(prot.data,by=c("Gene","Sample"))%>%
  tidyr::replace_na(list(numMuts=0,LogFoldChange=0,transcriptCounts=0))
#now we need to merge them into a single table and store them
  synTableStore(all.mol,'Cell Line Molecular Data',parentId='syn22128879')
  synTableStore(drugDat,'Cell Line Drug Data',parentId='syn22128879')

}


##get aml
loadBeatAMLData()

syn=synapseLogin()

## get depmap
cl.mol.dat<-syn$tableQuery('select * from  syn23004114')$asDataFrame()%>%
  dplyr::rename(mRNALevels='transcriptCounts',proteinLevels='LogFoldChange',geneMutations='numMuts')

#create manual mapping from AML to cell Line
drug.mapping <- syn$tableQuery('select * from syn23193914 where PTRC is not NULL')$asDataFrame()

##now match the drugs between cl and tumors
cl.auc.dat<-syn$tableQuery('select * from syn23004543')$asDataFrame()%>%
  subset(Metric=='auc')%>%
  subset(Condition%in%drug.mapping$drugName)%>%
  mutate(drugName=unlist(Condition))%>%
  left_join(drug.mapping,by=c('drugName','source'))%>%
  dplyr::select(Sample,PTRC,Value,source)%>%
  subset(Sample%in%cl.mol.dat$Sample)


##for each drug -  build a model in AML, eval in sanger and in PTRC

##now train model on AML and eval on depmap data
reg.preds<-purrr::map_df(list(mRNA='mRNALevels',
                              protein='proteinLevels',
                              gene='geneMutations'),~ drugMolRegression(auc.dat,
                                                                        pat.data,
                                                                        .x,category='Condition'))



