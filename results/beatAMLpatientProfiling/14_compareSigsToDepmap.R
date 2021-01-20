

##eval sigs on depmap
#dm<-query(eh,'depmap')

##first get eh <- ExperimentHub()
#all the samples, filtering for AML 


library(amlresistancenetworks)
library(ggplot2)
library(dplyr)

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
  auc%>%#subset(ccle_primary_hist=='haematopoietic_neoplasm')%>%
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
 # library("depmap")
  library("ExperimentHub")
  eh <- ExperimentHub()
  
  #then get sdrug, gene, protein, x
  sampdata<-eh[['EH3086']]#%>%
    #subset(lineage=='leukemia')
  #moresamps<-eh[['EH2266']]
#Taking the latest one for now
#expr<-eh[['EH2264']]%>%
#  subset(depmap_id%in%sampdata$depmap_id)
#moreexpr<-eh[['EH2554']]%>%
#  subset(depmap_id%in%sampdata$depmap_id)
  red.samp<-sampdata%>%
    select(depmap_id,stripped_cell_line_name,lineage)%>%
    distinct()
evenmoreexp<-eh[['EH3084']]%>%
  subset(depmap_id%in%red.samp$depmap_id)%>%
  left_join(red.samp)%>%
  dplyr::select(Gene='gene_name',transcriptCounts='expression',
                fullsamp='stripped_cell_line_name',lineage)%>%
  tidyr::separate(fullsamp,'_',into=c("Sample","rest"))%>%
  dplyr::select(-c(rest))

#muts<-eh[['EH2265']]%>%
#  subset(depmap_id%in%sampdata$depmap_id)
#moremuts<-eh[['EH2555']]%>%
#  subset(depmap_id%in%sampdata$depmap_id)
evenmoremuts<-eh[['EH3085']]%>%
  dplyr::filter(depmap_id%in%red.samp$depmap_id)%>%left_join(red.samp)%>%
  dplyr::select(Gene='gene_name',Sample='stripped_cell_line_name',lineage,var_annotation)%>%
  group_by(Gene,Sample)%>%
  dplyr::mutate(numMuts=n())%>%dplyr::select(-var_annotation)


#now the proteomics
prots<-read.csv(syn$get('syn22975106')$path)
noncols<-c("Protein_Id",'Description','Group_ID','Uniprot','Uniprot_Acc','Gene_Symbol')
pepcols<-colnames(prots)[grep('Peptides',colnames(prots))]
allcols<-setdiff(colnames(prots),union(pepcols,noncols))
#hemcols<-colnames(prots)#[grep('HAEMA',colnames(prots))]
othercols<-c('Gene_Symbol')

prot.data<-prots%>%dplyr::select(c(othercols,allcols))%>%
  tidyr::pivot_longer(cols=allcols,values_to='LogFoldChange',names_to='fullsamp')%>%
  dplyr::rename(Gene='Gene_Symbol')%>%
  tidyr::separate(fullsamp,'_',into=c("Sample","rest"))%>%
  dplyr::select(-c(rest))


all.mol<-full_join(evenmoremuts,evenmoreexp,by=c('Gene','Sample','lineage'))%>%
  full_join(prot.data,by=c("Gene","Sample"))%>%
  tidyr::replace_na(list(numMuts=0,LogFoldChange=0,transcriptCounts=0))
#now we need to merge them into a single table and store them
breakpoint=15000000
mol1<-all.mol[1:breakpoint,]
mol2<-all.mol[breakpoint+1:nrow(all.mol),]
synTableStore(mol1,'All Cell Line Molecular Data',parentId='syn22128879')
synTableStore(mol2,'All Cell Line Molecular Data',parentId='syn22128879')
#  synTableStore(drugDat,'Full Cell Line Drug Data',parentId='syn22128879')

}


##get aml
if(!exists('data.loaded')||!data.loaded){
  loadBeatAMLData()
  data.loaded=TRUE
}
syn=synapseLogin()


selectDataMatAndPlotCTRP<-function(compound,method,Molecular){
  print(paste(c(compound,method,Molecular)))
  if(Molecular=='proteomicNetworkDistance')
    dat.mat<-rename(cell.net.df,value=Molecular)
  else
    dat.mat<-rename(cl.mol.dat,value=Molecular)
  amlresistancenetworks::clusterSingleDrugEfficacy(drugName=compound,
                                    meth=method,
                                      data=Molecular,
                                      auc.dat=select(ctrp.cl.auc,-c(Value,source)),
                                      auc.thresh=100,
                                        new.results=rename(ctrp.full.res,var='compound'),
                                      data.mat=dat.mat,
                                    prefix='CTRP')

}


selectDataMatAndPlotSanger<-function(compound,method,Molecular){
  print(paste(c(compound,method,Molecular)))
  if(Molecular=='proteomicNetworkDistance')
    dat.mat<-rename(cell.net.df,value=Molecular)
  else
    dat.mat<-rename(cl.mol.dat,value=Molecular)
  
  amlresistancenetworks::clusterSingleDrugEfficacy(drugName=compound,
                                                   meth=method,
                                                   data=Molecular,
                                                   auc.dat=select(sanger.cl.auc,-c(Value,source)),
                                                   auc.thresh=100,
                                                   new.results=rename(sanger.full.res,var='compound'),
                                                   data.mat=dat.mat,
                                                   prefix='Sanger')
  
}


runAnalysis<-function(){

## get depmap
cells<-syn$tableQuery('select distinct Sample from  syn23004114 where ( LogFoldChange<>0 AND transcriptCounts<>0)')$asDataFrame()$Sample

cl.mol.dat<-syn$tableQuery('select * from syn23004114')$asDataFrame()%>%
  subset(Sample%in%cells)%>%
  dplyr::rename(mRNALevels='transcriptCounts',proteinLevels='LogFoldChange',
                geneMutations='numMuts')%>%mutate(Sample=unlist(Sample),Gene=unlist(Gene))

#create manual mapping from AML to cell Line
drug.mapping <- syn$tableQuery('select * from syn23193914 where PTRC is not NULL')$asDataFrame()

##now match the drugs between cl and tumors
cl.auc.dat<-syn$tableQuery('select * from syn23004543')$asDataFrame()%>%
  subset(Metric=='auc')%>%
  subset(Condition%in%drug.mapping$drugName)%>%
  mutate(drugName=unlist(Condition),Sample=unlist(Sample))%>%
  left_join(drug.mapping,by=c('drugName','source'))%>%
  dplyr::select(Sample,Condition='PTRC',Value,source)%>%
  subset(Sample%in%cl.mol.dat$Sample)%>%
  group_by(Sample,Condition,source)%>%summarize(Value=mean(Value))%>%unggroup()

sanger.cl.auc<-cl.auc.dat%>%
  subset(source=='Sanger')%>%
  mutate(AUC=(Value*100+50))

ctrp.cl.auc<-cl.auc.dat%>%
  subset(source=='CTRP')%>%
  mutate(AUC=Value*10)


###let's summarize the data here
mol.sum<- cl.mol.dat%>%group_by(Sample)%>%summarize(mutations=length(which(geneMutations>0)),proteomics=length(which(proteinLevels!=0)),mRNAseq=length(which(mRNALevels!=0)))
drug.sum<-cl.auc.dat%>%
  group_by(Sample)%>%
  summarize(numDrugs=n_distinct(Condition),numSources=n_distinct(source))
full.sum<-mol.sum%>%full_join(drug.sum)

cell.net <-querySynapseTable("syn23481350")
  filtered = cell.net%>% 
    dplyr::filter(net2_type=='community')%>%
    dplyr::filter(hyp1=='patients')%>%
    dplyr::filter(hyp2=='panCan')

cell.net.df<<-filtered%>%
    select(Community='net2',distance,Sample='net1')%>%
    select(proteomicNetworkDistance='distance',Gene='Community',Sample)

prot.net.df<-prot.nets%>%select(proteomicNetworkDistance='distance',Gene='Community',`AML sample`)


pn.reg.results<-do.call(rbind,drugMolRegressionEval(auc.dat,prot.net.df,'proteomicNetworkDistance', 
                                      ctrp.cl.auc,cell.net.df))%>%as.data.frame()
pn.lr.results<-do.call(rbind,drugMolLogRegEval(auc.dat,prot.net.df,'proteomicNetworkDistance',
                                 ctrp.cl.auc,cell.net.df))%>%as.data.frame()%>%
  mutate(testMSE=unlist(testMSE)*10000)



##for each drug -  build a model in AML, eval in sanger and in PTRC

##now train model on AML and eval on depmap data
reg.preds<-purrr::map_df(list(mRNA='mRNALevels',
                              protein='proteinLevels',
                              gene='geneMutations'),~ drugMolRegressionEval(auc.dat,
                                                                        pat.data,
                                                                        .x,
                                                                        ctrp.cl.auc,
                                                                        cl.mol.dat,
                                                                        category='Condition'))


log.reg.preds<-purrr::map_df(list(protein='proteinLevels',
                                  mRNA='mRNALevels',
                              gene='geneMutations'),~ drugMolLogRegEval(auc.dat,
                                                                         pat.data,
                                                                         .x,
                                                                         ctrp.cl.auc,
                                                                         cl.mol.dat,
                                                                         category='Condition'))%>%
  mutate(testMSE=unlist(testMSE)*10000)
##i think we have all the results now, we can join them

ctrp.full.res<-rbind(mutate(log.reg.preds,method='LogisticRegression'), 
                mutate(reg.preds,method='LassoRegression'),
                mutate(pn.reg.results,method='LassoRegression'),
                mutate(pn.lr.results,method='LogisticRegression'))




##now do the same for Sanger
##now train model on AML and eval on depmap data
reg.preds<-purrr::map_df(list(mRNA='mRNALevels',
                              protein='proteinLevels',
                              gene='geneMutations'),~ drugMolRegressionEval(auc.dat,
                                                                            pat.data,
                                                                            .x,
                                                                            sanger.cl.auc,
                                                                            cl.mol.dat,
                                                                            category='Condition'))

pn.reg.results<-do.call(rbind,drugMolRegressionEval(auc.dat,prot.net.df,'proteomicNetworkDistance', 
                                      sanger.cl.auc, cell.net.df))%>%as.data.frame()
pn.lr.results<-do.call(rbind,drugMolLogRegEval(auc.dat,prot.net.df,'proteomicNetworkDistance',
                                 sanger.cl.auc, cell.net.df))%>%as.data.frame()%>%
  mutate(testMSE=unlist(testMSE)*10000)


log.reg.preds<-purrr::map_df(list(protein='proteinLevels',
                                  mRNA='mRNALevels',
                                  gene='geneMutations'),~ amlresistancenetworksdrugMolLogRegEval(auc.dat,
                                                                            pat.data,
                                                                            .x,
                                                                            sanger.cl.auc,
                                                                            cl.mol.dat,
                                                                            category='Condition'))%>%
  mutate(testMSE=unlist(testMSE)*10000)
##i think we have all the results now, we can join them
sanger.full.res<-rbind(mutate(log.reg.preds,method='LogisticRegression'), 
                     mutate(reg.preds,method='LassoRegression'),
                          mutate(pn.reg.results,method='LassoRegression'),
                mutate(pn.lr.results,method='LogisticRegression'))%>%
  mutate(compound=as.character(compound),Molecular=as.character(Molecular))
}



plotAnalysis<-function(){
  
  
ctrp.full.res%>%
  subset(numFeatures>1)%>%
  rowwise()%>%mutate(selectDataMatAndPlotCTRP(compound,method,Molecular))

sanger.full.res%>%
  subset(numFeatures>1)%>%
  rowwise()%>%mutate(selectDataMatAndPlotSanger(compound,method,Molecular))

new.results<-rbind(mutate(sanger.full.res,source='Sanger'),
                   mutate(ctrp.full.res,source='CTRP'))%>%
  mutate(MSE=unlist(MSE),testMSE=unlist(testMSE),resCor=unlist(resCor),numFeatures=unlist(numFeatures),
         numSamples=unlist(numSamples),compound=unlist(compound),Molecular=unlist(Molecular))


new.results<-new.results%>%subset(numFeatures!=0)

p1<-ggplot(new.results,aes(x=numFeatures,y=testMSE,col=Molecular,shape=method,
                           size=numSamples,alpha=0.7))+geom_point()+facet_grid(~source)+scale_y_log10()+scale_color_viridis_d()
ggsave('cellLinepredictorSummary.png',p1,width=10)

p2<-ggplot(new.results,aes(x=compound,y=testMSE,col=Molecular,size=numSamples,shape=method))+
  geom_point()+facet_grid(source~.)+theme(axis.text.x=element_text(angle=90, hjust=1))+scale_y_log10()+scale_color_viridis_d()
ggsave('cellLinepredictorDotPlot.png',p2,width=10,height=10)

p3<-ggplot(new.results,aes(x=method,y=testMSE,fill=Molecular))+geom_boxplot()+
  facet_grid(source~.)+scale_y_log10()+scale_color_viridis_d()
ggsave('cellLinedataComparison.png',p3,width=10)
}
