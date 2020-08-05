## single drug experience. 

library(amlresistancenetworks)

if(!exists('dataLoaded')){
  source('beatAMLdata.R')
  dataLoaded=TRUE
}
#' How do we assess the efficacy of genes/transcripts/proteins in a single drug? 
#' 
#' 
#'First get the drugs
#'
#'then  cluster by top predictors

getFeaturesFromString<-function(string,prefix='mRNALevels'){
  stringr::str_remove_all(string,prefix)%>%stringr::str_split(';')%>%unlist()
}



getAllPreds<-function(){
  
  print("Getting phospho preds")
  substrate.dat<-pat.phos%>%
    dplyr::select(`AML sample`='Sample',Gene='site',Phosphosite='LogFoldChange')
  phospho.reg.results<-drugMolRegression(auc.dat,substrate.dat,'Phosphosite')
  phospho.rf.results<-drugMolRandomForest(auc.dat,substrate.dat,'Phosphosite')
  
  
  
  reg.preds<-purrr::map_df(list(mRNA='mRNALevels',
                                protein='proteinLevels',
                                gene='geneMutations'),~ drugMolRegression(auc.dat,
                                                                          pat.data,
                                                                          .x,category='Condition'))
  
  # cor.preds<-purrr::map_df(list(mRNA='mRNALevels',
  #                               protein='proteinLevels',
  #                               gene='geneMutations'),~ computeAUCCorVals(drug.dat,pat.data,.x))
  
  rf.preds<-purrr::map_df(list(mRNA='mRNALevels',
                               protein='proteinLevels',
                               gene='geneMutations'),~ drugMolRandomForest(auc.dat,
                                                                           pat.data,
                                                                           .x,category='Condition'))
  
  lv.df<-lv.df%>%rename(Gene='Latent_Variable',`AML sample`='AML_sample',`Latent Variable`='Loading')
  lv.reg.results<-drugMolRegression(auc.dat,lv.df,'Latent Variable')
  lv.rf.results<-drugMolRandomForest(auc.dat,lv.df,'Latent Variable')
  kinase.dat<-pat.kin%>%
    rename(`AML sample`='Sample',Gene='Kinase',KinaseExpr='meanLFC')
  kin.reg.results=drugMolRegression(auc.dat,kinase.dat,'KinaseExpr')
  kin.rf.results=drugMolRandomForest(auc.dat,kinase.dat,'KinaseExpr')
  
  
  #now bind them together with MSE and num Features to plot
  full.results<-rbind(reg.preds,lv.reg.results,kin.reg.results,phospho.reg.results)%>%
    mutate(method='LASSO')
  rf.results<-rbind(rf.preds,lv.rf.results,kin.rf.results,phospho.rf.results)%>%mutate(method='RandomForest')
  
  full.results<-rbind(full.results,rf.results)
  saveRDS(full.results,'mostlyCompletePredictions.rds')
 
  return(full.results)
  
}

plotAllPreds<-function(allPreds){
  library(ggplot2)
  #plot dotplot of MSE vs. num genes
  #shape is type of predictor
  #color is 
}

plotMostVarByDrug<-function(drugName,data,mostVar=50){
  library(pheatmap)
  
  fname=paste0(gsub(' ','',drugName),'_',mostVar,'mostVariable_',data,'vals.pdf')
  print(fname)
  
  ##ASSSUMES AUC.DAT is global
  
  if(data%in%c('proteinLevels','mRNALevels','geneMutations'))
    data.mat<-pat.data%>%rename(value=data)
  else if(data=='Latent Variable')
    data.mat<-lv.df%>%rename(Gene='Latent_Variable',`AML sample`='AML_sample',value='Loading')
  else if(data=='Kinase')
    data.mat<-pat.kin%>%rename(Gene='Kinase',`AML sample`='Sample',value='meanLFC')
  else
    data.mat<-pat.phos%>%dplyr::select(Gene='site',`AML sample`='Sample',value='LogFoldChange')
  
  drug.dat<-subset(auc.dat,Condition==drugName)%>%
    select(-c(Condition,medAUC,percAUC,overallSurvival,ageAtDiagnosis))%>%
    mutate(sensitive=if_else(AUC<100,'Sensitive','Resistant'))%>%
    distinct()%>%
    tibble::column_to_rownames('AML sample')
 
  pat.df<-data.mat%>%select('AML sample','Gene','value')%>%distinct()
  
  var.vals<-pat.df%>%ungroup()%>%group_by(Gene)%>%
    summarize(var=var(value))%>%
    arrange(desc(var))%>%
    select(Gene)
  print(var.vals$Gene[1:mostVar])
  
  pat.mat<-data.mat%>%select('AML sample','Gene','value')%>%
    subset(Gene%in%var.vals$Gene[1:mostVar])%>%
    subset(`AML sample`%in%rownames(drug.dat))%>%
    pivot_wider(values_from='value',
                names_from='AML sample', 
                values_fill =list(value=0.0),
                values_fn=list(value=mean))%>%
    tibble::column_to_rownames('Gene')%>%as.matrix()
  
  zvals<-which(apply(pat.mat,2,var)==0)
  if(length(zvals>0))
    pat.mat<-pat.mat[,-zvals]
  #print(pat.mat)
  #cluster all selections
  if(data=='mRNALevels')
    pat.mat<-log10(pat.mat+0.01)
  try(pheatmap::pheatmap(pat.mat,cellwidth = 10,cellheight=10,annotation_col = drug.dat,
                         clustering_distance_cols = 'correlation',
                         clustering_method = 'ward.D2',filename=fname))
  
}

#' selects genes from 
clusterSingleDrugEfficacy<-function(drugName='Doramapimod (BIRB 796)',
                                    meth='RandomForest',
                                    data='proteinLevels'){
  ##get gene,transcript, protein predictor
  library(pheatmap)
  
  fname=paste0(gsub(' ','',drugName),'_',meth,'selected_',data,'preds.pdf')
  print(fname)
  
  drug.dat<-NULL

    ##ASSSUMES AUC.DAT is global
    drug.dat<-subset(auc.dat,Condition==drugName)%>%
    select(-c(Condition,medAUC,percAUC,overallSurvival,ageAtDiagnosis))%>%
      mutate(sensitive=if_else(AUC<100,'Sensitive','Resistant'))%>%
      distinct()%>%
    tibble::column_to_rownames('AML sample')

  ##ASSUMES full.results is global
  oneRow=subset(new.results,var==drugName)%>%
    subset(method==meth)%>%
    subset(Molecular==data)%>%distinct()
 
  genes=getFeaturesFromString(oneRow$genes,data)
  
  print(genes)
 ##get dimensionality-reduced samples
  if(data%in%c('proteinLevels','mRNALevels','geneMutations'))
    data.mat<-pat.data%>%rename(value=data)
  else if(data=='Latent Variable')
    data.mat<-lv.df%>%rename(Gene='Latent_Variable',`AML sample`='AML_sample',value='Loading')
  else if(data=='KinaseExpr')
    data.mat<-pat.kin%>%rename(Gene='Kinase',`AML sample`='Sample',value='meanLFC')
  else
    data.mat<-pat.phos%>%dplyr::select(Gene='site',`AML sample`='Sample',value='LogFoldChange')
  
  pat.mat<-data.mat%>%select('AML sample','Gene','value')%>%
    subset(Gene%in%genes)%>%
    subset(`AML sample`%in%rownames(drug.dat))%>%
    pivot_wider(values_from='value',
                names_from='AML sample', 
                values_fill =list(value=0.0),
                values_fn=list(value=mean))%>%
    tibble::column_to_rownames('Gene')%>%as.matrix()
  
  zvals<-which(apply(pat.mat,2,var)==0)
  if(length(zvals>0))
    pat.mat<-pat.mat[,-zvals]
  #print(pat.mat)
  #cluster all selections
  if(data=='mRNALevels')
    pat.mat<-log10(pat.mat+0.01)
  
  try(pheatmap::pheatmap(pat.mat,cellwidth = 10,cellheight=10,annotation_col = drug.dat,
                         clustering_distance_cols = 'correlation',
                     clustering_method = 'ward.D2',filename=fname))
  
  plotMostVarByDrug(drugName,data)
  return(fname)
  
}

getPreds<-function(){
  library(ggplot2)
  full.results<-getAllPreds()
  new.results<-full.results%>%mutate(reducedData=Molecular%in%c('Latent Variable','KinaseExpr'))
  new.results<-subset(new.results,numFeatures>0)
  p1<-ggplot(new.results,aes(x=numFeatures,y=MSE,col=Molecular,shape=method,size=numSamples,alpha=0.7))+geom_point()+facet_grid(~reducedData)+scale_x_log10()
  ggsave('predictorSummary.png',p1,width=10)
  
  p2<-ggplot(new.results,aes(x=var,y=MSE,col=Molecular,size=numSamples,shape=method))+
    geom_point()+theme(axis.text.x=element_text(angle=90, hjust=1))
  ggsave('predictorDotPlot.png',p2,width=10)
  
  p3<-ggplot(new.results,aes(x=method,y=MSE,fill=Molecular))+geom_boxplot()
  ggsave('dataComparison.png',p3)
  return(new.results)
#plot all preds
#arbitrary filter
}
plotAllAUCs(auc.dat,'AUC')

#new.results<-getPreds()

subset(new.results,var%in%auc.dat$Condition)%>%
  subset(numFeatures>1)%>%
  rowwise()%>%mutate(clusterSingleDrugEfficacy(var,method,Molecular))


#subset set of preds