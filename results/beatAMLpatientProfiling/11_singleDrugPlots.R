## single drug experience. 

library(amlresistancenetworks)
source('beatAMLdata.R')
#' How do we assess the efficacy of genes/transcripts/proteins in a single drug? 
#' 
#' 
#'First get the drugs
#'
#'then  cluster by top predictors

getFeaturesFromString<-function(string,prefix='mRNALevels'){
  stringr::str_remove_all(string,prefix)%>%stringr::str_split(',')%>%unlist()
}

plotDrugResponse<-function(pred.df,drugName){
  #get drug data for annotaitons
  
}

getAllPreds<-function(){
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
    rename(`AML sample`='Sample',Gene='Kinase',SubstrateExpr='meanLFC')
  kin.reg.results=drugMolRegression(auc.dat,kinase.dat,'SubstrateExpr')
  kin.rf.results=drugMolRandomForest(auc.dat,kinase.dat,'SubstrateExpr')
  
  #now bind them together with MSE and num Features to plot
  full.results<-rbind(reg.preds,lv.reg.results,kinase.reg.results)%>%
    mutate(method='LASSO')
  rf.results<-rbind(rf.preds,lv.rf.results,kin.rf.results)%>%mutate(method='RandomForest')
  
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

clusterSingleDrugEfficacy<-function(drugName='Doramapimod (BIRB 796)',meth='RandomForest',data='proteinLevels'){
  ##get gene,transcript, protein predictor
  library(pheatmap)
  
  fname=paste0(gsub(' ','',drugName),'_',meth,'selected_',data,'preds.pdf')
  print(fname)
  
  drug.dat<-NULL
  #try(
    drug.dat<-subset(auc.dat,Condition==drugName)%>%
    select(-c(Condition,medAUC,percAUC))%>%
      distinct()%>%
    tibble::column_to_rownames('AML sample')
   # )
#if(is.null(drug.dat))
#    return(fname)
  
  oneRow=subset(full.results,var==drugName)%>%
    subset(method==meth)%>%
    subset(Molecular==data)%>%distinct()
 
 # print(oneRow)
  genes=getFeaturesFromString(oneRow$genes,data)
  
  print(genes)
  ##get dimensionality-reduced samples
  pat.mat<-pat.data%>%select('AML sample','Gene',data)%>%
    subset(Gene%in%genes)%>%
    subset(`AML sample`%in%rownames(drug.dat))%>%
    rename(value=data)%>%
    pivot_wider(values_from='value',
                names_from='AML sample', 
                values_fill =list(value=0.0),
                values_fn=list(value=mean))%>%
    tibble::column_to_rownames('Gene')%>%as.matrix()
  
  #print(dim(pat.mat))
  #cluster all selections
  try(pheatmap::pheatmap(pat.mat,cellwidth = 10,cellheight=10,annotation_col = drug.dat,
                     clustering_distance_cols = 'correlation',clustering_method = 'ward.D2',filename=fname))
  return(fname)
  
}

full.results<-getAllPreds()
ggplot(full.results,aes(x=numFeatures,y=MSE,col=Molecular,shape=method,size=numSamples))+geom_point()+facet_grid(~method)
ggsave('predictorSummary.png')
#plot all preds

#subset set of preds