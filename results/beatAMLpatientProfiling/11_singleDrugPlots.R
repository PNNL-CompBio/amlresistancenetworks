## single drug experience. 

library(amlresistancenetworks)

if(!exists('dataLoaded')){
  amlresistancenetworks::loadBeatAMLData()
  dataLoaded=TRUE
}
#' How do we assess the efficacy of genes/transcripts/proteins in a single drug? 
#' 
#' 
#'First get the drugs
#'
#'then  cluster by top predictors

#'getAllPreds as the name gets all types of predictors
#'based on all types of data/models
#'
getAllPreds<-function(doExtra=FALSE){
  
#  print('getting network preds')
#  mut.net.df<-mut.nets%>%select(mutationNetworkDistance='distance',Gene='Community',`AML sample`)
#  mn.reg.results<-drugMolRegression(auc.dat,mut.net.df,'mutationNetworkDistance')
#  mn.rf.results<-drugMolRandomForest(auc.dat,mut.net.df,'mutationNetworkDistance')
#  mn.lr.results<-drugMolLogReg(auc.dat,mut.net.df,'mutationNetworkDistance')
  extra.reg.results<-NULL
  extra.lr.results<-NULL
  if(doExtra){
    print('Getting proteomic network preds')
    prot.net.df<-prot.nets%>%dplyr::select(proteomicNetworkDistance='distance',Gene='Community',`AML sample`)
    pn.reg.results<-drugMolRegression(auc.dat,prot.net.df,'proteomicNetworkDistance')
    pn.rf.results<-drugMolRandomForest(auc.dat,prot.net.df,'proteomicNetworkDistance')
    pn.lr.results<-drugMolLogReg(auc.dat,prot.net.df,'proteomicNetworkDistance')
  
    print('getting LV preds')
    lv.df<-lv.df%>%rename(Gene='Latent_Variable',`AML sample`='AML_sample',`Latent Variable`='Loading')
    lv.reg.results<-drugMolRegression(auc.dat,lv.df,'Latent Variable')
  #  lv.rf.results<-drugMolRandomForest(auc.dat,lv.df,'Latent Variable')
    lv.lr.results<-drugMolLogReg(auc.dat,lv.df,'Latent Variable')
                                                                               
   print('getting kinase preds') 
    kinase.dat<-pat.kin%>%
      rename(`AML sample`='Sample',Gene='Kinase',KinaseExpr='meanLFC')
    kin.reg.results=drugMolRegression(auc.dat,kinase.dat,'KinaseExpr')
#    kin.rf.results=drugMolRandomForest(auc.dat,kinase.dat,'KinaseExpr')
    kin.lr.results<-drugMolLogReg(auc.dat,kinase.dat,'KinaseExpr')
  
    extra.reg.results<-rbind(pn.reg.results,lv.reg.results,kin.reg.results)%>%
      mutate(method='LASSO')
    extra.lr.results<-rbind(pn.lr.results,lv.lr.results,kin.lr.results)%>%
      mutate(method='LogisticReg')
  }
  
  print("Getting phospho preds")
  
  substrate.dat<-pat.phos%>%
    dplyr::select(`AML sample`='Sample',Gene='site',Phosphosite='LogFoldChange')
  phospho.reg.results<-drugMolRegression(auc.dat,substrate.dat,'Phosphosite')
#  phospho.rf.results<-drugMolRandomForest(auc.dat,substrate.dat,'Phosphosite')
  phospho.lr.results<-drugMolLogReg(auc.dat,substrate.dat,'Phosphosite')
  
  print('getting full preds')
  logr.preds<-purrr::map_df(list(mRNA='mRNALevels',
                                 protein='proteinLevels',
                                 gene='geneMutations'),~drugMolLogReg(auc.dat,pat.data,.x,category='Condition'))
  
  reg.preds<-purrr::map_df(list(mRNA='mRNALevels',
                                protein='proteinLevels',
                                gene='geneMutations'),~ drugMolRegression(auc.dat,
                                                                          pat.data,
                                                                          .x,category='Condition'))
#  rf.preds<-purrr::map_df(list(mRNA='mRNALevels',
#                               protein='proteinLevels',
#                               gene='geneMutations'),~ drugMolRandomForest(auc.dat,
#                                                                           pat.data,
#                                                                           .x,category='Condition'))
   #now bind them together with MSE and num Features to plot
  full.results<-rbind(reg.preds,phospho.reg.results)%>%
    mutate(method='LASSO')%>%rbind(extra.reg.results)

#  rf.results<-rbind(rf.preds,lv.rf.results,kin.rf.results,phospho.rf.results,mn.rf.results,pn.rf.results)%>%
#    mutate(method='RandomForest')
  
  lr.results<-rbind(logr.preds,phospho.lr.results)%>%
    mutate(method='LogisticReg')%>%
    rbind(extra.lr.results)%>%
    mutate(MSE=MSE*10000)
  
  full.results<-rbind(full.results,lr.results)
  saveRDS(full.results,'mostlyCompletePredictions.rds')
 
  return(full.results)
  
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
    dplyr::select(-c(Condition,medAUC,percAUC,overallSurvival,RNA,phosphoSites,mutations,proteins,ageAtDiagnosis))%>%
    mutate(sensitive=if_else(AUC<100,'Sensitive','Resistant'))%>%
    distinct()%>%
    tibble::column_to_rownames('AML sample')
 
  pat.df<-data.mat%>%dplyr::select('AML sample','Gene','value')%>%distinct()
  
  var.vals<-pat.df%>%ungroup()%>%group_by(Gene)%>%
    summarize(var=var(value))%>%
    arrange(desc(var))%>%
    dplyr::select(Gene)
  
  print(var.vals$Gene[1:mostVar])
  
  pat.mat<-data.mat%>%dplyr::select('AML sample','Gene','value')%>%
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

#'
#'selectDataMatAndPlot matches the output of a predictor (new.results) to the AUC data
#'and original molecular data to create heatmap
#'and combines ith with the specified data type and passes along a method to be used
#'@param drugName
#'@param meth
#'@param data
#'@return output of prediction
selectDataMatAndPlot<-function(drugName,meth,data,doEnrich=FALSE){
  #get dimensionality-reduced samples
  if(data%in%c('proteinLevels','mRNALevels','geneMutations') &&!is.null(pat.data))
    data.mat<-pat.data%>%dplyr::rename(value=data,Sample='AML sample')
  else if(data=='Latent Variable'&&!is.null(lv.df))
    data.mat<-lv.df%>%dplyr::rename(Gene='Latent_Variable',Sample='AML_sample',value='Loading')
  else if(data=='KinaseExpr'&&!is.null(pat.kin))
    data.mat<-pat.kin%>%dplyr::select(Gene='Kinase',Sample,value='meanLFC')
  else if(data=='Phosphosite'&&!is.null(pat.phos)){
    data.mat<-pat.phos%>%dplyr::select(Gene='site',Sample,value='LogFoldChange')
  }else if(data=='proteomicNetworkDistance'&&!is.null(prot.nets))
    data.mat<-prot.nets%>%dplyr::select(Gene='Community',value='distance',Sample=`AML sample`)
#  else if(data=='mutationNetworkDistances' && !is.null(mut.nets))
#    data.mat<-mut.nets%>%dplyr::select(Gene='Community',value='distance',Sample=`AML sample`)
  else{
    print(paste("Do not have data for",data))
    return(NULL)
  }
  
  auc.d<-auc.dat%>%
    dplyr::select(-c(medAUC,percAUC,overallSurvival,ageAtDiagnosis))%>%
    dplyr::rename(Sample='AML sample')

    clusterSingleDrugEfficacy(drugName,meth,data,doEnrich=doEnrich,auc.dat=auc.d,auc.thresh=100,
                            new.results,data.mat)  
}



getPreds<-function(){
  library(ggplot2)
  if(!file.exists('mostlyCompletePredictions.rds'))
    full.results<-getAllPreds()
  else
     full.results<-readRDS('mostlyCompletePredictions.rds')
    
  new.results<-full.results%>%
    mutate(reducedData=Molecular%in%c('proteomicNetworkDistance',
                                      'mutationNetworkDistance','Latent Variable',
                                      'KinaseExpr'))%>%
    subset(method!='RandomForest')
  
  new.results<-subset(new.results,numFeatures>0)
  p1<-ggplot(new.results,aes(x=numFeatures,y=MSE,col=Molecular,shape=method,
                             size=numSamples,alpha=0.7))+geom_point()+facet_grid(~reducedData)+scale_x_log10()+scale_color_viridis_d()
  ggsave('predictorSummary.png',p1,width=10)
  
  p2<-ggplot(new.results,aes(x=var,y=MSE,col=Molecular,size=numSamples,shape=method))+
    geom_point()+theme(axis.text.x=element_text(angle=90, hjust=1))+scale_color_viridis_d()
  ggsave('predictorDotPlot.png',p2,width=10,height=10)
  
  p3<-ggplot(new.results,aes(x=method,y=MSE,fill=Molecular))+geom_boxplot()+scale_fill_viridis_d()
  ggsave('dataComparison.png',p3)
#  newer.result<-new.results%>%mutate(perSampleError=MSE/numSamples)
#  p4<-ggplot(newer.result,aes(x=method,y=perSampleError,fill=Molecular))+geom_boxplot()
#  ggsave('normalizedDataComparison.png',p4,width=10)
 return(new.results)
#plot all preds
#arbitrary filter
}

new.results<-getPreds()

with.class<-drug.class%>%dplyr::rename(var="Condition")%>%
  inner_join(new.results)

drug.nums<-with.class%>%subset(!is.na(MSE))%>%group_by(family)%>%
  summarize(numDrugs=n_distinct(var))%>%
  subset(numDrugs>1)

#pclass <- with.class%>%
#  subset(family%in%drug.nums$family)%>%
#  subset(method%in%c('LogisticReg','LASSO'))%>%
#  ggplot(aes(x=family,y=MSE,fill=Molecular))+geom_boxplot()+facet_grid(method~.)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#ggsave('predictorPerfByClass.png',pclass,width=12)

#pclass2 <- with.class%>%
#  subset(family%in%drug.nums$family)%>%
 # subset(!reducedData)%>%
  # subset(method%in%c('LogisticReg','LASSO'))%>%
#  ggplot(aes(x=family,y=MSE,fill=Molecular))+geom_boxplot()+facet_grid(method~.)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#ggsave('predictorPerfByClassOrigData.png',pclass2,width=12)
doPlot=TRUE

if(doPlot){
  e.res<-subset(new.results,var%in%auc.dat$Condition)%>%
    subset(numFeatures>1)%>%
    mutate(doEnrich=TRUE)%>%
    rowwise()%>%mutate(enrich=list(selectDataMatAndPlot(var,method,Molecular,doEnrich)))

}
#subset set of preds