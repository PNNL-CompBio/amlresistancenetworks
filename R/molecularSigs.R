



#' drugMolRandomForest
#' builds random forest predictor of molecular data
#' @param clin.data
#' @param mol.data
#' @param mol.feature
#' @param category
#' @export
drugMolRandomForest<-function(clin.data,
                              mol.data,
                              mol.feature,
                              category='Condition'){
  
  
  if(length(mol.feature)==1){
    drug.mol<-clin.data%>%
      dplyr::select(`AML sample`,var=category,AUC)%>%
      group_by(`AML sample`,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(dplyr::select(mol.data,c(Gene,`AML sample`,!!mol.feature)),
                by='AML sample')
    
    
    reg.res<-drug.mol%>%group_by(var)%>%
      group_modify(~ miniForest(.x,mol.feature),.keep=T)%>%
      mutate(Molecular=mol.feature)
  }else{
    drug.mol<-clin.data%>%
      dplyr::select(`AML sample`,var=category,AUC)%>%
      group_by(`AML sample`,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(dplyr::select(mol.data,c(Gene,`AML sample`,mol.feature)),
                by='AML sample')
    reg.res<-drug.mol%>%group_by(var)%>%
      group_modify(~ combForest(.x,mol.feature),.keep=T)%>%
      mutate(Molecular=paste(mol.feature,collapse='_'))
  }
  return(reg.res)
  
}


#' working on this function still
#' The goal is to use a basic elastic net regression to identify how 
#' well each molecular feature predicts outcome as well as how many features
#' are selected
#' @import glmnet
#' @param clin.data is tidied clinical data
#' @param mol.data is tidied table of molecular data
#' @param mol.feature is name of column to select from mol.data
#' @param category can be either Condition or family
#' @export 
drugMolRegression<-function(clin.data,
                            mol.data,
                            mol.feature,
                            category='Condition'){
  
   
  if(length(mol.feature)==1){
    drug.mol<-clin.data%>%
      dplyr::select(`AML sample`,var=category,AUC)%>%
       group_by(`AML sample`,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(dplyr::select(mol.data,c(Gene,`AML sample`,!!mol.feature)),
                by='AML sample')
  
  
    reg.res<-drug.mol%>%group_by(var)%>%
      group_modify(~ miniReg(.x,mol.feature),.keep=T)%>%
      mutate(Molecular=mol.feature)
  }else{
    drug.mol<-clin.data%>%
      dplyr::select(`AML sample`,var=category,AUC)%>%
      group_by(`AML sample`,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(dplyr::select(mol.data,c(Gene,`AML sample`,mol.feature)),
                by='AML sample')
    reg.res<-drug.mol%>%group_by(var)%>%
      group_modify(~ combReg(.x,mol.feature),.keep=T)%>%
      mutate(Molecular=paste(mol.feature,collapse='_'))
  }
  return(reg.res)
  
}

#' drugMolLogReg
#' Computes logistic regression based on AUC threshold of 100
#' @param tab
#' @param feature
#' @param aucThresh
#' @export
drugMolLogReg<-function(clin.data, 
                        mol.data,
                        mol.feature,
                        category='Condition',
                        aucThresh=100){
  if(length(mol.feature)==1){
    drug.mol<-clin.data%>%
      dplyr::select(`AML sample`,var=category,AUC)%>%
      group_by(`AML sample`,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(dplyr::select(mol.data,c(Gene,`AML sample`,!!mol.feature)),
                by='AML sample')%>%
      mutate(sensitive=meanVal<aucThresh)
    
    
    reg.res<-drug.mol%>%group_by(var)%>%
      group_modify(~ miniLogR(.x,mol.feature),.keep=T)%>%
      mutate(Molecular=mol.feature)
  }else{
    drug.mol<-clin.data%>%
      dplyr::select(`AML sample`,var=category,AUC)%>%
      group_by(`AML sample`,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(dplyr::select(mol.data,c(Gene,`AML sample`,mol.feature)),
                by='AML sample')
    reg.res<-drug.mol%>%group_by(var)%>%
      group_modify(~ combDE(.x,mol.feature),.keep=T)%>%
      mutate(Molecular=paste(mol.feature,collapse='_'))
  }
  return(reg.res)
  
}

#' miniDE
#' Carries out differential expressiona nalysis using an AUC of 100
#' @param tab
#' @param mol.features
#' 
miniLogR<-function(tab,mol.feature){
#  irst build our feature matrix
  library(glmnet)
  set.seed(10101)

 mat<-buildFeatureMatrix(tab,mol.feature)
 
  #print(mat)
  if(is.null(dim(mat)))
    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=0,corVal=0))
  
 ##remove uninformative features, make sure we have enough at the end of it (at least 5)  
  cm<-apply(mat,1,mean)
  vm<-apply(mat,1,var)
  zvals<-union(which(cm==0),which(vm==0))
  if(length(zvals)>0)
    mat<-mat[-zvals,]
    
  zcols<-apply(mat,2,var)
  zvals<-which(zcols==0)
   #sprint(zvals)
   if(length(zvals)>0)
    mat<-mat[,-zvals]
    
   if(ncol(mat)<5 || nrow(mat)<5)
      return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=nrow(mat)),corVal=0)
  
  print(paste("Found",length(zvals),'features with no',mol.feature,'data across',
              ncol(mat),'features'))    

  #now collect our y output variable - AUC
  tmp<-tab%>%
     dplyr::select(sensitive,`AML sample`)%>%
     distinct()
  yvar<-tmp$sensitive

  names(yvar)<-tmp$`AML sample`
  yvar<-unlist(yvar[rownames(mat)])
  
  #use CV to get minimum error
  cv.res<-NULL
  try(cv.res<-cv.glmnet(x=mat,y=yvar,family='binomial',
                        type.measure='mse',nfolds=length(yvar)))
  if(is.null(cv.res))
    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=nrow(mat)))
  
  
  best.res<-data.frame(lambda=cv.res$lambda,MSE=cv.res$cvm)%>%
    subset(MSE==min(MSE))
  
  #then select how many elements
  full.res<-NULL
  try(full.res<-glmnet(x=mat,y=yvar,family='binomial',type.measure='mse'))
  
  if(is.null(full.res))
    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=nrow(mat)),corVal=0)
  preds <- predict.glmnet(full.res,newx=mat)[,which(full.res$lambda==best.res$lambda)]
  
  genes=names(which(full.res$beta[,which(full.res$lambda==best.res$lambda)]!=0))
  genelist<-paste(genes,collapse=';')
  cv = cor(preds,yvar,use='pairwise.complete.obs',method='spearman')
 # print(cv)
  return(data.frame(MSE=best.res$MSE,numFeatures=length(genes),
                    genes=genelist,
                    numSamples=length(yvar),
         corVal=cv))
}

#'combForest
#'Runs random forest on combination of feature types
#'@param feature.list
#'@param tab
#'@export
#'@return a data frame with 3 values
combForest<-function(tab,feature.list=c('proteinLevels','mRNAlevels','geneMutations')){
  comb.mat<-do.call('cbind',lapply(feature.list,function(x) buildFeatureMatrix(tab,x)))
  
  if(ncol(comb.mat)<5 || nrow(comb.mat)<5)
    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=nrow(comb.mat)))
  
  #now collect our y output variable
  tmp<-tab%>%
    dplyr::select(meanVal,`AML sample`)%>%
    distinct()
  yvar<-tmp$meanVal
  names(yvar)<-tmp$`AML sample`
  yvar<-unlist(yvar[rownames(comb.mat)])
               
  rf<-randomForest(comb.mat,yvar)
  # rf.pred<-predict(rf,mat)
  # mse=mean((rf.pred-yvar)^2)
  return(data.frame(MSE=min(rf$mse),
                    numFeatures=length(which(rf$importance!=0)),
                    genes=paste(names(rf$importance)[which(rf$importance!=0)],collapse=';',),
                    numSamples=length(yvar)))
  
}

#' combReg
#' Runs lasso regression on a combination of feature types
#' @param tab
#' @export
#' @param feature.list
#' @return a data frame with three values/columns
combReg<-function(tab,feature.list=c('proteinLevels','mRNALevels','geneMutations')){
  
   comb.mat<-do.call('cbind',lapply(feature.list,function(x) buildFeatureMatrix(tab,x)))
   
   
  # cm<-apply(comb.mat,1,mean)
  # zvals<-which(cm==0)
  # if(length(zvals)>0)
  #   comb.mat<-comb.mat[-zvals,]
   
   
  if(ncol(comb.mat)<5 || nrow(comb.mat)<5)
    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=nrow(comb.mat)))
  
  #now collect our y output variable
  tmp<-tab%>%
    dplyr::select(meanVal,`AML sample`)%>%
    distinct()
  yvar<-tmp$meanVal
  names(yvar)<-tmp$`AML sample`
  yvar<-unlist(yvar[rownames(comb.mat)])
  
  #use CV to get maximum AUC
  cv.res=cv.glmnet(x=comb.mat,y=yvar,type.measure='mse')
  best.res<-data.frame(lambda=cv.res$lambda,MSE=cv.res$cvm)%>%
    subset(MSE==min(MSE))
  
  #then select how many elements
  full.res<-NULL
  try(  full.res<-glmnet(x=comb.mat,y=yvar,type.measure='mse'))
  if(is.null(full.res))
    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=nrow(comb.mat)))
  genes=names(which(full.res$beta[,which(full.res$lambda==best.res$lambda)]!=0))
  genelist<-paste(genes,collapse=';')
  #print(paste(best.res$MSE,":",genelist))
  return(data.frame(MSE=best.res$MSE,numFeatures=length(genes),
                    genes=genelist,numSamples=length(yvar)))
  
  
}

#' buildFeatureMatrix
#' Builds a matrix for regression with rows as patients and columns as gene
#' @param tab
#' @param mol.feature
#' @param sampname
#' @return matrix
buildFeatureMatrix<-function(tab,mol.feature,sampname='AML sample'){
 # print(mol.feature)
  vfn=list(0.0)
  names(vfn)=mol.feature
  
  vfc<-list(mean)
  names(vfc)=mol.feature

  mat<-tab%>%
  dplyr::select(!!sampname,Gene,!!mol.feature)%>%
    subset(Gene!="")%>%
    tidyr::pivot_wider(names_from=Gene,values_from=mol.feature,
                     values_fill=vfn,values_fn = vfc,names_prefix=mol.feature)%>%
    tibble::column_to_rownames(sampname)

  tmat<-apply(mat,2,unlist)
  colnames(tmat)<-colnames(mat)
  rownames(tmat)<-rownames(mat)
  return(tmat)
}


#' miniForest
#' Runds random forest on the table and molecular feature of interest
#' @import randomForest
#' @export 
#' @return list
miniForest<-function(tab,mol.feature,quant=0.995){
  library(randomForest)
  
  #first build our feature matrix
  mat<-buildFeatureMatrix(tab,mol.feature)
  
  #rint(dim(mat))
  if(is.null(dim(mat)))
    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=0))
  
  cm<-apply(mat,1,mean)
  zvals<-which(cm==0.0)
  print(paste("Found",length(zvals),'features with no',mol.feature,'data across',ncol(mat),'features'))
  if(length(zvals)>0)
    mat<-mat[-zvals,]
  
  if(ncol(mat)<5 || nrow(mat)<5)
    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=nrow(mat)))
  
  #now collect our y output variable
  tmp<-tab%>%
    dplyr::select(meanVal,`AML sample`)%>%
    distinct()
  yvar<-tmp$meanVal
  names(yvar)<-tmp$`AML sample`
  yvar<-unlist(yvar[rownames(mat)])
  
  rf<-randomForest(mat,yvar)

  ##let's parse through the importance
 
  top5=quantile(rf$importance,quant)
  
  
  return(data.frame(MSE=min(rf$mse),
                    numFeatures=length(which(rf$importance>top5)),
                    genes=paste(rownames(rf$importance)[which(rf$importance>top5)],collapse=';'),
                    numSamples=length(yvar)))
         
  
}

#' miniReg
#' Runs lasso regression on a single feature from tabular data
#' @param tab with column names `AML sample`,meanVal,Gene, and whatever the value of 'mol.feature' is.
#' @param enet.alpha numeric vector specifying the alpha values to use when running glmnet
#' @export 
#' @return a data.frame with three values/columns: MSE, numFeatures, and Genes
miniReg<-function(tab,mol.feature, enet.alpha = c(1)){
  library(glmnet)
#  set.seed(1010101)
  set.seed(101010101)

  #first build our feature matrix
 mat<-buildFeatureMatrix(tab,mol.feature)
 #print(mat)
 if(is.null(dim(mat)))
   return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=0))
 
 ##remove uninformative values that can mess up regression
 cm<-apply(mat,1,mean)
 vm<-apply(mat,1,var)
 zvals<-union(which(cm==0),which(vm==0))
 if(length(zvals)>0)
   mat<-mat[-zvals,]
 
 
 zcols<-apply(mat,2,var)
 zvals<-which(zcols==0)
# print(zvals)
 if(length(zvals)>0)
   mat<-mat[,-zvals]
 
  ##let's assume we have 5 values shared!!!
  if(ncol(mat)<5 || nrow(mat)<5)
      return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=nrow(mat)))
 
  print(paste("Found",length(zvals),'features with no',mol.feature,'data across',ncol(mat),'features'))

  
  #now collect our y output variable - AUC
  tmp<-tab%>%
     dplyr::select(meanVal,`AML sample`)%>%
      distinct()
  yvar<-tmp$meanVal
  names(yvar)<-tmp$`AML sample`
  yvar<-unlist(yvar[rownames(mat)])
  
  best.res <- data.frame(lambda = numeric(0), 
                         MSE = numeric(0),
                         alpha = numeric(0))
  
  models <- list()
  ## Run glmnet for each alpha, saving the best lambda value every time
  for (alpha in enet.alpha) {
    model <- cv.glmnet(x = mat, y = yvar, alpha = alpha, 
                       type.measure = 'mse')
    best <- data.frame(lambda = model$lambda, MSE = model$cvm) %>%
      subset(MSE == min(MSE)) %>%
      mutate(alpha = alpha)
    best.res <- rbind(best.res, best)
    models[[as.character(alpha)]] <- model
  }
  
  ## Picking optimal (according to MSE) lambda and alpha
  best.res  <- best.res %>%
    subset(MSE == min(MSE))
  alpha = best.res$alpha %>%
    as.character()
  lambda = best.res$lambda
  full.res <- models[[alpha]]
  
  #then select how many elements
  preds <- predict(full.res,newx=mat, s = lambda)
  coefs <- coef(full.res, s = lambda)
  genes <- names(coefs[coefs[,1] != 0, ])[-1]
  genelist<-paste(genes,collapse=';')
  #print(paste(best.res$MSE,":",genelist))
  cv=cor(preds,yvar,use='pairwise.complete.obs',method='spearman')
 # print(cv)
  return(data.frame(alpha=alpha,MSE=best.res$MSE,numFeatures=length(genes),genes=genelist,
                    numSamples=length(yvar),corVal=cv))
}


#' how do we visualize the correlations of each drug/gene pair?
#' 
#' Computes drug by element correlation and plots them in various ways
#' @param cor.res
#' @param cor.thresh 
#' @import ggplot2
#' @import dplyr
#' @import cowplot
#' @import ggridges
#' @export
plotCorrelationsByDrug<-function(cor.res,cor.thresh){
  ##for each drug class - what is the distribution of correlations broken down by data type
  library(ggplot2)
  library(dplyr)
  do.p<-function(dat,cor.thresh){
    print(head(dat))
    fam=dat$family[1]
    #fam=dat%>%dplyr::select(family)%>%unlist()
    #fam=fam[1]
    fname=paste0(fam,'_correlations.png')
    p1<-ggplot(dat,aes(y=Condition,x=drugCor))+
      geom_density_ridges_gradient(aes(fill=feature,alpha=0.5))+
      scale_fill_viridis_d()+ggtitle(paste('Correlation with',fam))
    ##for each drug, how many genes have a corelation over threshold
    p2<-subset(dat,abs(drugCor)>cor.thresh)%>%
      ungroup()%>%
      group_by(Condition,feature,family)%>%
      summarize(CorVals=n_distinct(Gene))%>%ggplot(aes(x=Condition,y=CorVals,fill=feature))+
      geom_bar(stat='identity',position='dodge')+
      scale_fill_viridis_d()+ggtitle(paste("Correlation >",cor.thresh))+ 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    cowplot::plot_grid(p1,p2,nrow=2) 
    
    ggsave(fname)
    return(fname)
  }
  
  famplots<-cor.res%>%split(cor.res$family)%>%purrr::map(do.p,cor.thresh)                                                  
  
  lapply(famplots,synapseStore,'syn22130776')
  
}


#'Nots ure if this is used yet
#'Currently deprecated
computeDiffExByDrug<-function(sens.data){
  #samps<-assignSensResSamps(sens.data,'AUC',sens.val,res.val)
  
  result<-sens.data%>%
    dplyr::select('Sample',Gene,cellLine='Drug',value="LogFoldChange",treatment='Status')%>%
    distinct()%>%
    subset(!is.na(cellLine))%>%
    group_by(cellLine)%>%
    group_modify(~ computeFoldChangePvals(.x,control=NA,conditions =c("Sensitive","Resistant")),.keep=TRUE)%>%
    rename(Drug='cellLine')
  
  table(result%>%group_by(Drug)%>%subset(p_adj<0.1)%>%summarize(sigProts=n()))
  
  result
}


#' what molecules are
#' @import dplyr
#' @param clin.data
#' @param mol.data
#' @param mol.feature
#' @export  
computeAUCCorVals<-function(clin.data,mol.data,mol.feature){
  tdat<-mol.data%>%
    dplyr::select(Gene,`AML sample`,Mol=mol.feature)%>%
    subset(!is.na(Mol))%>%
    inner_join(clin.data,by='AML sample')
  
  print('here')
  dcors<-tdat%>%select(Gene,Mol,Condition,AUC)%>%
    distinct()%>%
    group_by(Gene,Condition)%>%
    mutate(numSamps=n(),drugCor=cor(Mol,AUC))%>%
    dplyr::select(Gene,Condition,numSamps,drugCor)%>%distinct()%>%
    arrange(desc(drugCor))%>%
    subset(numSamps>10)%>%
    mutate(feature=mol.feature)
  
  return(dcors)
  
}

#' getFeaturesFromString - separates out comma-delimited model features
#' @param string
#' @param prefix
#' @export
#' @return list of features
getFeaturesFromString<-function(string,prefix='mRNALevels'){
  stringr::str_remove_all(string,prefix)%>%stringr::str_split(';')%>%unlist()
}


#' selects AUC data and molecular data nand builds heatmap
#' of single drug and the results of the prediction
#' @param drugName Name of drug to select from prediction results
#' @param meth Method used to select predictors
#' @param data Type of data to be plotted
#' @param auc.dat AUC data from drug/samples
#' @param auc.thresh Threshold to determine if sample is sensitive
#' @param new.results Prediction data frame
#' @param data.mat Data frame of data to plot
#' @return filename of pdf
#' @export
clusterSingleDrugEfficacy<-function(drugName='Doramapimod (BIRB 796)',
                                    meth='RandomForest',
                                    data='proteinLevels',
                                    doEnrich=FALSE,
                                    auc.dat=auc.dat,
                                    auc.thresh=100,
                                    genes,
                                    data.mat=NULL,
                                    prefix=''){
  
  ##get gene,transcript, protein predictor
  library(pheatmap)
  library(wesanderson)
  pal<-wes_palette('Darjeeling1',100,type='continuous')
  
  fname=paste0(prefix,gsub(' ','',drugName),'_',meth,'selected_',data,'preds.pdf')
  print(fname)

  drug.dat<-subset(auc.dat,Condition==drugName)%>%
    mutate(SampleResponse=if_else(AUC<auc.thresh,'Sensitive','Resistant'))%>%
    dplyr::select(AUC,SampleResponse,Sample)%>%
    distinct()%>%
    tibble::column_to_rownames('Sample')

  pat.mat<-data.mat%>%select('Sample','Gene','value')%>%
    subset(Gene%in%genes)%>%
    subset(`Sample`%in%rownames(drug.dat))%>%
    pivot_wider(values_from='value',
                names_from='Sample', 
                values_fill =list(value=0.0),
                values_fn=list(value=mean))%>%
    tibble::column_to_rownames('Gene')%>%as.matrix()

  annote.colors<-list(SampleResponse=c(Sensitive='darkgrey',Resistant='white'))
 # names(annote.colors)<-setdiff(names(pat.vars),'overallSurvival')

  
  res=data.frame(ID='',Description='',pvalue=1.0,p.adjust=1.0) 
  #zvals<-which(apply(pat.mat,2,var)==0)
  #if(length(zvals>0))
  #  pat.mat<-pat.mat[,-zvals]
  if(nrow(pat.mat)<3)
    return(res)
  #print(pat.mat)
  #cluster all selections
  if(data=='mRNALevels')
    pat.mat<-log10(pat.mat+0.01)
  
  try(pheatmap::pheatmap(pat.mat,cellwidth = 10,cellheight=10,annotation_col = drug.dat,
                         clustering_distance_cols = 'euclidean', color=pal,annotation_colors=annote.colors,
                         clustering_method = 'ward.D2',filename=fname))

  if(doEnrich && length(rownames(pat.mat))>2){
    if(data=='Phosphosite')
      try(res<-doRegularKin(rownames(pat.mat)))
    else
      try(res<-doRegularGo(rownames(pat.mat)))
  }
  
  # plotMostVarByDrug(drugName,data)
  return(res)
  
}

