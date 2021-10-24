




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
drugMolRegressionEval<-function(clin.data,
                            mol.data,
                            mol.feature,
                            test.clin,
                            test.mol,
                            category='Condition'){
  
  ##first check to evaluate drug overlap
   drugs<-unlist(intersect(select(clin.data,cat=category)$cat,
                           select(test.clin,cat=category)$cat))
   print(paste('Found',length(drugs),'conditions that overlap between training and testing'))
   
    drug.mol<-clin.data%>%
      dplyr::select(`AML sample`,var=category,AUC)%>%
      subset(var%in%drugs)%>%
       group_by(`AML sample`,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(select(mol.data,c(Gene,`AML sample`,!!mol.feature)),
                by='AML sample')
  
    drug.test<-test.clin%>%
      dplyr::select(Sample,var=category,AUC)%>%
      subset(var%in%drugs)%>%
      group_by(Sample,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(select(test.mol,c(Gene,Sample,!!mol.feature)),by='Sample')
    
    reg.res<-lapply(unique(drug.mol$var),function(x){
      print(x)
      data.frame(miniRegEval(subset(drug.mol,var==x),subset(drug.test,var==x),mol.feature),
        compound=x,Molecular=mol.feature)})
  
  return(reg.res)
  
}

miniRegMod<-function(trainTab,mol.feature){
    #first build our feature matrix
  mat<-buildFeatureMatrix(trainTab,mol.feature)
  #print(mat)
  if(is.null(dim(mat)))
    return(ret.df)
  
  cm<-apply(mat,1,mean)
  vm<-apply(mat,1,var)
  zvals<-union(which(cm==0),which(vm==0))
  if(length(zvals)>0)
    mat<-mat[-zvals,]
  
  print(paste("Found",length(zvals),'patients with no',
              mol.feature,'data across',ncol(mat),'features'))
  
  zcols<-apply(mat,2,var)
  zvals<-which(zcols==0)
  # print(zvals)
  if(length(zvals)>0)
    mat<-mat[,-zvals]
  
  if(ncol(mat)<5 || nrow(mat)<5)
    return(ret.df)
  
  #now collect our y output variable for training
  tmp<-trainTab%>%
    dplyr::select(meanVal,`AML sample`)%>%
    distinct()
  yvar<-tmp$meanVal
  names(yvar)<-tmp$`AML sample`
  yvar<-unlist(yvar[rownames(mat)])
  cv.res=cv.glmnet(x=mat,y=yvar,type.measure='mse')
  
  best.res<-data.frame(lambda=cv.res$lambda,MSE=cv.res$cvm)%>%
    subset(MSE==min(MSE))
  
  #then select how many elements
  full.res<-glmnet(x=mat,y=yvar,type.measure='mse')
  
  genes=NULL
  try(genes<-names(which(full.res$beta[,which(full.res$lambda==best.res$lambda)]!=0)))
  
}

#' miniRegEval
#' Runs lasso regression on a single feature from tabular data and evaluates on a second set of data.
#' This is a little more complicated than LOO
#' analysis because we have to check for shared feature names
#' @param trainTab with column names `AML sample`,meanVal,Gene, and whatever the value of `mol.feature` is.
#' @param testTab with column names `Sample`, meanVal, Gene, and whatever the value of `mol.feature` is
#' @param mol.feature The molecular feature to be evaluated
#' @export 
#' @return a data.frame with three values/columns: MSE, numFeatures, and Genes
miniRegEval<-function(trainTab,testTab,mol.feature){
  library(glmnet)
  set.seed(10101)
  
  ret.df<-data.frame(MSE=0,testMSE=0,corVal=0,numFeatures=0,genes='',numSamples=0)
  
  tmat=NULL
  mat<-NULL
  
  try(mat<-buildFeatureMatrix(trainTab,mol.feature))
  
  try(tmat<-buildFeatureMatrix(testTab,mol.feature,'Sample'))
  
  
  if(is.null(mat)||is.null(tmat)||is.null(dim(mat)))
    return(ret.df)
  
  ##check for non-informative features
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
  
  if(ncol(mat)<5 || nrow(mat)<5)
    return(ret.df)
  
  print(paste("Found",length(zvals),'patients with no',
              mol.feature,'data across',ncol(mat),'features'))

  #now collect our y output variable for training
  tmp<-trainTab%>%
    dplyr::select(meanVal,`AML sample`)%>%
    distinct()
  yvar<-tmp$meanVal
  names(yvar)<-tmp$`AML sample`
  yvar<-unlist(yvar[rownames(mat)])
  
  #now get the y output for test
  ttmp<-testTab%>%
    dplyr::select(meanVal,Sample)%>%
    distinct()
  tyvar<-ttmp$meanVal
  names(tyvar)<-ttmp$Sample
  tyvar<-unlist(tyvar[rownames(tmat)])
  
  
  #use CV to get minimum MSE error
  cv.res=cv.glmnet(x=mat,y=yvar,type.measure='mse')
  best.res<-data.frame(lambda=cv.res$lambda,MSE=cv.res$cvm)%>%
    subset(MSE==min(MSE))
  
  ##now reduce test matrix to only those features in original model
  shared<-intersect(colnames(tmat),colnames(mat))
  missing<-setdiff(colnames(mat),colnames(tmat))
 
  print(paste('missing',length(missing),'genes and using',length(shared)))
  if(length(shared)==0){
    print("None of the genes are shared?")
    print(paste(colnames(tmat)[1:10],collapse=';'))
    return(ret.df)
  }
  
  ##remove features that aren't shared
  tmat<-tmat[,shared]
  if(length(missing)>0){
   #docker print(missing)
    newmat<-matrix(nrow=nrow(tmat),ncol=length(missing),data=0)
    colnames(newmat)<-missing
    tmat<-cbind(tmat,newmat)
  }
   
  #then build the full model
  full.res<-glmnet(x=mat,y=yvar,type.measure='mse')
  
  ##these are the genes that minimize the MSE
  genes=NULL
  try(genes<-names(which(full.res$beta[,which(full.res$lambda==best.res$lambda)]!=0)))

  if(is.null(genes))
    return(ret.df)
  genelist<-paste(genes,collapse=';')  

  ## get the prediction onto the new matrix (to assess correlation)
  t.res<-predict(full.res,newx=tmat,s=best.res$lambda)

  #use the assess function to get a new MSE
  res=assess.glmnet(full.res,newx=tmat,newy=tyvar,s=best.res$lambda)$mse
  
  res.cor=cor(t.res[,1],tyvar,method='spearman',use='pairwise.complete.obs')
  print(paste(best.res$MSE,":",res,':',res.cor))
  return(data.frame(MSE=best.res$MSE,testMSE=res,corVal=res.cor,numFeatures=length(genes),genes=as.character(genelist),
                    numSamples=length(yvar)))
}

#' drugMolLogReg
#' Computes logistic regression based on AUC threshold of 100
#' @param tab
#' @param feature
#' @param aucThresh
#' @export
drugMolLogRegEval<-function(clin.data, 
                        mol.data,
                        mol.feature,
                        test.clin,
                        test.mol,
                        category='Condition',
                        aucThresh=100){
  
  drugs<-unlist(intersect(select(clin.data,cat=category)$cat,
                          select(test.clin,cat=category)$cat))
  
  print(paste('Found',length(drugs),'conditions that overlap between training and testing'))
  
    drug.mol<-clin.data%>%
      dplyr::select(`AML sample`,var=category,AUC)%>%
      group_by(`AML sample`,var)%>%
      subset(var%in%drugs)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(select(mol.data,c(Gene,`AML sample`,!!mol.feature)),
                by='AML sample')%>%
      mutate(sensitive=meanVal<aucThresh)
    
    drug.test<-test.clin%>%
      dplyr::select(Sample,var=category,AUC)%>%
      subset(var%in%drugs)%>%
      group_by(Sample,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      left_join(select(test.mol,c(Gene,Sample,!!mol.feature)),by='Sample')%>%
      mutate(sensitive=meanVal<aucThresh)
    
    reg.res<-lapply(unique(drug.mol$var),function(x){
      data.frame(miniLogREval(subset(drug.mol,var==x),subset(drug.test,var==x),mol.feature),
        compound=x, Molecular=mol.feature)})
  
  return(reg.res)
  
}

#' miniLogReval: Assesses the performance of building a logistic regression model on one test of 
#' data `trainTab` and assessing it on a second `testTab`
#' @param trainTab - a tabular form of data with three columns
#' @param testTab - a tabular form of data with three columns 
#' @param mol.features
#' 
miniLogREval<-function(trainTab,testTab,mol.feature){
#  first build our feature matrix
  library(glmnet)
    set.seed(101010101)
    #empty data frame
    ret.df<-data.frame(MSE=0,testMSE=0,corVal=0,numFeatures=0,genes='',numSamples=0)
    
  tmat=NULL
  mat<-NULL

  try(mat<-buildFeatureMatrix(trainTab,mol.feature))
 
 try(tmat<-buildFeatureMatrix(testTab,mol.feature,'Sample'))

 
  if(is.null(mat)||is.null(tmat)||is.null(dim(mat)))
    return(ret.df)
    
 #remove uninformiatve features
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
      return(ret.df)
   
      
  print(paste("Found",length(zvals),'patients with no',mol.feature,
              'data across',ncol(mat),'features')) 
    #now collect our y output variables
  tmp<-trainTab%>%
     dplyr::select(sensitive,`AML sample`)%>%
     distinct()
  yvar<-tmp$sensitive
  #print(yvar)
  names(yvar)<-tmp$`AML sample`
  yvar<-unlist(yvar[rownames(mat)])
  
  
  #now get the y output for test
  ttmp<-testTab%>%
    dplyr::select(sensitive,Sample)%>%
    distinct()
  tyvar<-ttmp$sensitive
  names(tyvar)<-ttmp$Sample
  tyvar<-unlist(tyvar[rownames(tmat)])
  #use CV to get maximum AUC
  
  ##now reduce test matrix to only those features in original model
  shared<-intersect(colnames(tmat),colnames(mat))
  missing<-setdiff(colnames(mat),colnames(tmat))
  
  print(paste('missing',length(missing),'genes and using',length(shared)))
  
  if(length(shared)==0){
    print("None of the genes are shared?")
    print(paste(colnames(tmat)[1:10],collapse=';'))
    return(ret.df)
  }
  
  tmat<-tmat[,shared]
  if(length(missing)>0){
    newmat<-matrix(nrow=nrow(tmat),ncol=length(missing),data=0)
    colnames(newmat)<-missing
    tmat<-cbind(tmat,newmat)
  }

  #use CV to get minimum MSE
  cv.res<-NULL
  try(cv.res<-cv.glmnet(x=mat,y=yvar,family='binomial',
                        type.measure='mse',nfolds=length(yvar)))
  if(is.null(cv.res))
    return(ret.df)
  
  best.res<-data.frame(lambda=cv.res$lambda,MSE=cv.res$cvm)%>%
    subset(MSE==min(MSE))
  
  #then select how many elements
  full.res<-glmnet(x=mat,y=yvar,family='binomial',type.measure='mse')

  genes=NULL
  try(genes<-names(which(full.res$beta[,which(full.res$lambda==best.res$lambda)]!=0)))
  
  if(is.null(genes))
    return(ret.df)

  genelist<-paste(genes,collapse=';')
  t.res<-predict(full.res,newx=tmat,family='binomial',s=best.res$lambda)

  res=0
  try(res<-assess.glmnet(full.res,newx=tmat,newy=tyvar,s=best.res$lambda)$mse)
  #res.cor=cor(t.res[,1],tyvar,method='spearman',use='pairwise.complete.obs')
  res.cor<-0
  try(res.cor<-cor(t.res[,1],tyvar,method='spearman',use='pairwise.complete.obs'))
  print(paste(best.res$MSE,":",res,':',res.cor))
  
  return(data.frame(MSE=best.res$MSE,testMSE=res,corVal=res.cor,numFeatures=length(genes),genes=as.character(genelist),
                    numSamples=length(yvar)))
}

#'combForest
#'Runs random forest on combination of feature types
#'@param feature.list
#'@param tab
#'@export
#'@return a data frame with 3 values
combForestEval<-function(tab,feature.list=c('proteinLevels','mRNAlevels','geneMutations')){
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
combRegEval<-function(tab,feature.list=c('proteinLevels','mRNALevels','geneMutations')){
  
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
  
  #use CV to get maximum AUC
  cv.res=cv.glmnet(x=comb.mat,y=yvar,type.measure='mse')
  best.res<-data.frame(lambda=cv.res$lambda,MSE=cv.res$cvm)%>%
    subset(MSE==min(MSE))
  
  #then select how many elements
  full.res<-glmnet(x=comb.mat,y=yvar,type.measure='mse')
  genes=names(which(full.res$beta[,which(full.res$lambda==best.res$lambda)]!=0))
  genelist<-paste(genes,collapse=';')
  #print(paste(best.res$MSE,":",genelist))
  return(data.frame(MSE=best.res$MSE,numFeatures=length(genes),genes=genelist,numSamples=length(yvar)))
  
  
}


#' miniForest
#' Runds random forest on the table and molecular feature of interest
#' @import randomForest
#' @export 
#' @return list
miniForestEval<-function(tab,mol.feature,quant=0.995){
  library(randomForest)
  
  #first build our feature matrix
  mat<-buildFeatureMatrix(tab,mol.feature)
  
  #rint(dim(mat))
  if(is.null(dim(mat)))
    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=0))
  
  cm<-apply(mat,1,mean)
  zvals<-which(cm==0.0)
  print(paste("Found",length(zvals),'patients with no',mol.feature,'data across',ncol(mat),'features'))
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


