




#' drugMolRegressionEval
#' This is the primary function that runs continuous regression (enet) on separate train and test sets
#' of molecular features and drug responses. It will iterate through every drug in the table to create indivdiual models of each
#' 
#' @import glmnet
#' @param clin.data is tidied clinical data
#' @param mol.data is tidied table of molecular data
#' @param mol.feature is name of column to select from mol.data, e.g. Gene
#' @param mol.feature.name the name of teh value to select from mol.data
#' @param test.clin tidied clinical data to test
#' @param test.mol tidieid molecular data to test
#' @param category can be either Condition or family
#' @param doEnet run elastic net instead of LASSO
#' @export 
drugMolRegressionEval<-function(clin.data,
                            mol.data,
                            mol.feature,
                            mol.feature.name,
                            test.clin,
                            test.mol,
                            category='Condition',
                               doEnet=FALSE){
  
  ##first check to evaluate drug overlap
   drugs<-unlist(intersect(select(clin.data,cat=category)$cat,
                           select(test.clin,cat=category)$cat))
   
   message(paste('Found',length(drugs),'conditions that overlap between training and testing'))

  mol.feature<-unlist(mol.feature)
  mol.feature.name<-unlist(mol.feature.name)
#   print(mol.feature)
#   print(mol.feature.name)
  
     drug.mol<-clin.data%>%
     dplyr::select(`AML sample`,var=category,AUC)%>%
     group_by(`AML sample`,var)%>%
     subset(var%in%drugs)%>%
     summarize(meanVal=mean(AUC,na.rm=T))%>%
     inner_join(select(mol.data,c(unique(mol.feature),'AML sample',unique(mol.feature.name))),
                by='AML sample')#%>%
     #mutate(sensitive=meanVal<aucThresh)
   
   drug.test<-test.clin%>%
     dplyr::select(Sample,var=category,AUC)%>%
     subset(var%in%drugs)%>%
     group_by(Sample,var)%>%
     summarize(meanVal=mean(AUC,na.rm=T))%>%
     inner_join(select(test.mol,c(unique(mol.feature),'Sample',unique(mol.feature.name))),by='Sample')#%>%
 
     alpha=1.0
    if(doEnet)
      alpha=seq(0.1, 0.9, 0.1)
    #mol.feature=paste(mol.feature,collapse=';')
    
    reg.res<-lapply(unique(drug.mol$var),function(x){
      message(x)
      data.frame(miniRegEval(subset(drug.mol,var==x),subset(drug.test,var==x),
                             mol.feature,mol.feature.name,enet.alpha=alpha),
        compound=x,Molecular=paste(mol.feature.name,collapse=';'))})
    
  return(reg.res)
  
}


#' miniRegEval
#' Runs lasso regression on a single feature from tabular data and evaluates on a second set of data.
#' This is a little more complicated than LOO
#' analysis because we have to check for shared feature names
#' @param trainTab with column names `AML sample`,meanVal,Gene, and whatever the value of `mol.feature` is.
#' @param testTab with column names `Sample`, meanVal, Gene, and whatever the value of `mol.feature` is
#' @param mol.feature list of molecular features to be evaluated
#' @export 
#' @return a data.frame with three values/columns: MSE, numFeatures, and Genes
miniRegEval<-function(trainTab,testTab,mol.feature='Gene',feature.val='mRNALevels', enet.alpha = c(1)){
  library(glmnet)
  set.seed(10101)
  
  ret.df<-data.frame(MSE=0,testMSE=0,corVal=0,numFeatures=0,genes='',numSamples=0)
  
  tmat=NULL
  mat<-NULL


  #mol.feature=list(mol.feature)
  names(mol.feature)<-feature.val
 # print(mol.feature)
  #fnames=names(mol.feature)
  
  if(length(mol.feature)>1){
    try(mat<-do.call('cbind',lapply(feature.val,function(x) buildFeatureMatrix(trainTab,mol.feature=mol.feature[[x]],feature.val=x))))
    try(tmat<-do.call('cbind',lapply(feature.val,function(x) buildFeatureMatrix(testTab,mol.feature=mol.feature[[x]],feature.val=x,'Sample'))))
    
  }else{
    try(mat<-buildFeatureMatrix(trainTab,feature.val=feature.val[[1]],mol.feature=mol.feature[[1]]))
    
    try(tmat<-buildFeatureMatrix(testTab,feature.val=feature.val[[1]],mol.feature=mol.feature[[1]],sampname='Sample'))
  }
#  print(tmat[1:10,1:10])
  
     
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
  

  mol.feature<-paste(feature.val,collapse=';')
  message(paste("Found",length(zvals),'features with no',
              mol.feature,'information, keeping',ncol(mat),'features'))

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
  
  models <- list()
  
  best.res<-NULL
  ## Run glmnet for each alpha, saving the best lambda value every time
  #for loops in R are no good!
  #for (alpha in enet.alpha) {
  models<-lapply(enet.alpha,function(alpha){
    model <- cv.glmnet(x = mat, y = yvar, alpha = alpha, 
                       type.measure = 'mse')
    best <- data.frame(lambda = model$lambda, MSE = model$cvm) %>%
      subset(MSE == min(MSE)) %>%
      mutate(alpha = alpha)
    best.res <<- rbind(best.res, best)
    #models[[as.character(alpha)]] <- model
    model
  })
  names(models)<-enet.alpha
  
  ## Picking optimal (according to MSE) lambda and alpha
  best.res  <- best.res %>%
    subset(MSE == min(MSE))
  alpha = best.res$alpha %>%
    as.character()
  lambda = best.res$lambda
  full.res <- models[[alpha]]
  
  ##now reduce test matrix to only those features in original model
  shared<-intersect(colnames(tmat),colnames(mat))
  missing<-setdiff(colnames(mat),colnames(tmat))
 
  message(paste('missing',length(missing),'genes between train and test and using',length(shared)))
  if(length(shared)==0){
    warning("None of the genes are shared?")
    warning(paste(colnames(tmat)[1:10],collapse=';'))
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
  
  ##these are the genes that minimize the MSE
  genes=NULL
  coefs <- coef(full.res, s = lambda)
  try(genes<-names(coefs[coefs[,1] != 0, ])[-1])

  if(is.null(genes))
    return(ret.df)
  genelist<-paste(genes,collapse=';')  

  ## get the prediction onto the new matrix (to assess correlation)
  t.res<-predict(full.res,newx=tmat,s=lambda)

  #use the assess function to get a new MSE
  res=assess.glmnet(full.res,newx=tmat,newy=tyvar,s=lambda)$mse[[1]]
  
  res.cor=cor(t.res[,1],tyvar,method='spearman',use='pairwise.complete.obs')


  message(paste(best.res$MSE,":",res,':',res.cor))
  return(data.frame(MSE=best.res$MSE,testMSE=res,corVal=res.cor,numFeatures=length(genes),
                    genes=as.character(genelist),
                    numSamples=length(yvar)))
}

#' drugMolLogRegEval
#' Computes logistic regression and evaluates on test dataset for a particular aucThreshold
#' @param clin.data AUC data for training
#' @param mol.data molecular data table
#' @param mol.feature list of molecular features to evaluate
#' @param mol.feature.name list of molecular feature names
#' @param test.clin AUC data for testing
#' @param test.mol molecular data for testing
#' @param category Column in AUC data to use
#' @param aucThresh Threshold to call sample sensitive vs. resistant
#' @return 
#' @export
drugMolLogRegEval<-function(clin.data, 
                        mol.data,
                        mol.feature='Gene',
                        mol.feature.name,
                        test.clin,
                        test.mol,
                        category='Condition',
                        aucThresh=100){
  
  
  #names(mol.feature)<-mol.feature.name
  #message(mol.feature)
 # message(names(mol.feature.name))
 
  mol.feature<-unlist(mol.feature)
  mol.feature.name<-unlist(mol.feature.name)
  #print(mol.feature)
  #print(mol.feature.name)
  
  drugs<-unlist(intersect(select(clin.data,cat=category)$cat,
                          select(test.clin,cat=category)$cat))
  
  message(paste('Found',length(drugs),'conditions that overlap between training and testing'))
  
  drug.mol<-clin.data%>%
      dplyr::select(`AML sample`,var=category,AUC)%>%
      group_by(`AML sample`,var)%>%
      subset(var%in%drugs)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
      inner_join(select(mol.data,c(unique(mol.feature),'AML sample',unique(mol.feature.name))),
                by='AML sample')%>%
      mutate(sensitive=meanVal<aucThresh)
    #print(drug.mol)
    
    drug.test<-test.clin%>%
      dplyr::select(Sample,var=category,AUC)%>%
      subset(var%in%drugs)%>%
      group_by(Sample,var)%>%
      summarize(meanVal=mean(AUC,na.rm=T))%>%
    inner_join(select(test.mol,c(unique(mol.feature),'Sample',unique(mol.feature.name))),by='Sample')%>%
      mutate(sensitive=meanVal<aucThresh)
  #  print(drug.test)
    
    reg.res<-lapply(unique(drug.mol$var),function(x){
      message(x)
      data.frame(miniLogREval(subset(drug.mol,var==x),
                              subset(drug.test,var==x),mol.feature,mol.feature.name),
        compound=x, Molecular=paste(mol.feature.name,collapse=';'))})
  
  return(reg.res)
  
}

#' miniLogReval: Assesses the performance of building a logistic regression model on one test of 
#' data `trainTab` and assessing it on a second `testTab`
#' @param trainTab - a tabular form of data with three columns
#' @param testTab - a tabular form of data with three columns 
#' @param mol.feature is a list of column features, with the names of the column features representing the measurements. 
#' e.g. mRNALevels='Gene', or Phosphosite='site'
#' 
miniLogREval<-function(trainTab,testTab,mol.feature='Gene',feature.val='mRNALevels'){
#  first build our feature matrix
  library(glmnet)
    set.seed(101010101)
    #empty data frame
    ret.df<-data.frame(MSE=0,testMSE=0,corVal=0,numFeatures=0,genes='',numSamples=0)
    
  tmat=NULL
  mat<-NULL
  feature.val<-unlist(feature.val)
  
  names(mol.feature)<-feature.val
  
  ##make this into a list with names 
  #mol.feature=list(mol.feature)
  #names(mol.feature)<-feature.val
  #print(mol.feature)
  
  #fnames=names(mol.feature)
 # print(feature.val)
#  print(mol.feature)
  if(length(mol.feature)>1){
    try(mat<-do.call('cbind',lapply(feature.val,function(x) buildFeatureMatrix(trainTab,mol.feature=mol.feature[[x]],feature.val=x))))
    try(tmat<-do.call('cbind',lapply(feature.val,function(x) buildFeatureMatrix(testTab,mol.feature=mol.feature[[x]],feature.val=x,sampname='Sample'))))
    
  }else{
    try(mat<-buildFeatureMatrix(trainTab,feature.val=feature.val[[1]],mol.feature=mol.feature[[1]]))
    try(tmat<-buildFeatureMatrix(testTab,feature.val=feature.val[[1]],mol.feature=mol.feature[[1]],sampname='Sample'))
  }
 # print(tmat[1:10,1:10])
 
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
   

  feature.val<-paste(feature.val,collapse=';')

  message(paste("Found",length(zvals),'features with no',feature.val,
              'information across',ncol(mat),'features')) 
    #now collect our y output variables
  tmp<-trainTab%>%
     dplyr::select(sensitive,`AML sample`)%>%
     distinct()
  yvar<-tmp$sensitive
  
  names(yvar)<-tmp$`AML sample`
  yvar<-unlist(yvar[rownames(mat)])
  #print(yvar)
  
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
  
  message(paste('missing',length(missing),'genes and using',length(shared)))
  
  if(length(shared)==0){
    warning("None of the genes are shared?")
    warning(paste(colnames(tmat)[1:10],collapse=';'))
    return(ret.df)
  }
  
  tmat<-tmat[,shared]
  if(length(missing)>0){
    newmat<-matrix(nrow=nrow(tmat),ncol=length(missing),data=0)
    colnames(newmat)<-missing
    tmat<-cbind(tmat,newmat)
  }

  #print(mat[1:10,1:10])
  #print(yvar)
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
  message(paste(best.res$MSE,":",res,':',res.cor))
  
  return(data.frame(MSE=best.res$MSE,testMSE=res,corVal=res.cor,numFeatures=length(genes),
                    genes=as.character(genelist),
                    numSamples=length(yvar)))
}

#' #'combForest
#' #'Runs random forest on combination of feature types
#' #'@param feature.list
#' #'@param tab
#' #'@export
#' #'@return a data frame with 3 values
#' combForestEval<-function(trainTab,testTab,
#'                          feature.list=c('proteinLevels','mRNAlevels','geneMutations')){
#'   
#'   tr.comb.mat<-NULL
#'   te.comb.mat<-NULL
#'   try(tr.comb.mat<-do.call('cbind',lapply(feature.list,function(x) buildFeatureMatrix(trainTab,x))))
#'   try(te.comb.mat<-do.call('cbind',lapply(feature.list,function(x) buildFeatureMatrix(trainTab,x))))
#'   
#'   #f(ncol(tr.comb.mat)<5 || nrow(tr.comb.mat)<5)
#'   if(is.null(tr.comb.mat)||is.null(te.comb.mat)) 
#'    return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=0))
#'   
#'   #now collect our y output variable
#'   tmp<-tab%>%
#'     dplyr::select(meanVal,`AML sample`)%>%
#'     distinct()
#'   yvar<-tmp$meanVal
#'   names(yvar)<-tmp$`AML sample`
#'   yvar<-unlist(yvar[rownames(comb.mat)])
#'                
#'   rf<-randomForest(comb.mat,yvar)
#' 
#'   return(data.frame(MSE=min(rf$mse),
#'                     numFeatures=length(which(rf$importance!=0)),
#'                     genes=paste(names(rf$importance)[which(rf$importance!=0)],collapse=';',),
#'                     numSamples=length(yvar)))
#'   
#' }
#' 
#' #' combReg
#' #' Runs lasso regression on a combination of feature types
#' #' @param tab
#' #' @export
#' #' @param feature.list
#' #' @return a data frame with three values/columns
#' combRegEval<-function(trainTab,testTab,feature.list=c('proteinLevels','mRNALevels','geneMutations')){
#'   
#'   
#'   tr.comb.mat<-NULL
#'   te.comb.mat<-NULL
#'   try(tr.comb.mat<-do.call('cbind',lapply(feature.list,function(x) buildFeatureMatrix(trainTab,x))))
#'   try(te.comb.mat<-do.call('cbind',lapply(feature.list,function(x) buildFeatureMatrix(testTab,x))))
#'   
#'   #f(ncol(tr.comb.mat)<5 || nrow(tr.comb.mat)<5)
#'   if(is.null(tr.comb.mat)||is.null(te.comb.mat)) 
#'     return(data.frame(MSE=0,numFeatures=0,genes='',numSamples=0))
#'   
#'   
#'   
#'   #now collect our y output variable
#'   tmp<-tab%>%
#'     dplyr::select(meanVal,`AML sample`)%>%
#'     distinct()
#'   yvar<-tmp$meanVal
#'   names(yvar)<-tmp$`AML sample`
#'   yvar<-unlist(yvar[rownames(comb.mat)])
#'   
#'   #use CV to get maximum AUC
#'   cv.res=cv.glmnet(x=comb.mat,y=yvar,type.measure='mse')
#'   best.res<-data.frame(lambda=cv.res$lambda,MSE=cv.res$cvm)%>%
#'     subset(MSE==min(MSE))
#'   
#'   #then select how many elements
#'   full.res<-glmnet(x=comb.mat,y=yvar,type.measure='mse')
#'   genes=names(which(full.res$beta[,which(full.res$lambda==best.res$lambda)]!=0))
#'   genelist<-paste(genes,collapse=';')
#'   #print(paste(best.res$MSE,":",genelist))
#'   return(data.frame(MSE=best.res$MSE,numFeatures=length(genes),genes=genelist,numSamples=length(yvar)))
#'   
#'   
#' }


#' miniForestEval
#' Runds random forest on the table and molecular feature of interest
#' @title miniForestEval
#' @import randomForest
#' @param tab Data table of merged info
#' @param mol.feature
#' @param quant
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
  message(paste("Found",length(zvals),'patients with no',mol.feature,'data across',ncol(mat),'features'))
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



# miniRegMod<-function(trainTab,mol.feature){
#     #first build our feature matrix
#   mat<-buildFeatureMatrix(trainTab,mol.feature)
#   #print(mat)
#   if(is.null(dim(mat)))
#     return(ret.df)
#   
#   cm<-apply(mat,1,mean)
#   vm<-apply(mat,1,var)
#   zvals<-union(which(cm==0),which(vm==0))
#   if(length(zvals)>0)
#     mat<-mat[-zvals,]
#   
#   message(paste("Found",length(zvals),'patients with no',
#               mol.feature,'data across',ncol(mat),'features'))
#   
#   zcols<-apply(mat,2,var)
#   zvals<-which(zcols==0)
#   # print(zvals)
#   if(length(zvals)>0)
#     mat<-mat[,-zvals]
#   
#   if(ncol(mat)<5 || nrow(mat)<5)
#     return(ret.df)
#   
#   #now collect our y output variable for training
#   tmp<-trainTab%>%
#     dplyr::select(meanVal,`AML sample`)%>%
#     distinct()
#   yvar<-tmp$meanVal
#   names(yvar)<-tmp$`AML sample`
#   yvar<-unlist(yvar[rownames(mat)])
#   cv.res=cv.glmnet(x=mat,y=yvar,type.measure='mse')
#   
#   best.res<-data.frame(lambda=cv.res$lambda,MSE=cv.res$cvm)%>%
#     subset(MSE==min(MSE))
#   
#   #then select how many elements
#   full.res<-glmnet(x=mat,y=yvar,type.measure='mse')
#   
#   genes=NULL
#   try(genes<-names(which(full.res$beta[,which(full.res$lambda==best.res$lambda)]!=0)))
#   
# }



#' The goal is to pre filter features using limma, with the hope that
#' the smaller set of features consist of features which are more likely
#' to produce good models by glmnet.
#' @importFrom data.table rbindlist
#' @import limma
#' @param mat matrix samples as rows, features as columns
#' @param tmat see mat
#' @param yvar (auc) response. This is what we eventually want to predict, but here we use it to define categorical variables and train a model.
#' @param tyvar see yvar
#' @param p.val cutoff when filtering
#' @export 
findCandidates <- function(mat, tmat, yvar, tyvar, p.val = 0.001){
  
  sample.names <- c(rownames(mat), rownames(tmat))
  xx <- rbindlist(list(as.data.frame(mat), as.data.frame(tmat)), fill = TRUE, use.names = TRUE) %>%
    as.data.frame()
  rownames(xx) <- sample.names
  
  auc <- c(yvar, tyvar)
  
  ## Making sure names match
  if (!all(names(auc) == sample.names)){
    stop("Error, sample name mismatch")
  }
  
  ## making AUC factor from AUC.
  meta <- data.frame(Sample = sample.names) %>%
    mutate(auc.group = cut(auc, breaks = c(0,100,200,300)))
  auc.mod <- model.matrix(~0 + auc.group, meta)
  coefs <- colnames(auc.mod)
  xx <- t(as.matrix(xx))
  fit <- lmFit(xx, auc.mod)
  fit.smooth <- eBayes(fit)
  sig <- topTable(fit.smooth, number = nrow(xx), sort.by = "none", 
                  coef = coefs) %>%
    filter(P.Value < p.val)
  candidates <- rownames(sig)
  
  return(candidates)
}
