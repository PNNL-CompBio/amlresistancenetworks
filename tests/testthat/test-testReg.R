test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


##create gene by sample by protein table

test_that('multi-data regression', {
  library(dplyr)
  # dat<-load('testData')
  
  tr.samps<-unique(testData$mol.data$Sample)[1:40]
  te.samps<-setdiff(testData$mol.data$Sample,tr.samps)
  
  auc.dat<-testData$auc.data
  drugs<-unique(auc.dat$Condition)[1:5]
  
  
  tr.dat<-subset(testData$mol.data,Sample%in%tr.samps)
  te.dat<-subset(testData$mol.data,Sample%in%te.samps)
  
  auc.dat<-auc.dat%>%
    rename(AUC='Value')%>%
    subset(Condition%in%drugs)
  ##now train model on AML and eval on depmap data
  #eval.list<-list(transcriptCounts='Gene',LogFoldChange='Gene',numMuts='Gene')
  
  eval.list<-list(list(mf=c('Gene','Gene','Gene'),fn=c('transcriptCounts','LogFoldChange','numMuts')))
  
  ##now train model on AML and eval on depmap data
  reg.preds<-purrr::map_df(eval.list,
                           ~ amlresistancenetworks::drugMolRegressionEval(rename(auc.dat,`AML sample`='Sample'),
                                                                          rename(tr.dat,`AML sample`='Sample'),
                                                                          mol.feature=.x[1],
                                                                          mol.feature.name=.x[2],
                                                                          auc.dat,
                                                                          te.dat,
                                                                          category='Condition'))
  
  expect_equal(dim(reg.preds),c(5,8))
  
})

test_that('regression', {
  library(dplyr)
 # dat<-load('testData')
  
  tr.samps<-unique(testData$mol.data$Sample)[1:40]
  te.samps<-setdiff(testData$mol.data$Sample,tr.samps)
  
  auc.dat<-testData$auc.data
  drugs<-unique(auc.dat$Condition)[1:5]
  
  
  tr.dat<-subset(testData$mol.data,Sample%in%tr.samps)
  te.dat<-subset(testData$mol.data,Sample%in%te.samps)
  
  auc.dat<-auc.dat%>%
    rename(AUC='Value')%>%
    subset(Condition%in%drugs)
  ##now train model on AML and eval on depmap data
  
  eval.list<-list(c('Gene','transcriptCounts'),c('Gene','LogFoldChange'),c('Gene','numMuts'))
  
  ##now train model on AML and eval on depmap data
  reg.preds<-purrr::map_df(eval.list,
                           ~ amlresistancenetworks::drugMolRegressionEval(rename(auc.dat,`AML sample`='Sample'),
                                                                      rename(tr.dat,`AML sample`='Sample'),
                                                                      mol.feature=.x[1],
                                                                      mol.feature.name=.x[2],
                                                                      auc.dat,
                                                                      te.dat,
                                                                      category='Condition'))
  
  expect_equal(dim(reg.preds),c(15,8))
  
})

test_that('logistic regression',{
  
  library(dplyr)
 # dat<-load('testData')
  
  tr.samps<-unique(testData$mol.data$Sample)[1:40]
  te.samps<-setdiff(testData$mol.data$Sample,tr.samps)
  
  auc.dat<-testData$auc.data
  drugs<-unique(auc.dat$Condition)[1:5]
  
  
  tr.dat<-subset(testData$mol.data,Sample%in%tr.samps)
  te.dat<-subset(testData$mol.data,Sample%in%te.samps)
  
  auc.dat<-auc.dat%>%
    rename(AUC='Value')%>%
    subset(Condition%in%drugs)
  
  eval.list<-list(c('Gene','transcriptCounts'),c('Gene','LogFoldChange'),c('Gene','numMuts'))
  
  ##now train model on AML and eval on depmap data
  reg.preds<-purrr::map_df(eval.list,
                           ~ amlresistancenetworks::drugMolLogRegEval(rename(auc.dat,`AML sample`='Sample'),
                                                   rename(tr.dat,`AML sample`='Sample'),
                                                   mol.feature=.x[1],
                                                   mol.feature.name=.x[2],
                                                   auc.dat,
                                                   te.dat,
                                                   category='Condition',aucThresh=0.5))
  expect_equal(dim(reg.preds),c(15,8))
  
})

test_that('elastic net regression',{
  library(dplyr)
 # dat<-load('testData')
  
  tr.samps<-unique(testData$mol.data$Sample)[1:40]
  te.samps<-setdiff(testData$mol.data$Sample,tr.samps)
  
  auc.dat<-testData$auc.data
  drugs<-unique(auc.dat$Condition)[1:5]
  
  
  tr.dat<-subset(testData$mol.data,Sample%in%tr.samps)
  te.dat<-subset(testData$mol.data,Sample%in%te.samps)
  
  auc.dat<-auc.dat%>%
    rename(AUC='Value')%>%
    subset(Condition%in%drugs)
  ##now train model on AML and eval on depmap data
  eval.list<-list(c('Gene','transcriptCounts'),c('Gene','LogFoldChange'),c('Gene','numMuts'))
  
  ##now train model on AML and eval on depmap data
  reg.preds<-purrr::map_df(eval.list,
                           ~ amlresistancenetworks::drugMolRegressionEval(rename(auc.dat,`AML sample`='Sample'),
                                                                          rename(tr.dat,`AML sample`='Sample'),
                                                                          mol.feature=.x[1],
                                                                          mol.feature.name=.x[2],
                                                                          auc.dat,
                                                                          te.dat,
                                                                          category='Condition',doEnet=TRUE))
  expect_equal(dim(reg.preds),c(15,8))
  
})