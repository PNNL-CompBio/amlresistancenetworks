test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


##create gene by sample by protein table


test_that('regression', {
  library(dplyr)
  dat<-readRDS('testData.rda')
  
  tr.samps<-unique(dat$mol.data$Sample)[1:40]
  te.samps<-setdiff(dat$mol.data$Sample,tr.samps)
  
  auc.dat<-dat$auc.data
  drugs<-unique(auc.dat$Condition)[1:5]
  
  
  tr.dat<-subset(dat$mol.data,Sample%in%tr.samps)
  te.dat<-subset(dat$mol.data,Sample%in%te.samps)
  
  auc.dat<-auc.dat%>%
    rename(AUC='Value')%>%
    subset(Condition%in%drugs)
  ##now train model on AML and eval on depmap data
  reg.preds<-purrr::map_df(list(mRNA='transcriptCounts',
                                protein='LogFoldChange',
                                #  gene='geneMutations',
                                binGene='numMuts'),
                           ~ drugMolRegressionEval(rename(auc.dat,`AML sample`='Sample'),
                                                   rename(tr.dat,`AML sample`='Sample'),
                                                   .x,
                                                   auc.dat,
                                                   te.dat,
                                                   category='Condition'))
  expect_equal(dim(reg.preds),c(15,8))
  
})

test_that('logistic regression',{
  
  library(dplyr)
  dat<-readRDS('testData.rda')
  
  tr.samps<-unique(dat$mol.data$Sample)[1:40]
  te.samps<-setdiff(dat$mol.data$Sample,tr.samps)
  
  auc.dat<-dat$auc.data
  drugs<-unique(auc.dat$Condition)[1:5]
  
  
  tr.dat<-subset(dat$mol.data,Sample%in%tr.samps)
  te.dat<-subset(dat$mol.data,Sample%in%te.samps)
  
  auc.dat<-auc.dat%>%
    rename(AUC='Value')%>%
    subset(Condition%in%drugs)
  ##now train model on AML and eval on depmap data
  reg.preds<-purrr::map_df(list(mRNA='transcriptCounts',
                                protein='LogFoldChange',
                                #  gene='geneMutations',
                                binGene='numMuts'),
                           ~ drugMolLogRegEval(rename(auc.dat,`AML sample`='Sample'),
                                                   rename(tr.dat,`AML sample`='Sample'),
                                                   .x,
                                                   auc.dat,
                                                   te.dat,
                                                   category='Condition'))
  expect_equal(dim(reg.preds),c(15,8))
  
})

test_that('elastic net regression',{
  library(dplyr)
  dat<-readRDS('testData.rda')
  
  tr.samps<-unique(dat$mol.data$Sample)[1:40]
  te.samps<-setdiff(dat$mol.data$Sample,tr.samps)
  
  auc.dat<-dat$auc.data
  drugs<-unique(auc.dat$Condition)[1:5]
  
  
  tr.dat<-subset(dat$mol.data,Sample%in%tr.samps)
  te.dat<-subset(dat$mol.data,Sample%in%te.samps)
  
  auc.dat<-auc.dat%>%
    rename(AUC='Value')%>%
    subset(Condition%in%drugs)
  ##now train model on AML and eval on depmap data
  reg.preds<-purrr::map_df(list(mRNA='transcriptCounts',
                                protein='LogFoldChange',
                                #  gene='geneMutations',
                                binGene='numMuts'),
                           ~ drugMolRegressionEval(rename(auc.dat,`AML sample`='Sample'),
                                                   rename(tr.dat,`AML sample`='Sample'),
                                                   .x,
                                                   auc.dat,
                                                   te.dat,
                                                   category='Condition',
                                                   doEnet=TRUE))
  expect_equal(dim(reg.preds),c(15,8))
  
})