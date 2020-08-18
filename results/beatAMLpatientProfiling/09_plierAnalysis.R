##run PLIER analysis on data

library(amlresistancenetworks)
library(PLIER)
library(dplyr)
library(tidyr)
#rc2path = "https://ndownloader.figshare.com/files/10881866"

##download these behemoth separately
#https://figshare.com/articles/recount_rpkm_RData/5716033/4

#then we want to looad this file
plier.results <- readRDS("recount_PLIER_model.RDS")

#'
#'Collado-Torres L, Nellore A, Kammers K, Ellis SE, Taub MA, Hansen KD, Jaffe AE, Langmead B and Leek JT (2017). "Reproducible RNA-seq analysis using recount2." Nature Biotechnology. doi: 10.1038/nbt.3838
#'And the PLIER preprint, if you use the PLIER model:
#'  Mao W, Harmann B, Sealfon SC, Zaslavsky E, and Chikina M (2017). "Pathway-Level Information ExtractoR (PLIER) for gene expression data." bioRxiv. doi: 10.1101/116061
#'  
#
if(!exists('dataLoaded')){
  source('beatAMLdata.R')
  dataLoaded=TRUE
}

exprs.mat<-pat.data%>%
  select(`AML sample`,'Gene','mRNALevels')%>%
  pivot_wider(values_from=mRNALevels, names_from=`AML sample`,
              values_fn=list(mRNALevels=mean),
              values_fill=list(mRNALevels=0.01))%>% #no zeroes allowed!
  tibble::column_to_rownames('Gene')%>%as.matrix()

source("../../../multi-plier/util/plier_util.R")
pat.recount.b <- GetNewDataB(exprs.mat = exprs.mat,
                             plier.model = plier.results)


###now that all the data are loaded we can compute the correlation, regression, and random forest values

lv.df<-pat.recount.b%>%
    as.data.frame()%>%
    tibble::rownames_to_column("Latent Variable")%>%
    pivot_longer(-`Latent Variable`,names_to='AML sample',values_to="Loading")

#doesnt work
#synTableStore(as.data.frame(lv.df),tabname='BeatAML Patient Latent Variables',parentId='syn22128879')

lv.df<-rename(lv.df,Gene='Latent Variable')
reg.results<-drugMolRegression(auc.dat,lv.df,'Loading')
rf.results<-drugMolRandomForest(auc.dat,lv.df,'Loading')
cor.results<-computeAUCCorVals(auc.dat,lv.df,'Loading')


