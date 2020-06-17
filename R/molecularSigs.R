

#' working on this function still
#' The goal is to use a basic elastic net regression to identify how 
#' well each molecular feature predicts outcome as well as how many features
#' are selected
#' @import glmnet
#' @param clin.data is tidied clinical data
#' @param mol.data is tidied table of molecular data
#' @param mol.feature is name of column to select from mol.data
#' @param category can be either 'Condition' or 'family' which is how to group AUC values
#' @export 
drugMolRegression<-function(clin.data,
                            mol.data,
                            mol.feature,
                            category='Condition'){
  
  #' miniReg
  #' Runs lasso regression from tabular data
  #' @param tab with column names `AML sample`,Gene, and whatever the value of 'mol.feature' is.
  #' @return a data.frame with three values/colulmns: MSE, NumGenes, and Genes
  miniReg<-function(tab){
    library(glmnet)
    
    vfn=list(0.0)
    names(vfn)=mol.feature
    vfc<-list(mean)
    names(vfc)=mol.feature
    
    mat<-tab%>%
      dplyr::select(`AML sample`,Gene,!!mol.feature)%>%
      tidyr::pivot_wider(names_from="AML sample",values_from=mol.feature,
                         values_fill=vfn,values_fn = vfc)%>%
      subset(Gene!="")%>%
      tibble::column_to_rownames('Gene')
    
    cm<-apply(mat,2,mean)
    zvals<-which(cm==0)
    if(length(zvals)>0)
      mat<-mat[,-zvals]
    
    if(ncol(mat)<5)
      return(data.frame(MSE=0,numGenes=0,genes=''))
    
    tmp<-tab%>%
      dplyr::select(meanVal,`AML sample`)%>%
      distinct()
    yvar<-tmp$meanVal
    names(yvar)<-tmp$`AML sample`
    yvar<-yvar[colnames(mat)]
    
    #use CV to get maximum AUC
    cv.res=cv.glmnet(x=t(mat),y=yvar,type.measure='mse')
    best.res<-data.frame(lambda=cv.res$lambda,MSE=cv.res$cvm)%>%
      subset(MSE==min(MSE))
    
    #then select how many elements
    full.res<-glmnet(x=t(mat),y=yvar,type.measure='mse')
    genes=names(which(full.res$beta[,which(full.res$lambda==best.res$lambda)]!=0))
    genelist<-paste(genes,collapse=',')
    #print(paste(best.res$MSE,":",genelist))
    return(data.frame(MSE=best.res$MSE,numGenes=length(genes),genes=genelist))
  }
  
  
  drug.mol<-clin.data%>%
    dplyr::select(`AML sample`,var=category,percAUC)%>%
    #subset(clin.data,Metric%in%c("auc","AUC"))%>%
    #  dplyr::select(`AML sample`,var=category,Value)%>%
    group_by(`AML sample`,var)%>%
    summarize(meanVal=mean(percAUC,na.rm=T))%>%
    left_join(select(mol.data,c(Gene,`AML sample`,!!mol.feature)),by='AML sample')
  
  
  reg.res<-drug.mol%>%group_by(var)%>%
    group_modify(~ miniReg(.x),keep=T)%>%
    mutate(Molecular=mol.feature)
  
  return(reg.res)
  
}

#' how do we visualize the correlations of each drug/gene pair?
#' @param cor.res
#' @param cor.thresh 
#' Computes drug by element correlaction and plots them in various ways
#' @import ggplot2
#' @import dplyr
#' @import cowplot
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
  
  famplots<-all.cors%>%split(all.cors$family)%>%purrr::map(do.p,cor.thresh)                                                  
  
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
    group_modify(~ computeFoldChangePvals(.x,control=NA,conditions =c("Sensitive","Resistant")),keep=TRUE)%>%
    rename(Drug='cellLine')
  
  
  table(result%>%group_by(Drug)%>%subset(p_adj<0.1)%>%summarize(sigProts=n()))
  
  result
}


#' what molecules are
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
  dcors<-tdat%>%select(Gene,Mol,Condition,percAUC,family)%>%
    distinct()%>%
    group_by(Gene,Condition)%>%
    mutate(numSamps=n(),drugCor=cor(Mol,percAUC))%>%
    dplyr::select(Gene,Condition,numSamps,drugCor)%>%distinct()%>%
    arrange(desc(drugCor))%>%
    subset(numSamps>10)%>%
    mutate(feature=mol.feature)
  
  return(dcors)
  
}