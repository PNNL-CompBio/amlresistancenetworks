##proteomic stats





#'
#'plot single protein
#'@import dplyr
#'@import ggplot2
#'@export
#'@param gene
#'@param df tidied data fraome
#'@param prefix
plotSingleProtein<-function(gene,df,prefix=''){
    subset(df,Gene==gene)%>%
    ggplot2::ggplot(aes(x=treatment,y=value))+
      ggplot2::geom_jitter(aes(col=ligand),position=position_dodge(0.8)) +
      ggplot2::geom_boxplot(aes(fill=ligand,alpha=0.5),outlier.shape=NA,position=position_dodge(0.8))+
      ggplot2::facet_grid(~cellLine)+
      ggplot2::ggtitle(paste(gene,'Expression',prefix))
    ggsave(paste0(prefix,gene,'Expression.png'))
  
}



#'
#'limmaTwoFactorDEAnalysis
#'uses Osama's code to compute de from limma
#'@author Osama
#'@import limma
#'@param data matrix
#'@param group1 ids
#'@param group2 ids
limmaTwoFactorDEAnalysis <- function(dat, sampleIDs.group1, sampleIDs.group2) {
  # Conduct DE expression analysis using limma from the expression matrix dat (group2 vs group1, group1 is reference)
  #
  # Args:
  #   dat: Expression data matrix, rows are genes, columns are samples
  #   sampleIDs.group1: Vector with ids of samples in reference group (eg. normal samples)
  #   sampleIDs.group2: Vector with ids of samples in interest group (eg. tumor samples) 
  #
  # Returns:
  #   limma Differential Expression results.
  #
  #http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/R_guide.html
  #http://genomicsclass.github.io/book/pages/using_limma.html
  #https://wiki.bits.vib.be/index.php/Tutorial:_Testing_for_differential_expression_I
  library(limma)
  fac <- factor(rep(c(2,1), c(length(sampleIDs.group2), length(sampleIDs.group1))))
  design <- model.matrix(~fac)
  fit <- lmFit(dat[,c(sampleIDs.group2, sampleIDs.group1)], design)
  fit <- eBayes(fit)
  print(topTable(fit, coef=2))
  res <- topTable(fit, coef=2, number=Inf, sort.by="none")
  res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  return(res)
}


#' Reads in data frame and computes fold change
#' 
#' @export
#' @param data.frame 
#' @param treatmentCond
#' @param ligands
#' @return tidied data frame with p-values and fold change
#' @import dplyr
#' @import stringr
#' @import tidyr
#' @import purrr
#' 
computeFoldChangePvals<-function(g.data,
                                 control=c('None'),
                                 conditions=c("FLT3","FGF2")){
  
  print(unique(g.data$cellLine))
  
  data<-g.data%>%
    dplyr::select(cellLine,Gene,value,treatment,Sample)%>%
    subset(!is.na(value))%>%
    subset(!is.na(Gene))
  
  
  #remove those data that have only 1 replicate  
  reps<-data%>%
    group_by(cellLine,Gene,treatment)%>%
    dplyr::mutate(reps=length(value))%>%
    subset(reps>1)%>%
    ungroup()        
  
  
  lig.datasets<-purrr::map_df(conditions,function(lig){
    print(lig)
    g1<-subset(reps,treatment%in%control)%>%select(Sample)%>%distinct()
    g2<-subset(reps,treatment==lig)%>%select(Sample)%>%distinct()
    matp=reps%>%
      subset(treatment%in%c(control,lig))%>% #get only those valuese that are 'nOne' or the name of the ligand
      dplyr::select(Gene,Sample,value)%>%distinct()%>%
      tidyr::pivot_wider(values_from=value,names_from=Sample,values_fn=list(value=mean))%>%
      tibble::column_to_rownames("Gene")
    limmaTwoFactorDEAnalysis(matp,unlist(g1),unlist(g2))%>%
      dplyr::mutate(Condition=lig,Control=paste(control,collapse=','))%>%
      dplyr::select(Gene=featureID,Condition,Control,condition_to_control=logFC,p_val=P.Value,p_adj=adj.P.Val)
    })
 
  return(lig.datasets) 
}

##drug sensitivity statistics

#' Get correlated proteins
#' Computes protein correlations with dose response metric, either AUC or IC50
#' returns list of drugs, proteins and correlation values
#' @param data
#' @param metric
#' @return data.frame
#' @require ggplot2
#' @require dplyr
#' @export
getCorrelatedProteins<-function(data,metric="AUC"){
  cor.stats<-subset(data,Metric==metric)%>%
    group_by(Condition,Gene)%>%
    dplyr::select(`AML sample`,LogFoldChange,Value)%>%distinct()%>%
    mutate(correlation=cor(LogFoldChange,Value))%>%
    dplyr::select(Condition,correlation)%>%
    distinct()
  ggplot2::ggplot(cor.stats)+
    ggplot2::geom_histogram(mapping=ggplot2::aes(x=correlation,fill=Condition),
                            position='dodge')+
    ggplot2::ggtitle(paste('Correlation between protein change and drug',metric))
  ggplot2::ggsave(paste0('drugProtein',metric,'correlation.png'))
  
  univ=unique(cor.stats$Gene)
  sapply(unique(cor.stats$Condition),function(cond){
    subset(cor.stats,Condition==cond)%>%
      dplyr::select(Gene,value='correlation')%>%
      amlresistancenetworks::computeGSEA(prot.univ = univ,prefix=paste(cond,'correlation with',metric))
  })
  
  cor.stats
  
}


#' getCoClusteredProteins
#' Designed to evaluate time-course data
#' Gets groups of proteins that are clustered across time points in a condition of interest
#' @param timeCourseData
#' @return dataFrame
#' @require dplyr
#' @export 
getCoClusteredProteins<-function(timeCourseData){
  
  
}
