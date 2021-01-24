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
#'@import BiocManager
#'@export
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
  if(!require('limma')){
    BiocManager::install('limma')
    
    library(limma)
}
  fac <- factor(rep(c(2,1), c(length(sampleIDs.group2), length(sampleIDs.group1))))
  design <- model.matrix(~fac)
  fit <- lmFit(dat[,c(sampleIDs.group2, sampleIDs.group1)], design)
  fit <- eBayes(fit)
  print(topTable(fit, coef=2))
  res <- topTable(fit, coef=2, number=Inf, sort.by="none")
  res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  return(arrange(res,P.Value))
}

#' compute manual lfc
#' @export
#' @param data matrix
#' @param group1 ids
#' @param group2 ids
#' @return data frame with fold change, p-value
manualDEAnalysis<-function(data,group1,group2){
  
  if(length(group1)==1)
    control=data[,group1]
  else
    control=rowMeans(data[,group1])
  if(length(group2)==1)
    condition=data[,group2]
  else
    condition=rowMeans(data[,group2])
  res = data.frame(featureID=rownames(data),logFC=condition - control)
  res$absFc = abs(res$logFC)
  res<-subset(res,!is.na(logFC))
  res$P.Value = 1.0-(rank(res$absFc)/nrow(res))
  res$adj.P.Val = p.adjust(res$P.Value)
  return(res)
}


#' Reads in data frame and computes fold change
#' 
#' @export
#' @param data.frame 
#' @param control
#' @param doManual binary set to calculate fc manual
#' @param conditions
#' @return tidied data frame with p-values and fold change
#' @import dplyr
#' @import stringr
#' @import tidyr
#' @import purrr
#' 
computeFoldChangePvals<-function(g.data,
                                 control=c('None'),
                                 conditions=c("FLT3","FGF2"),
                                 doManual=FALSE){
  
  print(unique(g.data$cellLine))
  
  data<-g.data%>%
    dplyr::select(cellLine,Gene,value,treatment,Sample)%>%
    subset(!is.na(value))%>%
    subset(!is.na(Gene))
  
  #remove those data that have only 1 replicate  
  if(doManual)
    reps=data
  else
    reps<-data%>%
      group_by(cellLine,Gene,treatment)%>%
      dplyr::mutate(reps=length(value))%>%
      subset(reps>1)%>%
      ungroup()        
  
  if(doManual)
    myfun=manualDEAnalysis
  else
    myfun=limmaTwoFactorDEAnalysis
  
  lig.datasets<-purrr::map_df(conditions,function(lig){
    print(lig)
    g1<-subset(reps,treatment%in%control)%>%dplyr::select(Sample)%>%distinct()
    g2<-subset(reps,treatment==lig)%>%dplyr::select(Sample)%>%distinct()
    matp=reps%>%
      subset(treatment%in%c(control,lig))%>% #get only those valuese that are 'nOne' or the name of the ligand
      dplyr::select(Gene,Sample,value)%>%distinct()%>%
      tidyr::pivot_wider(values_from=value,names_from=Sample,values_fn=list(value=mean))%>%
      tibble::column_to_rownames("Gene")
     myfun(matp,unlist(g1),unlist(g2))%>%
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
#' @require amap
#' @export 
getCoClusteredProteins<-function(prot.data){
  
  library(amap)
  # print(cellLine,treatment)
  cellLine=unique(prot.data$cellLine)
  treatment=unique(prot.data$treatment)
  # print(head(prot.data))
  mat<-prot.data%>%
    dplyr::select(-c(cellLine,treatment))%>%
    mutate(hours=if_else(timePoint=='16 hr',16,if_else(timePoint=='3 hr',3,.5)))%>%
    mutate(val=as.numeric(as.character(LogFoldChange)))%>%
    dplyr::select(-c(LogFoldChange,timePoint))%>%
    tidyr::pivot_wider(values_from=val,names_from='hours',values_fn=list(val=mean))
  
  mv<-which(apply(mat,1,function(x) any(is.na(x))))
  if(length(mv)>0)
    mat<-mat[-mv,]
  # wss <-purrr::map_dbl(3:20, ~{Kmeans(dplyr::select(mat,-Gene),.,nstart=5,iter.max=30,method='pearson')$tot.withinss})
  # ggplot(data.frame(k=3:20,wss=wss))+geom_point(aes(x=k,y=wss))
  # ggsave(paste0(treatment,'_treated_',cellLine,'_cellsClusters.png'))
  #let's choose k =10
  clusters = Kmeans(dplyr::select(mat,-Gene),centers=10,nstart=20,iter.max=30,method='pearson')
  centers <- tibble::rownames_to_column(as.data.frame(clusters$centers), "cluster")%>%
    tidyr::pivot_longer(2:4,names_to='timepoint',values_to='logRatio')
  
  # ggplot(centers)+geom_path(aes(x=as.numeric(timepoint),y=logRatio,col=cluster,group=cluster))
  mat$cluster = clusters$cluster 
  
  res<-mat%>%tidyr::pivot_longer(2:4,values_to='logRatio',names_to='time')
  
  ggplot(res,aes(x=as.numeric(time),y=logRatio,col=as.factor(cluster),group=Gene))+
    geom_path()+facet_grid(.~cluster)+
    ggtitle(paste(treatment,'treated',cellLine,'cells'))
  ggsave(paste0(treatment,'_treated_',cellLine,'cells.png'))
  
  return(res)
}


#'get ordered phospho sites
#'@param phosData
#'@return data table of log fold change
#'@import dplyr
computePhosphoChanges<-function(phosData){
  dplyr::select(phosData,Gene,site,Peptide)%>%distinct()%>%
    dplyr::mutate(residue=stringr::str_replace(site,paste0(Gene,'-'),''))%>%
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([STY])", ";\\1"))%>%
    dplyr::mutate(residue=stringr::str_replace(residue,"^;", ""))%>%
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([sty])", ""))
  
}

#'mapPhosphoToKInase
#'reads in phospho data, looads of KSDB and then computes log fold change for each kinase
#'@param pat.phos
#'@return data table of merged data
#'@import dplyr
#'@export
mapPhosphoToKinase<-function(pat.phos){

  KSDB <- read.csv(system.file('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',
                               package='amlresistancenetworks'),stringsAsFactors = FALSE)

  clinvars<-setdiff(names(pat.phos),c('Protein','Gene','site','Peptide','LogRatio'))
  pat.phos$Gene<-unlist(pat.phos$Gene)
  gene.phos<-computePhosphoChanges(pat.phos)
  phos.with.subs<-gene.phos%>%
    inner_join(rename(KSDB,Gene='SUB_GENE',residue='SUB_MOD_RSD'),by=c('Gene','residue'))%>%
    inner_join(pat.phos,by=c('site','Gene'))

  pat.kin.scores<-phos.with.subs%>%
    dplyr::select('Sample','site',Gene,LogFoldChange,GENE,networkin_score)%>%distinct()%>%
    group_by(Sample,GENE)%>%
    summarize(meanLFC=mean(LogFoldChange),meanNKINscore=mean(networkin_score),numSubstr=n_distinct(Gene))%>%
    rename(Kinase='GENE')
  return(pat.kin.scores)
}

#'plotKinDat
#'@param kindat - output of `mapPhospoToKinase`
#'@param phosData - tidied phospho data
#'@param prefix - string for file
#'@param idcol- name of sample identifier
#'@param vars - sample-specific 
#'@export
#'@import pheatmap
plotKinDat<-function(kindat,phosData=phosData,prefix='all',
                     idcol='Sample',vars=c('Sample','cellLine','Ligand')){
  library(pheatmap)
  ##create matrix of kinase scores
  mat <-kindat%>%
    ungroup()%>%
    tidyr::pivot_wider(-c(meanNKINscore,numSubstr),
                                                values_from=meanLFC,
                                                names_from=Sample,
                                                values_fn=list(meanLFC=mean))%>%
    tibble::column_to_rownames('Kinase')
  kinAts<-kindat%>%
    ungroup()%>%
    dplyr::select(Kinase,numSubstr)%>%
    distinct()%>%
    group_by(Kinase)%>%
    summarize(substrates=mean(numSubstr))%>%
    tibble::remove_rownames()%>%
  tibble::column_to_rownames('Kinase')
  
  sampAts<-phosData%>%
    dplyr::select(vars)%>%
    distinct()%>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames(idcol)
  #sampAts$TimePoint=as.factor(sampAts$TimePoint)
  vars=names(sort(apply(mat,1,var),decreasing=T)[1:200])
 pheatmap(mat[vars,],cellwidth = 8,cellheight=8,clustering_distance_cols = 'correlation',
          clustering_distance_rows = 'correlation',
          annotation_row = kinAts,annotation_col=sampAts,
          file=paste0(prefix,'KinaseHeatmap.pdf'),height=25,width=8) 
}


