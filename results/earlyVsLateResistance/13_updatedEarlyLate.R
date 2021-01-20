##process cytokine data
library(amlresistancenetworks)
library(dplyr)

##first run gilteritinib data
protData<-querySynapseTable('syn22156807')%>%mutate(Gene=unlist(Gene))%>%
  dplyr::rename(sample='Sample')%>%
  dplyr::rename(LogRatio='value')
#readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
#prot.univ<-unique(gilt.data$Gene)
phosData<-querySynapseTable('syn22156809')%>%subset(!is.nan(value))%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))%>%
  dplyr::rename(LogRatio='value')%>%
  dplyr::rename(sample='Sample')


syn<-synapseLogin()
targ.data<-read.csv2(syn$get('syn22214245')$path,sep='\t',header=T)%>%
  dplyr::rename(Gene='T..Category')%>%
  tidyr::pivot_longer(-Gene,names_to='sample',values_to='value')%>%
  mutate(value=as.numeric(value))%>%
  rowwise()%>%
  mutate(samp2=stringr::str_replace(sample,'_3June20.+',''))%>%
  mutate(Patient=stringr::str_replace(samp2,'.*Ex16_',''))%>%
  dplyr::select(-c(samp2,sample))

norm.data<-targ.data%>%tidyr::pivot_wider(values_from=value,names_from=Patient)%>%
  tidyr::pivot_longer(-c(Gene,pool),values_to='value',names_to='Patient')%>%rowwise()%>%
  mutate(logRatio=log10(value)-log10(pool))%>%
  dplyr::select(-c(pool,value))%>%
  dplyr::rename(value='logRatio')

kindat<-mapPhosphoToKinase(dplyr::rename(phosData,LogFoldChange='LogRatio',Sample='sample'))

##
#' @param dat.table
plotAllData<-function(dat.table,vars=c('sample','cellLine','ligand','treatment')){
  library(ggfortify)
  met<-dat.table%>%dplyr::select(vars)%>%
    distinct()
  #%>%
  #  tibble::column_to_rownames('sample')
    
  mat<-dat.table%>%dplyr::select(Gene,LogRatio,sample)%>%
    distinct()%>%
    mutate(LogRatio=as.numeric(LogRatio))%>%
    tidyr::pivot_wider(names_from='sample',values_from='LogRatio',values_fn=list(LogRatio=function(x) mean(x,na.rm=T)),values_fill=list(LogRatio=0))%>%
    tibble::column_to_rownames('Gene')
  
  autoplot(prcomp(t(mat)),data=met,colour='treatment',shape='cellLine')
 
}

plotKinDat<-function(kindat,prefix='all',vars=c('sample','cellLine','ligand','treatment')){
  library(pheatmap)
  ##create matrix of kinase scores
  mat <-kindat%>%ungroup()%>%tidyr::pivot_wider(-c(meanNKINscore,numSubstr),
                                                values_from=meanLFC,
                                                names_from=Sample,
                                                values_fn=list(meanLFC=mean))%>%
    tibble::column_to_rownames('Kinase')
  kinAts<-kindat%>%ungroup()%>%dplyr::select(Kinase,numSubstr)%>%distinct()%>%
    group_by(Kinase)%>%summarize(substrates=mean(numSubstr))%>%
  tibble::column_to_rownames('Kinase')
  
  sampAts<-phosData%>%dplyr::select(vars)%>%
    distinct()%>%
    tibble::column_to_rownames('sample')
  #sampAts$TimePoint=as.factor(sampAts$TimePoint)
  vars=names(sort(apply(mat,1,var),decreasing=T)[1:150])
 pheatmap(mat[vars,],cellwidth = 8,cellheight=8,clustering_distance_cols = 'correlation',
          clustering_distance_rows = 'correlation',
          annotation_row = kinAts,annotation_col=sampAts,
          file=paste0(prefix,'KinaseHeatmap.pdf'),height=20,width=8) 
}



protMat<-protData%>%dplyr::select(sample,Gene,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('Gene')

phosMat<-phosData%>%dplyr::select(sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('site')


library(ggplot2)
library(ggalluvial)
library(ggridges)
#' plotConditionsInFlow
#' Tries to compare multiple logfc values
plotConditionsInFlow<-function(condList,title='',pvalThresh=0.05,upDown='logFC'){
  
  full.df<-purrr::map_df(names(condList),.f=function(clName){ 
    condList[[clName]]%>%
      tibble::rownames_to_column('Gene')%>%
      dplyr::select(Gene,logFC,adj.P.Val)%>%
      mutate(Condition=clName)%>%
      mutate(direction=ifelse(logFC>0,'UpReg','DownReg'))%>%
      mutate(signif=ifelse(adj.P.Val<pvalThresh,'Significant','Non-significant'))
  })
  
  p1<-ggplot(full.df,aes(x=logFC,y=Condition))+geom_density_ridges_gradient(aes(fill=signif))
  
  res.df<-full.df%>%
    subset(adj.P.Val<pvalThresh)
  if(nrow(res.df)==0)
    return(p1)
  
  p2<-res.df%>%#subset(Gene%in%changing$Gene)%>%
    ggplot(aes(x=Condition,stratum=direction,alluvium=Gene,fill=direction,label=Condition,alpha=0.5))+
    geom_flow(stat='alluvium',lode.guidance='frontback')+
    geom_stratum()+
    theme_minimal()+
    viridis::scale_fill_viridis(3,discrete=T)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p<-cowplot::plot_grid(p1,p2,nrow=1)+ggtitle(title)
  return(p)
}



#' plot all the GO 
#' @param condList
#' @return data frame
doAllGOplots<-function(condList){
  
  full.df<-purrr::map_df(names(condList),.f=function(clName){ 
    condList[[clName]]%>%
      tibble::rownames_to_column('Gene')%>%
      dplyr::select(Gene,value='logFC')%>%
      amlresistancenetworks::plotOldGSEA(.,prefix=clName,0.1)%>%
      as.data.frame()
  })
  
}

#' plot all the KSEA 
#' @param condList
#' @return data frame
doAllKSEAplots<-function(condList,pdat=phosData){
  
  gene.to.site<-dplyr::select(pdat,Gene,site,Peptide)%>%distinct()%>%
    dplyr::mutate(residue=stringr::str_replace(site,paste0(Gene,'-'),''))%>%
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([STY])", ";\\1"))%>%
    dplyr::mutate(residue=stringr::str_replace(residue,"^;", ""))%>%
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([sty])", ""))
  
  full.df<-purrr::map_df(names(condList),.f=function(clName){ 
    condList[[clName]]%>%
      tibble::rownames_to_column('site')%>%
      left_join(gene.to.site)%>%
      dplyr::select(Gene,Peptide,residue,value='logFC',p_adj='adj.P.Val')%>%
      amlresistancenetworks::computeKSEA(.,prefix=clName,0.05)%>%
      mutate(Condition=clName)%>%
      as.data.frame()
  })
  return(full.df)
  
}

#' build networks from data frame
#' @param data.res
#' @param gene.col
#' @param weight.col
#' @param condition.col
#' @return network list?
runNetworksFromDF<-function(data,gene.col='Kinase.Gene',
                              weight.col='aveSubstrateLog2FC',
                              condition.col='Condition',
                            extra.col=c('Substrate.Gene','Source','log2FC'),
                              signif=0.05){
  res = data%>%
   # dplyr::select(cond=condition.col,value=weight.col,Gene=gene.col,p.value)%>%
    mutate(signif=p.value<signif)%>%
      dplyr::select(c(condition.col,weight.col,gene.col,'signif',extra.col))%>%distinct()%>%
    dplyr::rename(cond=condition.col,value=weight.col,Gene=gene.col)%>%
    dplyr::select(c('cond','Gene','value',extra.col,'signif'))%>%
    group_by(cond)%>%
    group_map(~ amlresistancenetworks::computeProteinNetwork(.x),keep=TRUE)
  return(res)
}


doPlots=TRUE
if(doPlots){

  plotKinDat(kindat,'earlyLate')
  
  plots=list(plotAllData(protData),plotAllData(phosData))
  cowplot::plot_grid(plotlist=plots,labels=c("Bulk Proteomics",'Phosphoprotomics'),nrow=2)
  ggsave('pcaOfSamples.png')
  
  clinvars = c("sample","cellLine","treatment","ligand")
  ##what are we doing again?
  summary<-protData%>%dplyr::select(clinvars)%>%distinct()%>%rowwise()%>%
    mutate(Condition=stringr::str_c(cellLine,treatment,ligand,sep='_'))
  print(summary)
  
   earlyLateProt<-list(early_late_flt3 =limmaTwoFactorDEAnalysis(protMat,
                                                        filter(summary,Condition=='MOLM14_Early Gilteritinib_FLT3')$sample,
                                                        filter(summary,Condition=='MOLM14_Late Gilteritinib_FLT3')$sample),
                      early_late_fgf=limmaTwoFactorDEAnalysis(protMat,
                                                      filter(summary,Condition=='MOLM14_Early Gilteritinib_FGF2')$sample,
                                                      filter(summary,Condition=='MOLM14_Late Gilteritinib_FGF2')$sample),
                      early_late_combined=limmaTwoFactorDEAnalysis(protMat,
                                                           filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FGF2','MOLM14_Early Gilteritinib_FLT3'))$sample,
                                                           filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2','MOLM14_Late Gilteritinib_FLT3'))$sample))
  
  
earlyLatePhos<-list(early_late_flt3 =limmaTwoFactorDEAnalysis(phosMat,
                                                    filter(summary,Condition=='MOLM14_Early Gilteritinib_FLT3')$sample,
                                                    filter(summary,Condition=='MOLM14_Late Gilteritinib_FLT3')$sample),
                 early_late_fgf=limmaTwoFactorDEAnalysis(phosMat,
                                                  filter(summary,Condition=='MOLM14_Early Gilteritinib_FGF2')$sample,
                                                  filter(summary,Condition=='MOLM14_Late Gilteritinib_FGF2')$sample),
                 early_late_combined=limmaTwoFactorDEAnalysis(phosMat,
                                                  filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FGF2','MOLM14_Early Gilteritinib_FLT3'))$sample,
                                                  filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2','MOLM14_Late Gilteritinib_FLT3'))$sample))
earlyLateCombM13Phos<-list(early_both_m13_comb =limmaTwoFactorDEAnalysis(phosMat,
                                                                  filter(summary,Condition%in%c("MOLM14_None_None"))$sample,
                                                                  filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FLT3','MOLM14_Early Gilteritinib_FGF2'))$sample),
                        late_both_m13_comb=limmaTwoFactorDEAnalysis(phosMat,
                                                                filter(summary,Condition%in%c("MOLM14_None_None"))$sample,
                                                                filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3','MOLM14_Late Gilteritinib_FGF2'))$sample))
#early_late_both_combined=limmaTwoFactorDEAnalysis(phosMat,

earlyLateCombPhos<-list(early_both_comb =limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,Condition%in%c("MV411_None_None","MOLM14_None_None"))$sample,
                                                              filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FLT3','MV411_Early Gilteritinib_FLT3','MOLM14_Early Gilteritinib_FGF2','MV411_Early Gilteritinib_FGF2'))$sample),
                    late_both_comb=limmaTwoFactorDEAnalysis(phosMat,
                                                            filter(summary,Condition%in%c("MV411_None_None","MOLM14_None_None"))$sample,
                                                            filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3','MV411_Late Gilteritinib_FLT3','MOLM14_Late Gilteritinib_FGF2','MV411_Late Gilteritinib_FGF2'))$sample))
                    #early_late_both_combined=limmaTwoFactorDEAnalysis(phosMat,
                    #                                             filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FGF2','MOLM14_Early Gilteritinib_FLT3'))$sample,
                    #                                             filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2','MOLM14_Late Gilteritinib_FLT3'))$sample))


p1<-plotConditionsInFlow(earlyLateCombPhos,title='Combined results',0.05)
ggsave("earlyLateCombinedPhos.png",width=11,height=6)
doAllKSEAplots(earlyLateCombPhos)

p3<-plotConditionsInFlow(earlyLatePhos,title='Early Late Phos',0.05)
ggsave('earlyLatePhos.png',p3,width=11,height=6)
ph3<-doAllKSEAplots(earlyLatePhos)
#nets<-ph3%>%mutate(Condition=stringr::str_c(Condition,'_phos'))%>%
#xs  runNetworksFromDF()

ph5<-doAllKSEAplots(earlyLateCombM13Phos)


p4<-plotConditionsInFlow(earlyLateProt,title='Early Late Prot',0.05)
ggsave('earlyLateProt.png',p4,width=11,height=6)
ph4<-doAllGOplots(ph3)

resdf<-do.call(rbind,lapply(names(earlyLateProt),function(x) data.frame(earlyLateProt[[x]],Condition=x)))

#pnets<-resdf%>%mutate(Condition=stringr::str_c(Condition,'_prot'))%>%
#  dplyr::rename(p.value='adj.P.Val')%>%
#  runNetworksFromDF(.,gene.col='featureID',weight.col='logFC',condition.col='Condition',extra.col=c('AveExpr','t','B','P.Value'),signif=0.01)


#mcp1Resistnetworks<-runNetworksFromDF(ph3)


##plot single kinase/substrate expression of mapk3, mapk1, and mapk8
p5<-kindat%>%
  left_join(dplyr::rename(summary,Sample='sample'))%>%
  subset(Kinase%in%c('NRAS','AURKB'))%>%
  ggplot(aes(x=as.factor(ligand),y=meanLFC,fill=treatment))+
  geom_boxplot()+
  facet_grid(~cellLine+Kinase)+scale_fill_viridis_d()+
  ggtitle("Estimated Kinase Activity")
ggsave('estimatedAurbActivity.png',p5,width=10)
}