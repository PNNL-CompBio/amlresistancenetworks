##process cytokine data
library(amlresistancenetworks)
library(dplyr)
library(nationalparkcolors)
pal<-park_palette('Saguaro',5)
##first run quizartinib data
amlresistancenetworks::synapseLogin()
protData<-amlresistancenetworks::querySynapseTable('syn23595222')%>%
  subset(!is.nan(value))%>%
  mutate(Gene=unlist(Gene))%>%
  rename(LogRatio='value')
#  mutate(Molecular='proteinLevels')
 
phosData<-amlresistancenetworks::querySynapseTable('syn23595223')%>%
  subset(!is.nan(value))%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))%>%
  rename(LogRatio='value')
  #  mutate(Molecular='Phosphosite')

#full.quiz.prot<-quizPhosData%>%
#  dplyr::select(Gene='site',LogRatio='value',sample='Sample',CellType='cellLine',Molecular,Treatment='Ligand')%>%
#  rbind(dplyr::select(quizProtData,c(Gene,LogRatio='value',sample='Sample',CellType='cellLine',Molecular,Treatment='Ligand')))


kindat<-amlresistancenetworks::mapPhosphoToKinase(dplyr::rename(phosData,LogFoldChange='LogRatio'))

##
#' @param dat.table
plotAllData<-function(dat.table,vars=c('Sample','cellLine','Ligand')){
  library(ggfortify)
  met<-dat.table%>%dplyr::select(vars)%>%
    distinct()
  #%>%
  #  tibble::column_to_rownames('sample')
    
  mat<-dat.table%>%dplyr::select(Gene,LogRatio,Sample)%>%
    distinct()%>%
    mutate(LogRatio=as.numeric(LogRatio))%>%
    tidyr::pivot_wider(names_from='Sample',values_from='LogRatio',values_fn=list(LogRatio=function(x) mean(x,na.rm=T)),values_fill=list(LogRatio=0))%>%
    tibble::column_to_rownames('Gene')
  
  autoplot(prcomp(t(mat)),data=met,colour='Ligand',shape='cellLine')+scale_color_manual(values=pal)
 
}


protMat<-protData%>%dplyr::select(Sample,Gene,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=Sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('Gene')

phosMat<-phosData%>%dplyr::select(Sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=Sample,
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
  
  p1<-ggplot(full.df,aes(x=logFC,y=Condition))+geom_density_ridges_gradient(aes(fill=signif))+scale_color_manual(values=pal)
  
  res.df<-full.df%>%
    subset(adj.P.Val<pvalThresh)
  if(nrow(res.df)==0)
    return(p1)
  
  p2<-res.df%>%#subset(Gene%in%changing$Gene)%>%
    ggplot(aes(x=Condition,stratum=direction,alluvium=Gene,fill=direction,label=Condition,alpha=0.5))+
    geom_flow(stat='alluvium',lode.guidance='frontback')+
    geom_stratum()+
    theme_minimal()+
    scale_fill_manual(values=pal)+
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

  amlresistancenetworks::plotKinDat(kindat,phosData,
                                    idcol='Sample',prefix='quizResistance',
                                    vars=c("Sample","cellLine","Ligand"))
  
  plots=list(plotAllData(protData),plotAllData(phosData))
  cowplot::plot_grid(plotlist=plots,labels=c("Bulk Proteomics",'Phosphoprotomics'),nrow=2)
  ggsave('pcaOfQuizSamples.png')
  
  clinvars = c("Sample","cellLine","Ligand")
  ##what are we doing again?
  summary<-protData%>%dplyr::select(clinvars)%>%distinct()%>%rowwise()%>%
    mutate(Condition=stringr::str_c(cellLine,Ligand,sep='_'))
  print(summary)
  
   earlyLateProt<-list(quiz_late_vs_early_flt3 =limmaTwoFactorDEAnalysis(protMat,
                                                        filter(summary,Condition=="EARLY Quizartinib Resistance MOLM14_FLT3 ligand")$Sample,
                                                        filter(summary,Condition=="LATE Quizartinib Resistance MOLM14_FLT3 ligand")$Sample),
                      quiz_late_vs_early_fgf=limmaTwoFactorDEAnalysis(protMat,
                                                      filter(summary,Condition=='EARLY Quizartinib Resistance MOLM14_FGF2')$Sample,
                                                      filter(summary,Condition=='LATE Quizartinib Resistance MOLM14_FGF2')$Sample),
                      quiz_late_vs_early_combined=limmaTwoFactorDEAnalysis(protMat,
                                                           filter(summary,Condition%in%c('EARLY Quizartinib Resistance MOLM14_FGF22','EARLY Quizartinib Resistance MOLM14_FLT3 ligand'))$Sample,
                                                           filter(summary,Condition%in%c('LATE Quizartinib Resistance MOLM14_FGF2','LATE Quizartinib Resistance MOLM14_FLT3 ligand'))$Sample))
  
     earlyLatePhos<-list(quiz_late_vs_early_flt3 =limmaTwoFactorDEAnalysis(phosMat,
                                                        filter(summary,Condition=="EARLY Quizartinib Resistance MOLM14_FLT3 ligand")$Sample,
                                                        filter(summary,Condition=="LATE Quizartinib Resistance MOLM14_FLT3 ligand")$Sample),
                      quiz_late_vs_early_fgf=limmaTwoFactorDEAnalysis(phosMat,
                                                      filter(summary,Condition=='EARLY Quizartinib Resistance MOLM14_FGF2')$Sample,
                                                      filter(summary,Condition=='LATE Quizartinib Resistance MOLM14_FGF2')$Sample),
                      quiz_late_vs_early_combined=limmaTwoFactorDEAnalysis(phosMat,
                                                           filter(summary,Condition%in%c('EARLY Quizartinib Resistance MOLM14_FGF22','EARLY Quizartinib Resistance MOLM14_FLT3 ligand'))$Sample,
                                                           filter(summary,Condition%in%c('LATE Quizartinib Resistance MOLM14_FGF2','LATE Quizartinib Resistance MOLM14_FLT3 ligand'))$Sample))

  earlyParentalPhos<-list(quiz_early_flt3_vs_parental=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition=='MOLM14_None')$Sample,
                                                                    filter(summary,Condition=='EARLY Quizartinib Resistance MOLM14_FLT3 ligand')$Sample),
                          quiz_early_fgf2_vs_parental=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition=='MOLM14_None')$Sample,
                                                                    filter(summary,Condition=='EARLY Quizartinib Resistance MOLM14_FGF2')$Sample),
                           quiz_early_combined_vs_parental=limmaTwoFactorDEAnalysis(phosMat,
                                                                 filter(summary,Condition=='MOLM14_None')$Sample,
                                                                filter(summary,Condition%in%c('EARLY Quizartinib Resistance MOLM14_FGF2','EARLY Quizartinib Resistance MOLM14_FLT3 ligand'))$Sample))

                        
  lateParentalPhos<-list(late_flt3_parental=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition=='MOLM14_None')$Sample,
                                                                    filter(summary,Condition=='LATE Quizartinib Resistance MOLM14_FLT3 ligand')$Sample),
                         late_fgf2_parental=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition=='MOLM14_None')$Sample,
                                                                    filter(summary,Condition=='LATE Quizartinib Resistance MOLM14_FGF2')$Sample),
                           late_combined_parental=limmaTwoFactorDEAnalysis(phosMat,
                                                                 filter(summary,Condition=='MOLM14_None')$Sample,
                                                                filter(summary,Condition%in%c('LATE Quizartinib Resistance MOLM14_FGF2','LATE Quizartinib Resistance MOLM14_FLT3 ligand'))$Sample))

  resistParentalPhos<-list(resist_parental=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition=='MOLM14_None')$Sample,
                                                                    filter(summary,Condition=='RESISTANT MOLM14_None')$Sample))



p3<-plotConditionsInFlow(earlyLatePhos,title='Early Late Phos',0.05)
ggsave('earlyLatePhos.png',p3,width=11,height=6)
ph3<-doAllKSEAplots(earlyLatePhos)
#nets<-ph3%>%mutate(Condition=stringr::str_c(Condition,'_phos'))%>%
#  runNetworksFromDF()

p4<-plotConditionsInFlow(earlyLateProt,title='Early Late Prot',0.05)
ggsave('earlyLateProt.png',p4,width=11,height=6)
#ph4<-doAllGOplots(p4)

p5<-plotConditionsInFlow(earlyParentalPhos,title='Early Parental Phos',0.05)
ggsave('earlyParentalPhos.png',p5,width=11,height=6)
ph5<-doAllKSEAplots(earlyParentalPhos)
ph6<-doAllKSEAplots(lateParentalPhos)
ph7<-doAllKSEAplots(resistParentalPhos)


#resdf<-do.call(rbind,lapply(names(earlyLateProt),function(x) data.frame(earlyLateProt[[x]],Condition=x)))

#pnets<-resdf%>%mutate(Condition=stringr::str_c(Condition,'_prot'))%>%
#  dplyr::rename(p.value='adj.P.Val')%>%
#  runNetworksFromDF(.,gene.col='featureID',weight.col='logFC',condition.col='Condition',extra.col=c('AveExpr','t','B','P.Value'),signif=0.01)


#mcp1Resistnetworks<-runNetworksFromDF(ph3)


##plot single kinase/substrate expression of mapk3, mapk1, and mapk8
p5<-kindat%>%
  left_join(dplyr::rename(summary,Sample='Sample'))%>%
  subset(Kinase%in%c('NRAS','AURKB'))%>%
  ggplot(aes(x=as.factor(Ligand),y=meanLFC,fill=cellLine))+
  geom_boxplot()+
  facet_grid(~Kinase)+scale_fill_viridis_d()+
  ggtitle("Estimated Kinase Activity")
ggsave('estimatedQuizAurbActivity.png',p5,width=10)
}