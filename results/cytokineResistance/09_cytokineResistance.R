##process cytokine data
library(amlresistancenetworks)
library(dplyr)



###load all the data
protData<-querySynapseTable('syn22986326')%>%subset(!is.nan(LogRatio))%>%
  mutate(Gene=unlist(Gene))

phosData<-querySynapseTable('syn22986341')%>%subset(!is.nan(LogRatio))%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))

#otherData<-''
otherPhosData<-querySynapseTable('syn22255396')%>%
  subset(!is.nan(LogFoldChange))%>%
  subset(cellLine=='HL60')%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))%>%
  rowwise()%>%
  mutate(Condition=paste(treatment,timePoint,sep='_'))


clinvars<-phosData%>%
  dplyr::select(Sample='sample',CellType,TimePoint,Treatment)%>%
  distinct()

kindat<-mapPhosphoToKinase(dplyr::rename(phosData,Sample='sample', LogFoldChange='LogRatio'))

parental<-mapPhosphoToKinase(dplyr::rename(filter(phosData,CellType=='MOLM-13'),Sample='sample', LogFoldChange='LogRatio'))
                                    


##what are we doing again?
summary<-protData%>%
  dplyr::select(sample,CellType,TimePoint,Treatment)%>%
  distinct()%>%
  mutate(conditionName=stringr::str_c(CellType,TimePoint,Treatment,sep='_'))

print(summary)


protMat<-protData%>%dplyr::select(sample,Gene,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('Gene')

phosMat<-phosData%>%dplyr::select(sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('site')

otherPhosMat <-otherPhosData%>%ungroup()%>%dplyr::select(Sample,site,LogFoldChange)%>%
  tidyr::pivot_wider(values_from=LogFoldChange,names_from=Sample,
                     values_fn=list(LogFoldChange=mean),values_fill=list(LogFoldChange=0.0))%>%
  tibble::column_to_rownames('site')

tpm<-apply(otherPhosMat,2,as.numeric)
rownames(tpm)<-rownames(otherPhosMat)
otherPhosMat<-tpm

otherSummary<-otherPhosData%>%dplyr::select(Sample,Condition)%>%distinct()

##
#' @param dat.table
plotAllData<-function(dat.table){
  library(ggfortify)
  met<-dat.table%>%dplyr::select(sample,CellType,TimePoint,Treatment)%>%
    distinct()
  #%>%
  #  tibble::column_to_rownames('sample')
    
  mat<-dat.table%>%dplyr::select(Gene,LogRatio,sample)%>%
    distinct()%>%
    mutate(LogRatio=as.numeric(LogRatio))%>%
    tidyr::pivot_wider(names_from='sample',values_from='LogRatio',values_fn=list(LogRatio=function(x) mean(x,na.rm=T)),values_fill=list(LogRatio=0))%>%
    tibble::column_to_rownames('Gene')
  
  autoplot(prcomp(t(mat)),data=met,colour='Treatment',shape='CellType')
 
}

plotKinDat<-function(kindat,prefix='all'){
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
    tibble::column_to_rownames('Kinase')
  
  sampAts<-phosData%>%
    dplyr::select(sample,CellType,TimePoint,Treatment)%>%
    distinct()%>%
    tibble::column_to_rownames('sample')
  
  sampAts$TimePoint=as.factor(sampAts$TimePoint)
  
  vars=names(sort(apply(mat,1,var),decreasing=T)[1:150])
 
  pheatmap(mat[vars,],cellwidth = 8,cellheight=8,clustering_distance_cols = 'correlation',
          clustering_distance_rows = 'correlation',
          annotation_row = kinAts,annotation_col=sampAts,
          file=paste0(prefix,'cytokineKinaseHeatmap.pdf'),height=20,width=8) 
}

plotKinDat(kindat)
plotKinDat(parental,'molm13')

plots=list(plotAllData(protData),plotAllData(phosData))
cowplot::plot_grid(plotlist=plots,labels=c("Bulk Proteomics",'Phosphoprotomics'),nrow=2)
ggsave('pcaOfSamples.png')

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
  return(full.df)
  
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
                              condition.col='Condition',extra.col=c('Substrate.Gene','Source','log2FC'),
                              signif=0.05){
  res = data%>%
   # dplyr::select(cond=condition.col,value=weight.col,Gene=gene.col,p.value)%>%
    mutate(signif=p.value<signif)%>%
      dplyr::select(c(condition.col,weight.col,gene.col,'signif',extra.col))%>%distinct()%>%
    dplyr::rename(cond=condition.col,value=weight.col,Gene=gene.col)%>%
    group_by(cond)%>%
    dplyr::select(c('cond','Gene','value',extra.col,'signif'))%>%
    group_map(~ amlresistancenetworks::computeProteinNetwork(.x),keep=TRUE)
  return(res)
}

#####now do various comparisons
doLateComparisons<-function(){
  
  latePhos<-list(lateTram_vs_lateCombo=limmaTwoFactorDEAnalysis(phosMat,
                                                          filter(summary,conditionName=='Late MOLM-13_0_Trametinib')$sample,
                                                          filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample),
               m13_vs_lateTram=limmaTwoFactorDEAnalysis(phosMat,
                                                filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                filter(summary,conditionName=='Late MOLM-13_0_Trametinib')$sample),
               m13_vs_lateCombo=limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                                filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample))
  lateProt<-list(late_tram_vs_mcp1=limmaTwoFactorDEAnalysis(protMat,
                                                          filter(summary,conditionName=='Late MOLM-13_0_Trametinib')$sample,
                                                          filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample),
               m13_vs_lateTram=limmaTwoFactorDEAnalysis(protMat,
                                                        filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                        filter(summary,conditionName=='Late MOLM-13_0_Trametinib')$sample),
               m13_vs_lateCombo=limmaTwoFactorDEAnalysis(protMat,
                                                         filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                         filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample))

  earlyLatePhos<-list(early_60m_combo_vs_late_combo=limmaTwoFactorDEAnalysis(phosMat,
                                                                         filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample,
                                   filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample),
                      early_5m_combo_vs_late_combo=limmaTwoFactorDEAnalysis(phosMat,
                                                                         filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample,
                                   filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample),
                      early_both_combo_vs_late_combo=limmaTwoFactorDEAnalysis(phosMat,
                                                                               filter(summary,conditionName%in%c('MOLM-13_5_Trametinib+MCP-1','MOLM-13_60_Trametinib+MCP-1'))$sample,
                                                                              filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample))
                      
  earlyLateProt<-list(early_60m_combo_vs_late_combo=limmaTwoFactorDEAnalysis(protMat,
                                                                             filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample,
                                                                             filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample),
                      early_5m_combo_vs_late_combo=limmaTwoFactorDEAnalysis(protMat,
                                                                            filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample,
                                                                            filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample),
                      early_both_combo_vs_late_combo=limmaTwoFactorDEAnalysis(protMat,
                                                                              filter(summary,conditionName%in%c('MOLM-13_5_Trametinib+MCP-1','MOLM-13_60_Trametinib+MCP-1'))$sample,
                                                                              filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample))
  

lateP<-plotConditionsInFlow(lateProt,title='Bulk Proteomics in late',0.05)
ggsave('lateProt.png',lateP,width=11,height=6)
r2<-doAllGOplots(lateProt)              

latePh<-plotConditionsInFlow(latePhos,title='Phosphoproteomics in late',0.05)
ggsave('latePhos.png',latePh,width=11,height=6)

phresdf<-do.call(rbind,lapply(names(latePhos),function(x) data.frame(latePhos[[x]],Condition=x)))

ph3<-doAllKSEAplots(latePhos)


resdf<-do.call(rbind,lapply(names(lateProt),function(x) data.frame(lateProt[[x]],Condition=x)))

#pnets<-resdf%>%mutate(Condition=stringr::str_c(Condition,'_prot'))%>%
#  dplyr::rename(p.value='adj.P.Val')%>%
#  runNetworksFromDF(gene.col='featureID',weight.col='logFC',
#                    condition.col='Condition',extra.col=c('AveExpr','t','B','P.Value'),
#                    signif=0.005)


#lateNets<-runNetworksFromDF(ph3)

earlyLateP<-plotConditionsInFlow(earlyLateProt,title='Bulk Proteomics in early vs late',0.05)
ggsave('earlyLateProt.png',earlyLateP,width=11,height=6)
r2<-doAllGOplots(earlyLateProt)              

earlyLatePh<-plotConditionsInFlow(earlyLatePhos,title='Phosphoproteomics in late',0.05)
ggsave('earlyLatePhos.png',earlyLatePh,width=11,height=6)

earlyLatePhresdf<-do.call(rbind,lapply(names(earlyLatePhos),function(x) data.frame(earlyLatePhos[[x]],Condition=x)))

ph3<-doAllKSEAplots(earlyLatePhos)
lateNets<-runNetworksFromDF(ph3)


}


compareCellLineSig<-function(){
tramHLPhos<-list(hl60_tram_3_hour =manualDEAnalysis(otherPhosMat,
                                                    filter(otherSummary,Condition=='trametinib_3 hr')$Sample,
                                                    filter(otherSummary,Condition=='no treatment_3 hr')$Sample),
                 hl60_tram_16_hr=manualDEAnalysis(otherPhosMat,
                                                  filter(otherSummary,Condition=='trametinib_16 hr')$Sample,
                                                  filter(otherSummary,Condition=='no treatment_3 hr')$Sample),
                 hl60_tram_30min=manualDEAnalysis(otherPhosMat,
                                                  filter(otherSummary,Condition=='trametinib_30 min')$Sample,
                                                  filter(otherSummary,Condition=='no treatment_3 hr')$Sample))

p3<-plotConditionsInFlow(tramHLPhos,title='Effects of tram in hL60',0.05)
ggsave('tramHL60.png',p3,width=11,height=6)
ph3<-doAllKSEAplots(tramHLPhos,otherPhosData)
#mcp1Resistnetworks<-runNetworksFromDF(ph3)

}

compareResistantCells<-function(){

##now compute differences in conditions
t0Values<-list(molm13_vs_resistant=limmaTwoFactorDEAnalysis(protMat,
                                 filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                 filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample),
               resistant_vs_tramWithdrawn=limmaTwoFactorDEAnalysis(protMat,
                                  filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample,
                                  filter(summary,conditionName=='MOLM-13 Tr Resistant_0_Trametinib Withdrawn')$sample),
                molm13_vs_tramWithdrawn=limmaTwoFactorDEAnalysis(protMat,
                                                    filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                    filter(summary,conditionName=='MOLM-13 Tr Resistant_0_Trametinib Withdrawn')$sample))

p<-plotConditionsInFlow(t0Values,title='Effect of trametinib on resistance')    

ggsave("TrametinibResistanceConditions.png",p,width=11,height=6)
r1<-doAllGOplots(t0Values)

resdf<-do.call(rbind,lapply(names(t0Values),function(x) data.frame(t0Values[[x]],Condition=x)))

pnets<-resdf%>%mutate(Condition=stringr::str_c(Condition,'_prot'))%>%
  dplyr::rename(p.value='adj.P.Val')%>%
  runNetworksFromDF(gene.col='featureID',weight.col='logFC',condition.col='Condition',extra.col=c('AveExpr','t','B','P.Value'),signif=0.0001)


t0Phos=list(molm13_vs_resistant=limmaTwoFactorDEAnalysis(phosMat,
                                                                   filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                                   filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample),
                      resistant_vs_tramWithdrawn=limmaTwoFactorDEAnalysis(phosMat,
                                                                          filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample,
                                                                          filter(summary,conditionName=='MOLM-13 Tr Resistant_0_Trametinib Withdrawn')$sample),
                      molm13_vs_tramWithdrawn=limmaTwoFactorDEAnalysis(phosMat,
                                                                       filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                                       filter(summary,conditionName=='MOLM-13 Tr Resistant_0_Trametinib Withdrawn')$sample))

ph<-plotConditionsInFlow(t0Phos,title='Phoshpochanges - trametinib on resistance')   
ggsave('tramResistancePhos.png',ph,width=11,height=6)
ksea.res<-doAllKSEAplots(t0Phos)
#t0PhosNetworks=runNetworksFromDF(ksea.res)
}


compareTramTreatmentCombos<-function(){
  #' here we get various treatments of MCP1 and tram at different time points
  m13Values<-list(molm13_vs_tram_5min=limmaTwoFactorDEAnalysis(protMat,
                                                          filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                        filter(summary,conditionName=='MOLM-13_5_Trametinib')$sample),
                  molm13_vs_tram_60min=limmaTwoFactorDEAnalysis(protMat,
                                                               filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                               filter(summary,conditionName=='MOLM-13_60_Trametinib')$sample),
                  molm13_vs_MCP1_5min=limmaTwoFactorDEAnalysis(protMat,
                                                               filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                               filter(summary,conditionName=='MOLM-13_5_MCP-1')$sample),
                  molm13_vs_MCP1_60min=limmaTwoFactorDEAnalysis(protMat,
                                                              filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                              filter(summary,conditionName=='MOLM-13_60_MCP-1')$sample),
                  molm13_vs_MCP1_tram_5min=limmaTwoFactorDEAnalysis(protMat,
                                                             filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                             filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                molm13_vs_MCP1_tram_60min=limmaTwoFactorDEAnalysis(protMat,
                                                              filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                              filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample))

m13Phos<-list(molm13_vs_tram_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                             filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                             filter(summary,conditionName=='MOLM-13_5_Trametinib')$sample),
                molm13_vs_tram_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                              filter(summary,conditionName=='MOLM-13_60_Trametinib')$sample),
                molm13_vs_MCP1_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                             filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                             filter(summary,conditionName=='MOLM-13_5_MCP-1')$sample),
                molm13_vs_MCP1_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                              filter(summary,conditionName=='MOLM-13_60_MCP-1')$sample),
                molm13_vs_MCP1_tram_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                                  filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                                  filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                molm13_vs_MCP1_tram_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                                   filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                                   filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample))

p2<-plotConditionsInFlow(m13Values,title='Effects of tram with MCP1',0.05)
ggsave('Molm13Conditions.png',p2,width=11,height=6)
r2<-doAllGOplots(m13Values)              

p3<-doAllKSEAplots(m13Phos)

#m13phosNetworks=runNetworksFromDF(p3)

ph3<-plotConditionsInFlow(m13Phos,title='Phospho effects of tram with MCP1',0.05)
ggsave('molm13ConditionsPhos.png',ph3,width=11,height=6)
          


tramMCPValues<-list(tram_vs_mcp1tram_5min=limmaTwoFactorDEAnalysis(protMat,
                                                             filter(summary,conditionName=='MOLM-13_5_Trametinib')$sample,
                                                             filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                    tram_vs_mcp1tram_60min=limmaTwoFactorDEAnalysis(protMat,
                                                                   filter(summary,conditionName=='MOLM-13_60_Trametinib')$sample,
                                                                   filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample),
                    mcp1_vs_mcp1tram_5min=limmaTwoFactorDEAnalysis(protMat,
                                                                   filter(summary,conditionName=='MOLM-13_5_MCP-1')$sample,
                                                                   filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                    mcp1_vs_mcp1tram_60min=limmaTwoFactorDEAnalysis(protMat,
                                                                    filter(summary,conditionName=='MOLM-13_60_MCP-1')$sample,
                                                                    filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample))

p3<-plotConditionsInFlow(tramMCPValues,title='Effects of tram vs combo1',0.05)
ggsave('comboVsMCPTramIndiv.png',p3,width=11,height=6)

r3<-doAllGOplots(tramMCPValues)              


tramMCPPhos<-list(tram_vs_mcp1tram_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                                   filter(summary,conditionName=='MOLM-13_5_Trametinib')$sample,
                                                                   filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                    tram_vs_mcp1tram_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,conditionName=='MOLM-13_60_Trametinib')$sample,
                                                                    filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample),
                    mcp1_vs_mcp1tram_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                                   filter(summary,conditionName=='MOLM-13_5_MCP-1')$sample,
                                                                   filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                    mcp1_vs_mcp1tram_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,conditionName=='MOLM-13_60_MCP-1')$sample,
                                                                    filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample),
                  late_tram_vs_mcp1tram_5min=manualDEAnalysis(phosMat,
                                                                 filter(summary,conditionName=='Late MOLM-13_5_Trametinib')$sample,
                                                                 filter(summary,conditionName=='Late MOLM-13_5_Trametinib+MCP-1')$sample),
                  late_tram_vs_mcp1tram_60min=manualDEAnalysis(phosMat,
                                                                  filter(summary,conditionName=='Late MOLM-13_60_Trametinib')$sample,
                                                                  filter(summary,conditionName=='Late MOLM-13_60_Trametinib+MCP-1')$sample),
                  late_mcp1_vs_mcp1tram_5min=manualDEAnalysis(phosMat,
                                                                 filter(summary,conditionName=='Late MOLM-13_5_MCP-1')$sample,
                                                                 filter(summary,conditionName=='Late MOLM-13_5_Trametinib+MCP-1')$sample),
                  late_mcp1_vs_mcp1tram_60min=manualDEAnalysis(phosMat,
                                                                  filter(summary,conditionName=='Late MOLM-13_60_MCP-1')$sample,
                                                                  filter(summary,conditionName=='Late MOLM-13_60_Trametinib+MCP-1')$sample)
                  )

ph3<-doAllKSEAplots(tramMCPPhos)
tramMCP=runNetworksFromDF(ph3)

ph3<-plotConditionsInFlow(tramMCPPhos,title='Phospho effects of tram vs combo',0.05)
ggsave('phoscomboVsMCPTramIndiv.png',ph3,width=11,height=6)
}


compareTramInResistCells<-function(){
  #' how does trametinib affect resistant cells?
##des MCP1 loo like resistant cells?
mcp1Resist<-list(resist_vs_mcp1_5min=limmaTwoFactorDEAnalysis(protMat,
                                                              filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample,
                                                              filter(summary,conditionName=='MOLM-13 Tr Resistant_5_MCP-1')$sample),
                 resist_vs_mcp1_60min=limmaTwoFactorDEAnalysis(protMat,
                                                               filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample,
                                                               filter(summary,conditionName=='MOLM-13 Tr Resistant_60_MCP-1')$sample))

p3<-plotConditionsInFlow(mcp1Resist,title='Effects of mcp-1 in resistant',0.05)
ggsave('mcp1InResistantCells.png',p3,width=11,height=6)

r3<-doAllGOplots(mcp1Resist)              

mcp1ResistPhos<-list(resist_vs_mcp1_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample,
                                                              filter(summary,conditionName=='MOLM-13 Tr Resistant_5_MCP-1')$sample),
                 resist_vs_mcp1_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                               filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample,
                                                               filter(summary,conditionName=='MOLM-13 Tr Resistant_60_MCP-1')$sample))

p3<-plotConditionsInFlow(mcp1ResistPhos,title='Effects of tram/mcp-1 in resistant',0.05)
ggsave('mcp1PhosInResistantCells.png',p3,width=11,height=6)
ph3<-doAllKSEAplots(mcp1ResistPhos)
#mcp1Resistnetworks<-runNetworksFromDF(ph3)
}

##now pool time points for bulk differences

##plot single kinase/substrate expression of mapk3, mapk1, and mapk8
p5<-kindat%>%
  left_join(clinvars)%>%
  #subset(Treatment%in%c("MCP-1","none", 'Trametinib Withdrawn'))%>%
  subset(Kinase%in%c('MAPK1','MAPK3','MAPK8','PRKCD'))%>%
  ggplot(aes(x=as.factor(TimePoint),y=meanLFC,fill=Kinase))+
  geom_boxplot()+
  facet_grid(~CellType+Treatment)+scale_fill_viridis_d()+
  ggtitle("Estimated Kinase Activity")
ggsave('estimatedErkActivity.png',p5,width=10)