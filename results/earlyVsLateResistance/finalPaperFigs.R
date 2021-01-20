##process cytokine data
#library(amlresistancenetworks)
library(dplyr)
library(ggalluvial)
library(stringr)
library(reticulate)
library(KSEAapp)
library(readr)
library(ggfortify)
library(cowplot)
##put your own conda environment with synapseclient here
condaenv="C:\\Users\\gosl241\\OneDrive - PNNL\\Documents\\GitHub\\amlresistancenetworks\\renv\\python\\r-reticulate\\"

#' Logs into Synapse using local information
#' @import reticulate
#' @return Synapse login python entity
#' @export
synapseLogin<-function(){
  library(reticulate)
  reticulate::use_condaenv(condaenv)
  syn=reticulate::import('synapseclient')
  sync=syn$login()
  return(sync)
}

#' query synapse table
#' This is how you get data from the project
#' @param tableid
#' @export
querySynapseTable<-function(tableid){
  syn=synapseLogin()
  res<-syn$tableQuery(paste('select * from',tableid))$asDataFrame()
  if('Gene'%in%names(res))
    res$Gene<-unlist(res$Gene)
  if('site'%in%names(res))
    res$site<-unlist(res$site)
  return(res)
}


#' compute kinase substrate enrichment - osama's code wrapped in package.
#' @export
#' @import KSEAapp
#' @import readr
#' @import dplyr
#' @author Osama, Jason
#' @param genes.with.values of genes and difference values
#' @param prot.univ the space of all proteins we are considering
#' @return KSEA output type stuff
computeKSEA<-function(genes.with.values,ksea_FDR=0.05,prefix=''){
  library(KSEAapp)
  
  inputdfforKSEA <- data.frame(Protein=rep("NULL", nrow(genes.with.values)), 
                               Gene=genes.with.values$Gene,
                               Peptide=rep("NULL", nrow(genes.with.values)),
                               Residue.Both=genes.with.values$residue,
                               p=genes.with.values$p_adj,
                               FC=2^(genes.with.values$value), stringsAsFactors = F)
  
  #read kinase substrate database stored in data folder
  KSDB <- read.csv('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',stringsAsFactors = FALSE)
  
  #' * KSEA using not only the known substrates in PSP but also the predicted substrates in NetworKIN
  res<-KSEA.Complete(KSDB, inputdfforKSEA, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=ksea_FDR)
  file.remove("KSEA Bar Plot.tiff")#,paste0(prefix,'_KSEABarPlot.tiff'))
  subs <- read.csv("Kinase-Substrate Links.csv")
  #make own plot of KSEA results
  res_ksea <- readr::read_csv("KSEA Kinase Scores.csv")
  file.remove("Kinase-Substrate Links.csv")#,paste0(prefix,'_kinaseSubsLinks.csv'))
  
  plot_KSEA <- res_ksea %>% 
    #mutate(p.value = p.adjust(p.value)) %>% 
    filter(m >= 5) %>% 
    arrange(desc(z.score)) %>% 
    mutate(status = case_when(z.score > 0 & p.value <= ksea_FDR ~ "Up",
                              z.score < 0 & p.value <= ksea_FDR ~ "Down",
                              TRUE ~ "Not significant")) %>% 
    filter(status != "Not significant") %>%
    ggplot(aes(x=reorder(Kinase.Gene, z.score), y=z.score)) +
    geom_bar(stat='identity', aes(fill=status)) +
    scale_fill_manual(values = c("Down" = "blue", "Up" = "red", "Not significant" = "black")) +
    coord_flip() +
    theme(legend.position="none", 
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y=element_text(size = 14),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(y="Kinase z-score") #for some reason labs still works with orientation before cord flip so set y 
  file.remove("KSEA Kinase Scores.csv")#paste0(prefix,'kinaseScores.csv'))
  
  ##join results
  #kin_res
  res_ksea<-res_ksea%>%
    rename(`aveSubstrateLog2FC`='log2FC')%>%
    left_join(subs,by='Kinase.Gene')
  
  ggsave(paste0(prefix,"fig_KSEA.pdf"), plot_KSEA, height = 8.5, width = 11, units = "in")
  write.table(res_ksea,paste0(prefix,'_kseaRes.csv'),sep=',')
  return(res_ksea)
}

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
  tidyr::pivot_longer(-c(Gene,pool),values_to='value',names_to='Patient')%>%
  rowwise()%>%
  mutate(logRatio=log10(value)-log10(pool))%>%
  dplyr::select(-c(pool,value))%>%
  dplyr::rename(value='logRatio')


protMat<-protData%>%dplyr::select(sample,Gene,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('Gene')

phosMat<-phosData%>%dplyr::select(sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('site')

kindat<-mapPhosphoToKinase(dplyr::rename(phosData,LogFoldChange='LogRatio')%>%rename(Sample='sample'))

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
    tidyr::pivot_wider(names_from='sample',values_from='LogRatio',
                       values_fn=list(LogRatio=function(x) mean(x,na.rm=T)),values_fill=list(LogRatio=0))%>%
    tibble::column_to_rownames('Gene')
  
  autoplot(prcomp(t(mat)),data=met,colour='treatment',shape='cellLine')
 
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

  
  plots=list(plotAllData(protData),plotAllData(phosData))
  cowplot::plot_grid(plotlist=plots,labels=c("Bulk Proteomics",'Phosphoprotomics'),nrow=2)
  ggsave('pcaOfSamples.png')
  

    amlresistancenetworks::plotKinDat(kindat,phosData,
                                    idcol='sample',prefix='giltResistance',
                                    vars=c("sample","cellLine","ligand"))

    clinvars = c("sample","cellLine","treatment","ligand")
  ##what are we doing again?
  summary<-protData%>%dplyr::select(clinvars)%>%distinct()%>%rowwise()%>%
    mutate(Condition=stringr::str_c(cellLine,treatment,ligand,sep='_'))
  print(summary)
  
   plotKinDat(kindat,phosData,'giltResistance',idcol='sample',c('cellLine','treatment','ligand','sample'))
  
   earlyLateProt<-list(gilt_late_early_flt3 =limmaTwoFactorDEAnalysis(protMat,
                                                        filter(summary,Condition=='MOLM14_Early Gilteritinib_FLT3')$sample,
                                                        filter(summary,Condition=='MOLM14_Late Gilteritinib_FLT3')$sample),
                      gilt_late_early_fgf=limmaTwoFactorDEAnalysis(protMat,
                                                      filter(summary,Condition=='MOLM14_Early Gilteritinib_FGF2')$sample,
                                                      filter(summary,Condition=='MOLM14_Late Gilteritinib_FGF2')$sample),
                      gilt_late_early_combined=limmaTwoFactorDEAnalysis(protMat,
                                                           filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FGF2','MOLM14_Early Gilteritinib_FLT3',
                                                                                         'MV411_Early Gilteritinib_FLT3','MV411_Early Gilteritinib_FGF2'))$sample,
                                                           filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2','MOLM14_Late Gilteritinib_FLT3',
                                                                                         'MV411_Late Gilteritinib_FGF2','MV411_Late Gilteritinib_FLT3'))$sample))
  
  
earlyLateMolmPhos<-list(gilt_early_flt3_molm14 =limmaTwoFactorDEAnalysis(phosMat,
                                                    filter(summary,Condition=='MOLM14_None_None')$sample,
                                                    filter(summary,Condition=='MOLM14_Early Gilteritinib_FLT3')$sample),
                 gilt_early_fgf_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                  filter(summary,Condition=='MOLM14_None_None')$sample,
                                                  filter(summary,Condition=='MOLM14_Early Gilteritinib_FGF2')$sample),
                 gilt_late_flt3_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                  filter(summary,Condition%in%c('MOLM14_None_None'))$sample,
                                                  filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3'))$sample),
                  gilt_late_fgf_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                    filter(summary,Condition%in%c('MOLM14_None_None'))$sample,
                                                    filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2'))$sample),
                 gilt_late_vs_early_fgf_moml14=limmaTwoFactorDEAnalysis(phosMat,
                                                                        filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FGF2'))$sample,
                                                                        filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2'))$sample),
                  gilt_late_vs_early_flt3_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FLT3'))$sample,
                                                                          filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3'))$sample))


earlyLateMv411Phos<-list( gilt_early_flt3_mv411 =limmaTwoFactorDEAnalysis(phosMat,
                                                             filter(summary,Condition=='MV411_None_None')$sample,
                                                             filter(summary,Condition=='MV411_Early Gilteritinib_FLT3')$sample),
                         gilt_early_fgf_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                           filter(summary,Condition=='MV411_None_None')$sample,
                                                           filter(summary,Condition=='MV411_Early Gilteritinib_FGF2')$sample),
                         gilt_late_flt3_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                           filter(summary,Condition%in%c('MV411_None_None'))$sample,
                                                           filter(summary,Condition%in%c('MV411_Late Gilteritinib_FLT3'))$sample),
                         gilt_late_fgf_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                          filter(summary,Condition%in%c('MV411_None_None'))$sample,
                                                          filter(summary,Condition%in%c('MV411_Late Gilteritinib_FGF2'))$sample),
                         gilt_late_vs_early_fgf_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                                        filter(summary,Condition%in%c('MV411_Early Gilteritinib_FGF2'))$sample,
                                                                        filter(summary,Condition%in%c('MV411_Late Gilteritinib_FGF2'))$sample),
                         gilt_late_vs_early_flt3_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                                          filter(summary,Condition%in%c('MV411_Early Gilteritinib_FLT3'))$sample,
                                                                          filter(summary,Condition%in%c('MV411_Late Gilteritinib_FLT3'))$sample))

earlyLateCombPhos<-list( gilt_early_vs_parental_comb =limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,Condition%in%c("MV411_None_None","MOLM14_None_None"))$sample,
                                                              filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FLT3','MV411_Early Gilteritinib_FLT3',
                                                                                            'MOLM14_Early Gilteritinib_FGF2','MV411_Early Gilteritinib_FGF2'))$sample),
                     gilt_late_vs_parental_comb=limmaTwoFactorDEAnalysis(phosMat,
                                                            filter(summary,Condition%in%c("MV411_None_None","MOLM14_None_None"))$sample,
                                                            filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3','MV411_Late Gilteritinib_FLT3',
                                                                                          'MOLM14_Late Gilteritinib_FGF2','MV411_Late Gilteritinib_FGF2'))$sample),
                     gilt_late_vs_early_comb=limmaTwoFactorDEAnalysis(phosMat,
                                               filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FLT3','MV411_Early Gilteritinib_FLT3',
                                                                             'MOLM14_Early Gilteritinib_FGF2','MV411_Early Gilteritinib_FGF2'))$sample,
                                               filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3','MV411_Late Gilteritinib_FLT3',
                                                                             'MOLM14_Late Gilteritinib_FGF2','MV411_Late Gilteritinib_FGF2'))$sample))


  
ph2<-doAllKSEAplots(earlyLateCombPhos)

ph3<-doAllKSEAplots(earlyLateMv411Phos)

ph4<-doAllKSEAplots(earlyLateMolmPhos)

##plot single kinase/substrate expression of mapk3, mapk1, and mapk8
p5<-kindat%>%
  left_join(dplyr::rename(summary,Sample='sample'))%>%
  subset(Kinase%in%c('NRAS','AURKB'))%>%
  ggplot(aes(x=as.factor(ligand),y=meanLFC,fill=treatment))+
  geom_boxplot()+
  facet_grid(~Kinase)+scale_fill_viridis_d()+facet_grid(~cellLine)+
  ggtitle("Estimated Kinase Activity")
ggsave('estimatedGiltAurbActivity.png',p5,width=10)
