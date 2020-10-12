##process cytokine data
#library(amlresistancenetworks)
library(dplyr)
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


protMat<-protData%>%dplyr::select(sample,Gene,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('Gene')

phosMat<-phosData%>%dplyr::select(sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('site')


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
                                                           filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FGF2','MOLM14_Early Gilteritinib_FLT3',
                                                                                         'MV411_Early Gilteritinib_FLT3','MV411_Early Gilteritinib_FGF2'))$sample,
                                                           filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2','MOLM14_Late Gilteritinib_FLT3',
                                                                                         'MV411_Late Gilteritinib_FGF2','MV411_Late Gilteritinib_FLT3'))$sample))
  
  
earlyLateMolmPhos<-list(early_flt3_molm14 =limmaTwoFactorDEAnalysis(phosMat,
                                                    filter(summary,Condition=='MOLM14_None_None')$sample,
                                                    filter(summary,Condition=='MOLM14_Early Gilteritinib_FLT3')$sample),
                 early_fgf_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                  filter(summary,Condition=='MOLM14_None_None')$sample,
                                                  filter(summary,Condition=='MOLM14_Early Gilteritinib_FGF2')$sample),
                 late_flt3_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                  filter(summary,Condition%in%c('MOLM14_None_None'))$sample,
                                                  filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3'))$sample),
                 late_fgf_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                    filter(summary,Condition%in%c('MOLM14_None_None'))$sample,
                                                    filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2'))$sample))

earlyLateMv411Phos<-list(early_flt3_mv411 =limmaTwoFactorDEAnalysis(phosMat,
                                                             filter(summary,Condition=='MV411_None_None')$sample,
                                                             filter(summary,Condition=='MV411_Early Gilteritinib_FLT3')$sample),
                        early_fgf_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                           filter(summary,Condition=='MV411_None_None')$sample,
                                                           filter(summary,Condition=='MV411_Early Gilteritinib_FGF2')$sample),
                        late_flt3_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                           filter(summary,Condition%in%c('MV411_None_None'))$sample,
                                                           filter(summary,Condition%in%c('MV411_Late Gilteritinib_FLT3'))$sample),
                        late_fgf_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                          filter(summary,Condition%in%c('MV411_None_None'))$sample,
                                                          filter(summary,Condition%in%c('MV411_Late Gilteritinib_FGF2'))$sample))

earlyLateCombPhos<-list(early_vs_parental_comb =limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,Condition%in%c("MV411_None_None","MOLM14_None_None"))$sample,
                                                              filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FLT3','MV411_Early Gilteritinib_FLT3',
                                                                                            'MOLM14_Early Gilteritinib_FGF2','MV411_Early Gilteritinib_FGF2'))$sample),
                    late_vs_parental_comb=limmaTwoFactorDEAnalysis(phosMat,
                                                            filter(summary,Condition%in%c("MV411_None_None","MOLM14_None_None"))$sample,
                                                            filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3','MV411_Late Gilteritinib_FLT3',
                                                                                          'MOLM14_Late Gilteritinib_FGF2','MV411_Late Gilteritinib_FGF2'))$sample))
                    

  
p1<-plotConditionsInFlow(earlyLateCombPhos,title='MOLM14 and MV411',0.05)
ggsave("earlyLateCombinedPhos.png",width=11,height=6)
doAllKSEAplots(earlyLateCombPhos)

#p3<-plotConditionsInFlow(earlyLatePhos,title='Early Late Phos',0.05)
#ggsave('earlyLatePhos.png',p3,width=11,height=6)
ph3<-doAllKSEAplots(earlyLateMv411Phos)

ph3<-doAllKSEAplots(earlyLateMolmPhos)


##now do heatmap

#}