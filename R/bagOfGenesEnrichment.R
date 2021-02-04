
#' compute gene set enrichment - osama's code wrapped in package.
#' @export
#' @import WebGestaltR
#' @author Osama 
#' @param genes.with.values of genes and difference values
#' @param prot.univ the space of all proteins we are considering
#' @return gSEA output type stuff
computeGSEA<-function(genes.with.values,prefix,gsea_FDR=0.01){
  
  library(WebGestaltR)
  #library(ggplot2)
  inputdfforWebGestaltR <- genes.with.values%>%
    dplyr::rename(genes='Gene',scores='value')%>%
    dplyr::arrange(scores)
  
  
  #' * GSEA using gene ontology biological process gene sets
  
  go.bp.res.WebGestaltR <- WebGestaltR(enrichMethod = "GSEA", 
                                       organism="hsapiens", 
                                       enrichDatabase="geneontology_Biological_Process", 
                                       interestGene=inputdfforWebGestaltR, 
                                       interestGeneType="genesymbol", 
                                       collapseMethod="mean", perNum = 1000,
                                       fdrThr = gsea_FDR, nThreads = 2, isOutput = F)
  write.table(go.bp.res.WebGestaltR, paste0("proteomics_", prefix, "_gseaGO_result.txt"), sep="\t", row.names=FALSE, quote = F)

  all_gseaGO <- go.bp.res.WebGestaltR %>% 
     filter(FDR < gsea_FDR) %>% 
    dplyr::rename(pathway = description, NES = normalizedEnrichmentScore) %>% 
    arrange(NES) %>% 
    dplyr::mutate(status = case_when(NES > 0 ~ "Up",
                                     NES < 0 ~ "Down"),
                  status = factor(status, levels = c("Up", "Down"))) %>% 
    group_by(status) %>% 
    top_n(20, wt = abs(NES)) %>% 
    ungroup() %>% 
    ggplot2::ggplot(aes(x=reorder(pathway, NES), y=NES)) +
    ggplot2::geom_bar(stat='identity', aes(fill=status)) +
    ggplot2::scale_fill_manual(values = c("Up" = "darkred", "Down" = "dodgerblue4")) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y=element_text(size = 14),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") +
    ggplot2::labs(title = "", y="NES") +#for some reason labs still works with orientation before cord flip so set y
    ggplot2::ggtitle(paste('All',prefix))
  ggplot2::ggsave(paste0("allRegProts_", prefix,"_gseaGO_plot.pdf"), all_gseaGO, height = 8.5, width = 11, units = "in")
  
 
  return(go.bp.res.WebGestaltR) 
}



#' compute kinase substrate enrichment - osama's code wrapped in package.
#' @export
#' @import KSEAapp
#' @import readr
#' @import dplyr
#' @author Osama 
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
  KSDB <- read.csv(system.file('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',
                               package='amlresistancenetworks'),stringsAsFactors = FALSE)
  
  #' * KSEA using not only the known substrates in PSP but also the predicted substrates in NetworKIN
  res<-KSEA.Complete(KSDB, inputdfforKSEA, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5,
                     p.cutoff=ksea_FDR)
  file.rename("KSEA Bar Plot.tiff",paste0(prefix,'_KSEABarPlot.tiff'))
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
  res_ksea<-res_ksea%>%rename(`aveSubstrateLog2FC`='log2FC')%>%left_join(subs,by='Kinase.Gene')
  ggsave(paste0(prefix,"fig_KSEA.pdf"), plot_KSEA, height = 8.5, width = 11, units = "in")
  write.table(res_ksea,paste0(prefix,'_kseaRes.csv'),sep=',')
  return(res_ksea)
}



#' Old plot using clusterProfiler
#' @export 
#' @import BiocManager

plotOldGSEA<-function(genes.with.values,prefix,gsea_FDR=0.05){
  if(!require('org.Hs.eg.db')){
    BiocManager::install('org.Hs.eg.db')
    require(org.Hs.eg.db)
  }
  if(!require('clusterProfiler')){
    BiocManager::install('clusterProfiler')
    require('clusterProfiler')
  }

  genes.with.values<-arrange(genes.with.values,desc(value))
 # print(head(genes.with.values))
  genelist=genes.with.values$value
  names(genelist)=genes.with.values$Gene
  print(head(genelist))
  genelist<-sort(genelist,decreasing=TRUE)

  
  gr<-clusterProfiler::gseGO(unlist(genelist),ont="BP",keyType="SYMBOL",
                             OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH',pvalueCutoff = 0.5)#,eps=1e-10)
  #gr<-clusterProfiler::gseKEGG(genelist[!is.na(genelist)],organism='hsa',keyType="kegg",
  #OrgDb=org.Hs.eg.db,
  #                           pAdjustMethod = 'BH')#,eps=1e-10)
  
  # if(nrow(as.data.frame(gr))==0){
  #    gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="SYMBOL",
  #                             OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH',pvalueCutoff = 0.1)#,eps=1e-10)
  #  }
  
  res<-filter(as.data.frame(gr),p.adjust<gsea_FDR)
  if(nrow(res)==0)
    return(gr)
  
  all_gseaGO<-res %>% 
    dplyr::rename(pathway = 'Description') %>% 
    arrange(NES) %>% 
    dplyr::mutate(status = case_when(NES > 0 ~ "Up",
                                     NES < 0 ~ "Down"),
                  status = factor(status, levels = c("Up", "Down"))) %>% 
    group_by(status) %>% 
    top_n(20, wt = abs(NES)) %>% 
    ungroup() %>% 
    ggplot2::ggplot(aes(x=reorder(pathway, NES), y=NES)) +
    geom_bar(stat='identity', aes(fill=status)) +
    scale_fill_manual(values = c("Up" = "darkred", "Down" = "dodgerblue4")) +
    coord_flip() +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y=element_text(size = 14),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") +
    labs(title = "", y="NES") +#for some reason labs still works with orientation before cord flip so set y
    ggtitle(paste('All',prefix))
  ggsave(paste0("allRegProts_", prefix,"_gseaGO_plot.pdf"), all_gseaGO, height = 8.5, width = 11, units = "in")
  
  enrichplot::ridgeplot(gr,showCategory = 50,fill='pvalue')+ggplot2::ggtitle(paste0("GO Terms for ",prefix))
  ggplot2::ggsave(paste0(prefix,'_GO.pdf'),width=10,height=10)
  
  df<-as.data.frame(gr)%>%mutate(Condition=prefix)
  return(df)
}

#'Runs regular bag of genes enrichment
#'@export
#'@import BiocManager
doRegularGo<-function(genes,bg=NULL){
  if(!require('org.Hs.eg.db')){
    BiocManager::install('org.Hs.eg.db')
    require(org.Hs.eg.db)
  }
  if(!require('clusterProfiler')){
    BiocManager::install('clusterProfiler')
    require('clusterProfiler')
  }


  #genes<-unique(as.character(genes.df$Gene))
  mapping<-as.data.frame(org.Hs.egALIAS2EG)%>%
    dplyr::rename(Gene='alias_symbol')
  
  eg<-subset(mapping,Gene%in%genes)
  ret=data.frame(ID='',Description='',pvalue=1.0,p.adjust=1.0)
  
  try(res<-clusterProfiler::enrichGO(eg$gene_id,'org.Hs.eg.db',keyType='ENTREZID',ont='BP'))
    #sprint(res)
  ret<-as.data.frame(list(ID=NULL,Description=NULL,pvalue=NULL,p.adjust=NULL))
  try(ret<-as.data.frame(res)%>%dplyr::select(ID,Description,pvalue,p.adjust))
  return(ret)
  
  
}

#'Runs regular bag of genes enrichment
#'@export 
#'@import devtools

doRegularKin<-function(genes,bg=NULL){
  require(dplyr)
  if(!require('leapr')){
    devtools::install('biodataganache/leapr')
    require('leapr')
      }

   #read kinase substrate database stored in data folder
  KSDB <- read.csv(system.file('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',
                               package='amlresistancenetworks'),stringsAsFactors = FALSE)

  kdat<-KSDB%>%group_by(GENE)%>%select(SUB_GENE,SUB_MOD_RSD)%>%
    rowwise()%>%
    mutate(subval=paste(SUB_GENE,SUB_MOD_RSD,sep='-'))
  
  klist<-lapply(unique(kdat$GENE),function(x) unique(unlist(kdat[which(kdat$GENE==x),'subval'])))
  names(klist)<-unique(kdat$GENE)

  maxlen <- max(lengths(klist))
  kmat <- do.call(cbind,lapply(klist, function(lst) c(lst, rep(NA, maxlen - length(lst)))))
                  
  kslist<-list(names=names(klist),desc=paste(names(klist),'substrates'),
               sizes=unlist(lapply(klist,length)),matrix=t(kmat))


  ##now we need to fix the gene list, since it's not going to match    
  sgenes<-data.frame(genes=genes)%>%
    separate(genes, into=c('gene','mod'),sep='-')%>%
    mutate(modlist=strsplit(mod,split='s|t|y'))%>%
    apply(1,function(x) paste(x$gene,x$modlist,sep='-'))%>%
    unlist()%>%unique()
  
  
  print(paste("Found",length(sgenes),'substrates with known kinases'))
  
  ret<-as.data.frame(list(ID=NULL,Description=NULL,pvalue=NULL,p.adjust=NULL))

  #ret<-as.data.frame(list(Kinase=NULL,NumSubs=NULL,pvalue=NULL,p.adjust=NULL))
  if(length(sgenes)<2)
    return(ret)

  try(res <- leapR(geneset=kslist,
              enrichment_method='enrichment_in_sets',targets=sgenes))

  #print(head(res))
  #print(res)
  
  ret<-as.data.frame(res)%>%
    tibble::rownames_to_column('Kinase')%>%
    subset(ingroup_n>0)%>%
    dplyr::select(ID='Kinase',NumSubs='ingroup_n',pvalue,p.adjust='BH_pvalue')%>%rowwise()%>%
    mutate(Description=paste0(ID,': ',NumSubs,' targets'))%>%
    select(-NumSubs)
  
  
  return(ret)
  
}
  

