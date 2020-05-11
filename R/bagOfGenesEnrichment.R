
#' compute gene set enrichment - osama's code wrapped in package.
#' @export
#' @import WebGestaltR
#' @import ggplot2
#' @author Osama 
#' @param genes.with.values of genes and difference values
#' @param prot.univ the space of all proteins we are considering
#' @return gSEA output type stuff
computeGSEA<-function(genes.with.values,prefix,gsea_FDR=0.01){
  
  library(WebGestaltR)
  library(ggplot2)
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
  
  top_gseaGO <- go.bp.res.WebGestaltR %>% 
    filter(FDR < gsea_FDR) %>% 
    dplyr::rename(pathway = description, NES = normalizedEnrichmentScore) %>% 
    arrange(desc(NES)) %>% 
    dplyr::mutate(status = case_when(NES > 0 ~ "Up",
                                     NES < 0 ~ "Down"),
                  status = factor(status, levels = c("Up", "Down"))) %>% 
    #\group_by(status) %>% 
    top_n(30, wt = NES) %>% 
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
    ggtitle(paste('Up-regulated',prefix))
  ggsave(paste0("upRegProts_", prefix,"_gseaGO_plot.pdf"), top_gseaGO, height = 8.5, width = 11, units = "in")
  
  
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
  
  
  bot_gseaGO <- go.bp.res.WebGestaltR %>% 
    filter(FDR < gsea_FDR) %>% 
    dplyr::rename(pathway = description, NES = normalizedEnrichmentScore) %>% 
    arrange(NES) %>% 
    dplyr::mutate(status = case_when(NES > 0 ~ "Up",
                                     NES < 0 ~ "Down"),
                  status = factor(status, levels = c("Up", "Down"))) %>% 
    #group_by(status) %>% 
    top_n(40, wt = rev(NES)) %>% 
    ungroup() %>% 
    ggplot2::ggplot(aes(x=reorder(pathway, rev(NES)), y=NES)) +
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
    ggtitle(paste('Down-regulated',prefix))
  ggsave(paste0("downRegProts_", prefix,"_gseaGO_plot.pdf"), bot_gseaGO, height = 8.5, width = 11, units = "in")
  
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
  KSDB <- read.csv(system.file('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',package='amlresistancenetworks'),stringsAsFactors = FALSE)
  
  #' * KSEA using not only the known substrates in PSP but also the predicted substrates in NetworKIN
  res<-KSEA.Complete(KSDB, inputdfforKSEA, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=ksea_FDR)
  file.remove("KSEA Bar Plot.tiff")
  file.remove("Kinase-Substrate Links.csv")
  #make own plot of KSEA results
  res_ksea <- readr::read_csv("KSEA Kinase Scores.csv")
  
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
  file.remove("KSEA Kinase Scores.csv")
  ggsave(paste0(prefix,"fig_KSEA.pdf"), plot_KSEA, height = 8.5, width = 11, units = "in")
  
  return(res_ksea)
}



#' Old plot using clusterProfiler
#' @export 
#' @require org.Hs.eg.db
#' @import clusterProfiler
plotOldGSEA<-function(genes.with.values,prot.univ,prefix){
  require(org.Hs.eg.db)
  mapping<-as.data.frame(org.Hs.egALIAS2EG)%>%
    dplyr::rename(Gene='alias_symbol')
  
  genes.with.values<-genes.with.values%>%
    dplyr::left_join(mapping,by='Gene')%>%
    arrange(desc(value))
  
  genelist=genes.with.values$value
  names(genelist)=genes.with.values$gene_id
  
  # symbs<-names(genelist)[!is.na(genelist)]
  # xx <- as.list(org.Hs.egALIAS2EG)
  # ents<-unlist(sapply(intersect(names(xx),symbs), function(x) xx[[x]]))
  # print(ents)
  
  gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="SYMBOL",
                             OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH')#,eps=1e-10)
  #gr<-clusterProfiler::gseKEGG(genelist[!is.na(genelist)],organism='hsa',keyType="kegg",
  #OrgDb=org.Hs.eg.db,
  #                           pAdjustMethod = 'BH')#,eps=1e-10)
  
  # if(nrow(as.data.frame(gr))==0){
  #    gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="SYMBOL",
  #                             OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH',pvalueCutoff = 0.1)#,eps=1e-10)
  #  }
  
  enrichplot::ridgeplot(gr,showCategory = 50,fill='pvalue')+ggplot2::ggtitle(paste0("KEGG Terms for ",prefix))
  ggplot2::ggsave(paste0(prefix,'_KEGG.pdf'),width=10,height=10)
  
  
  return(gr)
}

#'Runs regular bag of
#'@export 
#'@require org.Hs.eg.db
#'@import clusterProfiler
doRegularGo<-function(genes,bg=NULL){
  require(org.Hs.eg.db)
  #genes<-unique(as.character(genes.df$Gene))
  mapping<-as.data.frame(org.Hs.egALIAS2EG)%>%
    dplyr::rename(Gene='alias_symbol')
  
  eg<-subset(mapping,Gene%in%genes)
  
  res<-clusterProfiler::enrichGO(eg$gene_id,'org.Hs.eg.db',keyType='ENTREZID',ont='BP')
    #sprint(res)
  ret=as.data.frame(res)%>%
    dplyr::select(ID,Description,pvalue,p.adjust)
  return(ret)
  
  
}
