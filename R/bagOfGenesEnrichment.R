
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
#' @import ggplot2
#' @import gridExtra
#' @import scales
#' @author Osama 
#' @param genes.with.values data frame containing a Gene column, as well as a Residue.Both column and FC column (see KSEA app).
#' @param prot.univ the space of all proteins we are considering
#' @return KSEA output type stuff
computeKSEA<-function(genes.with.values,ksea_FDR=0.05,prefix='', order_by = "z.score",
                      height = 8.5, width = 11){
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
  
  # Selecting only those Kinases with enough linked subtrates (m>5), as well as small enough false dicovery rate (p.adjust < KSEA_fdr)
  res_ksea <- res_ksea %>%
    mutate(p.adjust = p.adjust(p.value)) %>% 
    filter(m >= 5) %>% 
    arrange(desc(z.score)) %>% 
    mutate(status = case_when(z.score > 0 & p.adjust <= ksea_FDR ~ "Up",
                              z.score < 0 & p.adjust <= ksea_FDR ~ "Down",
                              TRUE ~ "Not significant")) %>% 
    filter(status != "Not significant")
  
  plot_KSEA <- res_ksea %>%
    ggplot(aes(x = z.score,y = reorder(Kinase.Gene, get(order_by)))) +
    geom_bar(stat='identity', aes(fill=status)) +
    scale_fill_manual(values = c("Down" = "dodgerblue3", "Up" = "firebrick2", "Not significant" = "black")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          legend.position="none", 
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y=element_text(size = 14),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Kinase z-score") + 
    ggtitle("Kinase z-score")
  
  plot_sig <- res_ksea %>%
    ggplot(aes(x = p.adjust, y = reorder(Kinase.Gene, get(order_by)))) +
    geom_bar(stat='identity') +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18), 
          legend.position="none", 
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.ticks.y = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.y = element_blank()) +
    scale_x_continuous(trans = reverselog_trans(10)) +
    labs(x = "Adjusted p-value") +
    ggtitle("Significance")
  
  arrange_matrix <- t(as.matrix(c(1,1,2)))
  plot_both <- grid.arrange(plot_KSEA, plot_sig, layout_matrix = arrange_matrix)
  
  ggsave(paste0("sig-included", prefix,"-ksea-plot.png"), plot_both, 
         height = height, width = width, units = "in")
  file.remove("KSEA Kinase Scores.csv")#paste0(prefix,'kinaseScores.csv'))
  
  ##join results
  #kin_res
  res_ksea<-res_ksea%>%rename(`aveSubstrateLog2FC`='log2FC')%>%left_join(subs,by='Kinase.Gene')
  ggsave(paste0(prefix,"fig_KSEA.pdf"), plot_KSEA, height = height, width = width, units = "in")
  write.table(res_ksea,paste0(prefix,'_kseaRes.csv'),sep=',')
  return(res_ksea)
}



#' Plot using clusterProfiler. 3 GSEA plots are saved to the working directory, 
#' two of which are custom to the amlresistancenetworks package.
#' @export 
#' @import BiocManager
#' @import ggplot2
#' @import gridExtra
#' @import scales
#' @import dplyr
#' @param genes.with.values A data frame with genes as row names, along with a column named "value". Usually this value column consists of log fold changes between two groups.
#' @param prefix string, used for naming the saved plots.
#' @param order.by This determines how the GO terms are sorted. Default is normalized enrichment score "NES", but can also use "p.adjust" to sort by significance of the terms.
plotOldGSEA<-function(genes.with.values, prefix, gsea_FDR=0.05, 
                      order.by = "NES", height = 8.5, width = 11){
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
  if(nrow(res)==0){
    return(gr)
  }
  
  all.gseaGO<-res %>% 
    dplyr::rename(pathway = 'Description') %>% 
    arrange(NES) %>% 
    dplyr::mutate(status = case_when(NES > 0 ~ "Up", NES < 0 ~ "Down"),
                  status = factor(status, levels = c("Up", "Down"))) %>% 
    group_by(status) %>% 
    top_n(20, wt = abs(NES)) %>% 
    ungroup()
  
  p.NES <- ggplot(all.gseaGO, aes(x = NES, y = reorder(pathway, get(order.by)))) +
    geom_bar(stat='identity', aes(fill=status)) +
    scale_fill_manual(values = c(Up = "firebrick2", Down = "dodgerblue3")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") + 
    labs(x = "NES") +
    ggtitle("Normalized Enrichment Score")
  
  ggsave(paste0("allRegProts_", prefix,"_gseaGO_plot.pdf"), p.NES, 
         height = 8.5, width = 11, units = "in")
  
  p.Pval <- ggplot(all.gseaGO, aes(x = p.adjust, y = reorder(pathway, get(order.by)))) +
    scale_x_continuous(trans = reverselog_trans(10)) +  
    theme_minimal() +
    geom_bar(stat = "identity") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18), 
          axis.title.x = element_text(size = 16), 
          axis.text.x = element_text(size = 14), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.line.y = element_line(color = "black"),
          axis.ticks.y = element_blank(), 
          legend.position = "none") + 
    labs(x = "Adjusted p-value") + 
    ggtitle("Significance")
  
  arrange_matrix <- t(as.matrix(c(1,1,1,2)))
  p.both <- grid.arrange(p.NES,p.Pval, layout_matrix = arrange_matrix)
  
  ggsave(paste0("sig-included", prefix,"-gseaGO-plot.png"), p.both, 
         height = height, width = width, units = "in")
  
  enrichplot::ridgeplot(gr,showCategory = 50,fill='pvalue')+ggplot2::ggtitle(paste0("GO Terms for ",prefix))
  ggplot2::ggsave(paste0(prefix,'_GO.pdf'),width=10,height=10)
  
  df<-as.data.frame(gr)%>%mutate(Condition=prefix)
  return(df)
}

#' Used to make reversed logarithmic scales
#' @import scales
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
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
  #ret<-as.data.frame(list(ID=NULL,Description=NULL,pvalue=NULL,p.adjust=NULL))
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
  
  #ret<-as.data.frame(list(ID=NULL,Description=NULL,pvalue=NULL,p.adjust=NULL))
  ret=data.frame(ID='',Description='',pvalue=1.0,p.adjust=1.0)
  
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
  

#' Plot using correlation enrichment from leapR package. 
#' A single plot is saved to the working directory
#' @export 
#' @import ggplot2
#' @import gridExtra
#' @import scales
#' @import dplyr
#' @import leapr
#' @param exprs A matrix of intensities with accessions as row names, along with samples in the columns.
#' @param prefix string, used for naming the saved plots.
#' @param order.by This determines how the pathways are sorted. Default is pathway correlation of "Ingroup mean", but can also use "BH_pvalue" to sort by significance of the pathways.
#' @param geneset Pathway/Kinase database, eg ncipid, msigdb, both of which are included in leapr.
#' @param clean.names Boolean, if TRUE removes the "_pathway" ending in pathway names, making the plot easier to read.
plotCorrelationEnrichment <- function(exprs, geneset, fdr.cutoff = 0.05, 
                                      corr.cutoff = 0.1, prefix, width = 11, 
                                      height = 8.5, order.by = "Ingroup mean", 
                                      clean.names = FALSE, ...) {
  
  corr.enrichment <- leapR(geneset, 
                           enrichment_method = "correlation_enrichment",
                           datamatrix = exprs) 
  corr.enrichment <- corr.enrichment %>%
    mutate(Pathway = rownames(.)) %>%
    rename(`Ingroup mean` = ingroup_mean,
           `Outgroup mean` = outgroup_mean) %>%
    mutate(Status = case_when(`Ingroup mean` > 0 ~ "Positively Correlated", 
                              `Ingroup mean` < 0 ~ "Negatively Correlated")) %>%
    select(Pathway, `Ingroup mean`, `Outgroup mean`, 
           ingroup_n, outgroup_n, pvalue, BH_pvalue, Status)
  
  corr.enrichment.filtered <- corr.enrichment %>%
    filter(BH_pvalue < fdr.cutoff & abs(`Ingroup mean`) > corr.cutoff) %>%
    mutate(BH_pvalue = case_when(BH_pvalue > 1e-10 ~ BH_pvalue,
                                 BH_pvalue < 1e-10 ~ 1e-10))
  
  if (clean.names) {
    corr.enrichment.filtered$Pathway <- sub("_pathway$", "", 
                                            corr.enrichment.filtered$Pathway)
  }
  
  p.corr <- ggplot(corr.enrichment.filtered, aes(x = `Ingroup mean`, 
                                                 y = reorder(Pathway, get(order.by)))) +
    geom_bar(stat='identity', aes(fill = Status)) +
    scale_fill_manual(values = c("Positively Correlated" = "dodgerblue3", 
                                 "Negatively Correlated" = "firebrick2")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") + 
    labs(x = "Average correlation") +
    ggtitle("Correlation Enrichment")
  
  p.pval <- ggplot(corr.enrichment.filtered, aes(x = BH_pvalue, 
                                                 y = reorder(Pathway, get(order.by)))) +
    geom_bar(stat='identity') +
    scale_x_continuous(trans = reverselog_trans(10)) + 
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_line(color = "black"),
          legend.position = "none") + 
    labs(x = "Adjusted p-value") +
    ggtitle("Significance")
  
  arrange_matrix <- t(as.matrix(c(1,1,1,2)))
  p.both <- grid.arrange(p.corr,p.pval, layout_matrix = arrange_matrix)
  
  ggsave(paste0("sig-included-", prefix,"-correlation-enrichment-plot.png"), p.both, 
         height = height, width = width, units = "in")
  return(corr.enrichment)
}





































