##proteomic stats


#' takes multiple lists of proteins and creates networks
#' then merges them to find a single community
#' @import multinet
#' @param protLists lists of proteins by condition
#' @return shared graph
createCommunityGraph<-function(protLists,nrand,beta){
  
}


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



#' @import PCSF
#' @export
#' 
computeProteinNetwork<-function(sig.vals,all.vals,phos.vals=NULL,nrand=100,beta=1000/nrow(sig.vals)){
  require(PCSF)

  data("STRING")
  
  ppi <- construct_interactome(STRING)
  
  
  terms<-sig.vals$value
 names(terms)<-sig.vals$Gene
  print(paste('Building subnetwork with',length(terms),'terminals'))
  
  subnet <- PCSF_rand(ppi,abs(terms), n=nrand, r=0.3,w = 4, b = beta, mu = 0.0005)
  #now add back in values from terminals as attributes
  
  lfcs<-all.vals$value[match(names(V(subnet)),all.vals$Gene)]
  lfcs[is.na(lfcs)]<-0.0
  subnet<-igraph::set.vertex.attribute(subnet,'logFoldChange',index=V(subnet),value=lfcs)
  
  if(!is.null(phos.vals)){
    lpc<-phos.vals$value[match(names(V(subnet)),all.vals$Gene)]
    lpc[is.na(lpc)]<-0.0
    subnet<-igraph::set.vertex.attribute(subnet,'phosphoLogFoldChange',index=V(subnet),value=lpc)
    
  }
  subnet  
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
                                 control='None',
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
    g1<-subset(reps,treatment==control)%>%select(Sample)%>%distinct()
    g2<-subset(reps,treatment==lig)%>%select(Sample)%>%distinct()
    matp=reps%>%
      subset(treatment%in%c(control,lig))%>% #get only those valuese that are 'nOne' or the name of the ligand
      dplyr::select(Gene,Sample,value)%>%distinct()%>%
      tidyr::pivot_wider(values_from=value,names_from=Sample,values_fn=list(value=mean))%>%
      tibble::column_to_rownames("Gene")
    limmaTwoFactorDEAnalysis(matp,unlist(g1),unlist(g2))%>%
      dplyr::mutate(Condition=lig,Control=control)%>%
      dplyr::select(Gene=featureID,Condition,Control,condition_to_control=logFC,p_val=P.Value,p_adj=adj.P.Val)
    })
 
  return(lig.datasets) 
}


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
#' 
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

