##proteomic stats


#' compute protein FC
#' @export
#' @param tidied data frame
#' @param condition of interest
#' @return data frame with 3 columns - genes with LFC values and significance values
computeProteinFC<-function(tidied.df){
  
}

#' @import PCSF
#' @export
computeProteinNetwork<-function(tidied.df){
  require(PCSF)
}


#' compute differences between gene lists
#' @export
#' @param protein FC data frame
#' @return genes with distance values
computeFCDistances<-function(proteins.with.lfc){
  
}


#' compute differences between networks
#' @export
#' @param network 1
#' @param network 2
#' @return genes and list of distance values
computeNetworkDifferences<-function(network1,network2){
  
}

#' compute gene set enrichment
#' @export
#' @importFrom clusterProfiler GSEA
#' @importFrom msigdbr msigdbr
#' @import org.Hs.eg.db
#' @param genes.with.values of genes and difference values
#' @param prot.univ the space of all proteins we are considering
#' @return gSEA output type stuff
computeGSEA<-function(genes.with.values,prot.univ){

    require(org.Hs.eg.db)
  genes.with.values<-genes.with.values%>%
      rename(Gene="SYMBOL")%>%
    left_join(AnnotationDbi::select(org.Hs.eg.db,keys=keys(org.Hs.eg.db,'ENTREZID'),
                                    columns=c("ENTREZID","SYMBOL"),by='SYMBOL'))
  
    genes.with.values<-genes.with.values%>%arrange(desc(value))
    
  genelist=genes.with.values$value
 # names(genelist)=genes.with.values$SYMBOL
   names(genelist)=genes.with.values$ENTREZID
  
 # t2g= msigdbr::msigdbr(species="Homo sapiens",category="C2")%>%dplyr::select(gs_id,human_gene_symbol)
    t2g= msigdbr::msigdbr(species="Homo sapiens",category="C5")%>%
      subset(human_gene_symbol%in%prot.univ)%>%
               dplyr::select(gs_name,entrez_gene)
  #  t2n =msigdbr::msigdbr(species="Homo sapiens",category="C2")%>%dplyr::select(gs_id,gs_name)
  
  gr<-clusterProfiler::GSEA(genelist[!is.na(genelist)],TERM2GENE = t2g,pAdjustMethod = 'fdr')
  if(nrow(as.data.frame(gr))==0){
    gr<-clusterProfiler::GSEA(genelist[!is.na(genelist)],TERM2GENE = t2g,pAdjustMethod = 'fdr',pvalueCutoff = 0.1)
  }
    return(gr)
  }

