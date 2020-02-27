##proteomic stats


#' compute protein FC
#' @export
#' @param tidied data frame
#' @param condition of interest
#' @return data frame with 3 columns - genes with LFC values and significance values
computeProteinFC<-function(tidied.df){
  
}

#' @require PCSF
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
#' @requires clusterProfiler
#' @param output of genes and difference values
#' @return gSEA output type stuff
computeGSEA<-function(genes.with.values){
  require(clusterProfiler)
}

