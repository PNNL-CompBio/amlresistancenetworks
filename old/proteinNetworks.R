##move network-related analyses here



#' @import PCSF
#' @export
#' 
computeProteinNetwork<-function(sig.vals,all.vals,phos.vals=NULL,nrand=100,beta=500/nrow(sig.vals)){
  require(PCSF)
  
  data("STRING")
  
  ppi <- construct_interactome(STRING)
  
  
  terms<-sig.vals$value
  names(terms)<-sig.vals$Gene
  print(paste('Building subnetwork with',length(terms),'terminals'))
  
  subnet <- PCSF_rand(ppi,abs(terms), n=nrand, r=0.2,w = 4, b = beta, mu = 0.0005)
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
