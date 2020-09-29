##move network-related analyses here


#' @import PCSF
#' @export
#' 
computeProteinNetwork<-function(all.vals,nrand=100){
  require(PCSF)
  
  data("STRING")
  print(head(all.vals))
  ppi <- construct_interactome(STRING)
  
  condname=unique(all.vals$cond)
  print(condname)
  gene.vals<-dplyr::select(all.vals,Gene,value,signif)%>%distinct()
  
  sig.vals<-subset(gene.vals,signif==TRUE)
  beta=500/nrow(sig.vals)
  terms<-sig.vals$value
  names(terms)<-sig.vals$Gene
  
  
  print(paste('Building subnetwork with',length(terms),'terminals'))
  
  subnet<-NULL
  try(
    subnet <- PCSF_rand(ppi,abs(terms), n=nrand, r=0.2,w = 4, b = beta, mu = 0.0005)
  )
  
  if(is.null(subnet))
    return("")
  #now add back in values from terminals as attributes
  
  #if(!is.null(all.vals)){
  lfcs<-gene.vals$value[match(names(V(subnet)),gene.vals$Gene)]
  lfcs[is.na(lfcs)]<-0.0
  
  print(lfcs)
  #sogs<-all.vals$signif[match(names(V(subnet)),all.vals$signif)]
  #sogs[is.na(sogs)]<-FALSE
  #print(sogs)
  subnet<-igraph::set.vertex.attribute(subnet,'logFoldChange',value=lfcs)
  subnet<-igraph::set.edge.attribute(subnet,'interactionType',value='protein-protein interaction')
  
  if('Substrate.Gene'%in%colnames(all.vals)){
  ksi <-all.vals%>%
    dplyr::select(Gene,Substrate.Gene,Source,log2FC)%>%distinct()
  edgelist<-c()
  ksi$Substrate.Gene<-as.character(ksi$Substrate.Gene)
  ksi<-subset(ksi,Gene%in%names(V(subnet)))
  for(i in 1:nrow(ksi))
    edgelist<-c(edgelist,ksi[i,c('Gene','Substrate.Gene')])
  newnodes<-setdiff(unlist(edgelist),names(V(subnet)))
  subnet<-igraph::add_vertices(subnet,nv=length(newnodes),name=unlist(newnodes),
                               type='Substrate',logFoldChange=0.0,prize=0)
  subnet <-igraph::add_edges(subnet,unlist(edgelist),attr=list(interactionType='kinase-substrate interaction',
                                                      source=as.character(ksi$Source),
                                                      weight=abs(as.numeric(ksi$log2FC))))
  }
 # subset<-igraph::set.vertex.attribute(subnet,'Significant',value=sogs)
  
  # if(!is.null(phos.vals)){
  #   lpc<-phos.vals$value[match(names(V(subnet)),all.vals$Gene)]
  #   lpc[is.na(lpc)]<-0.0
  #   subnet<-igraph::set.vertex.attribute(subnet,'phosphoLogFoldChange',index=V(subnet),value=lpc)
  #   
  # }
  
  write_graph(subnet,format='gml',file=paste0(condname,'.gml'))
  return(paste0(condname,'.gml'))
}
