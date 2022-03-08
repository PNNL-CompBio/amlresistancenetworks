##move network-related analyses here



#' computePhosphoNetwork
#' maps kinase activity to proteins in a sum
#' @param phos.vals named list of phosphosite values
#' @param prot.vals named list of protein values
#' @param gene.vals named list of gene mutations
#' @param rna.vals named list of transcripts
#' @param nrand number of randomizations
#' @param beta beta value
#' @param fname
#' @import remotes
#' @export
computePhosphoNetwork<-function(phos.vals=c(),prot.vals=c(),
                                gene.vals=c(),rna.vals=c(),
                                nrand=100,beta=2,fname){
  
  require(dplyr)
    #read kinase substrate database stored in data folder
  KSDB <- read.csv(system.file('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',
                               package='amlresistancenetworks'),stringsAsFactors = FALSE)
  
  kdat<-KSDB%>%group_by(GENE)%>%select(SUB_GENE,SUB_MOD_RSD)%>%
    rowwise()%>%
    mutate(subval=paste(SUB_GENE,SUB_MOD_RSD,sep='-'))
  
  allvals<-unique(kdat$subval)
  
  if(!require('PCSF')){
    remotes::install_github('sgosline/PCSF')
    require('PCSF')
    }
   data("STRING")
   #print(phos.vals)
   ##add phospho data to STRING
   if(length(phos.vals)>0){
     mval<-mean(STRING$cost)
     adf<-apply(kdat,1,function(x)
       #for each substrate interaction, add a link from the kinase gene -> substreate -> substrate gene
       data.frame(from=c(x[['GENE']],x[['subval']]),to=c(x[['subval']],x[['SUB_GENE']]),
        cost=c(mval*1.5,mval/4)))%>%  ##arbitrary costs based on mean cost of edges around network
       do.call(rbind,.)
   }else{
     adf<-data.frame()
   } 
   
   ##now we can do the same for mRNA prot interactions
   enrichr <- read.table(system.file('gene_attribute_edges.txt',package='amlresistancenetworks'),
                                     header=T)[-1,]
   if(length(rna.vals)>0){
     names(rna.vals)<-paste0(names(rna.vals),'_RNA')
     mval<-mean(STRING$cost)
     mdf = enrichr%>%
       rename(from='target')%>%
       select(from,source)%>%
       rowwise()%>%
       mutate(to=paste0(source,'_RNA'),cost=mval/4)%>%
       select(-source)
     #mdf<-apply(enrichr,1,function(x)
    #   #for each substrate interaction, add a link from the kinase gene -> substreate -> substrate gene
     #  data.frame(from=c(x[['source']]),to=c(x[['target']]),
    #              cost=c(mval*1.5,mval/4)))%>%  ##arbitrary costs based on mean cost of edges around network
    #   do.call(rbind,.)
     
   }else{
     mdf <- data.frame()
   }
   
 #  print(adf)
   ppi <- construct_interactome(rbind(STRING,adf,mdf))
   
   ##now run the code
   terms=c(phos.vals,prot.vals,gene.vals,rna.vals)
   subnet<-NULL
   try(
     subnet <- PCSF_rand(ppi,abs(terms), n=nrand, r=0.2,w = 4, b = beta, mu = 0.0005)
   )
   
   if(is.null(subnet))
     return("")
  
  lfcs<-terms[match(names(V(subnet)),names(terms))]
  lfcs[is.na(lfcs)]<-0.0
  
  ##assign proteins first
  types<-rep('proteins',length(names(V(subnet))))
  
  names(types)<-names(V(subnet))
  
  ##then assign phosphosites
  types[intersect(names(V(subnet)),allvals)]<-'phosphosite'
  
  ##assign gene mutations
  types[intersect(names(V(subnet)),names(gene.vals))]<-'mutation'

  ##assign mRNA
  mrna_nodes =intersect(names(V(subnet)),names(rna.vals))
  types[mrna_nodes]<-'mRNA'
   
  ##in case there are overlapping hits
  types[intersect(names(V(subnet)),intersect(names(gene.vals),names(prot.vals)))] <- 'mutation+proteomic'
 
  #then rename
  subnet <- igraph::set.vertex.attribute(subnet, 
                                         "name", 
                                         value=gsub('_RNA','',names(V(subnet))))
  
  subnet<-igraph::set.vertex.attribute(subnet,'logFoldChange',value=lfcs)
  subnet<-igraph::set.vertex.attribute(subnet,'nodeType',value=types)
  subnet<-igraph::set.edge.attribute(subnet,'interactionType',value='protein-protein interaction')

  
  write_graph(subnet,format='gml',file=paste0(fname,'.gml'))
  return(list(graph=subnet,fname=paste0(fname,'.gml')))
  #subnet
}

#' computeProteinNetwork
#' @import remotes
#' @param data.frame all.vals with required values, columsn should be: `Gene`,`value`,`signif`, and `cond`
#' @param nrand number of randomizations
#' @return file name of graph in gml format
#' @export
#' 
computeProteinNetwork<-function(all.vals,nrand=100){
  if(!require('PCSF')){
    remotes::install_github('sgosline/PCSF')
    require('PCSF')
  }
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

  
  write_graph(subnet,format='gml',file=paste0(condname,'.gml'))
  return(paste0(condname,'.gml'))
}
