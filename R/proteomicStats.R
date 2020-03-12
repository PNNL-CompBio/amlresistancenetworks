##proteomic stats


#' @import PCSF
#' @export
#' 
computeProteinNetwork<-function(sig.vals,all.vals,nrand=100){
  require(PCSF)

  data("STRING")
  
  ppi <- construct_interactome(STRING)
  
  
  terms<-sig.vals$value
 names(terms)<-sig.vals$Gene
  print(paste('Building subnetwork with',length(terms),'terminals'))
  
  subnet <- PCSF_rand(ppi,abs(terms), n=nrand, r=0.3,w = 4, b = 50, mu = 0.0005)
  #now add back in values from terminals as attributes
  
  lfcs<-all.vals$value[match(names(V(subnet)),all.vals$Gene)]
  lfcs[is.na(lfcs)]<-0.0
  subnet<-igraph::set.vertex.attribute(subnet,'logFoldChange',index=V(subnet),value=lfcs)
  subnet  
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
                                 conditions=c("FL","FGF2")){
  
  data<-g.data%>%
    dplyr::select(CellLine,Gene,value,treatment)%>%
    subset(!is.na(value))%>%
    subset(!is.na(Gene))
  
  
  #remove those data that have only 1 replicate  
  reps<-data%>%
    group_by(CellLine,Gene,treatment)%>%
    dplyr::mutate(reps=length(value))%>%
    subset(reps>1)%>%
    ungroup()        
  
  ##this function iterates through every ligand of interest and compares it to 'None'
  lig.signif<-purrr::map_df(conditions,function(lig){
    #compute the significance for FL-modulaetd changes in both cell lines
    reps%>%
      subset(treatment%in%c(control,lig))%>% #get only those valuese that are 'nOne' or the name of the ligand
      dplyr::group_by(Gene)%>% #group by cell line and gene to get per-gene values
      dplyr::mutate(treatment=as.factor(treatment))%>% #change to factor for t-test
      dplyr::mutate(conds=n_distinct(treatment))%>% subset(conds>1)%>% #only get those for which there are sufficient reps
      dplyr::group_modify(~infer::t_test(.,value~treatment,order=c(lig,control)),keep=TRUE)%>% #calculate p-value
      dplyr::select(c(Gene,p_value))%>%  #select the valuese of interest for us
      ungroup()%>%
    #  dplyr::group_by(CellLine)%>% #regroup so that we can correct the p-values
      dplyr::mutate(p_adj=p.adjust(p_value),Condition=lig,Control=control)%>% #adjust
      ungroup()
  })
  
  lig.fc<-purrr::map_df(conditions,function(lig){
    reps%>%
      subset(treatment%in%c(control,lig))%>% #get only those valuese that are 'nOne' or the name of the ligand
      dplyr::mutate(treatment=stringr::str_replace(treatment,lig,'condition'))%>% #rename so we can spread
      dplyr::mutate(treatment=stringr::str_replace(treatment,control,'conVal'))%>% #rename so we can spread
      ungroup()%>%dplyr::group_by(Gene,treatment)%>%
      dplyr::summarize(meanVal=mean(value))%>% #calculate mean across cell line/gene/treatment
      tidyr::pivot_wider(names_from=treatment,values_from = meanVal)%>% #spread to compute difference
      dplyr::mutate(condition_to_control=condition-conVal,Condition=lig,Control=control)%>% #compute difference
      dplyr::select(Gene,Condition,Control,condition_to_control)
  })
  lig.datasets<-
    lig.signif%>%full_join(lig.fc,by=c('Gene','Condition','Control'))
  
  return(lig.datasets) 
}


#' compute gene set enrichment
#' @export
#' @import clusterProfiler
#' @importFrom msigdbr msigdbr
#' @import org.Hs.eg.db
#' @param genes.with.values of genes and difference values
#' @param prot.univ the space of all proteins we are considering
#' @return gSEA output type stuff
computeGSEA<-function(genes.with.values,prot.univ,prefix){

    require(org.Hs.eg.db)
    genes.with.values<-genes.with.values%>%arrange(desc(value))
    
  genelist=genes.with.values$value
   names(genelist)=genes.with.values$Gene
  
  gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH')#,eps=1e-10)
  if(nrow(as.data.frame(gr))==0){
    gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH',pvalueCutoff = 0.1)#,eps=1e-10)
  }
  
  enrichplot::ridgeplot(gr,showCategory = 50,fill='pvalue')+ggplot2::ggtitle(paste0("GO Terms for ",prefix))
  ggplot2::ggsave(paste0(prefix,'_GO.png'),width=16,height=16)
  
  
  return(gr)
  }

