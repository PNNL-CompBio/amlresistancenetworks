##run through drug resistnace analysis


library(amlresistancenetworks)
library(dplyr)

##first run gilteritinib data
gilt.data<-readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(gilt.data$Gene)


early.data<-gilt.data%>%
  subset(treatment%in%(c('None','Early Gilteritinib')))%>%
  dplyr::select(Gene,Sample,cellLine,treatment,value,ligand)

late.data<-gilt.data%>%
  subset(treatment%in%(c('Early Gilteritinib','Late Gilteritinib')))%>%
  dplyr::select(Gene,Sample,cellLine,treatment,value,ligand)



#'this is our main analysis function  
#'it computes differential expression betwee control and various conditions
#'runs GSEA, plots, and creates a network
plotDataByCondition<-function(data,control='None',condition=c('Early Gilteritinib'),
                              cellLine='MOLM14',doNetworks=TRUE,prefix=''){
  
  total.mean.diffs<-data%>%
    dplyr::filter(cellLine==!!cellLine)%>%
    amlresistancenetworks::computeFoldChangePvals(control,condition)
  
  ##iterate over all conditions
  genes.with.values=purrr::map_df(condition,function(cond){
      diff.res<-total.mean.diffs%>%
          ungroup()%>%
          subset(Condition==cond)%>%
          dplyr::select(Gene,value=condition_to_control,p_adj,Condition)
  
      total.lig=amlresistancenetworks::computeGSEA(diff.res,
                                                     prefix=paste(prefix,cellLine,gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
  
    if(doNetworks){
    
      lig.network<-computeProteinNetwork(sig.vals=subset(diff.res,p_adj<0.05),
                                       all.vals=diff.res,nrand=1000)
      RCy3::createNetworkFromIgraph(lig.network,title=paste(prefix,cellLine,
                                                            gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
      }
    return(total.lig%>%mutate(Condition=cond))      
      })
  
  return(genes.with.values%>%mutate(CellLine=cellLine,Control=control))
}




plotHeatmapOfGenes<-function(genes,title,fname,gdat){
  print(title)
  print(fname)
#  print(genes)
  mat<-gdat%>%
    dplyr::select(Sample,Gene,value)%>%
    mutate(value=ifelse(is.na(value),0,value))%>%
    tidyr::pivot_wider(names_from=Sample,values_from=value,values_fn=list(value=mean,na.rm=T),values_fill=list(value=0))%>%
    subset(!is.na(Gene))%>%
    tibble::column_to_rownames('Gene')%>%
    as.matrix()
  
  metadat<-gdat%>%
    dplyr::select(Sample,cellLine,treatment,ligand)%>%
    distinct()%>%
    tibble::column_to_rownames('Sample')
  
  dmat=mat[intersect(genes,rownames(mat)),]
  print(dmat)
  pheatmap::pheatmap(dmat,annotation_col=metadat,main=title,cellwidth = 10,cellheight=10,filename=gsub(':','_',fname),)
dmat
}

plotLeadingEdgeGenes<-function(gdat,pathname='',respath='',prefix=''){
  tab<-read.csv(respath,sep='\t',header=T,stringsAsFactors = FALSE)%>%
    dplyr::select(geneSet,description,userId)
  matches<-grep(pathname,tab$description)
  print(paste('found',length(matches),'pathways with the term',pathname))
  tres<-tab[matches,]%>%
    rowwise()%>%
    mutate(title=description,id=geneSet,genes=stringr::str_split(userId,pattern=';'))%>%
    dplyr::select(title,id,genes)%>%
    mutate(fname=paste0(id,prefix,'.png'))
  
 # print(tres)
  apply(tres,1,function(x){
    plotHeatmapOfGenes(genes=x[['genes']],title=x[['title']],fname=x[['fname']],gdat=gdat)})
#  pr=group_map(tres, plotHeatmapOfGenes(genes,title,fname,gdat=gdat))
    
}

files=c('proteomics__MOLM14_None_vs_EarlyGilteritinib_gseaGO_result.txt','proteomics__MV411_None_vs_EarlyGilteritinib_gseaGO_result.txt')




rerunPrev=FALSE

if(rerunPrev){
  amlresistancenetworks::plotSingleProtein(gene='CPT1A',gilt.data,'abundance')
#  amlresistancenetworks::plotSingleProtein(gene='AURKB',gilt.data,'abundance')
 # amlresistancenetworks::plotSingleProtein(gene='CUL5',gilt.data,'abundance')
  
  dn=FALSE #no networks for one run
  
  
  #####here we run the results
  #quiz.data<-readRDS()
  
  #molm14
  m.vs.parental<-plotDataByCondition(gilt.data,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MOLM14',doNetwork=dn)
  
  m.late.early<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MOLM14',doNetwork=dn)
  
  
  #mv411
  v.vs.parental<-plotDataByCondition(gilt.data,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MV411',doNetwork=dn)
  
  v.late.early<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MV411',doNetwork=dn)
  
   
   sapply(c("FGF2","FLT3"),function(lig){
     m.early.net<-plotDataByCondition(subset(gilt.data,ligand%in%c('None',lig)),
                                    control='None',
                                      condition=c('Early Gilteritinib','Late Gilteritinib'),
                                      cellLine='MOLM14',doNetwork=dn,prefix=lig)
     v.early.net<-plotDataByCondition(subset(gilt.data,ligand%in%c('None',lig)),
                                      control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),
                                      cellLine='MV411',doNetwork=dn,prefix=lig)
   })
  
}

plotLeadingEdgeGenes(gdat=subset(gilt.data,cellLine=='MOLM14'),pathname='metab',respath=files[1],prefix='MOLM14_EarlyRes')
plotLeadingEdgeGenes(gdat=subset(gilt.data,cellLine=='MOLM14'),pathname='lipid',respath=files[1],prefix='MOLM14_EarlyRes')
plotLeadingEdgeGenes(gdat=subset(gilt.data,cellLine=='MV411'),pathname='metab',respath=files[2],prefix='MV411_EarlyRes')
plotLeadingEdgeGenes(gdat=subset(gilt.data,cellLine=='MV411'),pathname='lipid',respath=files[2],prefix='MV411_EarlyRes')
