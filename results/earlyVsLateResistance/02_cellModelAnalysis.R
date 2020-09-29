##run through drug resistance analysis


library(amlresistancenetworks)
library(dplyr)

##first run gilteritinib data
gilt.data<-querySynapseTable('syn22156807')
  #readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(gilt.data$Gene)

syn<-synapseLogin()
targ.data<-read.csv2(syn$get('syn22214245')$path,sep='\t',header=T)%>%
  dplyr::rename(Gene='T..Category')%>%
  tidyr::pivot_longer(-Gene,names_to='sample',values_to='value')%>%
  mutate(value=as.numeric(value))%>%
  rowwise()%>%
  mutate(samp2=stringr::str_replace(sample,'_3June20.+',''))%>%
  mutate(Patient=stringr::str_replace(samp2,'.*Ex16_',''))%>%
  dplyr::select(-c(samp2,sample))

norm.data<-targ.data%>%tidyr::pivot_wider(values_from=value,names_from=Patient)%>%tidyr::pivot_longer(-c(Gene,pool),values_to='value',names_to='Patient')%>%rowwise()%>%mutate(logRatio=log10(value)-log10(pool))%>%
  dplyr::select(-c(pool,value))%>%
  rename(value='logRatio')

add.cols<-norm.data%>%
  mutate(condition=stringr::str_replace(Patient,'[0-9]+',''))%>%
  rename(Sample='Patient')%>%
  mutate(ligand='none')%>%
  mutate(cellLine='Patient')

comb.data<-gilt.data%>%rename(condition='treatment')%>%
  dplyr::select(Gene,value,Sample,condition,ligand,cellLine)%>%
  subset(Gene%in%targ.data$Gene)%>%
  rbind(add.cols)%>%
  mutate(value=as.numeric(value))
                                             
ggplot2::ggplot(comb.data,ggplot2::aes(x=condition,fill=cellLine,y=value))+ggplot2::geom_boxplot()

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
                              cellLine='MOLM14',doNetworks=FALSE,doGSEA=FALSE,prefix=''){
  
  total.mean.diffs<-data%>%
    dplyr::filter(cellLine==!!cellLine)%>%
    amlresistancenetworks::computeFoldChangePvals(control,condition)
  
  ##iterate over all conditions
  genes.with.values=purrr::map_df(condition,function(cond){
      diff.res<-total.mean.diffs%>%
          ungroup()%>%
          subset(Condition==cond)%>%
          dplyr::select(Gene,value=condition_to_control,p_adj,Condition)
      
      if(doGSEA){
      total.lig=amlresistancenetworks::computeGSEA(diff.res,
                                                     prefix=paste(prefix,cellLine,gsub(' ','',control),'vs',gsub(' ','',cond),sep='_'))
      }else{
        total.lig <- diff.res
      }
    if(doNetworks){
    
      lig.network<-computeProteinNetwork(sig.vals=subset(diff.res,p_adj<0.01),
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


dn=TRUE
m.late.early<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='Early Gilteritinib',
                                  cellLine='MOLM14',doGSEA=TRUE)#doNetwork=dn)


rerunPrev=FALSE
if(rerunPrev){
  amlresistancenetworks::plotSingleProtein(gene='CPT1A',gilt.data,'abundance')
#  amlresistancenetworks::plotSingleProtein(gene='AURKB',gilt.data,'bundance')
 # amlresistancenetworks::plotSingleProtein(gene='CUL5',gilt.data,'abundance')
  
  dn=TRUE #no networks for one run
  
  
  #####here we run the results
  #quiz.data<-readRDS()
  
  #molm14
  m.vs.parental<-plotDataByCondition(gilt.data,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MOLM14',doNetwork=dn)


   v.vs.parental<-plotDataByCondition(gilt.data,control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),cellLine='MV411',doNetwork=dn)
  
   v.late.early<-plotDataByCondition(gilt.data,condition='Late Gilteritinib',control='Early Gilteritinib',cellLine='MV411',doNetwork=dn)
    #mv411
 
   
   sapply(c("FGF2","FLT3"),function(lig){
     m.early.net<-plotDataByCondition(subset(gilt.data,ligand%in%c('None',lig)),
                                    control='None',
                                      condition=c('Early Gilteritinib','Late Gilteritinib'),
                                      cellLine='MOLM14',doNetwork=dn,prefix=lig)
     v.early.net<-plotDataByCondition(subset(gilt.data,ligand%in%c('None',lig)),
                                      control='None',condition=c('Early Gilteritinib','Late Gilteritinib'),
                                      cellLine='MV411',doNetwork=dn,prefix=lig)
   })
  


plotLeadingEdgeGenes(gdat=subset(gilt.data,cellLine=='MOLM14'),pathname='metab',respath=files[1],prefix='MOLM14_EarlyRes')
plotLeadingEdgeGenes(gdat=subset(gilt.data,cellLine=='MOLM14'),pathname='lipid',respath=files[1],prefix='MOLM14_EarlyRes')
plotLeadingEdgeGenes(gdat=subset(gilt.data,cellLine=='MV411'),pathname='metab',respath=files[2],prefix='MV411_EarlyRes')
plotLeadingEdgeGenes(gdat=subset(gilt.data,cellLine=='MV411'),pathname='lipid',respath=files[2],prefix='MV411_EarlyRes')

}