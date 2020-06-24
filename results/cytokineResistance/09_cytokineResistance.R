##process cytokine data
library(amlresistancenetworks)

protData<-querySynapseTable('syn22213460')%>%subset(!is.nan(LogRatio))
phosData<-querySynapseTable('syn22213462')%>%subset(!is.nan(LogRatio))


##

plotAllData<-function(dat.table){
  library(ggfortify)
  met<-dat.table%>%dplyr::select(sample,CellType,TimePoint,Treatment)%>%
    distinct()
  #%>%
  #  tibble::column_to_rownames('sample')
    
  mat<-dat.table%>%dplyr::select(Gene,LogRatio,sample)%>%
    distinct()%>%
    mutate(LogRatio=as.numeric(LogRatio))%>%
    tidyr::pivot_wider(names_from='sample',values_from='LogRatio',values_fn=list(LogRatio=function(x) mean(x,na.rm=T)),values_fill=list(LogRatio=0))%>%
    tibble::column_to_rownames('Gene')
  
  autoplot(prcomp(t(mat)),data=met,colour='Treatment',shape='CellType')
 
}

plots=list(plotAllData(protData),plotAllData(phosData))
cowplot::plot_grid(plotlist=plots,labels=c("Bulk Proteomics",'Phosphoprotomics'),nrow=2)
ggsave('pcaOfSamples.png')


