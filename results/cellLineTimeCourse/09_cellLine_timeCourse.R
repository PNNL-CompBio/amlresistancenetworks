##time course analysis

library(amlresistancenetworks)

library(dplyr)

protSample=readRDS(system.file('timeCourseData.Rds',package='amlresistancenetworks'))
phosSample=readRDS(system.file('timeCoursePhosphoData.Rds',package='amlresistancenetworks'))


clusterProtsAcrossTreatments<-function(prot.data){
  ##group data by treatment, cell line
  #cluster genes together, how many clusters do we get? 
  
  #return prots and clusters
  
}


doEnrichmentOFClusters<-function(){
  
}