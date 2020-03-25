##tests differential expression of proteins 

library(amlresistancenetworks)
library(dplyr)

sens.data<-readRDS(system.file('gilteritinibSensitivityData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(sens.data$Gene)