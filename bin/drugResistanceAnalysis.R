##run through drug resistnace analysis


library(amlresistancenetworks)
library(dplyr)

gilt.data<-readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
prot.univ<-unique(gilt.data$Gene)



quiz.data<-readRDS()