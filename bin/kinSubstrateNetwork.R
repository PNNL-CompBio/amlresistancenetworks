##rebuild kinase-substrate interactome

##get networkin data
networkin_address='https://networkin.info/download/networkin_human_predictions_3.1.tsv.xz'
##get phosphositeplus

library(dplyr)
library(tidyr)

networkin <- download.file(networkin_address,'networkin.tsv.xz')
tab <- read.csv('networkin.tsv.xz',sep='\t')
##need to get kinase gene name

##need to get SUB_MOD_RSD value

psdb <-read.csv('inst/Kinase_Substrate_Dataset',sep='\t',skip=3)%>%
  subset(KIN_ORGANISM='human')
