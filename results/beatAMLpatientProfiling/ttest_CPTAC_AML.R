if(!require(dplyr))
  install.packages('dplyr')
if(!require(MASS))
  install.packages("MASS")
if(!require(stringr))
  install.packages("stringr")
if(!require(ggplot2))
  install.packages('ggplot2')
if(!require(purrr))
  install.packages("purrr")
library(stringr)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)


syn=amlresistancenetworks::synapseLogin()
df = amlresistancenetworks::querySynapseTable("syn22172602")

newdf <- df[df$Tumor.VAF >0 & (df$Gene == "NRAS" | df$Gene =="FLT3") ,]
# Compare number of rows
nrow(df)
nrow(newdf)##this has no rows, i don't see where you put in the tumor data
#try this:
mut_status <- subset(df,Gene%in%c("NRAS","FLT3"))%>%
  dplyr::mutate(mutated=(`Tumor VAF`>0))%>%
  dplyr::select(c("AML sample",Gene,"mutated"))

all_muts <- mut_status%>%tidyr::pivot_wider(values_from='mutated',
                                            names_from=Gene,values_fn=list(mutated=any))


AML = amlresistancenetworks::querySynapseTable("syn22288960")
community <- AML %>%dplyr::rename(`AML sample`='net1')%>%
  dplyr::left_join(all_muts)%>%
  subset(net2_type=='community')%>%
  dplyr::mutate(same=(hyp1==hyp2))%>%subset(same)

##here use a join to add in the mutational status

#community = AML %>% filter(net2_type == 'community')%>%
#  dplyr::mutate(same=(hyp1==hyp2))

do_ttest<-function(df,gn='NRAS'){
  res<-df%>%dplyr::select(c("AML sample","distance",gn))%>%
    tidyr::pivot_wider(values_from=distance,names_from=gn)
  t.test(res$`FALSE`,res$`TRUE`,na.rm=TRUE)$p.value
}


nras.pvals=community%>%
  group_by(hyp1,net2)%>%
  do(tidy(t.test(distance~NRAS,data=.)))%>%
  dplyr::select(hyp=hyp1,community=net2,PVal=p.value)%>%mutate(Gene='NRAS')

flt3.pvals=community%>%
  group_by(hyp1,net2)%>%
  do(tidy(t.test(distance~FLT3,data=.)))%>%
  dplyr::select(hyp=hyp1,community=net2,PVal=p.value)%>%mutate(Gene='FLT3')

all.pvals<-rbind(nras.pvals,flt3.pvals)
