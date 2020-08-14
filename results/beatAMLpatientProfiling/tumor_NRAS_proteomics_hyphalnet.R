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


syn=amlresistancenetworks::synapseLogin()
df = amlresistancenetworks::querySynapseTable("syn22172602")

#filter for the samples that have the NRAS and FLT3 mutations with a Tumor VAF>0

newdf <- df[df$Tumor.VAF >0 & (df$Gene == "NRAS" | df$Gene =="FLT3") ,]

AML = amlresistancenetworks::querySynapseTable("syn22288960")
#create two new columns for NRAS and FLT3 in AML dataframe, with "YES" in column depicting sample has mutation for that gene and "NO" depicting no mutation for that gene

AML$NRAS = NA
AML$FLT3 = NA
n = length(AML$net1)
m = length(newdf$AML.sample)
for(i in 1:n){
  for(j in 1:m){
    if(AML$net1[i] == newdf$AML.sample[j]){
      if(newdf$Gene[j] == "NRAS"){
        AML$NRAS[i] = "YES"
        AML$FLT3[i] = "NO"
      }
      if(newdf$Gene[j] == "FLT3"){
        AML$FLT3[i] = "YES"
        AML$NRAS[i] = "NO"
      }
    }
  }
}
AML[is.na(AML)] = "NO"

#filter for hyp1 and hyp2 to be the same, and filter for just communities
hyp_type = AML %>% filter(AML$hyp1 == "proteomics" & AML$hyp2 == "proteomics")
community = hyp_type %>% filter(hyp_type$net2_type == "community")

community$category = rep('NRAS Mutation',nrow(community))
community$category[grep('NO',community$NRAS)]<-'No NRAS Mutation'

#perform ttests for communities, based on grouping whether samples have mutation or not
ttest_comm = function(x){
  data_NRAS = community %>% filter(net2 == x)
  normal_NRAS = data_NRAS %>% filter(str_detect(NRAS,'NO'))
  mutation_NRAS = data_NRAS %>% filter(str_detect(NRAS,'YES'))
  model = t.test(mutation_NRAS$distance, normal_NRAS$distance, paired = FALSE)
  return (model)
}


community_list = unique(community$net2)
#print out ttests for each community
for (val in community_list){
  print(paste("Community number", val))
  model = ttest_comm(val)
  print(model)
}

pvals<-community_list%>%
  purrr::map(~ttest_comm(.)$p.value)%>%as.numeric()%>%unlist()

#now adding larger data frame to plot
res.df<-data.frame(Community=unlist(community_list),PValue=pvals)
community <- community%>%rename(Community=net2)%>%left_join(res.df)

ggplot(subset(community,PValue<0.3),aes(x=Community,y=distance,fill=category))+geom_boxplot()


