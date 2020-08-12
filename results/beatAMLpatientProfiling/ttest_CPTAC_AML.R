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

newdf <- df[df$Tumor.VAF >0 & (df$Gene == "NRAS" | df$Gene =="FLT3") ,]
# Compare number of rows
nrow(df)
nrow(newdf)##this has no rows, i don't see where you put in the tumor data
#try this:
mut_status <- subset(df,Gene%in%c("NRAS","FLT3"))%>%
  mutate(mutated=(`Tumor VAF`>0))%>%
  dplyr::select(c("AML sample",Gene,"mutated"))

AML = amlresistancenetworks::querySynapseTable("syn22288960")

##here use a join to add in the mutational status

community = AML %>% filter(net2_type == 'community')

select_data = function(hype,hype_type,community_number){
  data = community %>% filter(net2 == community_number)
  if (hype == 1){
    data = data %>% filter(hyp1 == hype_type)
  }
  else{
    data = data %>% filter(hyp2 == hype_type)
  }
  return (data)
}

ttest_comm = function(test_data,sample_number1,sample_number2){
  sample1 = data %>% filter(net1 == sample_number1)
  sample2 = data %>% filter(net1 == sample_number2)
  model = t.test(sample1$distance, sample2$distance, paired = FALSE)
  return (model$p.value)
}

community_list = unique(community$net2)
#Hyp value, 
#change this value to 1 for hyp1, or 2 for hyp2
hype = 1
#hyp type value
#change this to mutation to test for mutations or proteomics
hype_type = "mutations"

result <- data.frame("Communities" = NA, "Sample_1" = NA, "Sample_2" = NA,"Pvalue" = NA)
for (com in community_list){
  data = select_data(hype,hype_type,com)
  for(x in data$net1){
    for (y in data$net1){
      if (x != y){
        print(paste("Community number", com))
        print(paste("sample 1", x))
        print(paste("sample 2", y))
        pvalue = try(ttest_comm(data,x,y), silent = T)
        temp = try(data.frame("Communities" = com, "Sample_1" = x, "Sample_2" = y,"Pvalue" = pvalue),silent = T)
        result = rbind(result,temp)
        
      }
    }
  }
  
}

clean_result = result %>% filter(str_detect(result$Communities,"Error",negate = TRUE))






