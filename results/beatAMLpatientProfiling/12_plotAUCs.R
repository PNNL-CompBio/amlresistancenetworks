###plot AUC values by mutation status

#source("beatAMLdata.R")

plotAUCByMutation<-function(geneName){
  library(ggplot2)
  library(ggpubr)
  mut.scores<-subset(pat.data,Gene==geneName)%>%
    mutate(isMutated=(geneMutations>0))%>%
    dplyr::select(`AML sample`,isMutated)%>%distinct()
  
  df<-auc.dat%>%
    dplyr::select(`AML sample`,Condition,AUC)%>%
    left_join(mut.scores)%>%
   subset(!Condition%in%c('Staurosporine'))
  
  sigs<-NULL  
  try(sigs<-  df%>%group_by(Condition)%>%t_test(AUC ~ isMutated)%>% 
        adjust_pvalue(method = "bonferroni") %>%
           add_significance("p.adj")%>%
      add_xy_position(x = "Condition",dodge=0.8))
    
    p<-ggboxplot(data=df,x="Condition",y="AUC",fill="isMutated",title=paste(geneName,'mutations compared to drug response'),
                 add.params = list(group='isMutated'),palette='npg')+
        theme(axis.text.x=element_text(angle=90, hjust=1)) 

  p2<-p
  if(!is.null(sigs))
    p2<- p+ stat_pvalue_manual(
      sigs, label = "p", x='Condition',tip.length = 0.01,
      bracket.nudge.y = -2) 
  ggsave(paste0(geneName,'mutationsAcrossDrugs.png'),p2)
  
}

for(gs in c('FLT3','NRAS','KRAS','NPM1','SRSF2','AQR','DNMT3A'))
  plotAUCByMutation(gs)