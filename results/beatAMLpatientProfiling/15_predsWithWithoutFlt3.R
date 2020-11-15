

library(amlresistancenetworks)
library(ggplot2)
library(dplyr)


##get aml
if(!exists('data.loaded')||!data.loaded){
  loadBeatAMLData()
  data.loaded=TRUE
}
syn=synapseLogin()


##filter df prot.net and auc.dat by FLT3 status...

prot.net.df<-prot.nets%>%select(proteomicNetworkDistance='distance',Gene='Community',`AML sample`)


pn.reg.results<-do.call(rbind,drugMolRegressionEval(auc.dat,prot.net.df,'proteomicNetworkDistance', 
                                      ctrp.cl.auc,cell.net.df))%>%as.data.frame()
pn.lr.results<-do.call(rbind,drugMolLogRegEval(auc.dat,prot.net.df,'proteomicNetworkDistance',
                                 ctrp.cl.auc,cell.net.df))%>%as.data.frame()%>%
  mutate(testMSE=unlist(testMSE)*10000)



##for each drug -  build a model in AML, eval in sanger and in PTRC

##now train model on AML and eval on depmap data
reg.preds<-purrr::map_df(list(mRNA='mRNALevels',
                              protein='proteinLevels',
                              gene='geneMutations'),~ drugMolRegressionEval(auc.dat,
                                                                        pat.data,
                                                                        .x,
                                                                        ctrp.cl.auc,
                                                                        cl.mol.dat,
                                                                        category='Condition'))


log.reg.preds<-purrr::map_df(list(protein='proteinLevels',
                                  mRNA='mRNALevels',
                              gene='geneMutations'),~ drugMolLogRegEval(auc.dat,
                                                                         pat.data,
                                                                         .x,
                                                                         ctrp.cl.auc,
                                                                         cl.mol.dat,
                                                                         category='Condition'))%>%
  mutate(testMSE=unlist(testMSE)*10000)
##i think we have all the results now, we can join them
ctrp.full.res<-rbind(mutate(log.reg.preds,method='LogisticRegression'), 
                mutate(reg.preds,method='LassoRegression'),
                mutate(pn.reg.results,method='LassoRegression'),
                mutate(pn.lr.results,method='LogisticRegression'))%>%

selectDataMatAndPlotCTRP<-function(compound,method,Molecular){
  print(paste(c(compound,method,Molecular)))
  if(Molecular=='proteomicNetworkDistance')
    dat.mat<-rename(cell.net.df,value=Molecular)
  else
    dat.mat<-rename(cl.mol.dat,value=Molecular)
  amlresistancenetworks::clusterSingleDrugEfficacy(drugName=compound,
                                    meth=method,
                                      data=Molecular,
                                      auc.dat=select(ctrp.cl.auc,-c(Value,source)),
                                      auc.thresh=100,
                                        new.results=rename(ctrp.full.res,var='compound'),
                                      data.mat=dat.mat,
                                    prefix='CTRP')

}

ctrp.full.res%>%
  subset(numFeatures>1)%>%
  rowwise()%>%mutate(selectDataMatAndPlotCTRP(compound,method,Molecular))

##now do the same for Sanger
##now train model on AML and eval on depmap data
reg.preds<-purrr::map_df(list(mRNA='mRNALevels',
                              protein='proteinLevels',
                              gene='geneMutations'),~ drugMolRegressionEval(auc.dat,
                                                                            pat.data,
                                                                            .x,
                                                                            sanger.cl.auc,
                                                                            cl.mol.dat,
                                                                            category='Condition'))

pn.reg.results<-do.call(rbind,drugMolRegressionEval(auc.dat,prot.net.df,'proteomicNetworkDistance', 
                                      sanger.cl.auc, cell.net.df))%>%as.data.frame()
pn.lr.results<-do.call(rbind,drugMolLogRegEval(auc.dat,prot.net.df,'proteomicNetworkDistance',
                                 sanger.cl.auc, cell.net.df))%>%as.data.frame()%>%
  mutate(testMSE=unlist(testMSE)*10000)


log.reg.preds<-purrr::map_df(list(protein='proteinLevels',
                                  mRNA='mRNALevels',
                                  gene='geneMutations'),~ amlresistancenetworksdrugMolLogRegEval(auc.dat,
                                                                            pat.data,
                                                                            .x,
                                                                            sanger.cl.auc,
                                                                            cl.mol.dat,
                                                                            category='Condition'))%>%
  mutate(testMSE=unlist(testMSE)*10000)

##i think we have all the results now, we can join them
sanger.full.res<-rbind(mutate(log.reg.preds,method='LogisticRegression'), 
                     mutate(reg.preds,method='LassoRegression'),
                          mutate(pn.reg.results,method='LassoRegression'),
                mutate(pn.lr.results,method='LogisticRegression'))%>%
  mutate(compound=as.character(compound),Molecular=as.character(Molecular))

selectDataMatAndPlotSanger<-function(compound,method,Molecular){
  print(paste(c(compound,method,Molecular)))
  if(Molecular=='proteomicNetworkDistance')
    dat.mat<-rename(cell.net.df,value=Molecular)
  else
    dat.mat<-rename(cl.mol.dat,value=Molecular)
  
  amlresistancenetworks::clusterSingleDrugEfficacy(drugName=compound,
                                                   meth=method,
                                                   data=Molecular,
                                                   auc.dat=select(sanger.cl.auc,-c(Value,source)),
                                                   auc.thresh=100,
                                                   new.results=rename(sanger.full.res,var='compound'),
                                                   data.mat=dat.mat,
                                                   prefix='Sanger')
  
}

sanger.full.res%>%
  subset(numFeatures>1)%>%
  rowwise()%>%mutate(selectDataMatAndPlotSanger(compound,method,Molecular))

new.results<-rbind(mutate(sanger.full.res,source='Sanger'),
                   mutate(ctrp.full.res,source='CTRP'))%>%
  mutate(MSE=unlist(MSE),testMSE=unlist(testMSE),resCor=unlist(resCor),numFeatures=unlist(numFeatures),
         numSamples=unlist(numSamples),compound=unlist(compound),Molecular=unlist(Molecular))


new.results<-new.results%>%subset(numFeatures!=0)

p1<-ggplot(new.results,aes(x=numFeatures,y=testMSE,col=Molecular,shape=method,
                           size=numSamples,alpha=0.7))+geom_point()+facet_grid(~source)+scale_y_log10()+scale_color_viridis_d()
ggsave('FLT3predictorSummary.png',p1,width=10)

p2<-ggplot(new.results,aes(x=compound,y=testMSE,col=Molecular,size=numSamples,shape=method))+
  geom_point()+facet_grid(source~.)+theme(axis.text.x=element_text(angle=90, hjust=1))+scale_y_log10()+scale_color_viridis_d()
ggsave('FLT3predictorDotPlot.png',p2,width=10,height=10)

p3<-ggplot(new.results,aes(x=method,y=testMSE,fill=Molecular))+geom_boxplot()+
  facet_grid(source~.)+scale_y_log10()+scale_color_viridis_d()
ggsave('FLT3dataComparison.png',p3,width=10)

