library(amlresistancenetworks)
library(dplyr)

#get patient phospho/prot data

#get patient drug resposne for trametinib
source("../beatAMLpatientProfiling/beatAMLData.R")

#get trametinib treated cell lines

protData<-querySynapseTable('syn22217037')%>%subset(!is.nan(LogRatio))%>%
  mutate(Gene=unlist(Gene))

phosData<-querySynapseTable('syn22217040')%>%subset(!is.nan(LogRatio))%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))

#otherData<-''
otherPhosData<-querySynapseTable('syn22255396')%>%
  subset(!is.nan(LogFoldChange))%>%
  subset(cellLine=='HL60')%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))%>%
  rowwise()%>%
  mutate(Condition=paste(treatment,timePoint,sep='_'))

clinvars<-phosData%>%dplyr::select(Sample='sample',CellType,TimePoint,Treatment)%>%distinct()

kindat<-mapPhosphoToKinase(dplyr::rename(phosData,Sample='sample', LogFoldChange='LogRatio'))

protMat<-protData%>%dplyr::select(sample,Gene,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('Gene')

phosMat<-phosData%>%dplyr::select(sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('site')

otherPhosMat <-otherPhosData%>%ungroup()%>%dplyr::select(Sample,site,LogFoldChange)%>%
  tidyr::pivot_wider(values_from=LogFoldChange,names_from=Sample,
                     values_fn=list(LogFoldChange=mean),values_fill=list(LogFoldChange=0.0))%>%
  tibble::column_to_rownames('site')


kinMat<-kindat%>%dplyr::select(Sample,Kinase,meanLFC)%>%distinct()%>%
  tidyr::pivot_wider(values_from=meanLFC,names_from=Sample,values_fn=list(meanLFC=mean),
                     values_fill=list(meanLFC=0))%>%
  tibble::column_to_rownames('Kinase')

tpm<-apply(otherPhosMat,2,as.numeric)
rownames(tpm)<-rownames(otherPhosMat)
otherPhosMat<-tpm
#######
summary<-protData%>%dplyr::select(sample,CellType,TimePoint,Treatment)%>%distinct()%>%
  mutate(conditionName=stringr::str_c(CellType,TimePoint,Treatment,sep='_'))
print(summary)


##proteins, phosphosites, and kinases between resistant and sensitive cells
t0Comps=list(molm13_vs_resistant_phos=limmaTwoFactorDEAnalysis(phosMat,
                                                         dplyr::filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                         dplyr::filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample),
             molm13_vs_resistant_prot=limmaTwoFactorDEAnalysis(protMat,
                                                          dplyr::filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                          dplyr::filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample),
             molm13_vs_resistant_kin=limmaTwoFactorDEAnalysis(kinMat,
                                                          dplyr::filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                          dplyr::filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample))
          
##proteins, phosphopsites, and kinases between TRAM sens and resist patients
patPhosMat<-pat.phos%>%
  dplyr::select(Sample,site,LogFoldChange)%>%distinct()%>%
  tidyr::pivot_wider(values_from=LogFoldChange,names_from=Sample,values_fn=list(LogFoldChange=mean),
                     values_fill=list(LogFoldChange=0))%>%
  tibble::column_to_rownames('site')

pat.summmary<-auc.dat%>%subset(Condition=='Trametinib (GSK1120212)')%>%
  dplyr::select(`AML sample`,AUC)%>%mutate(Sensitive=AUC<100)%>%
  subset(`AML sample`%in%colnames(patPhosMat))%>%distinct()

patProtMat<-pat.data%>%
  subset(`AML sample`%in%pat.summmary$`AML sample`)%>%
  dplyr::select(Gene,`AML sample`,proteinLevels)%>%distinct()%>%
  tidyr::pivot_wider(values_from=proteinLevels,names_from=`AML sample`, values_fn=list(proteinLevels=mean),
                     values_fill=list(proteinLevels=0))%>%
  tibble::column_to_rownames('Gene')


patKinMat<-pat.kin%>%dplyr::select(Sample,Kinase,meanLFC)%>%distinct()%>%
  tidyr::pivot_wider(values_from=meanLFC,names_from=Sample,values_fn=list(meanLFC=mean),
                     values_fill=list(meanLFC=0))%>%
  tibble::column_to_rownames('Kinase')



patComps=list(tramSens_vs_resistant_phos=limmaTwoFactorDEAnalysis(patPhosMat,
                                                                  dplyr::filter(pat.summmary,!Sensitive)$`AML sample`,
                                                               dplyr::filter(pat.summmary,Sensitive)$`AML sample`),
             tramSens_vs_resistant_prot=limmaTwoFactorDEAnalysis(patProtMat,
                                                                 dplyr::filter(pat.summmary,!Sensitive)$`AML sample`,
                                                                 dplyr::filter(pat.summmary,Sensitive)$`AML sample`),
             tramSens_vs_resistant_kin=limmaTwoFactorDEAnalysis(patKinMat,
                                                                dplyr::filter(pat.summmary,!Sensitive)$`AML sample`,
                                                                dplyr::filter(pat.summmary,Sensitive)$`AML sample`))

plotGenesInMat<-function(limmaRes,resMat,pvalThresh=0.05, annotes,title){
  library(pheatmap)
  genes=c()
  try(genes<-subset(limmaRes,adj.P.Val<pvalThresh)$featureID)
  print(paste("Found",length(genes),'genes at corrected threshold of',pvalThresh))
  if(length(genes)==0){
    genes<-subset(limmaRes,P.Value<pvalThresh)$featureID
    print(paste('Now have',length(genes),'genes'))
  }
  red.mat<-resMat[intersect(genes,rownames(resMat)),]
 # print(red.mat)
  fname=paste0(title,'.pdf')
  pheatmap(red.mat,cellwidth = 10,cellheight=10,annotation_col = annotes,filename=fname,height=10)
  
}

##first plot patient data in cell lines
cellAnnotes<-summary%>%
  dplyr::select(CellType,TimePoint,Treatment,sample)%>%
  tibble::column_to_rownames('sample')

plotGenesInMat(patComps$tramSens_vs_resistant_phos,phosMat,0.005,cellAnnotes,'patientPhosSigInCellLines')
plotGenesInMat(patComps$tramSens_vs_resistant_prot,protMat,0.005,cellAnnotes,'patientProtSigInCellLines')
plotGenesInMat(patComps$tramSens_vs_resistant_kin,kinMat,0.05,cellAnnotes,'patientKinaseSigInCellLines')

##then plot cell line data in patients
patAnnotes<-pat.summmary%>%mutate(Sensitive=as.factor(Sensitive))%>%
  tibble::column_to_rownames('AML sample')
plotGenesInMat(patComps$tramSens_vs_resistant_prot,patProtMat,0.005,patAnnotes,'patientProtSigInPatient')

plotGenesInMat(t0Comps$molm13_vs_resistant_phos,patPhosMat,0.01,patAnnotes,'cellLinePhosSigInPatients')
plotGenesInMat(t0Comps$molm13_vs_resistant_prot,patProtMat,0.0001,patAnnotes,'cellLineProtSigInPatients')
plotGenesInMat(t0Comps$molm13_vs_resistant_kin,patKinMat,0.05,patAnnotes,'cellLineKinaseSigInPatients')
