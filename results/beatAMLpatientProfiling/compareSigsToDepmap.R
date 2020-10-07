

##eval sigs on depmap

library("depmap")
library("ExperimentHub")
eh <- ExperimentHub()

##first get eh <- ExperimentHub()
#all the samples, filtering for AML 
samps<-depmap::depmap_metadata()



#then get drug, gene, protein, x
drug<-depmap::depmap_drug_sensitivity()
proteomic<-depmap::depmap_proteomic()
rna<-depmap::depmap_TPM()
genes<-depmap::depmap_mutationCalls()


##load drug results
all.preds<-read.rds('mostlyCompletePredictions.rds')