#do gilteritinib analysis

library(amlresistancenetworks)

#lods the gilterinib data to analyze similarities
readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))

#step 1, compute LFC of FGF2, FL in early gilterinitib samples for both cell lines