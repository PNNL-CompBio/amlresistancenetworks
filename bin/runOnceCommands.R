##prep data, run once

library(amlresistancenetworks)

readAndTidyQuizProtMeasures()
readAndTidyQuizPhosphoProtMeasures()

##this data is for the time course study
getTimeCoursePhosphoData()
getTimeCourseData()

##this data is for the early/late resistance study
readAndTidyProtMeasures()
readAndTidyPhosphoProtMeasures()

#readAndTidyQuizProtMeasures()

##this data is for combination measurements
readAndTidySensProtMeasure()
readAndTidySensPhosMeasures()


#this data is for the patient data
getPatientPhosphoBaselines()
#getPatientBaselines()

#new patient sampling
getPatientMetadata()
getPatientMolecularData()

##add in agarwal data
 