##prep data, run once

library(amlresistancenetworks)

##this data is for the time course stud
getTimeCoursePhosphoData()
getTimeCourseData()

##only run this once
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