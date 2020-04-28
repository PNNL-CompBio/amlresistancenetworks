##prep data, run once

library(amlresistancenetworks)

##only run this once
readAndTidyProtMeasures()

readAndTidyPhosphoProtMeasures()

#readAndTidyQuizProtMeasures()

readAndTidySensProtMeasure()

readAndTidySensPhosMeasures()

getPatientPhosphoBaselines()

getPatientBaselines()