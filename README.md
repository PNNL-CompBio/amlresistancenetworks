# AML Resistance Networks

This package supports a series of projects that analyze various aspects of AML progression, drug sensitivity and drug resistance. 

<!-- badges: start -->
[![R build status](https://github.com/sgosline/amlresistancenetworks/workflows/R-CMD-check/badge.svg)](https://github.com/sgosline/amlresistancenetworks/actions)
<!-- badges: end -->

## PTRC Projects

The individual project files now listed in separate repositories to contribute to the analysis of the data and subsequent publication. These are supported through the NCI CPTAC PTRC work shared by OHSU and PNNL. 

All data for this project is stored on [Synapse](http://synapse.org/ptrc). To gain access to the data you must request access until it is made available for publication.

### BEAT AML Proteogenomic Pilot project.
This project is the primary driver of methods development, as it uses diverse omics measurements to ascertain the impact that proteomics can have in drug response prognosis. THis project is located at the [BEAT AML Proteomics Site](http://github.com/sgosline/beatAMLproteomics). This work will support two publications. The first paper is a pilot project studying the analysis of a subset of the patient data, the second will be a broader analysis of the full cohort.

Some of the data ingest code can be found in the [./bin/] directory in this repository. 

### Studying the effect of drug combinations
This study will focus on the the immediate effect of drugs on AML cell lines to determine how they could act in concert. The data ingeset for this project as well as the subsequent analysis will be at the [AML Drug Combination]() site. 


### Comaring effects of Early vs. Late Resistance in AML
This work compares the effects of gilteritinib and quizartinib in models of early and late resistance in AML. Most of the data ingest is in the [./bin] directory here but most of the analysis is in the [Early Late Resistance](https://github.com/sgosline/earlyLateAMLresistance) repository.

### Effects of Cytokines on drug resistance
This project led by the Agarwal lab is currently in process. Most of the data ingest can be found in the [./bin] directory but the analysis has been moved to its own repository. The [Cytokine drug resistance](https://github.com/sgosline/cytokineDrugResistance).



