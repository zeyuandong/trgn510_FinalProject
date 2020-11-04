# TRGN510_Final Project
## Title
**RNA sequence analysis was performed in caucasians aged 41-50 years and 51-60 years who died of bronchial and lung cancer**

## Author
**Zeyuan Dong**

**zeyuando@usc.edu**

## Overview of project

**1. Overview**
* RNA-seq is an important technique used to study gene expression. For the RNA-SEQ data, edgeR (Robinson, McCarthy, and Smyth 2010) and limma packets from the Bioconductor project (Huber et al. 2015) provide a complete statistical set method for dealing with this problem.Bioconductor's advantage is that it can quickly and efficiently analyze RNA sequencing data. After obtaining the RNA-Seq gene expression matrix, we need to preprocess the data and then conduct the difference analysis.

**2. Objectives**
* I will use the Bioconductor method to analyze RNA-seq from White people aged 20 to 80 years who had died of Bronchus and Lung cancer, and analyze the relationship between age and this cancer

**3. References/Links to Vignettes**
* https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow_CHN.html


## Data
* I will download the data from this website:https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.demographic.race%22%2C%22value%22%3A%5B%22white%22%5D%7D%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22cases.diagnoses.age_at_diagnosis%22%2C%22value%22%3A%5B7305%5D%7D%2C%22op%22%3A%22%3E%3D%22%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22cases.diagnoses.age_at_diagnosis%22%2C%22value%22%3A%5B29584%5D%7D%2C%22op%22%3A%22%3C%3D%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.disease_type%22%2C%22value%22%3A%5B%22squamous%20cell%20neoplasms%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22bronchus%20and%20lung%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-LUSC%22%5D%7D%7D%5D%7D&searchTableTab=files



## Milestone 1(November 3rd)
* I will download the data I need from the above website, obtain the RNA-SEQ gene expression matrix, and then preprocess the data in Rstudio.

## Milestone 2(November 12th)
* I want to analyze the correlation between different ages and Bronchus and Lung cancer through the analysis of differences, and I want to show them through the data graph.

## Deliverable(November 17th)
* R markdown
