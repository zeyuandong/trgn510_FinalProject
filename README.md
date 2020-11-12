# TRGN510_Final Project
## Title
**Analysis of RNA-seq in bladder cancer in White people and Asians**

## Author
**Zeyuan Dong**

**zeyuando@usc.edu**

## Overview of project

**1. Overview**
* RNA-seq is an important technique used to study gene expression. For the RNA-SEQ data, edgeR (Robinson, McCarthy, and Smyth 2010) and limma packets from the Bioconductor project (Huber et al. 2015) provide a complete statistical set method for dealing with this problem.Bioconductor's advantage is that it can quickly and efficiently analyze RNA sequencing data. After obtaining the RNA-Seq gene expression matrix, we need to preprocess the data and then conduct the difference analysis.

**2. Objectives**
* I will use https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow_CHN.html outlined in the working process of the analysis of data.I will analyze RNA-seq in bladder cancer in white and Asian people and analyze the relationship between race and this cancer. 

**3. References/Links to Vignettes**
* https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow_CHN.html

## Data
* From the national cancer institute GDC data portal (https://portal.gdc.cancer.gov/), I will download the publicly available data sets. I divided the data into two groups, one white and one Asian. I selected 34 data points for each group, including 30 men and four women (only 4 Asian women developed bladder cancer after screening).

## Milestone 1(November 3rd)
* By filtering, I get the RNA-SEQ data from the GDC data portal and read it into RStudio. And according to the process of data preprocessing at https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html.

* Install the required R package and set up the work path. Unzip the data, rename the file, and transfer it to the working directory. The count data is read in through the readDGE function and processed. Deletion of low-expression genes and normalization of gene expression distribution. Unsupervised clustering of samples and MDS diagram drawing.

## Milestone 2(November 12th)
* Conduct differential expression analysis, use the camera for genome testing to establish a design matrix, make the comparison, eliminate heteroscedastic in counting data, detect the number of DE genes, detect single DE genes from top to bottom, and carry out visual processing of data.

* Use camera for genomic testing, search for appropriate Hallmark genomic collections and use camera functionality for competitive testing to assess whether genes in a given gene set rank high in differential expression relative to genes outside of that set.

## Deliverable(November 17th)
* R markdown
* Workflow released on https://rpubs.com/Zeyuan0311/689508
