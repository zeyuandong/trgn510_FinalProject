# TRGN510_Final Project
## Title
**Analysis of RNA-seq from White people aged 40 to 60 years who had died of bronchus and lung cancer**

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

## Progress Update
* First, I download the files I need from the website. The Workflow Type I selected was HTseq-Counts.
* I installed TCGAbiolinks in R and loaded them.
<pre>
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", version = "3.8")
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query = query)
setwd("/Users/dongzeyuan/Desktop/510final") 
</pre>

* Move all files to the same folder and name it SampleFiles.
<pre>
dir.create("SampleFiles")
filepath <- dir(path ="./gdc_download_20201102_215714.577827",full.names = TRUE) 
for(wd in filepath){
  files <-dir(path = wd,pattern="gz$") # View files under the path that meet the criteria
  fromfilepath <- paste(wd,"/",files,sep ="")
  tofilepath <- paste("./SampleFiles/",files,sep ="")
  file.copy(fromfilepath,tofilepath)
}
</pre>

* Unzip all files & delete original files
<pre>
setwd("./SampleFiles")
countsFiles <-dir(path = "./",pattern="gz$") # View files under the path that meet the criteria
library(R.utils)
sapply(countsFiles, gunzip) # To unzip the function gunzip, install the R.utils package
</pre>

* Process the JSON file and get the Counts matrix
<pre>
library(rjson)
metadata_json_File <- fromJSON(file="../metadata.cart.2020-11-02.json") 
json_File_Info <- data.frame(filesName = c(),TCGA_Barcode = c())
for(i in 1:length(metadata_json_File)){
  TCGA_Barcode <- metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name <- metadata_json_File[[i]][["file_name"]]
  json_File_Info <- rbind(json_File_Info,data.frame(filesName = file_name,TCGA_Barcode = TCGA_Barcode))
}
rownames(json_File_Info) <- json_File_Info[,1]
write.csv(json_File_Info,file = "../json_File_Info.CSV")
</pre>
<pre>
filesName_To_TCGA_BarcodeFile <- json_File_Info[-1]
countsFileNames<-dir(pattern="counts$") 
allSampleRawCounts <- data.frame()
for(txtFile in countsFileNames){
  #One file at a time
  SampleCounts <- read.table(txtFile,header =FALSE)
  rownames(SampleCounts) <- SampleCounts[,1]
  SampleCounts <- SampleCounts[-1]
  # Column names are named based on the file name in the filesName_To_TCGA_BarcodeFile corresponding to the barcode in the file name
  colnames(SampleCounts) <- filesName_To_TCGA_BarcodeFile[paste(txtFile,".gz",sep = ""),]
  if (dim(allSampleRawCounts)[1]== 0){
    allSampleRawCounts <- SampleCounts
  }
  else 
  {allSampleRawCounts<- cbind(allSampleRawCounts,SampleCounts)}
}
write.csv(allSampleRawCounts,file = "../allSampleRawCounts.CSV")
ensembl_id <- substr(row.names(allSampleRawCounts),1,15)
rownames(allSampleRawCounts) <- ensembl_id
#The rawcard.csv file differs from the allsamplerawcard.csv file in that the Ensembl of the line name has removed the version number
write.csv(allSampleRawCounts,file = "../RawCounts.CSV")
</pre>

* ID conversion
<pre>
# Add a column of Ensemble_ID to the left of the RawCounts data box
RawCounts <- allSampleRawCounts
Ensembl_ID <- data.frame(Ensembl_ID = row.names(RawCounts))
rownames(Ensembl_ID) <- Ensembl_ID[,1]
RawCounts <- cbind(Ensembl_ID,RawCounts)
# Define a function that gets the Ensemble_ID corresponding to the gene name from the gencode.v35.annotation.gtf file
get_map = function(input) {
  if (is.character(input)) {
    if(!file.exists(input)) stop("Bad input file.")
    message("Treat input as file")
    input = data.table::fread(input, header = FALSE)
  } else {
    data.table::setDT(input)
  }
  
  input = input[input[[3]] == "gene", ]
  
  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  
  
  gene_id = sub(pattern_id, "\\1", input[[9]])
  gene_name = sub(pattern_name, "\\1", input[[9]])
  
  Ensembl_ID_TO_Genename <- data.frame(gene_id = gene_id,
                                        gene_name = gene_name,
                                        stringsAsFactors = FALSE)
  return(Ensembl_ID_TO_Genename)
}


Ensembl_ID_TO_Genename <- get_map("../gencode.v35.annotation.gtf")  
gtf_Ensembl_ID <- substr(Ensembl_ID_TO_Genename[,1],1,15)
Ensembl_ID_TO_Genename <- data.frame(gtf_Ensembl_ID,Ensembl_ID_TO_Genename[,2])
colnames(Ensembl_ID_TO_Genename) <- c("Ensembl_ID","gene_id")
write.csv(Ensembl_ID_TO_Genename,file = "../Ensembl_ID_TO_Genename.csv")
</pre>

* Get the final matrix
<pre>
#Data fusion
mergeRawCounts <- merge(Ensembl_ID_TO_Genename,RawCounts,by="Ensembl_ID")
#Sort by gene_ID column
mergeRawCounts <- mergeRawCounts[order(mergeRawCounts[,"gene_id"]),]
#Index based on the gene_ID column
index<-duplicated(mergeRawCounts$gene_id)
#We want the behavior to be FALSE, so i.e., plus "!" . False means no repetition
mergeRawCounts <- mergeRawCounts[!index,]
#Use the gene name as the row name
rownames(mergeRawCounts) <-mergeRawCounts[,"gene_id"]
#Delete the first two columns, because we want to have Numbers in the matrix when we do the operation
LUAD_Counts_expMatrix <- mergeRawCounts[,-c(1:2)]
#Save the file
write.csv(LUAD_Counts_expMatrix,file = "../LUAD_Counts_expMatrix.csv")
</pre>


## Milestone 1(November 3rd)
* I will download the data I need from the above website, obtain the RNA-SEQ gene expression matrix, and then preprocess the data in Rstudio.

## Milestone 2(November 12th)
* I want to analyze the correlation between different ages and Bronchus and Lung cancer through the analysis of differences, and I want to show them through the data graph.

## Deliverable(November 17th)
* R markdown
