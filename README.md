# Heterochromatin-quantitative-epiallele-formation
Main scripts for Methylome data analysis in this study "Heterochromatin is a quantitative trait locus associated with spontaneous epiallele formation
" is listed here. For questions, contact Yinwen Zhang at yz46606@uga.edu.
## System requirements
• R version 4.0.2 with the following installed packages: data.table, qtl, dplyr, xlsx, stringr, ggplot2, reshape2

• Sufficient amount of RAM memory to run the QTL analysis (recommended minimum is 8 GB)

• The administrative privileges are required to install and run R‑Studio utilities.

• A network connection for data recovering over network.

## Installation
To install the required packages type the following command in the R console.
```
install.packages(c('data.table', 'qtl', 'dplyr', 'xlsx', 'stringr', 'ggplot2', 'reshape2'))
```
To load the provided functions use load command in the R console.
```
source('functions/qtl-functions.R')
```

## Demo
examplary_run.R file presents the demo of running the pipeline for toy datasets located in toy_dataset directory (containing phenotype.csv, genotype.csv, positions.csv and marker.csv files).

To run the demo file, use the R (RStudio is recommended), change the first row of the script by typing your working directory where you have placed all the files from this repository, eg:
```
setwd('/home/robert/Zhang-et-al-project')
```
and run the script in R console.

## Expected output

The expected output is located in toy_outputs and it consists of the following files:

• PERM-all-traits_log_nm.Rdata - Rbinary dataset with the results from permutation test performed by rqtl package;

• QTL-direction-data-all-traits_log_nm.Rdata - Rbinary dataset with the effect direction for a given association (epimarker and phenotype);

• QTL-mapping-all-traits_log_nm.Rdata - Rbinary dataset with the results from the mapping procedure performed by rqtl package;

• REF-data-all-traits_log_nm.Rdata - Rbinary dataset with the crossing file used as the input for mapping procedure by rqtl package;

• SCANPOSINFO-all-traits_log_nm.Rdata - Rbinary dataset with positions of the markers after mapping procedure;

• gc_genotype_log_nm.csv - comma-separated csv file with the genotype data after preprocessing used for running the mapping procedure by rqtl package;

• mp_traits_log_nm.csv - comma-separated csv file with the phenotype data after preprocessing used for running the mapping procedure by rqtl package;

• info_peaks_log_nm.csv - comma-separated csv file with the positions of the peaks from mapping outputs.

## Expected runtime
For 20 phenotypic traits, the whole pipeline takes around 1 minute, whereas:

• Running a mapping_qtl function could take 2 seconds per one trait (40 seconds for toy datasets);

• Running a getQTLpeaks_new function could take less than 1 second per one trait (10 seconds for toy datasets);

• Running a qtlpeaks.plot and transcis.plot functions could take no longer than 10 seconds in total.


## Instructions for use
• Required input datasets:
##### phenotype dataset - csv comma divided file with quantitative traits for each line and phenotype
##### genotype dataset - csv comma divided file with epigenotype profile (M/U) for each epimarker and line
##### marker information dataset - directory to csv comma divided file with marker chromosome and start-end positions
##### positions_dir - directory to csv comma divided file with phenotype (epiRIL positions)
