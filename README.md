# Heterochromatin-quantitative-epiallele-formation
Main scripts for Methylome data analysis in this study "Heterochromatin is a quantitative trait locus associated with spontaneous epiallele formation
" is listed here. For questions, contact Yinwen Zhang at yz46606@uga.edu.


## System requirements
• perl 5 with the following installed modules: Getopt, Data, Math, Statistics

• Sufficient amount of RAM memory  (recommended minimum is 20 GB)


## Installation
To install the required packages type the following command in the R console.
```
cpanm Module::Name
```
requried modules were already imported at the head of perl scripts, so it may report errors if modules were not installed correctly before running scripts.


## Input file examples
The demo input files is located in example_input and it consists of the following files:

• test_103R.tsv - the output file generated from methylpy (a python-based pipeline for bisulfate sequencing data). it is a tab-separated text file, each row corresponds to one cyosine's location and methylation information in the genome. here is an example file for an ddm1 epiRIl sample.

[go to github of methylpy for a detailed illustration of output format!](https://github.com/yupenghe/methylpy)

• test_WT-parent.tsv - same as above, but for an col-0 wide type parent sample.

• test_Athaliana_447_Araport11.gene.gff3 - example gene annotation files used for the perl scripts

• test.primaryTranscriptOnly.fa - example primary transcripts annotation files used for the perl scripts

• test_gene.status.txt - a text file that give the methylation status (UM/gbM/teM) for each gene.


## Output examples

The demo output files is located in example_output and it consists of the following files:

• allc_WT.tsv.1.raw.out - 

• allc_WT.tsv.2.pvalue.out - 

• allc_103R.bin.table -


## Functions overview and Instructions for use 

• 
```
cpanm Module::Name
```

• 
```
cpanm Module::Name
```

• 
```
cpanm Module::Name
```


## Expected runtime
Here gives an estimated run time for real datasets

• the runtime for define_gbm.genes.v2.pl mainly depend on the sizes of input files, for an arabidopsis sample, it may take around 40-60 mins.

• the runtime for binominal.test.use.selfratio.pl will take around 5 mins.

• the runtime for formetaplot.v3.pl mainly depend on the sizes of input files, for an arabidopsis sample, it may take around 40-60 mins.

