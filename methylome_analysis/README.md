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


## Input examples
The demo input files is located in example_input and it consists of the following files:

• test_103R.tsv - the output file generated from methylpy (a python-based pipeline for bisulfate sequencing data). it is a tab-separated text file, each row corresponds to one cyosine's location and methylation information in the genome. here is an example file for an ddm1 epiRIl sample.

[go to github of methylpy for a detailed illustration of output format!](https://github.com/yupenghe/methylpy)

• test_WT-parent.tsv - same as above, but for an col-0 wide type parent sample.

• test_Athaliana_447_Araport11.gene.gff3 - example gene annotation files used for the perl scripts

• test.primaryTranscriptOnly.fa - example primary transcripts annotation files used for the perl scripts

• test_gene.status.txt - a text file that give the methylation status (UM/gbM/teM) for each gene.


## Output examples

The demo output files is located in example_output and it consists of the following files:

• allc_WT.tsv.1.raw.out - output file generated from `define_gbm.genes.v2.pl`

• allc_WT.tsv.2.pvalue.out - output file generated from `binominal.test.use.selfratio.pl`

• allc_103R.bin.table - output file generated from `formetaplot.v3.pl`


## Functions overview and Instructions for use 

• perl scripts `define_gbm.genes.v2.pl` is used to call the number of methylated cytocine sites for each C-context (CG/CHG/CHH) within CDS regions of each gene.

```
perl scripts/define_gbm.genes.v2.pl -fa example_input/test.primaryTranscriptOnly.fa -gff example_input/test_Athaliana_447_Araport11.gene.gff3 -tsv example_input/test_WT-parent.tsv -o testdir
```

• perl scripts `binominal.test.use.selfratio.pl` is used to perform the binorminal test based on output file generated from `define_gbm.genes.v2.pl`. The output gives a adjust-pvalue for each C-context (CG/CHG/CHH) on genes to indicate significance of methylation. therefore, this file is used for classify gene's methylation status (gbM/UM/teM).

```
perl scripts/binominal.test.use.selfratio.pl -i example_output/allc_WT.tsv.1.raw.out -o testdir
```

• perl scripts `formetaplot.v3.pl` is used to exploring DNA methylation patterns on genes or repeats. The gene or repeat body was divided into 20 windows. Additionally, regions 1,000 bp upstream and downstream were each divided into 20 50-bp windows. Methylation levels were calculated for each window. The mean methylation level for each window was then calculated for all genes, respectively. 

```
perl scripts/formetaplot.v3.pl -gff example_input/test_Athaliana_447_Araport11.gene.gff3 -l example_input/test_gene.status.txt -allc example_input/test_103R.tsv -o testdir
```

the output bin table from the scripts can be read in R to generate metaplot. Here provide an example of making metaplot in R.
```
setwd("your work dir")
library("ggplot2")
library("cowplot")
library("reshape")


d=read.table("./example_output/allc_103R.bin.table",header = T)
head(d)

d3=reshape(d, 
           direction = "long",
           varying = c("CG","CHG","CHH"),
           v.names = "Value",
           timevar = "context",
           times = c("CG","CHG","CHH"))
head(d3)
pd=c("#D00000","#3F88C5","#FFBA08")

p1=ggplot()+
  geom_line(data=d3,aes(x=binnum,y=Value,col=context),size=1)+
  facet_grid(.~Group)+
  scale_color_manual(values=pd)+
  #scale_colour_gradient(low = "yellow", high = "red", na.value = NA,guide = F)+
  theme_classic(base_size=25) +
  xlab("")+
  ylab("Methylation level (%)")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
p1

pdf('metaplot.example.pdf',
    width=12,
    height=6)

plot_grid(p1,nrow=1,align='hv',axis='l')

dev.off()

```


## Expected runtime
Here gives an estimated run time for real datasets

• the runtime for define_gbm.genes.v2.pl mainly depend on the sizes of input files, for an arabidopsis sample, it may take around 40-60 mins.

• the runtime for binominal.test.use.selfratio.pl will take around 5 mins.

• the runtime for formetaplot.v3.pl mainly depend on the sizes of input files, for an arabidopsis sample, it may take around 40-60 mins.

