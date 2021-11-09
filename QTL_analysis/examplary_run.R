setwd('set-your-working-dir')
rm(list=ls())

library(data.table)
library(qtl)
library(dplyr)
library(xlsx)
library(stringr)
library(ggplot2)
library(reshape2)
library(stringr)
source('functions/qtl-functions.R')

mp_traits_dir <- 'toy_dataset//phenotype.csv'
gc_genotype_dir <- 'toy_dataset//genotype.csv'
marker_dir  <- 'toy_dataset//markers.csv'
positions_dir <- 'toy_dataset//positions.csv'
output.dir <- 'examplary_outputs//'
input.dir <- output.dir 

mapping_qtl(mp_traits_dir, gc_genotype_dir, marker_dir, output.dir)  
getQTLpeaks_new(input.dir, marker_dir)

plot1 <- qtlpeaks.plot(input.dir, marker_dir)
plot2 <- transcis.plot(input.dir, marker_dir, positions_dir)
