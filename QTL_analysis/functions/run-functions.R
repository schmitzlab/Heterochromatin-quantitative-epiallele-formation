rm(list=ls())
setwd('set_your_working_dir')

#-----------------------------------------------------------------------------
# Step 1: Load the source code and required functions
library(data.table)
library(qtl)
library(dplyr)
library(xlsx)
library(stringr)
library(ggplot2)
library(reshape2)
library(stringr)
source('qtl-functions.R')
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Step 2: Assign directory pathways to phenotype (mp_traits_dir), genotype map (gc_genotype_dir) and markers map (marker_dir) and
#         output.dir where all of the codes are saved
mp_traits_dir <- 'phenotypes_sets/2.mchg.individual.gene.csv'
gc_genotype_dir <- 'mapping_marker_sets/epigenotypes_new.csv'
marker_dir  <- 'mapping_marker_sets/markers_new.csv'
output.dir <- ''
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Step 3. Run mapping_qtl function (preprocess the dataset, make a mapping and other dataset imporntant for the further steps of the analysis)
mapping_qtl(mp_traits_dir, gc_genotype_dir, marker_dir, output.dir)  
#-----------------------------------------------------------------------------  

#-----------------------------------------------------------------------------
# Step 4: Assign output.dir to input.dir (our directory where all of the datasets are saved in)
input.dir <- output.dir 
marker_dir <- 'mapping_marker_sets/markers_new.csv'
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Step 5: Perform getQTLpeaks_new functions to have the set of the QTL peaks with the given alpha levelfrom mapping for further steps
getQTLpeaks_new(input.dir = input.dir, marker_dir)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Step 6: Assign marker_dir and position_dir (locations of the phenotypes for transcis plot)
positions_dir <- 'phenotypes_positions_sets/2.location.gene.csv'
marker_dir <- 'mapping_marker_sets/markers_new.csv'
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Step 6: Plot QTL peaks and trans-cis plot
plot1 <- qtlpeaks.plot(input.dir, marker_dir)
plot2 <- transcis.plot(input.dir, marker_dir, positions_dir)
#-----------------------------------------------------------------------------









