library(data.table)
library(qtl)
library(dplyr)
library(xlsx)
library(stringr)
library(ggplot2)
library(reshape2)
library(stringr)

#-----------------------------------------------------------------------------
# Required input datasets
#-----------------------------------------------------------------------------
# mp_traits_dir - directory to csv comma divided file with quantitative traits for each line and phenotype
# (column names: phenotype/epiRILs; row names: lines)
#-----------------------------------------------------------------------------
# gc_genotype_dir - directory to csv comma divided file with epigenotype profile (M/U) for each epimarker and line
# (column names: epimarkers; row names: lines)
# first row: chromosome for a given marker
#-----------------------------------------------------------------------------
# marker_dir - directory to csv comma divided file with marker chromosome and start-end positions
# (column names: marker, chr, start, end)
#-----------------------------------------------------------------------------
# positions_dir - directory to csv comma divided file with phenotype (epiRIL positions)
# (column names: name, chr, start, end)
#-----------------------------------------------------------------------------

##### mapping_qtl function #####
# mapping_qtl - perform preprocessing of the data (based on the log transformation and preparing datasets for QTL analysis)
#             - perform a mappinh and permutations (QTL analysis)
#             - calculate genotype probabilites and map
#             - calculate effect direction (phenotype ~ epigenotype)
mapping_qtl <- function(mp_traits_dir, gc_genotype_dir, marker_dir, output.dir)
{
  # Step 1: Read required dataset for phenotypes (mp_traits), genotype map (gc_genotype) and marker map (markers_info)
  
  mp_traits <- fread(mp_traits_dir, header = TRUE)
  gc_genotype <- read.csv(gc_genotype_dir)
  markers_info <- fread(marker_dir)
  
  # Step 2: Organize the data by the unique IDs of the lines (it should match trait and genotype datasets) 
  
  IDs <- gc_genotype$ID[-1]
  for (i in 1:8)
  {
    IDs <- str_remove(IDs, paste("L", i, sep = ""))
  }
  IDs <- str_remove_all(IDs, "[R_merged]")
  IDs <- paste(IDs, "R", sep = "")
  
  IDs_traits <- mp_traits$sample
  for (i in 1:8)
  {
    IDs_traits <- str_remove(IDs_traits, paste("L", i, sep = ""))
  }
  IDs_traits <- str_remove_all(IDs_traits, "[R_merged]")
  IDs_traits <- paste(IDs_traits, "R", sep = "")
  mp_traits$sample <- IDs_traits
  
  gc_genotype <- gc_genotype[c(1,match(mp_traits$sample, IDs)+1),]
  gc_genotype$ID <- c("", IDs[match(mp_traits$sample, IDs)])
  colnames(mp_traits) <- c("ID", colnames(mp_traits)[-1])
  
  traits <- colnames(mp_traits)
  traits <- str_remove_all(traits, "-")
  colnames(mp_traits) <- traits
  mp_traits <- as.matrix(mp_traits)
  
  
  # Step 3: Make log transformation of the data and save if it is closer to normal distribution 
  #         (evaluation is done by checking the Shapiro Wilk's p-value - if it is higher after transformation, it is saved to the dataset)
  
  for (i in 2:dim(mp_traits)[2])
  {
    dat <- as.numeric(mp_traits[,i])
    pval_before <- shapiro.test(dat)$p.val
    dat <- dat + 1
    dat_log <- log(dat)
    pval_after <- shapiro.test(dat_log)$p.val
    if (pval_after > pval_before)
    {
      message("log trait: ", colnames(mp_traits)[i])
      mp_traits[,i] <- dat_log
    }
  }
  mp_traits <- as.data.frame(mp_traits)
  
  # Step 4: Save datasets for mapping in output.dir directory
  
  fwrite(mp_traits, paste(output.dir, 'mp_traits_log_nm.csv', sep = ""))
  fwrite(gc_genotype, paste(output.dir, 'gc_genotype_log_nm.csv', sep = ""))
  
  # Step 5: Set the vector of phenotypic traits
  traits <- setdiff(traits, "ID")
  
  # Step 6: Perform a crossing between phenotypes and genotypes used for the analysis
  qtl.data<-read.cross("csvs",
                       dir = "", 
                       genfile = paste(output.dir, 'gc_genotype_log_nm.csv', sep = ""), 
                       phefile = paste(output.dir, 'mp_traits_log_nm.csv', sep = ""),
                       genotypes=c("M","U"),
                       na.string=c("NA"), 
                       estimate.map=T)
  
  # Step 7: Calculate conditional genotype probabilities
  qtl.data<-calc.genoprob(qtl.data, step=2, off.end=0, error.prob=0,
                          map.function=c("haldane"), stepwidth=c("fixed"))
  
  # Step 8: Set the number of permutations and perform a permutation test
  n.perm <- 1000
  out.perm<-scanone(qtl.data, pheno.col=traits, model="normal",
                    method="hk", use="all.obs", addcovar=NULL, intcovar=NULL, weights=NULL,
                    upper=FALSE, ties.random=FALSE, start=NULL, maxit=4000,
                    tol=1e-4, n.perm=n.perm, perm.Xsp=FALSE, perm.strata=NULL, verbose=T,
                    batchsize=250, n.cluster=1)
  
  # Step 9: Perform a mapping procedure
  out.inter<-scanone(qtl.data, pheno.col=traits, model="normal",
                     method="hk", use="all.obs", addcovar=NULL, intcovar=NULL, weights=NULL,
                     upper=FALSE, ties.random=FALSE, start=NULL, maxit=4000,
                     tol=1e-4, perm.Xsp=FALSE, perm.strata=NULL, verbose=T,
                     batchsize=250, n.cluster=1)
  
  # Step 10: Make a genotype dataframe (AA -> 1, AB -> 2)
  prob.geno<-cbind(
    round(qtl.data$geno$`1`$prob[,,"AB"], 0) + 1,
    round(qtl.data$geno$`2`$prob[,,"AB"], 0) + 1,
    round(qtl.data$geno$`3`$prob[,,"AB"], 0) + 1,
    round(qtl.data$geno$`4`$prob[,,"AB"], 0) + 1,
    round(qtl.data$geno$`5`$prob[,,"AB"], 0) + 1)
  
  # Step 11: Pull out the phenotypes from crossing 
  pheno<-pull.pheno(qtl.data)
  
  # Step 12: Calculate effect direction based on the linear model between the vector of epigenotypes and phenotypes for
  #          a given phenotype - marker association (-1 - negative, +1 - positive)
  effect.direction<-matrix(,nrow=dim(out.inter)[1], ncol=dim(out.inter)[2]-2)
  
  for (a in 1:length(traits))
  {
    cat("Start: effect direction: ", traits[a], "\n")
    pheno.temp<-pheno[,which(colnames(pheno) == traits[a])]
    
    for (b in 1:ncol(prob.geno))
    {
      geno.temp<-prob.geno[,b]
      test.out<-summary(lm(pheno.temp~geno.temp))
      test.out<-test.out$coefficients[2]
      
      if (test.out > 0)
      {effect.direction[b,a]<-  1}
      if (test.out <=0)
      {effect.direction[b,a]<- -1}
    }
  }
  colnames(effect.direction)<-traits
  rownames(effect.direction)<-rownames(out.inter)
  
  # Step 13: Make a dataset with cM position for further analysis
  chr.numbers<-as.numeric(names(nmar(qtl.data)))
  fullmap.info<-NULL
  for (a in 1:length(chr.numbers))
  {
    temp1<-qtl.data$geno[[chr.numbers[a]]]
    temp2<-temp1$prob
    temp3<-attr(temp2,"map")
    chr.vec<-rep(chr.numbers[a],length(temp3))
    chr.all<-data.frame(names(temp3),chr.vec,temp3)
    colnames(chr.all)<-c("marker","chr","cM.pos")
    rownames(chr.all)<-NULL
    fullmap.info[a]<-list(chr.all)
  }
  names(fullmap.info)<-as.numeric(chr.numbers) 
  
  
  # Step 14: Save all of the datasets
  ## Reading out the effect direction
  dput(effect.direction, paste(output.dir, "QTL-direction-data-all-traits_log_nm.Rdata", sep = ''))
  ## Reading out the reference data
  dput(qtl.data, paste(output.dir, "REF-data-all-traits_log_nm.Rdata", sep = ''))
  ## Reading out the scan object
  dput(out.inter, paste(output.dir, "QTL-mapping-all-traits_log_nm.Rdata", sep = ''))
  ## Reading out the permutation results
  dput(out.perm, paste(output.dir, "PERM-all-traits_log_nm.Rdata", sep= ''))
  ## Reading out the full map info
  dput(fullmap.info, paste(output.dir, "SCANPOSINFO-all-traits_log_nm.Rdata", sep = ''))
}

##### getQTLpeaks_new function #####
# getQTLpeaks_new - select a significant peaks with a given alpha level
getQTLpeaks_new<-function(input.dir, marker_dir, alpha = 0.05)
{
  # Step 1: Read required dataset for reference genome, permutation data, mapping data, marker positions and effect direcion
  ref.data<-dget(paste(input.dir, "REF-data-all-traits_log_nm.Rdata", sep="")) #the traits were inputted incorrectly
  perm.data<-dget(paste(input.dir, "PERM-all-traits_log_nm.Rdata", sep="")) #permutation results dataframe
  qtl.data<-dget(paste(input.dir, "QTL-mapping-all-traits_log_nm.Rdata", sep="")) #mapping results dataframe
  qtl.direction<-dget(paste(input.dir, "QTL-direction-data-all-traits_log_nm.Rdata", sep="")) #effect direction dataset
  
  marker.pos<-read.csv(paste(marker_dir, sep = ""), header=T) #marker positions dataset
  
  # Step 2: Select permutation threshold for each trait separately
  perm.thresh<-summary(perm.data, alpha=alpha)
  # Step 3: Pull out map necessary datasets
  map.info<-pull.map(ref.data, as.table=T) #pull out genotype map (with cM positions)
  chr.in<-as.character(unique(qtl.data[,1])) #pull out chromosome unique set
  trait.in<-colnames(perm.thresh) #pull out phenotypic traits
  peaks.collect<-NULL
  peaks.new<-NULL
  # Step 4: For each phenotypic trait for each chromosome separately find a peak
  for (a in 1:length(trait.in))
  {
    perm.temp<-perm.thresh[, which(colnames(perm.thresh) == trait.in[a])]
    trait.collect<-NULL
    peaks.collect<-NULL
    
    for (b in 1:length(chr.in))
    {
      qtl.data.temp<-qtl.data[which(as.character(qtl.data$chr) == chr.in[b]), c(1,2, which(colnames(qtl.data) == trait.in[a]))]
      qtl.data.temp<-qtl.data.temp[which.max(qtl.data.temp[,3]), ]
      colnames(qtl.data.temp)[3]<-"LOD"
      peaks.collect<-rbind(peaks.collect, qtl.data.temp)
      trait.collect<-c(trait.collect, trait.in[a])
    }
    
    peaks.collect<-data.frame(peaks.collect, trait.collect)
    peaks.collect$marker<-rownames(peaks.collect)
    rownames(peaks.collect)<-NULL
    peaks.new<-rbind(peaks.new, peaks.collect[which(peaks.collect[,3] >= perm.temp),])
    
  }
  colnames(peaks.new)[4]<-"trait"
  colnames(peaks.new)[5]<-"peak.marker"
  flanking.marker<-NULL
  for (a in 1:nrow(peaks.new))
  {
    flanking.marker[a]<-find.marker(ref.data, peaks.new[a,1], peaks.new[a,2])
  }
  peaks.new$nearest.marker<-flanking.marker
  pos.bp<-NULL
  for (a in 1:nrow(peaks.new))
  {
    
    mp.temp<-marker.pos[which(as.character(marker.pos[,1]) == peaks.new[a, "nearest.marker"]),]
    pos.bp[a]<-0.5*(mp.temp[,3] + mp.temp[,4])
  }
  peaks.new$pos.bp<-pos.bp
  peaks.new<-peaks.new[, c(5, 6, 1, 2, 7, 3, 4)]
  for (e in 1:nrow(peaks.new))
  {
    temp<-qtl.direction[which(rownames(qtl.direction) == as.character(peaks.new$peak.marker[e])),]
    temp<-as.numeric(temp[which(names(temp) == peaks.new$trait[e])])
    peaks.new$increasing.epiallel[e]<-ifelse(temp == 1, "U", "M")
    
    int.temp<-lodint(qtl.data, chr=peaks.new$chr[e], qtl.index, drop=1, 
                     lodcolumn=which(colnames(qtl.data) == peaks.new$trait[e])-2, expandtomarkers=TRUE)
    peaks.new$CI.lower.marker[e]<-rownames(int.temp)[1]
    peaks.new$CI.upper.marker[e]<-rownames(int.temp)[3]
    peaks.new$CI.lower.cM[e]<-int.temp[1,2]
    peaks.new$CI.upper.cM[e]<-int.temp[3,2]
  }
  colnames(peaks.new)[4]<-"pos.cM"
  fwrite(peaks.new, file = paste(input.dir, 'info_peaks_log_nm.csv', sep = ""))
  return(peaks.new)
}

##### qtlpeaks.plot function #####
# qtlpeaks.plot - plot a LOD score thresholds for the peaks selected in getQTLpeaks_new
qtlpeaks.plot <- function(input.dir, marker_dir)
{
  # Step 1. Load required datasets
  out.inter <- dget(paste(input.dir, "QTL-mapping-all-traits_log_nm.Rdata", sep = "")) #data after scanone function with LOD values
  perm.data <- dget(paste(input.dir, "PERM-all-traits_log_nm.Rdata", sep="")) #permutation data used for selecting threshold for LOD values as significant
  marker_data <- fread(marker_dir) #marker bp positions
  #information about start-end positions for each chromosome
  chromosome_info <- data.frame(chr = seq(1,5),
                                start = rep(1,5),
                                end = c(30427671, 19698289, 23459830, 18585056, 26975502))
  #information about start-end pericentrometric regions for each chromosome
  pericentrometric_info <- data.frame(chr = seq(1,5),
                                      start = c(9698788, 0, 7298763, 0, 6999996),
                                      end = c(19897560, 9492918, 18289014, 9150003, 17332770))
  add_chromosome_bp <- rep(0,5) 
  for (i in 2:5) #start with 2nd chromosome
  {
    chromosome_info_i <- chromosome_info$end[1:i-1]
    add_chromosome_bp[i] <- sum(chromosome_info_i)
  }
  pericentrometric_info_bp <- data.frame(chr.trait = pericentrometric_info$chr,
                                         chr.marker = pericentrometric_info$chr,
                                         start = add_chromosome_bp + pericentrometric_info$start,
                                         end = add_chromosome_bp + pericentrometric_info$end)
  add_chromosome_bp <- data.frame(add.x = add_chromosome_bp,
                                  add.y = max(add_chromosome_bp),
                                  chr.marker = seq(1,5))
  chromosome_info_bp <- data.frame(chr.marker = chromosome_info$chr,
                                   x = 0)
  for (element in 1:5)
  {
    if (element - 1 == 0)
    {
      chromosome_info_bp_element <- chromosome_info$end[1]
    } else {
      chromosome_info_bp_element <- chromosome_info$end[1:element]
    }
    chromosome_info_bp$x[element] <- sum(chromosome_info_bp_element)
  }
  
  # Step 2. Data manipulation
  colnames(out.inter) <- c("chr", "position", colnames(out.inter)[-c(1,2)]) #name the columns (optional)
  out.inter$marker <- rownames(out.inter) #add information about marker as the column
  lods_info <- melt(out.inter, id.vars=c("chr", "position", "marker")) #extract from the scanone object the chromosome, position and marker
  
  # Step 3. Calculate the threshold using permutation dataset
  perm.thresh<-summary(perm.data, alpha=0.05) #calculate threshold
  perm.thresh <- data.frame(trait = colnames(out.inter)[-c(1,2, length(colnames(out.inter)))],
                            threshold = as.vector(perm.thresh)) #make a dataframe (optional)
  
  # Step 4. Calculate on the basis of the treshhold new LOD value
  lods_info$new_value <- NA
  for (element in 1:dim(lods_info)[1])
  {
    lods_info$new_value[element] <- lods_info$value[element] / perm.thresh$threshold[which(perm.thresh$trait == lods_info$variable[element])]
  }
  
  # Step 5. Modify pericentrometric file by adding info about y axes in the plot
  peri_info_bp <- data.frame(chr = pericentrometric_info$chr, 
                             xmin = pericentrometric_info$start, 
                             xmax = pericentrometric_info$end,
                             ymin = - Inf, ymax = +Inf)
  
  # Step 6. Delete from the dataset the information about pseudomarkers and
  #         add the information about real markers position in bp
  lods_info_nopseudo <-   lods_info[lods_info$marker %like% "MM",] #delete pseudomarkers
  marker_data$mid <- 0.5*(marker_data$start + marker_data$end) #calculate mid value position for markers
  lods_info_nopseudo$position_bp <- NA 
  for (i in 1:dim(lods_info_nopseudo)[1]) #add information about mid-marker positions
  {
    lods_info_nopseudo$position_bp[i] <- marker_data$mid[which(marker_data$marker == lods_info_nopseudo$marker[i])]
  }
  
  # Step 7. Add values from bp
  for (element in 1:dim(lods_info_nopseudo)[1])
  {
    lods_info_nopseudo$position_bp[element] <- lods_info_nopseudo$position_bp[element] + add_chromosome_bp$add.x[which(add_chromosome_bp$chr.marker == lods_info_nopseudo$chr[element])]
  }
  pericentrometric_info_bp$chr <- pericentrometric_info_bp$chr.marker
  chromosome_info_bp$chr <- chromosome_info_bp$chr.marker
  
  # Step 8. Make a final plot
  alog_nm <- ggplot(data=lods_info_nopseudo) +
    geom_line(aes(x=position_bp, y=new_value, color = variable), size = 0.6) +
    ggtitle('epiRILs phenotypes') +
    geom_point(data=chromosome_info_bp, aes(x=x, y=0), color = 'white') +
    geom_rect(data=pericentrometric_info_bp,
              aes(xmin=start, xmax=end, 
                  ymin=-Inf, ymax=+Inf), fill="grey", alpha=0.3) +
    facet_grid(~chr, scales = "free_x", space = 'free_x',
               labeller=label_wrap_gen(width = 10, multi_line = TRUE)) +
    theme_classic(base_size = 14) +
    xlab('') + ylab('LOD value') + 
    labs(color = "Feature") +
    geom_hline(yintercept=1, color = "red", linetype = "dashed") +
    scale_x_continuous(breaks = c(1*10^7, 10*10^7)) +
    theme(legend.position = "none") + 
    theme(strip.background = element_blank())
  return(alog_nm)
}

##### transcis.plot function #####
# transcis.plot - plot a LOD score thresholds for the peaks selected in getQTLpeaks_new
transcis.plot <- function(input.dir, marker_dir, positions_dir)
{
  # Step 1: Load required datasets
  peaks_NAD_data <- fread(paste(input.dir, 'info_peaks_log_nm.csv', sep = '')) #information about peaks from NAD regions
  marker_data <- fread(marker_dir) #marker bp positions
  phenotype_positions <- fread(positions_dir)
  positions_NAD_data <- unique(peaks_NAD_data$trait)
  positions_NAD_data <- data.frame(ID = unique(peaks_NAD_data$trait),
                                   seqnames = NA,
                                   position = NA)
  for (element in 1:dim(positions_NAD_data)[1])
  {
    ID_element <- positions_NAD_data$ID[element]
    positions_NAD_data$seqnames[element] <- phenotype_positions$chr[which(phenotype_positions$trait == ID_element)]
    positions_NAD_data$position[element] <- 0.5*(phenotype_positions$start[which(phenotype_positions$trait == ID_element)] + phenotype_positions$end[which(phenotype_positions$trait == ID_element)])
  }
  #information about start-end positions for each chromosome
  chromosome_info <- data.frame(chr = seq(1,5),
                                start = rep(1,5),
                                end = c(30427671, 19698289, 23459830, 18585056, 26975502))
  #information about start-end pericentrometric regions for each chromosome
  pericentrometric_info <- data.frame(chr = seq(1,5),
                                      start = c(9698788, 0, 7298763, 0, 6999996),
                                      end = c(19897560, 9492918, 18289014, 9150003, 17332770))
  
  # Step 2: Create dataframe corresponding to peaks data, whereas:
  # a) marker id, NAD region id
  # b) start, mid, end positions for markers
  # c) positions for NAD regions
  # will be saved for firther plotting
  peaks_NAD_data_new <- data.frame(matrix(vector(), ncol = 10, nrow = dim(peaks_NAD_data)[1]))
  colnames(peaks_NAD_data_new) <- c('chr.marker', 'chr.trait',
                                    'marker.id', 'trait.id', 
                                    'mid.marker.bp', 'mid.trait.bp',
                                    'start.marker.bp', 'start.trait.bp',
                                    'end.marker.bp', 'end.trait.bp')
  
  # Step 3: (OPTIONAL, only for the approach with overall positions, not for each chromosome)
  #         calculate the adding value for each chromosome position
  add_chromosome_bp <- rep(0,5) 
  for (i in 2:5) #start with 2nd chromosome
  {
    chromosome_info_i <- chromosome_info$end[1:i-1]
    add_chromosome_bp[i] <- sum(chromosome_info_i)
  }
  
  # Step 4: Fill the data for peaks_NAD_data
  for (ind in 1:dim(peaks_NAD_data)[1])
  {
    peaks_NAD_data_new$marker.id[ind] <- peaks_NAD_data$nearest.marker[ind] #marker ID
    peaks_NAD_data_new$trait.id[ind] <- peaks_NAD_data$trait[ind] #NAD ID
    peaks_NAD_data_new$chr.marker[ind] <- peaks_NAD_data$chr[ind] #marker chromosome
    peaks_NAD_data_new$chr.trait[ind] <- positions_NAD_data$seqnames[which(positions_NAD_data$ID == peaks_NAD_data_new$trait.id[ind])] #NAD chromosome
    
    # start-end marker positions
    peaks_NAD_data_new$start.marker.bp[ind] <- marker_data$start[which(marker_data$marker == peaks_NAD_data_new$marker.id[ind])]
    peaks_NAD_data_new$end.marker.bp[ind] <- marker_data$end[which(marker_data$marker == peaks_NAD_data_new$marker.id[ind])]
    # start-end NAD positions
    peaks_NAD_data_new$start.trait.bp[ind] <- positions_NAD_data$position[which(positions_NAD_data$ID == peaks_NAD_data_new$trait.id[ind])]
    peaks_NAD_data_new$end.trait.bp[ind] <- positions_NAD_data$position[which(positions_NAD_data$ID == peaks_NAD_data_new$trait.id[ind])]
    
    # add the corresponding bp for each chromosome to make it real position
    peaks_NAD_data_new$start.marker.bp[ind] <- add_chromosome_bp[peaks_NAD_data_new$chr.marker[ind]] + peaks_NAD_data_new$start.marker.bp[ind]
    peaks_NAD_data_new$end.marker.bp[ind] <- add_chromosome_bp[peaks_NAD_data_new$chr.marker[ind]] + peaks_NAD_data_new$end.marker.bp[ind]
    peaks_NAD_data_new$start.trait.bp[ind] <- add_chromosome_bp[peaks_NAD_data_new$chr.trait[ind]] + peaks_NAD_data_new$start.trait.bp[ind]
    peaks_NAD_data_new$end.trait.bp[ind] <- add_chromosome_bp[peaks_NAD_data_new$chr.trait[ind]] + peaks_NAD_data_new$end.trait.bp[ind]
  }
  # calculate mid value by adding start to end and dividing to it by 2
  peaks_NAD_data_new$mid.marker.bp <- 0.5*(peaks_NAD_data_new$start.marker.bp + peaks_NAD_data_new$end.marker.bp)
  peaks_NAD_data_new$mid.trait.bp <- 0.5*(peaks_NAD_data_new$start.trait.bp + peaks_NAD_data_new$end.trait.bp)
  
  # Step 5: Pericentrometric regions real value by adding the corresponding bp for each chromosome to make it real bp position
  pericentrometric_info_bp <- data.frame(chr.trait = pericentrometric_info$chr,
                                         chr.marker = pericentrometric_info$chr,
                                         start = add_chromosome_bp + pericentrometric_info$start,
                                         end = add_chromosome_bp + pericentrometric_info$end)
  add_chromosome_bp <- data.frame(add.x = add_chromosome_bp,
                                  add.y = max(add_chromosome_bp),
                                  chr.marker = seq(1,5))
  
  # Step 6: To fix the problems with y scale: make it in Mbp
  # x/y ticks for each chromosome
  scale_mbp_ticks <- c(0, 
                       30.42,
                       50.13, 
                       73.59, 
                       92.19)
  
  chromosome_info_bp <- data.frame(chr.marker = chromosome_info$chr,
                                   x = 0)
  for (element in 1:5)
  {
    if (element - 1 == 0)
    {
      chromosome_info_bp_element <- chromosome_info$end[1]
    } else {
      chromosome_info_bp_element <- chromosome_info$end[1:element]
    }
    chromosome_info_bp$x[element] <- sum(chromosome_info_bp_element)
  }
  
  # Step 7: Make the plot using ggplot2 envo
  tclog_nm <- ggplot(data = peaks_NAD_data_new) +
    geom_point(aes(x=mid.marker.bp, y=mid.trait.bp)) + #scatterplot of mid-val marker (bp) and trait (Mbp)
    geom_point(data = chromosome_info_bp, aes(x = x, y = 0), color = 'white') +
    geom_rect(data=pericentrometric_info_bp,
              aes(xmin=start, 
                  xmax=end, 
                  ymin=start, 
                  ymax=end), 
              fill="yellow", alpha=0.6) + #make pericentrometric regions in yellow
    geom_hline(yintercept = add_chromosome_bp$add.x[2:5], linetype="dashed", color = "red") +
    xlab('Marker position') + ylab('Trait position') + #set scale in Mbp
    facet_grid(~chr.marker, scale = 'free_x',
               labeller=label_wrap_gen(width = 10, multi_line = TRUE)) +
    scale_y_continuous(labels = paste0(scale_mbp_ticks, " "),
                       breaks = 10^6 * scale_mbp_ticks) + 
    scale_x_continuous(labels = paste0(scale_mbp_ticks, " "),
                       breaks = 10^6 * scale_mbp_ticks,
                       guide = guide_axis(check.overlap = TRUE)) +
    theme_classic(base_size = 12) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(strip.background = element_blank())
  return(tclog_nm)
}
