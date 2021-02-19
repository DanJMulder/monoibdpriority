### Monogenic IBD Variant Prioritization
##  Written by Daniel Mulder, February 2021
#   Part 1 Family member vcf pipeline

# TO USE THIS SCRIPT, the following user changes are needed:
  # 1. input the local path to the folder containing the family member VCFs to process (to replace "path/to/files" below)
    # the folder should contain only the VCF files to be processed
    # the folder should not contain the proband VCF files which are processed in a separate script
  # 2. change the "skip" argument in the read_delim() function below depending on how many header rows need to be skipped in the VCF


library(tidyverse)

files <- list.files("path/to/files")

fam_memb_processing <- function(input) {
  #load in vcf
  vcf <- read_delim(input, 
                    "\t", escape_double = FALSE, trim_ws = TRUE, 
                    skip = 296, col_types = cols('#CHROM' = col_factor()))
  names(vcf)[names(vcf) == "#CHROM"] <- "chr"
  vcf$chr <- as.factor(gsub('chr', '', vcf$chr))
  
  #creating address column to enable merge with proband vcf
  vcf$address <- paste0(vcf$chr, ":", vcf$POS)
  vcf$address <- paste0(vcf$address, ":", vcf$REF)
  vcf$address <- paste0(vcf$address, ":", vcf$ALT)
  
  #create het_or_hom from GENOTYPE column########################
  names(vcf)[10] <- "het_or_hom"
  het_or_hom <- vcf$het_or_hom
  vcf$het_or_hom <- substring(het_or_hom, 1, 3)
  vcf$het_or_hom <- as.factor(vcf$het_or_hom)
  
  vcf$het_or_hom <- gsub('0/1', 'het', vcf$het_or_hom)
  vcf$het_or_hom <- gsub('1/1', 'hom_alt', vcf$het_or_hom)
  vcf$het_or_hom <- gsub('1/2', 'hom_alt', vcf$het_or_hom)
  vcf$het_or_hom <- as.factor(vcf$het_or_hom)
  
  rm(het_or_hom)
  
  vcf <- subset(vcf, select = c(address, het_or_hom))
  
  write_csv(vcf, input)
  
  rm(vcf)
}
  
#loop through input files
for (i in 1:length(files)) {
  fam_memb_processing(files[i])
}
