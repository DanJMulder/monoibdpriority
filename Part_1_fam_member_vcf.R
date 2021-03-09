### Monogenic IBD Variant Prioritization
##  Written by Daniel Mulder, February 2021
#   Part 1 Family member vcf pipeline


# To use this script, the following user changes are needed:
# 1. change "/path/to/files" on the line below to the local path to the folder containing only the family member VCFs (not proband VCFs)
input_vcf_dir <- "path/to/files"
# 2. change the "number_of_header_rows" variable on the line below to how many header rows in your VCFs (that need to be skipped)
skip <- number_of_header_rows


library(tidyverse)
library(glue) #for interpreting string literals

files <- list.files(glue("{input_vcf_dir}"))

fam_memb_processing <- function(input, output, num_header) {
  #load in vcf
  vcf <- read_delim(input,
                    "\t", escape_double = FALSE, trim_ws = TRUE,
                    skip = num_header, col_types = cols('#CHROM' = col_factor()))
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

  write_csv(vcf, output)

  rm(vcf)
}

#loop through input files
for (i in seq_along(files)) {
  fam_memb_processing(files[i], files[i], skip)
}
