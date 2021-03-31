### Monogenic IBD Variant Prioritization
##  Written by Daniel Mulder, February 2021
#   Part 1 Family member vcf pipeline


# To use this script, the following user changes are needed:
# 1. change "/path/to/files" on the line below to the local path to the folder containing only the family member VCFs (not proband VCFs)
input_vcf_dir <- "path/to/files"


library(tidyverse)

files <- list.files(input_vcf_dir)

vcf_header_count <- function(in_file) {
  # Returns the number of lines that the header section contains
  # This function reads in data line-by-line, so it should remain fast and lightweight in most cases
  # Assumption: The given file is non-empty
  file_handle <- file(in_file, "r")
  header_size <- 0
  while (TRUE) {
    line <- readLines(file_handle, n = 1)
    if (grepl("^#", line)) {
      header_size <- header_size + 1
    } else {
      break
    }
  }
  return(header_size)
}

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
  vcf.now <- paste0(input_vcf_dir, files[i])
  skip <- vcf_header_count(vcf.now) - 1
  fam_memb_processing(vcf.now, paste0(vcf.now, ".txt"), skip)
}
