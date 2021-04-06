#!/usr/bin/env Rscript
### Monogenic IBD Variant Prioritization
##  Written by Daniel Mulder, February 2021

## Input files
# OPTION 1: Provide a samplesheet and VCF directory
vcf_dir <- ""
family_units <- ""
# OPTION 2: Specify trio
proband_fn <- ""
proband_sex <- ""
paternal_fn <- ""
paternal_affected <- ""
maternal_fn <- ""
maternal_affected <- ""

# Output
output_dir <- ""

# Options
num.cpu <- 8

# Tools
# 1. change to the folder containing the dbNSFP4.1a java applet (search_dbNSFP41a.class)
path_to_dbNSFP4.1a <- "/path/to/dbNSFP4.1a"
# 2. change the lines below to the folder containing the annovar libraries and db
path_to_annovar <- "/path/to/annovar"
path_to_annovar_db <- "/path/to/annovar/humandb"
# 3. change the line below to the local path to the folder containing the CADD scripts
path_to_cadd <- "/path/to/CADD-scripts-CADD1.6"
# 4. change the line below to the local path to the LOEUF table downloaded from gnomAD
path_to_LOEUF_table <- "/path/to/gnomad.v2.1.1.lof_metrics.by_gene.txt"
# 5. change the line below to the directory containing the 4 bed files for the repository
path_to_monoibdpriority <- "/path/to/monoibdpriority"

library(plyr)
library(tidyverse)
library(data.table)

# Helper function to determine if a parameter is truthy
is.specified <- function(value) {
  return(!is.na(value) &&
           (value != "NA" || value != "unknown") &&
           length(value) > 0)
}

# Helper function to help determine how many lines a VCF header contains, includes the column names
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
  close(file_handle)
  return(header_size)
}


fam_memb_processing <- function(input, output) {
  #load in vcf
  vcf <- read_delim(input,
                    "\t", escape_double = FALSE, trim_ws = TRUE,
                    skip = vcf_header_count(input) - 1, col_types = cols('#CHROM' = col_factor()))
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


proband_annotation <- function(input, output_dir) {
  # Part 2 Script Structure Overview

  # annotate each variant with:
  # 1. dbNSFP (v4.1a) using dbNSFP java applet
  # 2. MAF (from gnomAD) and RefSeq using annovar
  # 3. CADD score (phred, v 1.6)
  # 4. LOEUF (from gnomAD) using offline table
  # 5. Gene_lists (monoIBD, IBD_GWAS, PID, CDG) from local BED files
  # 6. Het or Homozygous from alignment call (in VCF from the beginning)
  print(paste("Processing proband:", input))
  # Step 1. Load in and prepare the VCF file ####

  # Create a string variable w just the input file name (as the patient ID)
  patient_id <- gsub('.vcf', '', basename(input))
  output_dir <- normalizePath(output_dir)

  # Load the vcf into R
  vcf <- read_delim(input, "\t", escape_double = FALSE, trim_ws = TRUE,
                    skip = vcf_header_count(input) - 1, col_types = cols('#CHROM' = col_factor()))
  names(vcf)[names(vcf) == "#CHROM"] <- "Chr"

  # Remove the 'chr' from the start of the chromosome column
  vcf$Chr <- gsub('chr', '', vcf$Chr)

  # Step 2. dbNSFP annotation ####

  # Select only the necessary columns for dbNSFP annotation
  vcf_for_dbnsfp <- subset(vcf, select = c(Chr, POS, REF, ALT))

  vcf_for_dbnsfp$Chr <- as.factor(vcf_for_dbnsfp$Chr)

  # Save smaller vcf dataframe to folder so dbNSFP java applet can use it
  input_in <- file.path(output_dir, paste0(patient_id, ".in"))
  dbnsfp_output <- file.path(output_dir, paste0(patient_id, ".out"))
  write_tsv(vcf_for_dbnsfp, input_in, na = ".", quote = FALSE)
  rm(vcf_for_dbnsfp)

  # Running dbNSFP java applet via command line to add dbNSFP annotations
  # Note, changing the directory here does not affect R's current working directory.
  # We need to use a slightly different command in order to run without any GUIs
  dbnsfp_phrase <- as.character(paste("module load java; cd", path_to_dbNSFP4.1a,
                                      "; java search_dbNSFP41a -i", input_in, "-o", dbnsfp_output, "-v hg38"))
  if (!file.exists(dbnsfp_output)) {
    print(paste0("Running command:\n", dbnsfp_phrase))
    system(dbnsfp_phrase)
  } else {
    print(paste("Skipping dbnsfp because", dbnsfp_output, "exists."))
  }

  # load the dbnsfp result into R and add an end column
  vcf_post_dbnsfp <- read_delim(dbnsfp_output,
                                "\t",
                                escape_double = FALSE,
                                col_names = TRUE,
                                trim_ws = TRUE,
                                na = ".",
                                guess_max = 200000,
                                col_types = cols('#chr' = col_factor(), 'hg19_chr' = col_factor(),
                                                 'hg18_chr' = col_factor()))

  # reorder the columns to match the original vcf order (to prevent confusion b/c lots of chr/start columns in dbnsf output for different ref genomes)
  vcf_post_dbnsfp <- vcf_post_dbnsfp[, c(1:7, 12:463)]

  # rename/clean up the "vcf" columns so they will match the "vcf_post_dbnsfp" columns
  names(vcf)[names(vcf) == "Chr"] <- "chr"
  names(vcf)[names(vcf) == "POS"] <- "start"
  vcf <- select(vcf, !ID)
  names(vcf)[names(vcf) == "REF"] <- "ref"
  names(vcf)[names(vcf) == "ALT"] <- "alt"
  names(vcf)[9] <- "genotype"

  # remove the vcf_post_dbnsfp columns not needed (steps i to iii)

  # clean up column names for ease of use
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "#chr"] <- "chr"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "pos(1-based)"] <- "start"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "M-CAP_pred"] <- "M_CAP_pred"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "M-CAP_rankscore"] <- "M_CAP_rankscore"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "LIST-S2_rankscore"] <- "LIST_S2_rankscore"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "LIST-S2_pred"] <- "LIST_S2_pred"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "fathmm-MKL_coding_rankscore"] <- "fathmm_MKL_coding_rankscore"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "fathmm-MKL_coding_pred"] <- "fathmm_MKL_coding_pred"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "fathmm-XF_coding_rankscore"] <- "fathmm_XF_coding_rankscore"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "fathmm-XF_coding_pred"] <- "fathmm_XF_coding_pred"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "Eigen-phred_coding"] <- "Eigen_phred_coding"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "Eigen-PC-phred_coding"] <- "Eigen_PC_phred_coding"
  names(vcf_post_dbnsfp)[names(vcf_post_dbnsfp) == "H1-hESC_fitCons_rankscore"] <- "H1_hESC_fitCons_rankscore"

  # make a list of the desired columns
  desired_columns <- c("chr",
                       "start",
                       "ref",
                       "alt",
                       "SIFT_converted_rankscore",
                       "SIFT_pred",
                       "SIFT4G_converted_rankscore",
                       "SIFT4G_pred",
                       "Polyphen2_HDIV_rankscore",
                       "Polyphen2_HDIV_pred",
                       "Polyphen2_HVAR_rankscore",
                       "Polyphen2_HVAR_pred",
                       "LRT_converted_rankscore",
                       "LRT_pred",
                       "MutationTaster_converted_rankscore",
                       "MutationTaster_pred",
                       "MutationAssessor_rankscore",
                       "MutationAssessor_pred",
                       "FATHMM_converted_rankscore",
                       "FATHMM_pred",
                       "PROVEAN_converted_rankscore",
                       "PROVEAN_pred",
                       "VEST4_rankscore",
                       "MetaSVM_rankscore",
                       "MetaSVM_pred",
                       "MetaLR_rankscore",
                       "MetaLR_pred",
                       "M_CAP_rankscore",
                       "M_CAP_pred",
                       "REVEL_rankscore",
                       "MutPred_rankscore",
                       "MVP_rankscore",
                       "MPC_rankscore",
                       "PrimateAI_rankscore",
                       "PrimateAI_pred",
                       "DEOGEN2_rankscore",
                       "DEOGEN2_pred",
                       "BayesDel_addAF_rankscore",
                       "BayesDel_addAF_pred",
                       "BayesDel_noAF_rankscore",
                       "BayesDel_noAF_pred",
                       "ClinPred_rankscore",
                       "ClinPred_pred",
                       "LIST_S2_rankscore",
                       "LIST_S2_pred",
                       "Aloft_pred",
                       "CADD_phred",
                       "DANN_rankscore",
                       "fathmm_MKL_coding_rankscore",
                       "fathmm_MKL_coding_pred",
                       "fathmm_XF_coding_rankscore",
                       "fathmm_XF_coding_pred",
                       "Eigen_phred_coding",
                       "Eigen_PC_phred_coding",
                       "GenoCanyon_rankscore",
                       "integrated_fitCons_rankscore",
                       "GM12878_fitCons_rankscore",
                       "H1_hESC_fitCons_rankscore",
                       "HUVEC_fitCons_rankscore",
                       "LINSIGHT_rankscore",
                       "GERP++_RS_rankscore",
                       "phyloP100way_vertebrate_rankscore",
                       "phyloP30way_mammalian_rankscore",
                       "phyloP17way_primate_rankscore",
                       "phastCons100way_vertebrate_rankscore",
                       "phastCons30way_mammalian_rankscore",
                       "phastCons17way_primate_rankscore",
                       "SiPhy_29way_logOdds_rankscore",
                       "bStatistic_converted_rankscore"
  )

  # get down to the desired (useful) columns
  vcf_post_dbnsfp <- subset(vcf_post_dbnsfp, select = desired_columns)

  # creating "address" column as a primary key for each data frame to enable merge
  vcf$address <- paste0(vcf$chr, ":", vcf$start)
  vcf$address <- paste0(vcf$address, ":", vcf$ref)
  vcf$address <- paste0(vcf$address, ":", vcf$alt)

  vcf_post_dbnsfp$address <- paste0(vcf_post_dbnsfp$chr, ":", vcf_post_dbnsfp$start)
  vcf_post_dbnsfp$address <- paste0(vcf_post_dbnsfp$address, ":", vcf_post_dbnsfp$ref)
  vcf_post_dbnsfp$address <- paste0(vcf_post_dbnsfp$address, ":", vcf_post_dbnsfp$alt)

  # merge the "vcf_post_dbnsfp" dataframe with the "vcf" dataframe (so that the many variants not called by dbnsfp and the het/homozygous status from "vcf" is preserved going forward)
  vcf <- merge(vcf, vcf_post_dbnsfp, by = "address", all.x = TRUE)
  rm(vcf_post_dbnsfp)

  # clean up the column names and remove redundant columns from the merge
  names(vcf)[names(vcf) == "chr.x"] <- "chr"
  names(vcf)[names(vcf) == "start.x"] <- "start"
  names(vcf)[names(vcf) == "ref.x"] <- "ref"
  names(vcf)[names(vcf) == "alt.x"] <- "alt"
  vcf <- subset(vcf, select = -c(chr.y, start.y, ref.y, alt.y))
  vcf <- vcf[, c(2:75, 1)]

  # Step 3. Add maf and RefSeq annotations via annovar ####

  # Calculate/add the "end" column so appropriate for annovar input
  vcf$variant_length <- nchar(vcf$alt) - 1
  vcf$end <- vcf$start + vcf$variant_length
  vcf <- vcf[, c(1, 2, 77, 3:75)]

  annovar_save <- file.path(output_dir, paste0(patient_id, ".txt"))
  annovar_out_prefix <- file.path(output_dir, patient_id)
  annovar_out <- paste0(annovar_out_prefix, ".hg38_multianno.txt")

  write_delim(vcf, annovar_save, delim = "\t", na = ".", col_names = TRUE, quote_escape = FALSE)

  annovar_script <- file.path(path_to_annovar, "table_annovar.pl")
  annovar_phrase <- as.character(paste("perl", annovar_script, annovar_save, path_to_annovar_db,
                                       "-buildver hg38 -out", annovar_out_prefix, "-remove",
                                       "-protocol gnomad211_exome,refGene -operation f,g -nastring .",
                                       "--thread", num.cpu))

  # Annotation w MAF (gnomAD 2.1.1) and refGene via annovar via terminal command
  if (!file.exists(annovar_out)) {
    print(paste("Running command:", annovar_phrase, sep = "\n"))
    system(annovar_phrase)
  } else {
    print(paste("Skipping annovar because", annovar_out, "exists."))
  }

  # load in the results of annovar (will also include dbNSFP annotations)
  post_annovar <- read_delim(annovar_out,
                             "\t",
                             escape_double = FALSE,
                             trim_ws = TRUE,
                             col_types = cols('Chr' = col_factor())
  )
  post_annovar <- post_annovar[post_annovar$Chr != "chr",]

  # rename columns
  names(post_annovar)[names(post_annovar) == "Chr"] <- "chr"
  names(post_annovar)[names(post_annovar) == "Start"] <- "start"
  names(post_annovar)[names(post_annovar) == "End"] <- "end"
  names(post_annovar)[names(post_annovar) == "Ref"] <- "ref"
  names(post_annovar)[names(post_annovar) == "Alt"] <- "alt"
  names(post_annovar)[names(post_annovar) == "Gene.refGene"] <- "gene"

  # create "address" column for post_annovar as a primary key to enable merge
  post_annovar$address <- paste0(post_annovar$chr, ":", post_annovar$start)
  post_annovar$address <- paste0(post_annovar$address, ":", post_annovar$ref)
  post_annovar$address <- paste0(post_annovar$address, ":", post_annovar$alt)

  # merge the "vcf_post_post_annovar" dataframe with the "vcf" dataframe (so that the 10000 variants not called by post_annovar and the het/homozygous status from "vcf" is preserved going forward)
  vcf <- merge(vcf, post_annovar, by = "address", all.x = TRUE)

  # clean up the column names and remove redundant columns
  names(vcf)[names(vcf) == "chr.x"] <- "#CHROM"  # CADD requires the capitalization and the hash sign.
  names(vcf)[names(vcf) == "start.x"] <- "start"
  names(vcf)[names(vcf) == "end.x"] <- "end"
  names(vcf)[names(vcf) == "ref.x"] <- "ref"
  names(vcf)[names(vcf) == "alt.x"] <- "alt"
  vcf <- subset(vcf, select = -c(chr.y, start.y, end.y, ref.y, alt.y, address))

  rm(post_annovar)

  # Step 4. Add CADD annotation ####

  cadd_save <- file.path(output_dir, paste0(patient_id, ".forcadd.vcf"))
  cadd_output <- file.path(output_dir, paste0(patient_id, ".forcadd.tsv.gz"))

  write_delim(vcf, cadd_save, delim = "\t", na = ".", col_names = TRUE, quote_escape = FALSE)

  cadd_phrase <- paste(file.path(path_to_cadd, "CADD.sh"), "-p -c", num.cpu, cadd_save)

  # terminal command to run offline local cadd command line program
  if (!file.exists(cadd_output)) {
    print(paste("Running command:", cadd_phrase, sep = "\n"))
    system(cadd_phrase)
  } else {
    print(paste("Skipping CADD because", cadd_output, "exists."))
  }

  # creating "address" column as a primary key for each data frame to enable merge
  vcf$address <- paste0(vcf[, "#CHROM"], ":", vcf$start)
  vcf$address <- paste0(vcf$address, ":", vcf$ref)
  vcf$address <- paste0(vcf$address, ":", vcf$alt)

  cadd_score <- read_delim(cadd_output, "\t", guess_max = 200000, escape_double = FALSE, trim_ws = TRUE, skip = 1,)

  # clean up the annotations
  names(cadd_score)[names(cadd_score) == "#Chrom"] <- "CHROM"  # CADD outputs header slightly differently
  names(cadd_score)[names(cadd_score) == "PHRED"] <- "CADD16_PHRED"

  # create a primary key for cadd data frame to enable merge
  cadd_score$address <- paste0(cadd_score$CHROM, ":", cadd_score$Pos)
  cadd_score$address <- paste0(cadd_score$address, ":", cadd_score$Ref)
  cadd_score$address <- paste0(cadd_score$address, ":", cadd_score$Alt)

  vcf$AF <- as.character(vcf$AF)

  vcf <- merge(vcf, cadd_score, by = "address", all.x = TRUE)

  rm(cadd_score)


  # Step 5. Summarize dbNSFP scores ####
  # Note: the code for summarizing and Simplifying Prediction Scores one by one is very long, about 700 lines

  #SIFT_pred to SIFT_pred_binary############################################################3

  #Create SIFT_pred binary; remove punctuation so can count letters in next step
  vcf$SIFT_pred_binary_prep <- gsub(';', '', vcf$SIFT_pred)
  vcf$SIFT_pred_binary_prep <- gsub("\\.", "", vcf$SIFT_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$SIFT_pred_Donly <- gsub('T', '', vcf$SIFT_pred_binary_prep)
  vcf$SIFT_pred_Tonly <- gsub('D', '', vcf$SIFT_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$SIFT_pred_Dcount <- nchar(vcf$SIFT_pred_Donly)
  vcf$SIFT_pred_Tcount <- nchar(vcf$SIFT_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$SIFT_pred_proportion <- (vcf$SIFT_pred_Dcount / (vcf$SIFT_pred_Dcount + vcf$SIFT_pred_Tcount))

  vcf$SIFT_pred_proportion[is.nan(vcf$SIFT_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$SIFT_pred_binary <- ifelse(vcf$SIFT_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(SIFT_pred, SIFT_pred_binary_prep, SIFT_pred_Donly, SIFT_pred_Tonly, SIFT_pred_Dcount, SIFT_pred_Tcount, SIFT_pred_proportion))


  #SIFT4G_pred to SIFT4G_pred_binary############################################################3

  #Create SIFT4G_pred binary; remove punctuation so can count letters in next step
  vcf$SIFT4G_pred_binary_prep <- gsub(';', '', vcf$SIFT4G_pred)
  vcf$SIFT4G_pred_binary_prep <- gsub("\\.", "", vcf$SIFT4G_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$SIFT4G_pred_Donly <- gsub('T', '', vcf$SIFT4G_pred_binary_prep)
  vcf$SIFT4G_pred_Tonly <- gsub('D', '', vcf$SIFT4G_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$SIFT4G_pred_Dcount <- nchar(vcf$SIFT4G_pred_Donly)
  vcf$SIFT4G_pred_Tcount <- nchar(vcf$SIFT4G_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$SIFT4G_pred_proportion <- (vcf$SIFT4G_pred_Dcount / (vcf$SIFT4G_pred_Dcount + vcf$SIFT4G_pred_Tcount))

  vcf$SIFT4G_pred_proportion[is.nan(vcf$SIFT4G_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$SIFT4G_pred_binary <- ifelse(vcf$SIFT4G_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(SIFT4G_pred, SIFT4G_pred_binary_prep, SIFT4G_pred_Donly, SIFT4G_pred_Tonly, SIFT4G_pred_Dcount, SIFT4G_pred_Tcount, SIFT4G_pred_proportion))


  #Polyphen2_HDIV_pred to Polyphen2_HDIV_pred_binary##################################################3

  #http://genetics.bwh.harvard.edu/pph2/dokuwiki/overview

  #levels:
  # B == benign
  # P == possibly damaging
  # D == probably damaging

  #Create Polyphen2_HDIV_pred binary; remove punctuation so can count letters in next step
  vcf$Polyphen2_HDIV_pred_binary_prep <- gsub(';', '', vcf$Polyphen2_HDIV_pred)
  vcf$Polyphen2_HDIV_pred_binary_prep <- gsub("\\.", "", vcf$Polyphen2_HDIV_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$Polyphen2_HDIV_pred_Bonly <- gsub('P', '', vcf$Polyphen2_HDIV_pred_binary_prep)
  vcf$Polyphen2_HDIV_pred_Bonly <- gsub('D', '', vcf$Polyphen2_HDIV_pred_Bonly)
  vcf$Polyphen2_HDIV_pred_PandDonly <- gsub('B', '', vcf$Polyphen2_HDIV_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$Polyphen2_HDIV_pred_PandDcount <- nchar(vcf$Polyphen2_HDIV_pred_PandDonly)
  vcf$Polyphen2_HDIV_pred_Bcount <- nchar(vcf$Polyphen2_HDIV_pred_Bonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$Polyphen2_HDIV_pred_proportion <- (vcf$Polyphen2_HDIV_pred_PandDcount / (vcf$Polyphen2_HDIV_pred_PandDcount + vcf$Polyphen2_HDIV_pred_Bcount))

  vcf$Polyphen2_HDIV_pred_proportion[is.nan(vcf$Polyphen2_HDIV_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$Polyphen2_HDIV_pred_binary <- ifelse(vcf$Polyphen2_HDIV_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(Polyphen2_HDIV_pred, Polyphen2_HDIV_pred_binary_prep, Polyphen2_HDIV_pred_PandDonly, Polyphen2_HDIV_pred_Bonly, Polyphen2_HDIV_pred_PandDcount, Polyphen2_HDIV_pred_Bcount, Polyphen2_HDIV_pred_proportion))


  #Polyphen2_HVAR_pred to Polyphen2_HVAR_pred_binary######################################3

  #http://genetics.bwh.harvard.edu/pph2/dokuwiki/overview

  #levels:
  # B == benign
  # P == possibly damaging
  # D == probably damaging

  #Create Polyphen2_HVAR_pred binary; remove punctuation so can count letters in next step
  vcf$Polyphen2_HVAR_pred_binary_prep <- gsub(';', '', vcf$Polyphen2_HVAR_pred)
  vcf$Polyphen2_HVAR_pred_binary_prep <- gsub("\\.", "", vcf$Polyphen2_HVAR_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$Polyphen2_HVAR_pred_Bonly <- gsub('P', '', vcf$Polyphen2_HVAR_pred_binary_prep)
  vcf$Polyphen2_HVAR_pred_Bonly <- gsub('D', '', vcf$Polyphen2_HVAR_pred_Bonly)
  vcf$Polyphen2_HVAR_pred_PandDonly <- gsub('B', '', vcf$Polyphen2_HVAR_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$Polyphen2_HVAR_pred_PandDcount <- nchar(vcf$Polyphen2_HVAR_pred_PandDonly)
  vcf$Polyphen2_HVAR_pred_Bcount <- nchar(vcf$Polyphen2_HVAR_pred_Bonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$Polyphen2_HVAR_pred_proportion <- (vcf$Polyphen2_HVAR_pred_PandDcount / (vcf$Polyphen2_HVAR_pred_PandDcount + vcf$Polyphen2_HVAR_pred_Bcount))

  vcf$Polyphen2_HVAR_pred_proportion[is.nan(vcf$Polyphen2_HVAR_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$Polyphen2_HVAR_pred_binary <- ifelse(vcf$Polyphen2_HVAR_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(Polyphen2_HVAR_pred, Polyphen2_HVAR_pred_binary_prep, Polyphen2_HVAR_pred_PandDonly, Polyphen2_HVAR_pred_Bonly, Polyphen2_HVAR_pred_PandDcount, Polyphen2_HVAR_pred_Bcount, Polyphen2_HVAR_pred_proportion))


  #LRT_pred to LRT_pred_binary############################################################3

  #https://genome.cshlp.org/content/19/9/1553.abstract

  #levels:
  # N == neutral
  # U == unknown
  # D == damaging

  #Create LRT_pred binary; remove punctuation so can count letters in next step
  vcf$LRT_pred_binary_prep <- gsub(';', '', vcf$LRT_pred)
  vcf$LRT_pred_binary_prep <- gsub("\\.", "", vcf$LRT_pred_binary_prep)

  #this score just gives one letter per row or "."
  #so I can just convert D to TRUE and N, U, and . to FALSE
  vcf$LRT_pred_binary <- ifelse(vcf$LRT_pred_binary_prep == "D", TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(LRT_pred, LRT_pred_binary_prep))


  #MutationTaster_pred to MutationTaster_pred_binary####################################3

  #levels:
  # P == polymorphism_automatic
  # N == polymorphism
  # D == disease_causing
  # A == disease_causing_automatic

  # Thus I want to count TRUE if D+A/((D+A) + (P+N)) (ie, P and N are negative; D and A are predicting damaging)

  #Create MutationTaster_pred binary; remove punctuation so can count letters in next step
  vcf$MutationTaster_pred_binary_prep <- gsub(';', '', vcf$MutationTaster_pred)
  vcf$MutationTaster_pred_binary_prep <- gsub("\\.", "", vcf$MutationTaster_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$MutationTaster_pred_PandN <- gsub('D', '', vcf$MutationTaster_pred_binary_prep)
  vcf$MutationTaster_pred_PandN <- gsub('A', '', vcf$MutationTaster_pred_PandN)
  vcf$MutationTaster_pred_DandA <- gsub('P', '', vcf$MutationTaster_pred_binary_prep)
  vcf$MutationTaster_pred_DandA <- gsub('N', '', vcf$MutationTaster_pred_DandA)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$MutationTaster_pred_PandNcount <- nchar(vcf$MutationTaster_pred_PandN)
  vcf$MutationTaster_pred_DandAcount <- nchar(vcf$MutationTaster_pred_DandA)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$MutationTaster_pred_proportion <- (vcf$MutationTaster_pred_DandAcount / (vcf$MutationTaster_pred_DandAcount + vcf$MutationTaster_pred_PandNcount))

  vcf$MutationTaster_pred_proportion[is.nan(vcf$MutationTaster_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$MutationTaster_pred_binary <- ifelse(vcf$MutationTaster_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(MutationTaster_pred, MutationTaster_pred_binary_prep, MutationTaster_pred_PandN, MutationTaster_pred_DandA, MutationTaster_pred_PandNcount, MutationTaster_pred_DandAcount, MutationTaster_pred_proportion))


  #MutationAssessor_pred to MutationAssessor_pred_binary##############################3

  #levels:
  # H == high
  # M == medium
  # L == low
  # N == neutral

  # Thus I want to count TRUE if H+M/((H+M) + (L+N)) (ie, L and N are negative; H and M are predicting damaging)

  #Create MutationAssessor_pred binary; remove punctuation so can count letters in next step
  vcf$MutationAssessor_pred_binary_prep <- gsub(';', '', vcf$MutationAssessor_pred)
  vcf$MutationAssessor_pred_binary_prep <- gsub("\\.", "", vcf$MutationAssessor_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$MutationAssessor_pred_LandN <- gsub('H', '', vcf$MutationAssessor_pred_binary_prep)
  vcf$MutationAssessor_pred_LandN <- gsub('M', '', vcf$MutationAssessor_pred_LandN)
  vcf$MutationAssessor_pred_MandH <- gsub('L', '', vcf$MutationAssessor_pred_binary_prep)
  vcf$MutationAssessor_pred_MandH <- gsub('N', '', vcf$MutationAssessor_pred_MandH)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$MutationAssessor_pred_LandNcount <- nchar(vcf$MutationAssessor_pred_LandN)
  vcf$MutationAssessor_pred_MandHcount <- nchar(vcf$MutationAssessor_pred_MandH)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$MutationAssessor_pred_proportion <- (vcf$MutationAssessor_pred_MandHcount / (vcf$MutationAssessor_pred_MandHcount + vcf$MutationAssessor_pred_LandNcount))

  vcf$MutationAssessor_pred_proportion[is.nan(vcf$MutationAssessor_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$MutationAssessor_pred_binary <- ifelse(vcf$MutationAssessor_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(MutationAssessor_pred, MutationAssessor_pred_binary_prep, MutationAssessor_pred_LandN, MutationAssessor_pred_MandH, MutationAssessor_pred_LandNcount, MutationAssessor_pred_MandHcount, MutationAssessor_pred_proportion))


  #FATHMM_pred to FATHMM_pred_binary############################################################3

  #Create FATHMM_pred binary; remove punctuation so can count letters in next step
  vcf$FATHMM_pred_binary_prep <- gsub(';', '', vcf$FATHMM_pred)
  vcf$FATHMM_pred_binary_prep <- gsub("\\.", "", vcf$FATHMM_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$FATHMM_pred_Donly <- gsub('T', '', vcf$FATHMM_pred_binary_prep)
  vcf$FATHMM_pred_Tonly <- gsub('D', '', vcf$FATHMM_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$FATHMM_pred_Dcount <- nchar(vcf$FATHMM_pred_Donly)
  vcf$FATHMM_pred_Tcount <- nchar(vcf$FATHMM_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$FATHMM_pred_proportion <- (vcf$FATHMM_pred_Dcount / (vcf$FATHMM_pred_Dcount + vcf$FATHMM_pred_Tcount))

  vcf$FATHMM_pred_proportion[is.nan(vcf$FATHMM_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$FATHMM_pred_binary <- ifelse(vcf$FATHMM_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(FATHMM_pred, FATHMM_pred_binary_prep, FATHMM_pred_Donly, FATHMM_pred_Tonly, FATHMM_pred_Dcount, FATHMM_pred_Tcount, FATHMM_pred_proportion))


  #PROVEAN_pred to PROVEAN_pred_binary############################################################3

  #Create PROVEAN_pred binary; remove punctuation so can count letters in next step
  vcf$PROVEAN_pred_binary_prep <- gsub(';', '', vcf$PROVEAN_pred)
  vcf$PROVEAN_pred_binary_prep <- gsub("\\.", "", vcf$PROVEAN_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$PROVEAN_pred_Donly <- gsub('N', '', vcf$PROVEAN_pred_binary_prep)
  vcf$PROVEAN_pred_Nonly <- gsub('D', '', vcf$PROVEAN_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$PROVEAN_pred_Dcount <- nchar(vcf$PROVEAN_pred_Donly)
  vcf$PROVEAN_pred_Ncount <- nchar(vcf$PROVEAN_pred_Nonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$PROVEAN_pred_proportion <- (vcf$PROVEAN_pred_Dcount / (vcf$PROVEAN_pred_Dcount + vcf$PROVEAN_pred_Ncount))

  vcf$PROVEAN_pred_proportion[is.nan(vcf$PROVEAN_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$PROVEAN_pred_binary <- ifelse(vcf$PROVEAN_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(PROVEAN_pred, PROVEAN_pred_binary_prep, PROVEAN_pred_Donly, PROVEAN_pred_Nonly, PROVEAN_pred_Dcount, PROVEAN_pred_Ncount, PROVEAN_pred_proportion))


  #MetaSVM_pred to MetaSVM_pred_binary############################################################3

  #Create MetaSVM_pred binary; remove punctuation so can count letters in next step
  vcf$MetaSVM_pred_binary_prep <- gsub(';', '', vcf$MetaSVM_pred)
  vcf$MetaSVM_pred_binary_prep <- gsub("\\.", "", vcf$MetaSVM_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$MetaSVM_pred_Donly <- gsub('T', '', vcf$MetaSVM_pred_binary_prep)
  vcf$MetaSVM_pred_Tonly <- gsub('D', '', vcf$MetaSVM_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$MetaSVM_pred_Dcount <- nchar(vcf$MetaSVM_pred_Donly)
  vcf$MetaSVM_pred_Tcount <- nchar(vcf$MetaSVM_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$MetaSVM_pred_proportion <- (vcf$MetaSVM_pred_Dcount / (vcf$MetaSVM_pred_Dcount + vcf$MetaSVM_pred_Tcount))

  vcf$MetaSVM_pred_proportion[is.nan(vcf$MetaSVM_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$MetaSVM_pred_binary <- ifelse(vcf$MetaSVM_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(MetaSVM_pred, MetaSVM_pred_binary_prep, MetaSVM_pred_Donly, MetaSVM_pred_Tonly, MetaSVM_pred_Dcount, MetaSVM_pred_Tcount, MetaSVM_pred_proportion))


  #MetaLR_pred to MetaLR_pred_binary############################################################3

  #Create MetaLR_pred binary; remove punctuation so can count letters in next step
  vcf$MetaLR_pred_binary_prep <- gsub(';', '', vcf$MetaLR_pred)
  vcf$MetaLR_pred_binary_prep <- gsub("\\.", "", vcf$MetaLR_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$MetaLR_pred_Donly <- gsub('T', '', vcf$MetaLR_pred_binary_prep)
  vcf$MetaLR_pred_Tonly <- gsub('D', '', vcf$MetaLR_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$MetaLR_pred_Dcount <- nchar(vcf$MetaLR_pred_Donly)
  vcf$MetaLR_pred_Tcount <- nchar(vcf$MetaLR_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$MetaLR_pred_proportion <- (vcf$MetaLR_pred_Dcount / (vcf$MetaLR_pred_Dcount + vcf$MetaLR_pred_Tcount))

  vcf$MetaLR_pred_proportion[is.nan(vcf$MetaLR_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$MetaLR_pred_binary <- ifelse(vcf$MetaLR_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(MetaLR_pred, MetaLR_pred_binary_prep, MetaLR_pred_Donly, MetaLR_pred_Tonly, MetaLR_pred_Dcount, MetaLR_pred_Tcount, MetaLR_pred_proportion))


  #M_CAP_pred to M_CAP_pred_binary############################################################3

  #Create M_CAP_pred binary; remove punctuation so can count letters in next step
  vcf$M_CAP_pred_binary_prep <- gsub(';', '', vcf$M_CAP_pred)
  vcf$M_CAP_pred_binary_prep <- gsub("\\.", "", vcf$M_CAP_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$M_CAP_pred_Donly <- gsub('T', '', vcf$M_CAP_pred_binary_prep)
  vcf$M_CAP_pred_Tonly <- gsub('D', '', vcf$M_CAP_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$M_CAP_pred_Dcount <- nchar(vcf$M_CAP_pred_Donly)
  vcf$M_CAP_pred_Tcount <- nchar(vcf$M_CAP_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$M_CAP_pred_proportion <- (vcf$M_CAP_pred_Dcount / (vcf$M_CAP_pred_Dcount + vcf$M_CAP_pred_Tcount))

  vcf$M_CAP_pred_proportion[is.nan(vcf$M_CAP_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$M_CAP_pred_binary <- ifelse(vcf$M_CAP_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(M_CAP_pred, M_CAP_pred_binary_prep, M_CAP_pred_Donly, M_CAP_pred_Tonly, M_CAP_pred_Dcount, M_CAP_pred_Tcount, M_CAP_pred_proportion))


  #PrimateAI_pred to PrimateAI_pred_binary############################################################3

  #Create PrimateAI_pred binary; remove punctuation so can count letters in next step
  vcf$PrimateAI_pred_binary_prep <- gsub(';', '', vcf$PrimateAI_pred)
  vcf$PrimateAI_pred_binary_prep <- gsub("\\.", "", vcf$PrimateAI_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$PrimateAI_pred_Donly <- gsub('T', '', vcf$PrimateAI_pred_binary_prep)
  vcf$PrimateAI_pred_Tonly <- gsub('D', '', vcf$PrimateAI_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$PrimateAI_pred_Dcount <- nchar(vcf$PrimateAI_pred_Donly)
  vcf$PrimateAI_pred_Tcount <- nchar(vcf$PrimateAI_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$PrimateAI_pred_proportion <- (vcf$PrimateAI_pred_Dcount / (vcf$PrimateAI_pred_Dcount + vcf$PrimateAI_pred_Tcount))

  vcf$PrimateAI_pred_proportion[is.nan(vcf$PrimateAI_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$PrimateAI_pred_binary <- ifelse(vcf$PrimateAI_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(PrimateAI_pred, PrimateAI_pred_binary_prep, PrimateAI_pred_Donly, PrimateAI_pred_Tonly, PrimateAI_pred_Dcount, PrimateAI_pred_Tcount, PrimateAI_pred_proportion))


  #DEOGEN2_pred to DEOGEN2_pred_binary############################################################3

  #Create DEOGEN2_pred binary; remove punctuation so can count letters in next step
  vcf$DEOGEN2_pred_binary_prep <- gsub(';', '', vcf$DEOGEN2_pred)
  vcf$DEOGEN2_pred_binary_prep <- gsub("\\.", "", vcf$DEOGEN2_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$DEOGEN2_pred_Donly <- gsub('T', '', vcf$DEOGEN2_pred_binary_prep)
  vcf$DEOGEN2_pred_Tonly <- gsub('D', '', vcf$DEOGEN2_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$DEOGEN2_pred_Dcount <- nchar(vcf$DEOGEN2_pred_Donly)
  vcf$DEOGEN2_pred_Tcount <- nchar(vcf$DEOGEN2_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$DEOGEN2_pred_proportion <- (vcf$DEOGEN2_pred_Dcount / (vcf$DEOGEN2_pred_Dcount + vcf$DEOGEN2_pred_Tcount))

  vcf$DEOGEN2_pred_proportion[is.nan(vcf$DEOGEN2_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$DEOGEN2_pred_binary <- ifelse(vcf$DEOGEN2_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(DEOGEN2_pred, DEOGEN2_pred_binary_prep, DEOGEN2_pred_Donly, DEOGEN2_pred_Tonly, DEOGEN2_pred_Dcount, DEOGEN2_pred_Tcount, DEOGEN2_pred_proportion))


  #BayesDel_addAF_pred to BayesDel_addAF_pred_binary###########################################3

  #Create BayesDel_addAF_pred binary; remove punctuation so can count letters in next step
  vcf$BayesDel_addAF_pred_binary_prep <- gsub(';', '', vcf$BayesDel_addAF_pred)
  vcf$BayesDel_addAF_pred_binary_prep <- gsub("\\.", "", vcf$BayesDel_addAF_pred_binary_prep)
  vcf$BayesDel_addAF_pred_binary_prep <- gsub("TRUE", "T", vcf$BayesDel_addAF_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$BayesDel_addAF_pred_Donly <- gsub('T', '', vcf$BayesDel_addAF_pred_binary_prep)
  vcf$BayesDel_addAF_pred_Tonly <- gsub('D', '', vcf$BayesDel_addAF_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$BayesDel_addAF_pred_Dcount <- nchar(vcf$BayesDel_addAF_pred_Donly)
  vcf$BayesDel_addAF_pred_Tcount <- nchar(vcf$BayesDel_addAF_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$BayesDel_addAF_pred_proportion <- (vcf$BayesDel_addAF_pred_Dcount / (vcf$BayesDel_addAF_pred_Dcount + vcf$BayesDel_addAF_pred_Tcount))

  vcf$BayesDel_addAF_pred_proportion[is.nan(vcf$BayesDel_addAF_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$BayesDel_addAF_pred_binary <- ifelse(vcf$BayesDel_addAF_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(BayesDel_addAF_pred, BayesDel_addAF_pred_binary_prep, BayesDel_addAF_pred_Donly, BayesDel_addAF_pred_Tonly, BayesDel_addAF_pred_Dcount, BayesDel_addAF_pred_Tcount, BayesDel_addAF_pred_proportion))


  #BayesDel_noAF_pred to BayesDel_noAF_pred_binary############################################3

  #Create BayesDel_noAF_pred binary; remove punctuation so can count letters in next step
  vcf$BayesDel_noAF_pred_binary_prep <- gsub(';', '', vcf$BayesDel_noAF_pred)
  vcf$BayesDel_noAF_pred_binary_prep <- gsub("\\.", "", vcf$BayesDel_noAF_pred_binary_prep)
  vcf$BayesDel_noAF_pred_binary_prep <- gsub("TRUE", "T", vcf$BayesDel_noAF_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$BayesDel_noAF_pred_Donly <- gsub('T', '', vcf$BayesDel_noAF_pred_binary_prep)
  vcf$BayesDel_noAF_pred_Tonly <- gsub('D', '', vcf$BayesDel_noAF_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$BayesDel_noAF_pred_Dcount <- nchar(vcf$BayesDel_noAF_pred_Donly)
  vcf$BayesDel_noAF_pred_Tcount <- nchar(vcf$BayesDel_noAF_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$BayesDel_noAF_pred_proportion <- (vcf$BayesDel_noAF_pred_Dcount / (vcf$BayesDel_noAF_pred_Dcount + vcf$BayesDel_noAF_pred_Tcount))

  vcf$BayesDel_noAF_pred_proportion[is.nan(vcf$BayesDel_noAF_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$BayesDel_noAF_pred_binary <- ifelse(vcf$BayesDel_noAF_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(BayesDel_noAF_pred, BayesDel_noAF_pred_binary_prep, BayesDel_noAF_pred_Donly, BayesDel_noAF_pred_Tonly, BayesDel_noAF_pred_Dcount, BayesDel_noAF_pred_Tcount, BayesDel_noAF_pred_proportion))


  #ClinPred_pred to ClinPred_pred_binary############################################################3

  #Create ClinPred_pred binary; remove punctuation so can count letters in next step
  vcf$ClinPred_pred_binary_prep <- gsub(';', '', vcf$ClinPred_pred)
  vcf$ClinPred_pred_binary_prep <- gsub("\\.", "", vcf$ClinPred_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$ClinPred_pred_Donly <- gsub('T', '', vcf$ClinPred_pred_binary_prep)
  vcf$ClinPred_pred_Tonly <- gsub('D', '', vcf$ClinPred_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$ClinPred_pred_Dcount <- nchar(vcf$ClinPred_pred_Donly)
  vcf$ClinPred_pred_Tcount <- nchar(vcf$ClinPred_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$ClinPred_pred_proportion <- (vcf$ClinPred_pred_Dcount / (vcf$ClinPred_pred_Dcount + vcf$ClinPred_pred_Tcount))

  vcf$ClinPred_pred_proportion[is.nan(vcf$ClinPred_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$ClinPred_pred_binary <- ifelse(vcf$ClinPred_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(ClinPred_pred, ClinPred_pred_binary_prep, ClinPred_pred_Donly, ClinPred_pred_Tonly, ClinPred_pred_Dcount, ClinPred_pred_Tcount, ClinPred_pred_proportion))


  #LIST_S2_pred to LIST_S2_pred_binary############################################################3

  #Create LIST_S2_pred binary; remove punctuation so can count letters in next step
  vcf$LIST_S2_pred_binary_prep <- gsub(';', '', vcf$LIST_S2_pred)
  vcf$LIST_S2_pred_binary_prep <- gsub("\\.", "", vcf$LIST_S2_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$LIST_S2_pred_Donly <- gsub('T', '', vcf$LIST_S2_pred_binary_prep)
  vcf$LIST_S2_pred_Tonly <- gsub('D', '', vcf$LIST_S2_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$LIST_S2_pred_Dcount <- nchar(vcf$LIST_S2_pred_Donly)
  vcf$LIST_S2_pred_Tcount <- nchar(vcf$LIST_S2_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$LIST_S2_pred_proportion <- (vcf$LIST_S2_pred_Dcount / (vcf$LIST_S2_pred_Dcount + vcf$LIST_S2_pred_Tcount))

  vcf$LIST_S2_pred_proportion[is.nan(vcf$LIST_S2_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$LIST_S2_pred_binary <- ifelse(vcf$LIST_S2_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(LIST_S2_pred, LIST_S2_pred_binary_prep, LIST_S2_pred_Donly, LIST_S2_pred_Tonly, LIST_S2_pred_Dcount, LIST_S2_pred_Tcount, LIST_S2_pred_proportion))


  #Aloft_pred to Aloft_pred_binary############################################################3

  #levels:
  #Recessive == Damaging w Recessive Inheritance
  #Dominant == Damaging w Dominant Inheritance
  #Tolerant == Tolerant

  #Create Aloft_pred binary; remove punctuation so can count letters in next step
  vcf$Aloft_pred_binary_prep <- gsub(';', '', vcf$Aloft_pred)
  vcf$Aloft_pred_binary_prep <- gsub("\\.", "", vcf$Aloft_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$Aloft_pred_Tonly <- gsub('Recessive', '', vcf$Aloft_pred_binary_prep)
  vcf$Aloft_pred_Tonly <- gsub('Dominant', '', vcf$Aloft_pred_Tonly)
  vcf$Aloft_pred_Tonly <- gsub('Tolerant', 'T', vcf$Aloft_pred_Tonly)
  vcf$Aloft_pred_RandDonly <- gsub('Tolerant', '', vcf$Aloft_pred_binary_prep)
  vcf$Aloft_pred_RandDonly <- gsub('Recessive', 'R', vcf$Aloft_pred_RandDonly)
  vcf$Aloft_pred_RandDonly <- gsub('Dominant', 'D', vcf$Aloft_pred_RandDonly)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$Aloft_pred_RandDcount <- nchar(vcf$Aloft_pred_RandDonly)
  vcf$Aloft_pred_Tcount <- nchar(vcf$Aloft_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$Aloft_pred_proportion <- (vcf$Aloft_pred_RandDcount / (vcf$Aloft_pred_RandDcount + vcf$Aloft_pred_Tcount))

  vcf$Aloft_pred_proportion[is.nan(vcf$Aloft_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$Aloft_pred_binary <- ifelse(vcf$Aloft_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(Aloft_pred, Aloft_pred_binary_prep, Aloft_pred_RandDonly, Aloft_pred_Tonly, Aloft_pred_RandDcount, Aloft_pred_Tcount, Aloft_pred_proportion))


  #fathmm_MKL_coding_pred to fathmm_MKL_coding_pred_binary##################################3

  #Create fathmm_MKL_coding_pred binary; remove punctuation so can count letters in next step
  vcf$fathmm_MKL_coding_pred_binary_prep <- gsub(';', '', vcf$fathmm_MKL_coding_pred)
  vcf$fathmm_MKL_coding_pred_binary_prep <- gsub("\\.", "", vcf$fathmm_MKL_coding_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$fathmm_MKL_coding_pred_Donly <- gsub('T', '', vcf$fathmm_MKL_coding_pred_binary_prep)
  vcf$fathmm_MKL_coding_pred_Tonly <- gsub('D', '', vcf$fathmm_MKL_coding_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$fathmm_MKL_coding_pred_Dcount <- nchar(vcf$fathmm_MKL_coding_pred_Donly)
  vcf$fathmm_MKL_coding_pred_Tcount <- nchar(vcf$fathmm_MKL_coding_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$fathmm_MKL_coding_pred_proportion <- (vcf$fathmm_MKL_coding_pred_Dcount / (vcf$fathmm_MKL_coding_pred_Dcount + vcf$fathmm_MKL_coding_pred_Tcount))

  vcf$fathmm_MKL_coding_pred_proportion[is.nan(vcf$fathmm_MKL_coding_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$fathmm_MKL_coding_pred_binary <- ifelse(vcf$fathmm_MKL_coding_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(fathmm_MKL_coding_pred, fathmm_MKL_coding_pred_binary_prep, fathmm_MKL_coding_pred_Donly, fathmm_MKL_coding_pred_Tonly, fathmm_MKL_coding_pred_Dcount, fathmm_MKL_coding_pred_Tcount, fathmm_MKL_coding_pred_proportion))


  #fathmm_XF_coding_pred to fathmm_XF_coding_pred_binary######################################3

  #Create fathmm_XF_coding_pred binary; remove punctuation so can count letters in next step
  vcf$fathmm_XF_coding_pred_binary_prep <- gsub(';', '', vcf$fathmm_XF_coding_pred)
  vcf$fathmm_XF_coding_pred_binary_prep <- gsub("\\.", "", vcf$fathmm_XF_coding_pred_binary_prep)

  #create columns that are just D or just T values so can be counted in the next step
  vcf$fathmm_XF_coding_pred_Donly <- gsub('T', '', vcf$fathmm_XF_coding_pred_binary_prep)
  vcf$fathmm_XF_coding_pred_Tonly <- gsub('D', '', vcf$fathmm_XF_coding_pred_binary_prep)

  #create columns that are counts of D or T so proportion (D/(D+T)) can be calculated in the next step
  vcf$fathmm_XF_coding_pred_Dcount <- nchar(vcf$fathmm_XF_coding_pred_Donly)
  vcf$fathmm_XF_coding_pred_Tcount <- nchar(vcf$fathmm_XF_coding_pred_Tonly)

  #Proportions of D out of total (D+T) so can get a simple binary output in the next step
  vcf$fathmm_XF_coding_pred_proportion <- (vcf$fathmm_XF_coding_pred_Dcount / (vcf$fathmm_XF_coding_pred_Dcount + vcf$fathmm_XF_coding_pred_Tcount))

  vcf$fathmm_XF_coding_pred_proportion[is.nan(vcf$fathmm_XF_coding_pred_proportion)] <- 0

  #function to count as TRUE if D/D+T > 0.5 or FALSE if (D/D+T <=0.5 | blank (if any blanks still remain))
  vcf$fathmm_XF_coding_pred_binary <- ifelse(vcf$fathmm_XF_coding_pred_proportion > 0.5, TRUE, FALSE)

  #remove columns used in the middle
  vcf <- subset(vcf, select = -c(fathmm_XF_coding_pred, fathmm_XF_coding_pred_binary_prep, fathmm_XF_coding_pred_Donly, fathmm_XF_coding_pred_Tonly, fathmm_XF_coding_pred_Dcount, fathmm_XF_coding_pred_Tcount, fathmm_XF_coding_pred_proportion))

  #Add a total_binary_score column###########################################################3

  #count the trues in the 21 binary columns:
  vcf_binary_only <- subset(vcf, select = c(
    SIFT_pred_binary,
    SIFT4G_pred_binary,
    Polyphen2_HDIV_pred_binary,
    Polyphen2_HVAR_pred_binary,
    LRT_pred_binary,
    MutationTaster_pred_binary,
    MutationAssessor_pred_binary,
    FATHMM_pred_binary,
    PROVEAN_pred_binary,
    MetaSVM_pred_binary,
    MetaLR_pred_binary,
    M_CAP_pred_binary,
    PrimateAI_pred_binary,
    DEOGEN2_pred_binary,
    BayesDel_addAF_pred_binary,
    BayesDel_noAF_pred_binary,
    ClinPred_pred_binary,
    LIST_S2_pred_binary,
    Aloft_pred_binary,
    fathmm_MKL_coding_pred_binary,
    fathmm_XF_coding_pred_binary
  ))

  vcf_binary_all_empty <- ifelse(rowSums(is.na(vcf_binary_only)) == ncol(vcf_binary_only), "empty_row", "not_empty_row")

  dbNSFP_count <- rowSums(vcf_binary_only, na.rm = TRUE)

  vcf <- cbind(vcf, dbNSFP_count, vcf_binary_all_empty)

  binary_just_empty <- subset(vcf, vcf_binary_all_empty == "empty_row")
  binary_not_empty <- subset(vcf, vcf_binary_all_empty == "not_empty_row")

  binary_just_empty$dbNSFP_count <- NA

  vcf <- rbind(binary_just_empty, binary_not_empty)

  rm(vcf_binary_all_empty, vcf_binary_only, dbNSFP_count, binary_just_empty, binary_not_empty)

  names(vcf)[names(vcf) == "#CHROM"] <- "CHROM"

  vcf <- select(vcf, address, CHROM, start, end, ref, alt, QUAL, FILTER, INFO, FORMAT, genotype, AF, Func.refGene,
                gene, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene, CADD16_PHRED, dbNSFP_count)

  #Step 6. LOEUF Annotation############################################################

  gnomad_constraints <- read_delim(path_to_LOEUF_table, "\t", escape_double =
    FALSE, trim_ws = TRUE)

  gnomad_LOEUF <- subset(gnomad_constraints, select = c(gene, oe_lof_upper))

  vcf <- merge(x = vcf, y = gnomad_LOEUF, by = "gene", all.x = TRUE)

  names(vcf)[names(vcf) == "oe_lof_upper"] <- "LOEUF"
  names(vcf)[names(vcf) == "CHROM"] <- "chr"

  rm(gnomad_constraints, gnomad_LOEUF)

  #Step 7. Gene Lists Annotation (denotes a logical if the variant is within the coordinates of the gene####

  #A. monogenic ibd 99 genelist annotation as 'monoibd99'
  vcf$chr <- as.factor(vcf$chr)
  monoibd99_bed <- read.csv(file.path(path_to_monoibdpriority, "monoibd99.bed"))

  #remove the 'chr' from the start of the chromosome column
  monoibd99_bed$chrom <- gsub('chr', '', monoibd99_bed$chrom)
  names(monoibd99_bed)[names(monoibd99_bed) == "chrom"] <- "chr"
  monoibd99_bed$chr <- as.factor(monoibd99_bed$chr)

  #get a list of monogenic genes so can filter on gene name
  monoibd99_list <- monoibd99_bed$gene

  #make a list out of the gene column from the vcf
  variant_gene_list <- vcf$gene

  #create a monoibd99 column that labels the variant "TRUE" if the gene name is found in the list of the 99 known monogenic genes
  vcf$monoibd99_by_name <- ifelse(variant_gene_list %in% monoibd99_list, TRUE, FALSE)
  #   sum(vcf$monoibd99_by_name)

  #create a monoibd99 column that labels the variant "TRUE" if the interval is within a monoibd99 bed interval
  #set up a data.table form the bed df so can use the foverlaps function (that is within the data.table package)
  monoibd99_bed_dt <- subset(monoibd99_bed, select = c(chr, start, end))
  monoibd99_bed_dt <- as.data.table(monoibd99_bed_dt)

  setkey(monoibd99_bed_dt)

  #make a data.table from the vcf that is also amenable to the foverlaps function
  vcf_simple <- subset(vcf, select = c(chr, start, end))
  vcf_simple$chr <- as.factor(vcf_simple$chr)
  vcf_simple <- as.data.table(vcf_simple)

  #find overlaps
  mono_by_interval <- foverlaps(vcf_simple, monoibd99_bed_dt, by.x = names(vcf_simple), type = "within", mult = "all", nomatch = 0L)

  #make a short_address in the by_interval datatable in order to be able to merge the data
  mono_by_interval$short_address <- paste0(mono_by_interval$chr, ":", mono_by_interval$i.start)

  #making a short_address secondary key in the vcf in order to be able to merge the data
  vcf$short_address <- paste0(vcf$chr, ":", vcf$start)

  vcf <- merge(vcf, mono_by_interval, by = "short_address", all.x = TRUE)

  vcf$monoibd99_via_bed <- ifelse(!is.na(vcf$i.start), TRUE, FALSE)

  #combine both mono criteria to make one column of monoibd99 TRUE/FALSE
  vcf$monoibd99 <- ifelse(((vcf$monoibd99_by_name == TRUE) | (vcf$monoibd99_via_bed == TRUE)), TRUE, FALSE)
  #sum(vcf$monoibd99)

  #clean up vcf object so can be run through next annotation
  vcf <- subset(vcf, select = -c(short_address, monoibd99_by_name, chr.y, start.y, end.y, i.start, i.end, monoibd99_via_bed))
  vcf <- rename(vcf, c("chr" = "chr.x", "start" = "start.x", "end" = "end.x"))
  rm(mono_by_interval, monoibd99_list, monoibd99_bed, monoibd99_bed_dt, vcf_simple)


  #B. IBD GWAS genelist annotation as 'GWAS440' - see section A for detailed comments on this process
  vcf$chrom <- as.factor(vcf$chr)
  GWAS440_bed <- read.csv(file.path(path_to_monoibdpriority, "ibdgwas440.bed"))

  GWAS440_bed$chrom <- gsub('chr', '', GWAS440_bed$chrom)
  names(GWAS440_bed)[names(GWAS440_bed) == "chrom"] <- "chr"
  GWAS440_bed$chr <- as.factor(GWAS440_bed$chr)

  GWAS440_list <- GWAS440_bed$gene

  variant_gene_list <- vcf$gene

  vcf$GWAS440_by_name <- ifelse(variant_gene_list %in% GWAS440_list, TRUE, FALSE)
  #sum(vcf$GWAS440_by_name)

  GWAS440_bed_dt <- subset(GWAS440_bed, select = c(chr, start, end))
  GWAS440_bed_dt <- as.data.table(GWAS440_bed_dt)
  setkey(GWAS440_bed_dt)

  vcf_simple <- subset(vcf, select = c(chr, start, end))
  vcf_simple$chr <- as.factor(vcf_simple$chr)
  vcf_simple <- as.data.table(vcf_simple)

  mono_by_interval <- foverlaps(vcf_simple, GWAS440_bed_dt, by.x = names(vcf_simple), type = "within", mult = "all", nomatch = 0L)

  mono_by_interval$short_address <- paste0(mono_by_interval$chr, ":", mono_by_interval$i.start)

  vcf$short_address <- paste0(vcf$chr, ":", vcf$start)

  vcf <- merge(vcf, mono_by_interval, by = "short_address", all.x = TRUE)

  vcf$GWAS440_via_bed <- ifelse(!is.na(vcf$i.start), TRUE, FALSE)

  vcf$GWAS440 <- ifelse(((vcf$GWAS440_by_name == TRUE) | (vcf$GWAS440_via_bed == TRUE)), TRUE, FALSE)
  # sum(vcf$GWAS440)

  vcf <- subset(vcf, select = -c(short_address, GWAS440_by_name, chr.y, start.y, end.y, i.start, i.end, GWAS440_via_bed))
  vcf <- rename(vcf, c("chr" = "chr.x", "start" = "start.x", "end" = "end.x"))
  rm(mono_by_interval, GWAS440_list, GWAS440_bed, GWAS440_bed_dt, vcf_simple)


  #C. pid400 genelist annotation as 'pid400' - see section A for detailed comments on this process
  vcf$chrom <- as.factor(vcf$chr)
  pid400_bed <- read.csv(file.path(path_to_monoibdpriority, "pid400.bed"))

  pid400_bed$chrom <- gsub('chr', '', pid400_bed$chrom)
  names(pid400_bed)[names(pid400_bed) == "chrom"] <- "chr"
  pid400_bed$chr <- as.factor(pid400_bed$chr)

  pid400_list <- pid400_bed$gene

  variant_gene_list <- vcf$gene

  vcf$pid400_by_name <- ifelse(variant_gene_list %in% pid400_list, TRUE, FALSE)
  # sum(vcf$pid400_by_name)

  pid400_bed_dt <- subset(pid400_bed, select = c(chr, start, end))
  pid400_bed_dt <- as.data.table(pid400_bed_dt)
  setkey(pid400_bed_dt)

  vcf_simple <- subset(vcf, select = c(chr, start, end))
  vcf_simple$chr <- as.factor(vcf_simple$chr)
  vcf_simple <- as.data.table(vcf_simple)

  mono_by_interval <- foverlaps(vcf_simple, pid400_bed_dt, by.x = names(vcf_simple), type = "within", mult = "all", nomatch = 0L)

  mono_by_interval$short_address <- paste0(mono_by_interval$chr, ":", mono_by_interval$i.start)

  vcf$short_address <- paste0(vcf$chr, ":", vcf$start)

  vcf <- merge(vcf, mono_by_interval, by = "short_address", all.x = TRUE)

  vcf$pid400_via_bed <- ifelse(!is.na(vcf$i.start), TRUE, FALSE)

  vcf$pid400 <- ifelse(((vcf$pid400_by_name == TRUE) | (vcf$pid400_via_bed == TRUE)), TRUE, FALSE)
  sum(vcf$pid400)

  vcf <- subset(vcf, select = -c(short_address, pid400_by_name, chr.y, start.y, end.y, i.start, i.end, pid400_via_bed))
  vcf <- rename(vcf, c("chr" = "chr.x", "start" = "start.x", "end" = "end.x"))
  rm(mono_by_interval, pid400_list, pid400_bed, pid400_bed_dt, vcf_simple)


  #D. cdg468 genelist annotation as 'cdg468 - see section A for detailed comments on this process
  vcf$chrom <- as.factor(vcf$chr)
  cdg468_bed <- read.csv(file.path(path_to_monoibdpriority, "cdg468.bed"))

  cdg468_bed$chrom <- gsub('chr', '', cdg468_bed$chrom)
  names(cdg468_bed)[names(cdg468_bed) == "chrom"] <- "chr"
  cdg468_bed$chr <- as.factor(cdg468_bed$chr)

  cdg468_list <- cdg468_bed$gene

  variant_gene_list <- vcf$gene

  vcf$cdg468_by_name <- ifelse(variant_gene_list %in% cdg468_list, TRUE, FALSE)
  #sum(vcf$cdg468_by_name)

  cdg468_bed_dt <- subset(cdg468_bed, select = c(chr, start, end))
  cdg468_bed_dt <- as.data.table(cdg468_bed_dt)
  setkey(cdg468_bed_dt)

  vcf_simple <- subset(vcf, select = c(chr, start, end))
  vcf_simple$chr <- as.factor(vcf_simple$chr)
  vcf_simple <- as.data.table(vcf_simple)

  mono_by_interval <- foverlaps(vcf_simple, cdg468_bed_dt, by.x = names(vcf_simple), type = "within", mult = "all", nomatch = 0L)

  mono_by_interval$short_address <- paste0(mono_by_interval$chr, ":", mono_by_interval$i.start)

  vcf$short_address <- paste0(vcf$chr, ":", vcf$start)

  vcf <- merge(vcf, mono_by_interval, by = "short_address", all.x = TRUE)

  vcf$cdg468_via_bed <- ifelse(!is.na(vcf$i.start), TRUE, FALSE)

  vcf$cdg468 <- ifelse(((vcf$cdg468_by_name == TRUE) | (vcf$cdg468_via_bed == TRUE)), TRUE, FALSE)
  #sum(vcf$cdg468)

  vcf <- subset(vcf, select = -c(short_address, cdg468_by_name, chr.y, start.y, end.y, i.start, i.end, cdg468_via_bed))
  vcf <- rename(vcf, c("chr" = "chr.x", "start" = "start.x", "end" = "end.x"))
  rm(mono_by_interval, cdg468_list, cdg468_bed, cdg468_bed_dt, vcf_simple)


  #Step 8. Annotate het_or_hom from GENOTYPE column in VCF and clean up data frame####

  het_or_hom <- vcf$genotype
  vcf$het_or_hom <- substring(het_or_hom, 1, 3)
  vcf$het_or_hom <- as.factor(vcf$het_or_hom)

  vcf$het_or_hom <- gsub('0/1', 'het', vcf$het_or_hom)
  vcf$het_or_hom <- gsub('1/1', 'hom_alt', vcf$het_or_hom)
  vcf$het_or_hom <- gsub('1/2', 'hom_alt', vcf$het_or_hom)
  vcf$het_or_hom <- as.factor(vcf$het_or_hom)

  rm(het_or_hom)

  vcf$AF[vcf$AF == "."] <- 0
  vcf[vcf == "."] <- NA

  vcf$AF <- as.numeric(vcf$AF)

  # vcf <- subset(vcf, select = -c(chrom))

  vcf <- select(vcf, gene,
                Func.refGene,
                ExonicFunc.refGene,
                het_or_hom,
                AF,
                CADD16_PHRED,
                dbNSFP_count,
                LOEUF,
                monoibd99,
                GWAS440,
                pid400,
                cdg468,
                GeneDetail.refGene,
                AAChange.refGene,
                chr,
                start,
                end,
                ref,
                alt,
  )

  #save the annotated VCF, ready for Part 3 (filtering)
  annotation_output_phrase <- as.character(file.path(output_dir, paste0(patient_id, "_annotated.csv")))

  write_csv(vcf, annotation_output_phrase)

}

filter.annotations <- function(out.prefix, proband.fn, proband.sex,
                               paternal.fn, paternal.affected,
                               maternal.fn, maternal.affected) {
  # Filtering the variants
  # Script overview:
  # 1. load in proband vcf (from Part 2 script) and family member vcfs (from Part 1 script) to annotate proband vcf
  # 2. filter variants
  # 3. reorder variants so they are prioritized (highest are most likely damaging)
  # 4. save output to a .csv file for further analysis
  # Part 1 - Loading VCFs ####

  # proband; create address column as primary key; add sex
  proband.vcf <- read_csv(proband.fn, col_types = cols(chr = col_character())) %>%
    mutate(address = paste(chr, start, ref, alt, sep = ":"), sex = proband.sex) %>%
    rename(proband_het_or_hom = het_or_hom)  # Rename column to avoid conflicts

  # paternal VCF
  if (!is.na(paternal.fn)) {
    pat.vcf <- read_csv(paternal.fn)
    proband.vcf <- mutate(proband.vcf, father_affected = paternal.affected) %>%
      left_join(pat.vcf, by = "address", na_matches = "never") %>%
      mutate(pat_het_or_hom = if_else(is.na(het_or_hom), "unknown", het_or_hom)) %>%
      select(-het_or_hom)
  } else {
    proband.vcf <- mutate(proband.vcf, father_affected = "unknown",
                          pat_het_or_hom = "unknown")
  }

  # maternal VCF
  if (!is.na(maternal.fn)) {
    mat.vcf <- read_csv(maternal.fn)
    proband.vcf <- mutate(proband.vcf, mother_affected = maternal.affected) %>%
      left_join(mat.vcf, by = "address", na_matches = "never") %>%
      mutate(mat_het_or_hom = if_else(is.na(het_or_hom), "unknown", het_or_hom)) %>%
      select(-het_or_hom)
  } else {
    proband.vcf <- mutate(proband.vcf, mother_affected = "unknown",
                          mat_het_or_hom = "unknown")
  }

  # Fix that name
  proband.vcf <- rename(proband.vcf, het_or_hom = proband_het_or_hom)

  # Part 2 - Inheritance Annotation ####

  #RECESSIVE INHERITANCE OF VARAINT (full penetrance not assumed, just inheritance at this stage, no affectation assumption)
  proband.vcf <- mutate(
    proband.vcf,
    recessive_inheritance = if_else(
      ((pat_het_or_hom == "het" | pat_het_or_hom == "hom_alt") &
        (mat_het_or_hom == "het" | mat_het_or_hom == "hom_alt") &
        het_or_hom == "hom_alt"),
      TRUE, FALSE))
  #DOMINANT INHERITANCE - the variant is het or hom and comes from an affected family member
  proband.vcf <- mutate(
    proband.vcf,
    dominant_inheritance = if_else(
      (((pat_het_or_hom == "het" | pat_het_or_hom == "hom_alt") & father_affected == "TRUE") |
        ((mat_het_or_hom == "het" | mat_het_or_hom == "hom_alt") & mother_affected == "TRUE") |
        het_or_hom == "het"),
      TRUE, FALSE))
  #X LINKED - the variant is het or hom_alt and is on X and patient is male
  proband.vcf <- mutate(
    proband.vcf,
    XL_inheritance = if_else(chr == "X" & sex == "male", TRUE, FALSE))

  # combine patterns
  proband.vcf <- mutate(
    proband.vcf,
    inheritance = case_when(recessive_inheritance == TRUE ~ "recessive",
                            dominant_inheritance == TRUE ~ "dominant",
                            XL_inheritance == TRUE ~ "x_linked",
                            TRUE ~ "other"))
  proband.vcf <- distinct(proband.vcf)

  #Part 3 - Filtering ####

  proband.vcf <- select(proband.vcf,
                        gene,
                        Func.refGene,
                        ExonicFunc.refGene,
                        het_or_hom,
                        AF,
                        CADD16_PHRED,
                        dbNSFP_count,
                        LOEUF,
                        monoibd99,
                        GWAS440,
                        pid400,
                        cdg468,
                        GeneDetail.refGene,
                        AAChange.refGene,
                        everything(),
                        -address,
                        -recessive_inheritance,
                        -dominant_inheritance,
                        -XL_inheritance
  )

  # save unfiltered vcf
  write_csv(proband.vcf, paste0(out.prefix, "_unfiltered.csv"))

  # filtering parameters (can easily be adjusted at this step)
  proband.vcf <- subset(proband.vcf,
                        (Func.refGene == "exonic" | is.na(Func.refGene)) &
                          ((AF < 0.003) &
                            ((CADD16_PHRED > 18) | is.na(CADD16_PHRED)) &
                            ((dbNSFP_count >= 3) | is.na(dbNSFP_count)) &
                            ((LOEUF < 1.50) | is.na(LOEUF)) &
                            ((monoibd99 == TRUE) |
                              (pid400 == TRUE))))

  # sort variants so most likely monogenic at top
  proband.vcf <- arrange(proband.vcf, desc(monoibd99), desc(ExonicFunc.refGene), AF, desc(CADD16_PHRED))

  # save filtered vcf
  write_csv(proband.vcf, paste0(out.prefix, "_filtered.csv"))
}

process_trio <- function(out.dir, proband.fn, proband.sex,
                         pat.fn = NA, pat.affect = "unknown",
                         mat.fn = NA, mat.affect = "unknown") {
  dir.create(out.dir, showWarnings = F)

  if (is.specified(pat.fn)) {
    fmted.pat.fn <- file.path(out.dir, gsub(".vcf", ".txt", basename(pat.fn), fixed = TRUE))
    if (!file.exists(fmted.pat.fn)) {
      fam_memb_processing(pat.fn, fmted.pat.fn)
    }
  } else {
    fmted.pat.fn <- NA
    pat.affect <- "unknown"
  }
  if (is.specified(mat.fn)) {
    fmted.mat.fn <- file.path(out.dir, gsub(".vcf", ".txt", basename(mat.fn), fixed = TRUE))
    if (!file.exists(fmted.mat.fn)) {
      fam_memb_processing(mat.fn, fmted.mat.fn)
    }
  } else {
    fmted.pat.fn <- NA
    mat.affect <- "unknown"
  }
  annotated.proband.fn <- file.path(out.dir, gsub(".vcf", "_annotated.csv",
                                                  basename(proband.fn), fixed = TRUE))
  if (!file.exists(annotated.proband.fn)) {
    proband_annotation(proband.fn, out.dir)
  }
  final_out.prefix <- file.path(out.dir, gsub(".vcf", "", basename(proband.fn), fixed = TRUE))
  filter.annotations(final_out.prefix, annotated.proband.fn, proband.sex,
                     fmted.pat.fn, pat.affect, fmted.mat.fn, mat.affect)

  return()
}

# Main function // Run the analysis with given inputs at the top of the script

if (is.specified(vcf_dir) && is.specified(family_units)) {
  # Option 1 used
  if (endsWith(family_units, ".xlsx") || endsWith(family_units, ".xls")) {
    library(readxl)
    family_data <- read_excel(family_units)
  } else {
    family_data <- read_csv(family_units)
  }
  pwalk(family_data, function(proband, sex, paternal, paternal_affected, maternal, maternal_affected,
                              vcf_suffix = ".hc.vcf") {
    proband_file <- file.path(vcf_dir, paste0(proband, vcf_suffix))
    paternal_file <- ifelse(paternal == "NA", NA, file.path(vcf_dir, paste0(paternal, vcf_suffix)))
    maternal_file <- ifelse(maternal == "NA", NA, file.path(vcf_dir, paste0(maternal, vcf_suffix)))
    process_trio(output_dir, proband_file, sex, paternal_file, paternal_affected,
                 maternal_file, maternal_affected)
  })
} else {
  # Option 2 used
  process_trio(output_dir, proband_fn, proband_sex,
               paternal_fn, paternal_affected, maternal_fn, maternal_affected)
}
