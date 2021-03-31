### Monogenic IBD Variant Prioritization
##  Written by Daniel Mulder, February 2021
#   Part 3 Filtering the variants

# This script needs to be run individually for each proband, so it is suggested that this TEMPLATE be copied and modified proband-by-proband
# If one or both parent VCFs are not available, the "paternal VCF" and/or "maternal VCF" sections should be removed, as should the "Part 2 - Inheritance Annotation" section

# To use this script, the following user changes are needed:
#1. change the line below to replace "proband_id" with the patient's ID in the line below (do not remove the quotation marks)
proband.annotated.snp.fn <- "proband_id"
#2. change the line below to replace "sex" with "male" or "female" as it pertains to the proband (do not remove the quotation marks)
proband_sex <- "sex"
#3. change the line below to replace "paternal_filename" with the proband's father's VCF filename (without the .vcf extension) in the line below (do not remove the quotation marks)
paternal.vcf.fn <- "paternal_filename"
#4. change the line below to replace "paternal_status" with the proband's father's affectation status as TRUE or FALSE (no quotation marks needed)
paternal_affectation <- paternal_status
#5. change the line below to replace "maternal_filename" with the proband's mother's VCF filename (without the .vcf extension) in the line below (do not remove the quotation marks)
maternal.vcf.fn <- "maternal_filename"
#6. change the line below to replace "maternal_status" with the proband's mother's affectation status as TRUE or FALSE (no quotation marks needed)
maternal_affectation <- maternal_status
#7. The family's VCFs (proband + available parents) should all be placed in the same folder, working directory should be changed below by changing "path/to/files" with the proper path
setwd("path/to/files")

# Script overview:
# 1. load in proband vcf (from Part 2 script) and family member vcfs (from Part 1 script) to annotate proband vcf
# 2. filter variants
# 3. reorder variants so they are prioritized (highest are most likely damaging)
# 4. save output to a .csv file for further analysis


library(tidyverse)
library(glue) #for interpreting string literals


filter.annotations <- function(out.prefix, proband.fn, proband.sex,
                               paternal.fn, paternal.affected,
                               maternal.fn, maternal.affected) {
  # Part 1 - Loading VCFs ####

  # proband; create address column as primary key; add sex
  proband.vcf <- read_csv(proband.fn, col_types = cols(chr = col_character())) %>%
    mutate(address = paste(chr, start, ref, alt, sep = ":"), sex = proband.sex) %>%
    rename(proband_het_or_hom = het_or_hom)  # Rename column to avoid conflicts

  # paternal VCF
  if (!is.null(paternal.fn) || length(paternal.fn) > 0) {
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
  if (!is.null(maternal.fn) | length(maternal.fn) > 0) {
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
    inheritance = case_when(pat_het_or_hom == "unknown" | mat_het_or_hom == "unknown" ~ "unknown",
                            recessive_inheritance == TRUE ~ "recessive",
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

