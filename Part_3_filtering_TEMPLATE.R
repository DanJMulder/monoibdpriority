### Monogenic IBD Variant Prioritization
##  Written by Daniel Mulder, February 2021
#   Part 3 Filtering the variants

# This script needs to be run individually for each proband, so it is suggested that this TEMPLATE be copied and modified proband-by-proband
# If one or both parent VCFs are not available, the "paternal VCF" and/or "maternal VCF" sections should be removed, as should the "Part 2 - Inheritance Annotation" section

# To use this script, the following user changes are needed:
  #1. change the line below to replace "proband_id" with the patient's ID in the line below (do not remove the quotation marks)
    patient_id <- "proband_id"
  #2. change the line below to replace "sex" with "male" or "female" as it pertains to the proband (do not remove the quotation marks)
    proband_sex <- "sex"
  #3. change the line below to replace "paternal_filename" with the proband's father's VCF filename (without the .vcf extension) in the line below (do not remove the quotation marks)
    paternal_id <- "paternal_filename"
  #4. change the line below to replace "paternal_status" with the proband's father's affectation status as TRUE or FALSE (no quotation marks needed)
    paternal_affectation <- paternal_status
  #5. change the line below to replace "maternal_filename" with the proband's mother's VCF filename (without the .vcf extension) in the line below (do not remove the quotation marks)
    maternal_id <- "maternal_filename"
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
  
  
# Part 1 - Loading VCFs ####
  
  # proband
  vcf <- read.csv(glue("{patient_id}_annotated.csv"))
  # creating "address" column as a primary key for each variant
  vcf$address <- paste0(vcf$chr, ":", vcf$start)
  vcf$address <- paste0(vcf$address, ":", vcf$ref)
  vcf$address <- paste0(vcf$address, ":", vcf$alt)
  #adding the proband's sex
  vcf$sex <- proband_sex
  
  # paternal VCF
  pat_vcf <- read.csv(glue("{paternal_id}.vcf"))
  #paternal affectiation status:
  vcf$father_affected <- paternal_affectation
  #annotating for inheritance of variant from father
  vcf$pat_het_or_hom <- pat_vcf$het_or_hom[match(vcf$address, pat_vcf$address)]
  
  #maternal VCF
  mat_vcf <- read.csv(glue("{maternal_id}.vcf"))
  #maternal affectiation status:
  vcf$mother_affected <- maternal_affectation
  #annotating for inheritance of variant from mother
  vcf$mat_het_or_hom <- mat_vcf$het_or_hom[match(vcf$address, mat_vcf$address)]
  
  
# Part 2 - Inheritance Annotation ####
  
  #RECESSIVE INHERITANCE OF VARAINT (full penetrance not assumed, just inheritance at this stage, no affectation assumption)
  vcf$recessive_inheritance <- ifelse(
    (((vcf$pat_het_or_hom == "het") | (vcf$pat_het_or_hom == "hom_alt")) &
       ((vcf$mat_het_or_hom == "het") | (vcf$mat_het_or_hom == "hom_alt")) &
       (vcf$het_or_hom == "hom_alt")),
    TRUE,
    FALSE)
  vcf$recessive_inheritance <- ifelse(is.na(vcf$recessive_inheritance), FALSE, vcf$recessive_inheritance)
  
  #DOMINANT INHERITANCE - the variant is het or hom and comes from an affected family member
  vcf$dominant_inheritance <- ifelse(
    ((((vcf$pat_het_or_hom == "het") | (vcf$pat_het_or_hom == "hom_alt")) & (vcf$father_affected == TRUE)) |
       (((vcf$mat_het_or_hom == "het") | (vcf$mat_het_or_hom == "hom_alt")) & (vcf$mother_affected == TRUE)) |
       (vcf$het_or_hom == "het")),
    TRUE,
    FALSE)
  
  #X LINKED - the variant is het or hom_alt and is on X and patient is male
  vcf$XL_inheritance <- ifelse(
    ((vcf$chr == "X") & (vcf$sex == "male")),
    TRUE,
    FALSE)
  
  #combine patterns
  
  vcf$inheritance <- "other"
  vcf$inheritance <- ifelse(vcf$recessive_inheritance == TRUE, "recessive",
                            ifelse(vcf$dominant_inheritance == TRUE, "dominant",
                                   ifelse(vcf$XL_inheritance == TRUE, "x_linked", "other")))
  
  vcf <- distinct(vcf)
  
#Part 3 - Filtering ####
  
  vcf <- select(vcf,
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
  write_csv(vcf, glue("{patient_id}_unfiltered.csv"))
  
  # filtering parameters (can easily be adjusted at this step)
  vcf <- subset(vcf,
                (Func.refGene == "exonic" | is.na(Func.refGene)) &
                  ((AF < 0.003) &
                     ((CADD16_PHRED > 18) | is.na(CADD16_PHRED)) &
                     ((dbNSFP_count >= 3) | is.na(dbNSFP_count)) &
                     ((LOEUF <1.50) | is.na(LOEUF)) &
                     ((monoibd99 == TRUE) |
                        (pid400 == TRUE))))
  
  # sort variants so most likely monogenic at top
  vcf <- arrange(vcf, desc(monoibd99), desc(ExonicFunc.refGene), AF, desc(CADD16_PHRED))
  
  # save filtered vcf
  write_csv(vcf, glue("{patient_id}_filtered.csv"))
