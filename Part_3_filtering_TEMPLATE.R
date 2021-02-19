### Monogenic IBD Variant Prioritization
##  Written by Daniel Mulder, February 2021
#   Part 3 Filtering the variants

# This script needs to be run individually for each proband, so it is suggested that this TEMPLATE be copied and modified proband-by-proband
# If one or both parent VCFs are not available, the "paternal VCF" and/or "maternal VCF" sections should be removed, as should the "Part 2 - Inheritance Annotation" section

# Script overview:
  # 1. load in proband vcf (from Part 2 script) and family member vcfs (from Part 1 script) to annotate proband vcf
  # 2. filter variants
  # 3. reorder variants so they are prioritized (highest are most likely damaging)
  # 4. save output to a .csv file for further analysis

# For each run of this script it needs to be customized as follows:
  #1. replace {patient_id} with the patient's ID in the line below (do not remove the quotation marks, but do remove the curly brackets):
    patient_id <- as.character(glue("{patient_id}")) 
  #2. input {sex} below with "male" or "female" as it pertains to the proband
  #3. replace {paternal_id} with the paternal VCF filename
  #4. input the father's affectation status
  #5. replace {maternal_id} with the maternal VCF filename
  #6. input the mother's affectation status

library(tidyverse)
    
    
# Part 1 - Loading VCFs ####

  # proband
    vcf <- read.csv("~/{patient_id}_annotated.csv")
    #adding the proband's sex
    vcf$sex <- "female"
  
  # paternal VCF
    pat_vcf <- read.csv("{paternal_id}.vcf")
    #paternal affectiation status:
    vcf$father_affected <- FALSE
    #annotating for inheritance of variant from father
    vcf$pat_het_or_hom <- pat_vcf$het_or_hom[match(vcf$address, pat_vcf$address)]
  
  #maternal VCF
    mat_vcf <- read.csv("{maternal VCF}.vcf")
    #maternal affectiation status:
    vcf$mother_affected <- FALSE
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
                  -QUAL,
                  -FILTER,
                  -INFO,
                  -FORMAT,
                  -recessive_inheritance,
                  -dominant_inheritance,
                  -XL_inheritance
  )
  
  # save unfiltered vcf
  write_csv(vcf, "{patient_id}_unfiltered.csv")
  
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
  write_csv(vcf, "{patient_id}_filtered.csv")