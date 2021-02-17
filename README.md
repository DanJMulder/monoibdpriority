# Variant Prioritization for Monogenic IBD Diagnosis

## Background
Genetic sequencing is increasingly used to diagnose monogenic disease. Genetic diagnosis can be highly beneficial to patients - informing prognosis, inheritance, and especially allowing for targeted therapy. Patients presenting with clinical features of very early onset inflammatory bowel disease (VEO-IBD) are "highly recommended" to undergo genetic sequencing analysis (Kelsen et al. JPGN. 2020;70(3);389-403).

Project Rationale
We have accrued a large cohort of pediatric patients with IBD (n = 1004 probands). Patients and available family members underwent whole exome sequencing (WES). Initially, we performed manual analysis for variants suspicous for causing monogenic disease. As we have previously reported, we were able to verify 38 disease causing monogenic variants (~3% of the cohort) through functional laboratory testing (Crowley et al. Gastroenterology. 2020;158(8):2208-2220).

In the study related to this repository, we used statistical analysis of our large cohort WES data to identify a small set of characteristics of in silico predictive features of a variant that were highly predictive of a disease causing variant.

Pipeline Overview
This pipeline is was designed and tested on WES VCF files aligned to both the GRCh37 and GRCh38 reference genomes. This pipeline has not been tested on whole genome sequencing or gene panel files. Raw WES VCFs should be filtered for quality separately prior to using this pipeline. The default is set to GRCh38. Dependencies for annotating VCFs are:
1. A local instance of dbNSFP (version 4.1a, https://sites.google.com/site/jpopgen/dbNSFP)
2. A local instance of annovar (2019Oct24 version, https://annovar.openbioinformatics.org/en/latest/user-guide/download/) with the gnomad211_exome and refGene databases for both GRCh37 and GRCh38 attached.
3. A local instance of CADD (v1.6, https://github.com/kircherlab/CADD-scripts/)
4. A local copy of the gene constraint metrics from gnomAD (v 2.1.1, https://gnomad.broadinstitute.org/downloads#v2-constraint)

This pipeline consists of 3 steps:
Step 1. Family member VCFs are processed into the appropriate format to allow inheritance modelling
Step 2. Proband VCF annotation including: a dbNSFP custom damaging score count, refSeq gene details, minor allele frequency (via gnomad versiona 2.1.1), CADD score, and LOEUF
Step 3. Inheritance modeling for each variant where family members are available, and variant filtering and prioritization based on features likely to identify a disease causing variant

This pipeline can accept a singleton VCF file, but including one or both parent VCFs can highly increase the diagnostic yield. In step 3, adding the affectation status of the parents (if diagnosed with IBD) and the sex of the proband will also add to to the likelihood of identifying a disease causing variant.

The pipeline produces a table of a short list of variants that meet the filtration criteria (output filename "{patient_id}_filtered.csv"). In our cohort, VCFs with 100,000-150,000 variant calls were filtered down to an average of 15 variant calls per patient without removing any of the known disease causing variants. In addition, the pipeline also produces a full table without any variants removed, but with full annotation and inheritance modeling as above, to allow for further well-informed manula curation and analysis of variants (output filename "{patient_id}_unfiltered.csv").
