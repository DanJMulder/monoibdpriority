# Variant Prioritization for Monogenic IBD Diagnosis

## Background
Genetic sequencing is increasingly used to diagnose monogenic disease. Genetic diagnosis can be highly beneficial to patients - informing prognosis, inheritance, and especially allowing for targeted therapy. Patients presenting with clinical features of very early onset inflammatory bowel disease (VEO-IBD) are "highly recommended" to undergo genetic sequencing analysis (Kelsen et al. JPGN. 2020;70(3);389-403).

## Project Design 
We have accrued a large cohort of pediatric patients with IBD (n = 1005 probands). Patients and available family members underwent whole exome sequencing (WES). Initially, we performed manual analysis for variants suspicous for causing monogenic disease. As we have previously reported, we were able to verify 38 disease causing monogenic variants (~3% of the cohort) through functional laboratory testing (Crowley et al. Gastroenterology. 2020;158(8):2208-2220).

In the study related to this repository, we used statistical analysis of our large cohort WES data to identify a small set of characteristics of in silico predictive features of a variant that are highly predictive of a disease causing variant.

## Pipeline Overview
Input: Input data can be specified in one of two ways. All files should be uncompressed.
1. Specifying a single set of data comprising the proband VCF +/- one or both parental VCFs
2. An Excel or CSV file specifying multiple datasets to operate as a batch and the directory containing all the VCFs.
   The file suffix should be specified (e.g. "hc.vcf").
    * The first line of the Excel or CSV file must contain the following column names as the header:
        1. proband
            * The ID of the proband. 
        2. sex
            * Values must be either male or female
        3. paternal
            * The ID of the paternal parent
        4. paternal_affected
            * Values must either be TRUE, FALSE, or unknown
        5. maternal
            * The ID of the maternal parent
        6. maternal_affected
            * Values must either be TRUE, FALSE, or unknown
    * VCF files are assumed to be located in a specified VCF directory and named as {id}{vcf_suffix}

Output: annotated VCF-style files as comma separated value files (both filtered+prioritizes and unfiltered versions for further manual curation/analysis)

This pipeline is was designed and tested on WES VCF files aligned to both the GRCh37 and GRCh38 reference genomes. This pipeline has not been tested on whole genome sequencing (WGS) or gene panel files. Raw WES VCFs should be filtered for quality separately prior to using this pipeline. The scripts assume the VCFs are aligned to GRCh38 but this can be manually adjusted for each annotation (dbNSFP, Annovar, CADD, gene lists).

Dependencies for annotating VCFs are:
1. A local instance of dbNSFP (version 4.1a, https://sites.google.com/site/jpopgen/dbNSFP)
2. A local instance of annovar (2019Oct24 version, https://annovar.openbioinformatics.org/en/latest/user-guide/download/) with the gnomad211_exome and refGene databases for the appropriate reference genome attached.
3. A local instance of CADD annotation tool (v1.6, https://github.com/kircherlab/CADD-scripts/)
4. A local copy of the gene constraint metrics from gnomAD (v 2.1.1, https://gnomad.broadinstitute.org/downloads#v2-constraint)
5. A local copy of the bed files containing the relevant gene co-ordinates

Required software packages are listed at the beginning of each script and can be installed using the `install.
packages()` function. You will need:
* tidyverse
* data.table

Running this pipeline consists of 3 steps, which are separated into 3 numbered R scripts:
1. Family member VCF Processing. Processed into the appropriate format to remove unused data and allow inheritance modelling.
2. Proband VCF annotation. Includes a dbNSFP custom damaging score count, refSeq gene details, minor allele frequency (via gnomad versiona 2.1.1), CADD score, and LOEUF.
3. Inheritance modeling. For each variant where family members are available, and variant filtering and prioritization based on features likely to identify a disease causing variant. Filtering and prioritization parameters can be adjusted by the user to improve identification of a monogenic variant.

This pipeline can accept a singleton VCF file, but including one or both parent VCFs can highly increase the diagnostic yield. In step 3, adding the affectation status of the parents (if diagnosed with IBD) and the sex of the proband will also increase the likelihood of identifying a disease causing variant.

The pipeline produces a table of a short list of variants that meet the filtration criteria (output filename `{patient_id}_filtered.csv`). In our cohort, VCFs with 100,000-150,000 variant calls were filtered down to an average of 15 variant calls per patient without removing any of the known disease causing variants. In addition, the pipeline also produces a full table without any variants removed, but with full annotation and inheritance modeling as above, to allow for further well-informed manual curation and analysis of variants (output filename `{patient_id}_unfiltered.csv`).
