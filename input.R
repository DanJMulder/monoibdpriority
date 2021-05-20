## Input files
# OPTION 1: Provide a samplesheet and VCF directory
vcf_dir <- ""
family_units <- ""  # Must have the following columns: proband, sex, paternal, paternal_affected, maternal,
# maternal_affected
vcf_suffix <- ".hc.vcf"  # VCFs will assume to be {vcf_dir}/{id}{vcf_suffix}
# OPTION 2: Specify trio
proband_fn <- ""
proband_sex <- ""  # One of "male" or "female"
paternal_fn <- ""
paternal_affected <- ""  # One of "TRUE", "FALSE", or "unknown"
maternal_fn <- ""
maternal_affected <- ""  # One of "TRUE", "FALSE", or "unknown"

# Output, will be created if not present
output_dir <- ""

# Options
num.cpu <- 8
