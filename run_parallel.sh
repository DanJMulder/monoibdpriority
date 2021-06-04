#!/usr/bin/env bash
#PBS -l walltime=23:59:59
#PBS -l nodes=1:ppn=12
#PBS -l mem=32g,vmem=32g
#PBS -N SNP_Prioritize
#PBS -j oe
#PBS -i 3
#PBS -o logs
set -x

cd $PBS_O_WORKDIR

# Normalize paths
vcf=$(readlink -f "$1")
manifest=$(readlink -f "$2")
output_dir=$(readlink -f "$3")
dependencies=$(readlink -f "$4")
r_script=$(readlink -f "$5")

# Simple error check
if [ ! -r "$vcf" ]; then
        echo "VCF file (${vcf}) is not readable.";
        exit 1;
elif [ ! -r "$manifest" ]; then
        echo "Manifest file (${manifest}) is not readable.";
        exit 1;
elif [ -z "$PBS_ARRAYID" ] && [ -z "$TEST" ]; then
        echo "This script must be submitted via qsub with -t";
        echo "Usage: qsub -F \"project_vcf manifest out_dir dependencies.R\" -t 1-100%12 run_prioritization.sh";
        exit 1;
elif [ ! -r "$dependencies" ]; then
        echo "The specified dependencies file ($dependencies) is not readable.";
        exit 1;
elif [ ! -r "$r_script" ]; then
        echo "The R script ($r_script) is not readable.";
        exit 1;
fi

# Ensure commands succeed after this point
set -eo pipefail

# Make the output directory, if missing
if [ ! -e "$output_dir" ]; then
        mkdir "${output_dir}";
elif [ ! -d "$output_dir" ] || [ ! -r "$output_dir" ]; then
        echo "${output_dir} exists, but is not a directory or readable}!";
        exit 2;
fi
mkdir -p "${output_dir}/filtered"
mkdir -p "${output_dir}/unfiltered"
mkdir -p "${output_dir}/dbnsfp"

# Get the line containing what we need to run
proband_id=$(sed "${PBS_ARRAYID}q;d" "${manifest}" | cut -f 1)
proband_name=$(sed "${PBS_ARRAYID}q;d" "${manifest}" | cut -f 2)
proband_gender=$(sed "${PBS_ARRAYID}q;d" "${manifest}" | cut -f 3)
maternal_name=$(sed "${PBS_ARRAYID}q;d" "${manifest}" | cut -f 4)
maternal_status=$(sed "${PBS_ARRAYID}q;d" "${manifest}" | cut -f 5)
paternal_name=$(sed "${PBS_ARRAYID}q;d" "${manifest}" | cut -f 6)
paternal_status=$(sed "${PBS_ARRAYID}q;d" "${manifest}" | cut -f 7)

# Go to working directory
# In a task array, the job ID looks like 1234[1], so the folder is unique
#workdir="/localhd/${PBS_JOBID}/"
workdir="$(readlink -f ./out/)/work/${proband_id}/"
mkdir -p "${workdir}" && cd "${workdir}"

# Load modules late
module load anaconda3 bcftools/1.11 java

source activate ibd_script

# Make individual VCFs from the project VCF file
# 1. Extract sample in question, drop orphaned alt alleles
# 2. Missing genotypes filled with reference
# 3. Genotype must have at least one alt allele
if [ ! -r proband.vcf ]; then
  bcftools view -s "${proband_name}" -m2 -M2 -O u "${vcf}" | \
    bcftools +setGT - -- -t . -n 0 | \
    bcftools view -O v -o proband.vcf -i 'GT="alt"' -;
fi
if [ ! -r maternal.vcf ] && [ -n "${maternal_name}" ]; then
  bcftools view -s "${maternal_name}" -m2 -M2 -O u "${vcf}" | \
    bcftools +setGT - -- -t . -n 0 | \
    bcftools view -O v -o maternal.vcf -i 'GT="alt"' -;
fi
if [ ! -r paternal.vcf ] && [ -n "${paternal_name}" ]; then
  bcftools view -s "${paternal_name}" -m2 -M2 -O u "${vcf}" | \
    bcftools +setGT - -- -t . -n 0 | \
    bcftools view -O v -o paternal.vcf -i 'GT="alt"' -;
fi

# Make the input file
cat > input.R <<EOF
vcf_dir <- ""
family_units <- ""

proband_fn <- "${workdir}/proband.vcf"
proband_sex <- "${proband_gender}"
paternal_fn <- ifelse(nchar("${paternal_name}"), "${workdir}/paternal.vcf", "")
paternal_affected <- "${paternal_status}"
maternal_fn <- ifelse(nchar("${maternal_name}"), "${workdir}/maternal.vcf", "")
maternal_affected <- "${maternal_status}"

output_dir <- "${workdir}"
num.cpu <- 12
EOF

Rscript "${r_script}" "input.R" "${dependencies}"

mv "proband_filtered.csv" "${output_dir}/filtered/${proband_name}_filtered.csv"
mv "proband_unfiltered.csv" "${output_dir}/unfiltered/${proband_name}_unfiltered.csv"
mv "proband_dbnsfp_full.tsv" "${output_dir}/dbnsfp/${proband_name}_dbnsfp.tsv"

exit 0;
