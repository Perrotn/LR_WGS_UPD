#!/bin/bash

# A4_TrioMix.sh
# Per-trio worker: resolves the proband/mother/father BAMs and invokes triomix.py.
# The rerun gate (check for existing plot.pdf) lives in A2; A4 always runs TrioMix.

# Default values
sample_data=""
virtual_env=""
patient_ID=""
output_dir=""
genome_fasta=""
nthreads=1
triomix_py="" # git clone https://github.com/cjyoon/triomix.git

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_data)
        sample_data="$2"
        shift 2
        ;;
        --patient_ID)
        patient_ID="$2"
        shift 2
        ;;
        --genome_fasta)
        genome_fasta="$2"
        shift 2
        ;;
        --virtual_env)
        virtual_env="$2"
        shift 2
        ;;
        --nthreads)
        nthreads="$2"
        shift 2
        ;;
        --output_dir)
        output_dir="$2"
        shift 2
        ;;
        --triomix_py)
        triomix_py="$2"
        shift 2
        ;;
        *)
            echo "Unknown parameter passed: $1"
            exit 1
            ;;
    esac
done

# Check for required arguments
if [[ -z "$sample_data" || -z "$patient_ID" || -z "$virtual_env" || -z "$output_dir" || -z "$genome_fasta" || -z "$triomix_py" ]]; then
    echo "Error: Missing one or more required arguments."
    echo "Usage: $0 --sample_data <file> --patient_ID <character> --genome_fasta <file> --virtual_env <dir> --output_dir <dir> --triomix_py <file> [--nthreads <int>]"
    exit 1
fi


module load python/3.13.2
module load samtools/1.22.1
module load r/4.5.0



# Activate the virtual environment
source "$virtual_env"

echo "Running script with the following parameters:"
echo "Sample Data: $sample_data"
echo "Patient ID: $patient_ID"
echo "Genome FASTA: $genome_fasta"
echo "Output Directory: $output_dir"
echo "Threads: $nthreads"
echo "TrioMix Script: $triomix_py"



echo "Processing: $patient_ID"

unset COL_NUM

declare -A COL_NUM
count=0

IFS=$'\t' read -r -a header < "$sample_data"
for name in "${header[@]}"; do
  ((count++))
  COL_NUM[$name]=$count
done


family=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $1=="'$patient_ID'"{print $('${COL_NUM["Trio"]}')}' $sample_data | sort -u)

# A2 normally guarantees a non-empty Trio ID; this guard catches the case where
# A4 is invoked directly with a proband whose CSV row has no family ID.
if [[ -z "$family" ]]; then
    echo "Error, empty Trio for $patient_ID"
    exit 1
fi

father_ID=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM["Trio"]}')=="'$family'" && $('${COL_NUM["formatted_role"]}')=="father" {print $1}' $sample_data | sort -u)
mother_ID=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM["Trio"]}')=="'$family'" && $('${COL_NUM["formatted_role"]}')=="mother" {print $1}' $sample_data | sort -u)

proband_bam=$(bash ${HOME}/projects/ctb-rallard/COMMUN/common_scripts/find_file_PacBioData.sh \
--base_dir="/home/${USER}/projects/ctb-rallard/COMMUN/PacBioData/S3-Storage" \
--dir_variations="/${family}/proband/,/${patient_ID}/" \
--suffix_dir="_NO_SUFFIX_" \
--file_variations="${patient_ID}.GRCh38.haplotagged.bam,${patient_ID}.haplotagged.bam")

mother_bam=$(bash ${HOME}/projects/ctb-rallard/COMMUN/common_scripts/find_file_PacBioData.sh \
--base_dir="/home/${USER}/projects/ctb-rallard/COMMUN/PacBioData/S3-Storage" \
--dir_variations="/${family}/mother/,/${mother_ID}/" \
--suffix_dir="_NO_SUFFIX_" \
--file_variations="${mother_ID}.GRCh38.haplotagged.bam,${mother_ID}.haplotagged.bam")

father_bam=$(bash ${HOME}/projects/ctb-rallard/COMMUN/common_scripts/find_file_PacBioData.sh \
--base_dir="/home/${USER}/projects/ctb-rallard/COMMUN/PacBioData/S3-Storage" \
--dir_variations="/${family}/father/,/${father_ID}/" \
--suffix_dir="_NO_SUFFIX_" \
--file_variations="${father_ID}.GRCh38.haplotagged.bam,${father_ID}.haplotagged.bam")


if [[ ! -s "$father_bam" ]] || [[ ! -s "$mother_bam" ]] || [[ ! -s "$proband_bam" ]]
then
    echo "Error, missing bam file, patientID:" $patient_ID
    rm -rf $output_dir
    exit 1
fi


rm -rf $output_dir
mkdir -p $output_dir

eval "python $triomix_py \
--father $father_bam \
--mother $mother_bam \
--child $proband_bam \
--reference $genome_fasta \
--thread $nthreads \
--snp "${HOME}/projects/ctb-rallard/COMMUN/common_software/triomix/common_snp/grch38_common_snp.bed.gz" \
--output_dir $output_dir"

echo "Triomix done"
