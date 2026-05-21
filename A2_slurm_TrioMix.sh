#!/bin/bash
#
#SBATCH --time=01:00:00
#SBATCH --account=def-rallard
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# A2_slurm_TrioMix.sh
# Orchestrator: converts sample CSV to TSV, filters to Trio probands, applies the
# update_only rerun gate, and submits A3 as a SLURM array (one task per trio).

sample_data=""
update_only="false"
version=""
my_virtualenv=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_data)
            sample_data="$2"
            shift 2
            ;;
        --update_only)
            update_only="$2"
            shift 2
            ;;
        --version)
            version="$2"
            shift 2
            ;;
        --my_virtualenv)
            my_virtualenv="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter passed: $1"
            exit 1
            ;;
    esac
done


version_json=${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_tertiary_analysis_report/versions/tertiary_analysis_report_v${version}.json

triomix_path=$(jq -r '.Report_Versions.UPD_Analysis.TrioMix.Path' "$version_json")
genome_fasta=$(jq -r '.Report_Versions.UPD_Analysis.TrioMix.Reference_Files.Reference_Genome.Path' "$version_json")


mkdir -p ${HOME}/scratch/triomix/
mkdir -p ${HOME}/scratch/temp/

clean_sample_data="${HOME}/scratch/triomix/clean_sample_data.tsv"
unique_sample_data="${HOME}/scratch/temp/temp_sample_data_UPD_unique_v${version}.tsv"
todo_sample_data="${HOME}/scratch/temp/temp_sample_data_UPD_todo_v${version}.tsv"

bash ${HOME}/projects/ctb-rallard/COMMUN/common_scripts/clean_convert_file.sh \
    --input "$sample_data" \
    --output "$clean_sample_data" \
    --temp-folder ${HOME}/scratch


unset COL_NUM

declare -A COL_NUM
count=0

IFS=$'\t' read -r -a header < "$clean_sample_data"
for name in "${header[@]}"; do
  ((count++))
  COL_NUM[$name]=$count
done


# Filter to Trio probands only
awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM["sing_trio_duo"]}')=="Trio" && $('${COL_NUM["formatted_role"]}')=="proband" {print $1}' "$clean_sample_data" | sort -u > "$unique_sample_data"

if [[ ! -s "$unique_sample_data" ]]; then
    echo "Error: no Trio probands found in $clean_sample_data"
    exit 1
fi


output_base="${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_UPD_report/patients_data_v${version}"

# When update_only=true, drop probands whose plot PDF already exists
if [[ "$update_only" == "true" ]]; then
    > "$todo_sample_data"
    while IFS= read -r pid; do
        plot_pdf="${output_base}/${pid}/${pid}.child.counts.plot.pdf"
        if [[ -s "$plot_pdf" ]]; then
            echo "Skipping (already done): $pid"
        else
            echo "$pid" >> "$todo_sample_data"
        fi
    done < "$unique_sample_data"
else
    cp "$unique_sample_data" "$todo_sample_data"
fi

num_tasks=$(wc -l < "$todo_sample_data")

if [[ "$num_tasks" -le 0 ]]; then
    echo "No tasks to run. All probands already processed."
    exit 0
fi

echo "Submitting array job with $num_tasks tasks."

sbatch --array=1-${num_tasks}%10 \
    --job-name="UPD_TrioMix_v${version//./_}" \
    "${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_UPD_report/A3_slurm_TrioMix.sh" \
    --todo_sample_data  "$todo_sample_data" \
    --clean_sample_data "$clean_sample_data" \
    --version           "$version" \
    --my_virtualenv     "$my_virtualenv" \
    --genome_fasta      "$genome_fasta" \
    --triomix_py        "$triomix_path"

echo "A2 orchestrator done."
