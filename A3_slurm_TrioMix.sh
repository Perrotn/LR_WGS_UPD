#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-rallard
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --output=slurm_TrioMix_%A_%a.out
#SBATCH --error=slurm_TrioMix_%A_%a.err

# A3_slurm_TrioMix.sh
# Array dispatcher: one task per trio proband. Reads its patient ID from the todo
# TSV using SLURM_ARRAY_TASK_ID as the line number, then calls A4_TrioMix.sh.

todo_sample_data=""
clean_sample_data=""
version=""
my_virtualenv=""
genome_fasta=""
triomix_py=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --todo_sample_data)
            todo_sample_data="$2"
            shift 2
            ;;
        --clean_sample_data)
            clean_sample_data="$2"
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
        --genome_fasta)
            genome_fasta="$2"
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

line_number=$((SLURM_ARRAY_TASK_ID))
patient_ID=$(sed -n "${line_number}p" "$todo_sample_data" | cut -f1)

echo "Processing trio proband: $patient_ID (Task ID: $SLURM_ARRAY_TASK_ID)"

output_dir="${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_UPD_report/patients_data_v${version}/${patient_ID}"

bash "${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_UPD_report/A4_TrioMix.sh" \
    --sample_data  "$clean_sample_data" \
    --patient_ID   "$patient_ID" \
    --genome_fasta "$genome_fasta" \
    --virtual_env  "$my_virtualenv" \
    --nthreads     "$SLURM_CPUS_PER_TASK" \
    --output_dir   "$output_dir" \
    --triomix_py   "$triomix_py"
