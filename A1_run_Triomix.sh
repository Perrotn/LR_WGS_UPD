#!/bin/bash
#
# A1_run_Triomix.sh
# Thin entry point. Submits one A2 SLURM orchestrator job per pipeline version.
#
# Optional flags:
#   --version <v>       Only submit for this version (e.g. 0.9.1). When omitted,
#                       submits for every version listed in default_versions below.
#   --update_only <t|f> Defaults to "true" (skip probands whose plot.pdf exists).
#                       Pass "false" to force a full rerun of every accepted trio.

sample_data="${HOME}/projects/ctb-rallard/COMMUN/Data_resources/Cohorts/Bioinfo-LR_SampleData.csv"
my_virtualenv="/home/pernic02/projects/ctb-rallard/pernic02/my_virtualenv/bin/activate"

default_versions=("0.9" "0.9.1")
version=""
update_only="true"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --version)     version="$2";     shift 2 ;;
        --update_only) update_only="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
done

if [[ -n "$version" ]]; then
    versions_to_run=("$version")
else
    versions_to_run=("${default_versions[@]}")
fi

for v in "${versions_to_run[@]}"; do
    rm -f "${HOME}/scratch/slurm_TrioMix_v${v}.out"
    sbatch --output="slurm_TrioMix_v${v}.out" \
        "${HOME}/projects/ctb-rallard/COMMUN/LR_WGS_UPD_report/A2_slurm_TrioMix.sh" \
        --sample_data    "$sample_data" \
        --version        "$v" \
        --my_virtualenv  "$my_virtualenv" \
        --update_only    "$update_only"
done

# grep "HSJ-059-03" /home/pernic02/scratch/slurm_TrioMix_v0.9.1.out
# grep "HSJ-059-03" /home/pernic02/scratch/slurm_TrioMix_v0.9.out
