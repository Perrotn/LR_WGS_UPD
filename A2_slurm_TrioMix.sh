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


skipped_sample_data="${HOME}/scratch/temp/temp_sample_data_UPD_skipped_v${version}.tsv"

export CLEAN_SAMPLE_DATA="$clean_sample_data"
export UNIQUE_SAMPLE_DATA="$unique_sample_data"
export SKIPPED_SAMPLE_DATA="$skipped_sample_data"

# Identify complete trios.
# Accepts: sing_trio_duo in {Trio, QUAD} AND formatted_role==proband AND
#          Trio cell is a real family ID (not empty / "#N/A" / "NA" / "N/A" / "#NA") AND
#          the same Trio ID has at least one row with formatted_role==father
#          AND at least one row with formatted_role==mother.
# Every field is .strip()ed before comparison. Pattern mirrors
# LR_WGS_STR/denovo1_run_STR_denovo.sh lines 52-92.
python3 <<PYEOF
import csv, os, sys

clean = os.environ["CLEAN_SAMPLE_DATA"]
unique_out = os.environ["UNIQUE_SAMPLE_DATA"]
skipped_out = os.environ["SKIPPED_SAMPLE_DATA"]

INVALID_TRIO_IDS = {"", "#N/A", "NA", "N/A", "#NA"}
ACCEPTED_FAMILY_TYPES = {"Trio", "QUAD"}

# Pass 1: read all rows with stripped fields, build by-Trio role map.
rows = []
# Sample CSV often contains Latin-1 accented characters in HPO/Commentaire fields;
# latin-1 decodes any byte stream losslessly, and the ASCII fields we filter on
# (PatientID/Trio/formatted_role/sing_trio_duo) are unaffected.
with open(clean, newline="", encoding="latin-1") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        rows.append({k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()})

by_trio_roles = {}
for r in rows:
    trio = r.get("Trio", "")
    role = r.get("formatted_role", "")
    pid = r.get("PatientID", "")
    if trio in INVALID_TRIO_IDS or not pid:
        continue
    by_trio_roles.setdefault(trio, set()).add(role)

accepted = set()
skipped_lonely = 0
skipped_incomplete = 0
skipped_empty_id = 0

with open(skipped_out, "w", newline="", encoding="utf-8") as fh_sk:
    sk = csv.writer(fh_sk, delimiter="\t", lineterminator="\n")
    sk.writerow(["PatientID", "reason"])
    for r in rows:
        sing = r.get("sing_trio_duo", "")
        role = r.get("formatted_role", "")
        if sing not in ACCEPTED_FAMILY_TYPES or role != "proband":
            continue  # Not a candidate at all; not a "skip" for diagnostic purposes.

        pid = r.get("PatientID", "")
        trio = r.get("Trio", "")

        if not pid:
            skipped_empty_id += 1
            sk.writerow(["", "empty_PatientID"])
            continue
        if trio in INVALID_TRIO_IDS:
            skipped_lonely += 1
            sk.writerow([pid, f"lonely:trio_id_invalid({trio or 'empty'})"])
            continue

        roles_in_trio = by_trio_roles.get(trio, set())
        missing = []
        if "father" not in roles_in_trio:
            missing.append("father")
        if "mother" not in roles_in_trio:
            missing.append("mother")
        if missing:
            skipped_incomplete += 1
            sk.writerow([pid, f"incomplete_trio:missing_{'+'.join(missing)}"])
            continue

        accepted.add(pid)

with open(unique_out, "w", encoding="utf-8") as fh_uq:
    for pid in sorted(accepted):
        fh_uq.write(pid + "\n")

print(f"Trio identification: accepted={len(accepted)} skipped_lonely={skipped_lonely} skipped_incomplete={skipped_incomplete} skipped_empty_id={skipped_empty_id}")
PYEOF
status=$?
if [[ $status -ne 0 ]]; then
    echo "Error: trio validator failed (exit $status)"
    exit $status
fi

if [[ ! -s "$unique_sample_data" ]]; then
    echo "Error: no complete trios found in $clean_sample_data (see $skipped_sample_data for reasons)"
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
