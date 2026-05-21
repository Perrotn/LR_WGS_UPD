# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Module Purpose

Per-patient UPD (Uniparental Disomy) detection from PacBio HiFi trio BAMs using **TrioMix**. Produces the `.child.counts.plot.pdf`, `.child.counts.summary.tsv`, and `.child.counts.upd.segments.tsv` consumed by the Shiny tertiary analysis dashboard (`LR_WGS_tertiary_analysis_report/`) under its UPD tab.

Runs on **Alliance Canada HPC** under SLURM account `def-rallard`. Only proband rows where `sing_trio_duo == "Trio"` are processed. Repository: <https://github.com/Perrotn/LR_WGS_UPD>.

---

## Key Commands

```bash
# Submit all trio probands for both pipeline versions (current default behaviour of A1).
# Each version submits a short orchestrator job that fans out a SLURM array
# (one task per trio, up to 10 concurrent).
bash [A1_run_Triomix.sh](A1_run_Triomix.sh)

# Submit a single version directly (skips probands whose plot.pdf already exists)
sbatch --output=slurm_TrioMix_v0.9.1.out \
    [A2_slurm_TrioMix.sh](A2_slurm_TrioMix.sh) \
    --sample_data    ~/projects/ctb-rallard/COMMUN/LR_WGS_STR_Report/Reference_files/Bioinfo-LR_SampleData.csv \
    --version        0.9.1 \
    --my_virtualenv  ~/projects/ctb-rallard/pernic02/my_virtualenv/bin/activate \
    --update_only    true

# Force a full rerun (every proband gets a fresh array task; per-patient output dir is wiped before TrioMix runs)
sbatch --output=slurm_TrioMix_v0.9.1.out \
    [A2_slurm_TrioMix.sh](A2_slurm_TrioMix.sh) \
    --sample_data <...> --version 0.9.1 \
    --my_virtualenv <...> --update_only false

# Run a single patient locally (no SLURM); useful for debugging one trio
bash [A4_TrioMix.sh](A4_TrioMix.sh) \
    --sample_data   ~/scratch/triomix/clean_sample_data.tsv \
    --patient_ID    HSJ-010-03 \
    --genome_fasta  ~/projects/ctb-rallard/COMMUN/Data_resources/hifi-wdl-resources-v2.0.0/GRCh38/human_GRCh38_no_alt_analysis_set.fasta \
    --virtual_env   ~/projects/ctb-rallard/pernic02/my_virtualenv/bin/activate \
    --triomix_py    ~/projects/ctb-rallard/COMMUN/common_software/triomix/triomix.py \
    --output_dir    ~/projects/ctb-rallard/COMMUN/LR_WGS_UPD_report/patients_data_v0.9.1/HSJ-010-03 \
    --nthreads      24
```

`A1_run_Triomix.sh` currently submits versions 0.9 and 0.9.1 in sequence; edit the `version=` lines to add/remove versions.

---

## Architecture

```
A1_run_Triomix.sh      → submits one A2 SLURM job per pipeline version
A2_slurm_TrioMix.sh    → orchestrator (1 h / 1 CPU / 2 GB): reads version JSON, cleans sample CSV,
                          filters Trio probands, applies update_only rerun gate, writes todo TSV,
                          submits A3 as --array=1-N%10
A3_slurm_TrioMix.sh    → array dispatcher (per-task 24 h / 24 CPU / 32 GB): reads its patient ID
                          via sed -n "${SLURM_ARRAY_TASK_ID}p" from the todo TSV, calls A4
A4_TrioMix.sh          → per-trio worker: resolves trio BAMs, runs triomix.py
```

A2 is now a lightweight orchestrator — the heavy work runs as a SLURM array of per-trio tasks with concurrency capped at `%10`. The rerun gate (skip probands whose `.child.counts.plot.pdf` already exists) lives in A2, so completed patients never occupy an array slot. A4 always runs TrioMix; if it was already done, A2 would not have enqueued it.

### Version JSON as single source of truth

Tool paths and reference files are read at runtime from:
```
~/projects/ctb-rallard/COMMUN/LR_WGS_tertiary_analysis_report/versions/
    tertiary_analysis_report_v{version}.json
```

A2 reads two keys with `jq`:
- `Report_Versions.UPD_Analysis.TrioMix.Path` → triomix Python entry point
- `Report_Versions.UPD_Analysis.TrioMix.Reference_Files.Reference_Genome.Path` → GRCh38 FASTA

Both paths use `${HOME}` literally — A4 passes them to `eval "python ..."`, which expands the variable. Never hard-code absolute paths; add new tools/files to the JSON instead.

The common SNP BED is currently hard-coded in [A4_TrioMix.sh](A4_TrioMix.sh) (the `--snp` line of the `python $triomix_py …` call): `~/projects/ctb-rallard/COMMUN/common_software/triomix/common_snp/grch38_common_snp.bed.gz`. If this needs to vary by version, promote it into the version JSON.

### Sample metadata flow

Starts from `Bioinfo-LR_SampleData.csv` (semicolon-delimited). A2 normalizes it to TSV via `common_scripts/clean_convert_file.sh` → `~/scratch/triomix/clean_sample_data.tsv`. Required columns:

| Column | Use |
|--------|-----|
| `PatientID` (col 1) | primary key; selects proband rows |
| `Trio` | family ID — used to find matching father/mother rows |
| `sing_trio_duo` | must equal `Trio` for the patient to be processed |
| `formatted_role` | `proband` / `father` / `mother` — drives BAM lookup |

A2 writes the filtered list of Trio probands to `~/scratch/temp/temp_sample_data_UPD_unique_v{version}.tsv` and the post-rerun-gate todo list to `~/scratch/temp/temp_sample_data_UPD_todo_v{version}.tsv`. The array reads line `${SLURM_ARRAY_TASK_ID}` of the todo TSV.

BAM resolution (in A4) uses `common_scripts/find_file_PacBioData.sh` against `~/projects/ctb-rallard/COMMUN/PacBioData/S3-Storage`, trying both `{family}/{role}/` and `{patient_ID}/` directory layouts and both `{ID}.GRCh38.haplotagged.bam` and `{ID}.haplotagged.bam` filename variants. If **any** of the three BAMs (proband/father/mother) is missing, A4 deletes the patient's output dir and exits with status 1 — this is intentional cleanup, not an error to suppress.

### Incremental rerun gate

A2 skips a proband when **all** of:
- `--update_only true`
- `{output_dir}/{patient_ID}.child.counts.plot.pdf` exists and is non-empty

Otherwise the proband ID is added to the todo TSV and gets its own array task. A4 itself always wipes the output dir and runs TrioMix fresh.

### Per-patient outputs

`patients_data_v{version}/{patient_ID}/` contains:
- `*.child.counts` — raw allele counts
- `*.child.counts.plot.pdf` — UPD visualization (the A2 rerun gate)
- `*.child.counts.summary.tsv` — joint-estimate contamination/UPD fractions
- `*.child.counts.upd.segments.tsv` — called UPD segments
- `*.x2a.depth.tsv` — chrX vs autosome depth
- `*_tmp/` — per-window mpileup intermediates (safe to delete after success)

### TrioMix software

Vendored fork of <https://github.com/cjyoon/triomix> at `~/projects/ctb-rallard/COMMUN/common_software/triomix/` — Python orchestrator (`triomix.py`) plus R MLE scripts (`mle.R`, `mle_parent.R`, `plot_variant*.R`). Loads `python/3.13.2`, `samtools/1.22.1`, `r/4.5.0` plus the project virtualenv (which provides Python deps that triomix imports).

### Slurm logs

- **A2 (orchestrator)**: A1 redirects its stdout to `slurm_TrioMix_v{version}.out` in the working directory; pre-existing copies are removed first.
- **A3 array (per-task)**: each task writes `slurm_TrioMix_{job_id}_{task_id}.out` and `.err` in the working directory (declared with `%A_%a` in A3's `#SBATCH --output`/`--error`).

Failed patients can be inspected by grepping `Error, missing bam file` or `patient_ID:` across the per-task `.out` files.
