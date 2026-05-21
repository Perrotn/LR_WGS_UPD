# LR_WGS_UPD — Uniparental Disomy detection pipeline

Per-trio UPD detection from PacBio HiFi long-read whole-genome sequencing data, using [**TrioMix**](https://github.com/cjyoon/triomix). Produces the plots, segment calls, and contamination/UPD-fraction TSVs consumed by the tertiary analysis Shiny app under its **UPD** tab.

Runs on **Alliance Canada HPC** under SLURM account `def-rallard`.

---

## Pipeline overview

```
Input: trio haplotagged BAMs (PacBio)  +  sample metadata CSV
              │
   ┌──────────▼──────────┐
   │      A-series        │  TrioMix UPD detection (one SLURM task per trio, ≤10 concurrent)
   │  A2 orchestrator     │  Cleans CSV, filters Trio probands, applies update_only
   │                      │  rerun gate, builds todo list, submits A3 array
   │  A3 array (per-pt)   │  24 h / 24 CPU / 32 GB per trio
   │  A4 TrioMix worker   │  bcftools mpileup → MLE in R → plot + segments + summary TSVs
   └──────────┬──────────┘
              │
        Shiny app reads {patient_ID}.child.counts.* outputs
```

---

## Quick start

```bash
# Submit every Trio proband for both pipeline versions (incremental — skips
# probands whose plot PDF already exists). This is the default entry point.
bash A1_run_Triomix.sh

# Submit a single version directly (incremental)
sbatch --output=slurm_TrioMix_v0.9.1.out \
    A2_slurm_TrioMix.sh \
    --sample_data    ~/projects/ctb-rallard/COMMUN/LR_WGS_STR_Report/Reference_files/Bioinfo-LR_SampleData.csv \
    --version        0.9.1 \
    --my_virtualenv  ~/projects/ctb-rallard/pernic02/my_virtualenv/bin/activate \
    --update_only    true

# Force a full rerun — every Trio proband becomes an array task; per-patient
# output dir is wiped before TrioMix runs.
sbatch --output=slurm_TrioMix_v0.9.1.out \
    A2_slurm_TrioMix.sh \
    --sample_data    <...> --version 0.9.1 \
    --my_virtualenv  <...> --update_only false

# Run one trio locally without SLURM (handy for debugging)
bash A4_TrioMix.sh \
    --sample_data   ~/scratch/triomix/clean_sample_data.tsv \
    --patient_ID    HSJ-010-03 \
    --genome_fasta  ~/projects/ctb-rallard/COMMUN/Data_resources/hifi-wdl-resources-v2.0.0/GRCh38/human_GRCh38_no_alt_analysis_set.fasta \
    --virtual_env   ~/projects/ctb-rallard/pernic02/my_virtualenv/bin/activate \
    --triomix_py    ~/projects/ctb-rallard/COMMUN/common_software/triomix/triomix.py \
    --output_dir    ~/projects/ctb-rallard/COMMUN/LR_WGS_UPD_report/patients_data_v0.9.1/HSJ-010-03 \
    --nthreads      24
```

Per-task SLURM logs land in the working directory as `slurm_TrioMix_{job_id}_{task_id}.out` / `.err`.

---

## Repository layout

| File | Role |
|---|---|
| `A1_run_Triomix.sh` | Entry point — submits one A2 SLURM job per pipeline version |
| `A2_slurm_TrioMix.sh` | Orchestrator — cleans the sample CSV, filters Trio probands, applies the `--update_only` rerun gate, writes the todo TSV, submits A3 as an array |
| `A3_slurm_TrioMix.sh` | Array dispatcher — one SLURM task per trio (24 h / 24 CPU / 32 GB); reads its patient ID via `sed -n "${SLURM_ARRAY_TASK_ID}p"` and calls A4 |
| `A4_TrioMix.sh` | Per-trio worker — resolves the three trio BAMs via `find_file_PacBioData.sh`, then invokes `triomix.py` |
| `CLAUDE.md` | Guidance for the Claude Code agent working in this repo |

---

## Sample metadata

Canonical input: a sample CSV (semicolon-delimited) such as `Bioinfo-LR_SampleData.csv`. A2 normalizes it to TSV via `common_scripts/clean_convert_file.sh` before parsing.

Columns the pipeline reads:

| Column | Use |
|---|---|
| `PatientID` (col 1) | Primary key. Selects proband rows. |
| `Trio` | Family ID, used to locate the matching father/mother rows. |
| `sing_trio_duo` | Must equal `Trio` for the patient to be processed. |
| `formatted_role` | `proband` / `father` / `mother` — selects the proband to run, then resolves the parent BAMs. |

Effective filter: `sing_trio_duo == "Trio" AND formatted_role == "proband"`.

Trio BAMs are looked up under `~/projects/ctb-rallard/COMMUN/PacBioData/S3-Storage/`, trying both `{family}/{role}/` and `{patient_ID}/` directory layouts and both `{ID}.GRCh38.haplotagged.bam` and `{ID}.haplotagged.bam` filename variants. If any of the three BAMs is missing, A4 deletes the patient's output dir and exits 1.

---

## Outputs

Per-patient outputs land in `patients_data_v{version}/{patient_ID}/`:

| File | Description |
|---|---|
| `{ID}.child.counts.plot.pdf` | UPD visualization. **This file is the rerun-gate sentinel** — A2 skips probands whose plot PDF already exists when `--update_only true`. |
| `{ID}.child.counts.summary.tsv` | Joint-estimate contamination / UPD fractions. |
| `{ID}.child.counts.upd.segments.tsv` | Called UPD segments. |
| `{ID}.child.counts` | Raw allele counts. |
| `{ID}.x2a.depth.tsv` | chrX vs autosome depth. |
| `{ID}_tmp/` | Per-window mpileup intermediates (safe to delete after a successful run). |

---

## Version JSON

Tool paths and reference-file paths are read at runtime from the sister project's version JSON:

```
~/projects/ctb-rallard/COMMUN/LR_WGS_tertiary_analysis_report/versions/
    tertiary_analysis_report_v{version}.json
```

A2 reads two keys with `jq`:

- `Report_Versions.UPD_Analysis.TrioMix.Path` — TrioMix Python entry point
- `Report_Versions.UPD_Analysis.TrioMix.Reference_Files.Reference_Genome.Path` — GRCh38 FASTA

Paths in the JSON use `${HOME}` literally; A4 lets the shell expand it via `eval`. **Never hard-code absolute paths in pipeline scripts — add them to the JSON instead.**

---

## Module + virtualenv requirements

A4 loads these modules on the compute node:

```bash
module load python/3.13.2
module load samtools/1.22.1
module load r/4.5.0
source ~/projects/ctb-rallard/pernic02/my_virtualenv/bin/activate
```

The virtualenv provides the Python dependencies that `triomix.py` imports.

---

## Upstream tool

TrioMix is vendored at `~/projects/ctb-rallard/COMMUN/common_software/triomix/` — a fork of <https://github.com/cjyoon/triomix>. Companion files used by A4:

- `triomix.py` — Python orchestrator
- `common_snp/grch38_common_snp.bed.gz` — common SNP BED passed via `--snp`
- `mle.R` / `mle_parent.R` / `plot_variant*.R` — maximum-likelihood and plotting helpers
