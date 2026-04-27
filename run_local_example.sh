#!/usr/bin/env bash
set -euo pipefail

# Review and edit this script before running. It intentionally does not execute
# automatically when the pipeline package is generated.

USER_NAME="${USER:-qisun}"
PROJECT_DIR="/workdir/${USER_NAME}/cutntag-pe-nextflow"
WORK_DIR="/workdir/${USER_NAME}/nextflow_work/cutntag_pe"
TMPDIR="/workdir/${USER_NAME}/tmp"

mkdir -p "${WORK_DIR}" "${TMPDIR}"
export TMPDIR

module load nextflow/25.4.3

nextflow run "${PROJECT_DIR}/main.nf" \
  -profile singularity,local \
  -work-dir "${WORK_DIR}" \
  --container "/workdir/${USER_NAME}/containers/cutntag-pe.sif" \
  --input_dir "/workdir/${USER_NAME}/data/fastq" \
  --outdir "/workdir/${USER_NAME}/results/cutntag" \
  --genome hg38 \
  --bowtie2_index "/workdir/${USER_NAME}/refs/hg38/bowtie2/hg38" \
  --chrom_sizes "/workdir/${USER_NAME}/refs/hg38/hg38.chrom.sizes" \
  --blacklist "/workdir/${USER_NAME}/refs/hg38/hg38.blacklist.bed" \
  --annotation_gtf "/workdir/${USER_NAME}/refs/hg38/gencode.gtf" \
  --tss_bed "/workdir/${USER_NAME}/refs/hg38/gencode.tss.bed" \
  --association_csv "/workdir/${USER_NAME}/metadata/association.csv" \
  --threads 8 \
  --normalize_using CPM \
  --dedup_mode picard \
  -resume
