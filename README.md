# Paired-end CUT&Tag Nextflow DSL2 Pipeline

Version: 0.1.0

This repository implements a production-oriented Nextflow DSL2 workflow for paired-end CUT&Tag data. It scans FASTQ pairs, validates sample metadata and treatment-control relationships, runs FastQC/MultiQC, trims adapters, aligns with Bowtie2, filters and deduplicates BAMs, creates BigWig/BedGraph tracks, calls MACS3 peaks per sample and merged group, builds consensus peaks, annotates peaks, calculates FRiP, computes deepTools/TSS enrichment QC, and produces integrated MultiQC reports.

## Quick start

On the BioHPC-style server, keep the project and all input/output paths under `/workdir/$USER`. Do not launch long Nextflow runs from this chat session; review `run_local_example.sh`, edit paths, then run it in a separate shell.

```bash
cd /workdir/$USER/cutntag-pe-nextflow
bash run_local_example.sh
```

The script uses:

```bash
module load nextflow/25.4.3
nextflow run main.nf -profile singularity,local -work-dir /workdir/$USER/nextflow_work/cutntag_pe ...
```

Singularity/Apptainer is recommended on HPC. Provide an all-in-one SIF with `--container /workdir/$USER/containers/cutntag-pe.sif`, or run with `-profile conda,local` after creating the conda environment. Docker is not recommended on this server; if you must build or inspect Docker images on this system, use the provided `docker1` wrapper, for example `docker1 images` or `docker1 build ...`.

## Required inputs

Required parameters:

```text
--input_dir        Directory containing paired-end FASTQ files
--outdir           Output directory
--genome           Genome identifier, e.g. hg38, mm10, rn6, custom
--bowtie2_index    Bowtie2 index basename
--chrom_sizes      Chromosome sizes file
--association_csv  Metadata/control/grouping CSV
```

Required for annotation and TSS enrichment:

```text
--annotation_gtf   GTF/GFF annotation file
--tss_bed          BED file of transcription start sites
```

Optional:

```text
--blacklist        BED file of blacklisted regions
--adapter_fwd      Forward adapter, default CTGTCTCTTATACACATCT
--adapter_rev      Reverse adapter, default CTGTCTCTTATACACATCT
--paired_pattern   Regex with named groups (?P<sample>...) and (?P<read>[12])
```

## Association CSV specification

The required `association_csv` must contain:

```csv
sample_id,species,genome,antibody,condition,replicate,group_id,is_control,control_group_id,merge_group_id,peak_calling_mode,notes
```

Column definitions:

- `sample_id`: Must exactly match the FASTQ-derived sample ID.
- `species`: `human`, `mouse`, `rat`, or `custom`.
- `genome`: Genome build such as `hg38`, `mm10`, `rn6`.
- `antibody`: Target antibody such as `H3K27ac`, `H3K4me3`, `IgG`, `Input`.
- `condition`: Biological condition, treatment, cell type, or experimental group.
- `replicate`: Biological replicate identifier.
- `group_id`: Unique biological sample or control group.
- `is_control`: `true` or `false`.
- `control_group_id`: For non-controls, references a control `group_id` or control `merge_group_id`.
- `merge_group_id`: Samples merged for group peak calling, usually species + antibody + condition.
- `peak_calling_mode`: `narrow`, `broad`, or `auto`. Auto uses broad for broad histone marks such as H3K27me3 and narrow otherwise.
- `notes`: Optional free text.

Validation fails early if FASTQ samples are missing from the CSV, extra rows are present without `--allow_extra_association_rows`, merge groups mix incompatible species/genome/antibody/condition values, or treatment groups lack a valid control unless `--allow_control_free_peak_calling true` is set.

See `examples_association.csv` for a minimal example.

## Output structure

```text
results/
  00_fastq_pairs/
  01_fastqc_raw/
  02_trimmed/
  03_fastqc_trimmed/
  04_alignment/
  05_filtering/
  06_bigwig/
  07_bedgraph/
  08_peaks_per_sample/
  09_peaks_merged_groups/
  10_peak_annotation/
  11_frip/
  12_deeptools/
  13_multiqc/
  logs/
  pipeline_info/
```

## Pipeline architecture

The workflow uses explicit validation-generated tables to keep treatment-control matching reproducible:

- `fastq_pairs.tsv`: scanner output with sample ID, R1, R2, input directory, detected pattern, and status.
- `samplesheet.validated.tsv`: FASTQ paths joined to association metadata.
- `sample_peak_jobs.tsv`: per-sample MACS3 treatment/control jobs.
- `group_members.tsv`: sample-to-merge-group map.
- `group_peak_jobs.tsv`: merged treatment/control peak-calling jobs.
- `consensus_jobs.tsv`: comparable group sets for union peaks and heatmaps.

This design avoids fragile implicit grouping in channel code and makes failures inspectable.

## Key defaults and recommendations

- Bowtie2 uses `--very-sensitive --dovetail --no-mixed --no-discordant`, a common paired-end CUT&Tag choice.
- Filtering keeps properly paired reads, removes low MAPQ alignments, removes mitochondrial reads, optionally removes blacklist overlaps, then deduplicates.
- Picard duplicate removal is the default. Use `--dedup_mode umi_tools --umi_enabled true` only when UMIs are present in read names or have been extracted upstream.
- BigWigs default to CPM. RPGC is available but requires `--effective_genome_size`.
- MACS3 uses BAMPE mode. Controls are strongly recommended for CUT&Tag, especially IgG or input controls.

## Installation

First check locally installed packages on the cluster software page:

```text
https://biohpc.cornell.edu/lab/userguide.aspx?a=software
```

If using conda/mamba:

```bash
mamba env create -f environment.yml
mamba activate cutntag-pe-nextflow
```

On the target server, prefer:

```bash
module load nextflow/25.4.3
nextflow run main.nf -profile singularity,local ...
```

## Example command

```bash
nextflow run /workdir/$USER/cutntag-pe-nextflow/main.nf \
  -profile singularity,local \
  -work-dir /workdir/$USER/nextflow_work/cutntag_pe \
  --input_dir /workdir/$USER/data/fastq \
  --outdir /workdir/$USER/results/cutntag \
  --genome hg38 \
  --bowtie2_index /workdir/$USER/refs/hg38/bowtie2/hg38 \
  --chrom_sizes /workdir/$USER/refs/hg38/hg38.chrom.sizes \
  --blacklist /workdir/$USER/refs/hg38/hg38.blacklist.bed \
  --annotation_gtf /workdir/$USER/refs/hg38/gencode.gtf \
  --tss_bed /workdir/$USER/refs/hg38/gencode.tss.bed \
  --association_csv /workdir/$USER/metadata/association.csv \
  --threads 8 \
  -resume
```

## MultiQC custom content

The pipeline writes MultiQC-compatible TSVs for read retention, peak counts, FRiP, annotation categories, and TSS enrichment. `assets/multiqc_config.yaml` registers these sections, while each `_mqc.tsv` file also includes MultiQC custom-content headers for portability.

## Testing strategy

For a minimal synthetic test, create a tiny chromosome, Bowtie2 index, chromosome sizes file, TSS BED, GTF, and paired FASTQs with a few thousand simulated fragments. Run the pipeline with `--allow_control_free_peak_calling true` first to validate mechanics, then add a small IgG/control group to test control-aware MACS3. Use `nextflow run ... -resume` after intentional restarts to verify resumability.

## Version history

- 0.1.0: Initial production-oriented DSL2 implementation for paired-end CUT&Tag.
