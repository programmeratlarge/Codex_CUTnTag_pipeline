#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SCAN_FASTQ_PAIRS }       from './modules/scan_pairs'
include { VALIDATE_ASSOCIATIONS }  from './modules/validate_associations'
include { FASTQC as FASTQC_RAW }   from './modules/fastqc'
include { FASTQC as FASTQC_TRIM }  from './modules/fastqc'
include { MULTIQC_RAW; MULTIQC_TRIM; MULTIQC_ALIGNMENT; MULTIQC_FINAL } from './modules/multiqc'
include { COUNT_FASTQ as COUNT_RAW }  from './modules/count_fastq'
include { COUNT_FASTQ as COUNT_TRIM } from './modules/count_fastq'
include { CUTADAPT_PE }            from './modules/cutadapt'
include { BOWTIE2_ALIGN }          from './modules/bowtie2_align'
include { ALIGNMENT_QC }           from './modules/alignment_qc'
include { FILTER_BAM }             from './modules/filter_bam'
include { DEDUP_BAM }              from './modules/dedup_bam'
include { READ_RETENTION }         from './modules/read_retention'
include { BAMCOVERAGE_BW; BAMCOVERAGE_BEDGRAPH } from './modules/deeptools_tracks'
include { MERGE_GROUP_BAMS }       from './modules/merge_group_bams'
include { MACS3_SAMPLE_PEAKS; MACS3_GROUP_PEAKS } from './modules/macs3_peaks'
include { CONSENSUS_PEAKS }        from './modules/consensus_peaks'
include { ANNOTATE_PEAKS as ANNOTATE_SAMPLE_PEAKS } from './modules/annotate_peaks'
include { ANNOTATE_PEAKS as ANNOTATE_GROUP_PEAKS }  from './modules/annotate_peaks'
include { ANNOTATE_PEAKS as ANNOTATE_CONSENSUS_PEAKS } from './modules/annotate_peaks'
include { FRIP }                   from './modules/frip'
include { DEEPTOOLS_TSS; DEEPTOOLS_PEAKS; TSS_ENRICHMENT } from './modules/deeptools_qc'
include { TOOL_VERSIONS }          from './modules/tool_versions'


def requireParam(String name) {
    if (!params[name]) {
        throw new IllegalArgumentException("Missing required parameter --${name}")
    }
}

def requirePath(String name, boolean directory=false) {
    requireParam(name)
    def p = file(params[name])
    if (!p.exists()) {
        throw new IllegalArgumentException("Parameter --${name} does not exist: ${params[name]}")
    }
    if (directory && !p.isDirectory()) {
        throw new IllegalArgumentException("Parameter --${name} must be a directory: ${params[name]}")
    }
}

def validateParams() {
    requirePath('input_dir', true)
    requireParam('outdir')
    requireParam('genome')
    requireParam('bowtie2_index')
    requirePath('chrom_sizes')
    requirePath('association_csv')

    if (params.annotation_gtf) {
        requirePath('annotation_gtf')
    }
    if (params.tss_bed) {
        requirePath('tss_bed')
    }
    if (params.blacklist) {
        requirePath('blacklist')
    }
    if (params.normalize_using == 'RPGC' && !params.effective_genome_size) {
        throw new IllegalArgumentException("--effective_genome_size is required when --normalize_using RPGC")
    }
    if (params.umi_enabled && !['umi_tools', 'none'].contains(params.dedup_mode)) {
        throw new IllegalArgumentException("--umi_enabled requires --dedup_mode umi_tools or none")
    }
}


workflow CUTNTAG_PE {
    validateParams()

    SCAN_FASTQ_PAIRS(params.input_dir, params.paired_pattern ?: '')
    VALIDATE_ASSOCIATIONS(
        SCAN_FASTQ_PAIRS.out.pairs,
        file(params.association_csv),
        params.allow_extra_association_rows,
        params.allow_control_free_peak_calling,
        params.call_control_peaks
    )

    samples_ch = VALIDATE_ASSOCIATIONS.out.samplesheet
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            def meta = [
                sample_id: row.sample_id,
                species: row.species,
                genome: row.genome,
                antibody: row.antibody,
                condition: row.condition,
                replicate: row.replicate,
                group_id: row.group_id,
                is_control: row.is_control.toBoolean(),
                control_group_id: row.control_group_id,
                resolved_control_merge_group_id: row.resolved_control_merge_group_id,
                merge_group_id: row.merge_group_id,
                peak_calling_mode: row.peak_calling_mode
            ]
            tuple(meta, file(row.r1), file(row.r2))
        }

    TOOL_VERSIONS()

    FASTQC_RAW(samples_ch)
    COUNT_RAW(samples_ch, 'raw')
    raw_fastqc_files = FASTQC_RAW.out.zip.map { meta, zips -> zips }.flatten()
    raw_count_files = COUNT_RAW.out.counts.map { meta, counts -> counts }
    MULTIQC_RAW(raw_fastqc_files.collect())

    CUTADAPT_PE(samples_ch)
    COUNT_TRIM(CUTADAPT_PE.out.fastq, 'trimmed')
    FASTQC_TRIM(CUTADAPT_PE.out.fastq)
    trimmed_fastqc_files = FASTQC_TRIM.out.zip.map { meta, zips -> zips }.flatten()
    trimmed_count_files = COUNT_TRIM.out.counts.map { meta, counts -> counts }
    cutadapt_report_files = CUTADAPT_PE.out.report.map { meta, report -> report }
    MULTIQC_TRIM(trimmed_fastqc_files.collect(), cutadapt_report_files.collect())

    BOWTIE2_ALIGN(CUTADAPT_PE.out.fastq)
    ALIGNMENT_QC(BOWTIE2_ALIGN.out.bam)

    FILTER_BAM(BOWTIE2_ALIGN.out.bam, params.blacklist ?: '')
    DEDUP_BAM(FILTER_BAM.out.bam)
    alignment_stats_files = ALIGNMENT_QC.out.stats.flatten()
    filter_count_files = FILTER_BAM.out.counts.map { meta, counts -> counts }
    dedup_count_files = DEDUP_BAM.out.counts.map { meta, counts -> counts }
    dedup_metric_files = DEDUP_BAM.out.metrics.map { meta, metrics -> metrics }

    READ_RETENTION(
        raw_count_files.collect(),
        trimmed_count_files.collect(),
        filter_count_files.collect(),
        dedup_count_files.collect()
    )

    BAMCOVERAGE_BW(DEDUP_BAM.out.bam)
    BAMCOVERAGE_BEDGRAPH(DEDUP_BAM.out.bam)

    final_bams_by_id = DEDUP_BAM.out.bam.map { meta, bam, bai -> tuple(meta.sample_id, meta, bam, bai) }

    group_member_ch = VALIDATE_ASSOCIATIONS.out.group_members
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            tuple(
                row.sample_id,
                row.merge_group_id,
                row.is_control,
                row.control_group_id,
                row.resolved_control_merge_group_id,
                row.peak_calling_mode,
                row.genome,
                row.species,
                row.antibody,
                row.condition,
                row.replicate
            )
        }

    group_inputs_ch = group_member_ch
        .join(final_bams_by_id, by: 0)
        .map { sample_id, merge_group_id, is_control, control_group_id, resolved_control_merge_group_id,
                peak_calling_mode, genome, species, antibody, condition, replicate, sample_meta, bam, bai ->
            def member = [
                sample_id: sample_id,
                genome: genome,
                species: species,
                antibody: antibody,
                condition: condition,
                replicate: replicate,
                is_control: is_control,
                control_group_id: control_group_id,
                resolved_control_merge_group_id: resolved_control_merge_group_id,
                peak_calling_mode: peak_calling_mode
            ]
            tuple(merge_group_id, member, bam, bai)
        }
        .groupTuple(by: 0)
        .map { merge_group_id, members, bams, bais ->
            def first = members[0]
            def meta = [
                group_id: merge_group_id,
                genome: first.genome,
                species: first.species,
                antibody: first.antibody,
                condition: first.condition,
                is_control: first.is_control.toBoolean(),
                control_group_id: first.control_group_id,
                resolved_control_merge_group_id: first.resolved_control_merge_group_id,
                peak_calling_mode: first.peak_calling_mode,
                sample_ids: members.collect { it.sample_id }.join(';')
            ]
            tuple(meta, bams, bais)
        }

    MERGE_GROUP_BAMS(group_inputs_ch)

    MACS3_SAMPLE_PEAKS(
        DEDUP_BAM.out.bam.map { meta, bam, bai -> bam }.collect(),
        DEDUP_BAM.out.bam.map { meta, bam, bai -> bai }.collect(),
        VALIDATE_ASSOCIATIONS.out.sample_peak_jobs
    )

    MACS3_GROUP_PEAKS(
        MERGE_GROUP_BAMS.out.bam.map { meta, bam, bai -> bam }.collect(),
        MERGE_GROUP_BAMS.out.bam.map { meta, bam, bai -> bai }.collect(),
        VALIDATE_ASSOCIATIONS.out.group_peak_jobs
    )

    CONSENSUS_PEAKS(MACS3_GROUP_PEAKS.out.peaks, VALIDATE_ASSOCIATIONS.out.consensus_jobs)

    annotation_mqc_ch = Channel.empty()
    if (params.annotation_gtf) {
        ANNOTATE_SAMPLE_PEAKS(MACS3_SAMPLE_PEAKS.out.peaks, file(params.annotation_gtf), params.tss_bed ?: '')
        ANNOTATE_GROUP_PEAKS(MACS3_GROUP_PEAKS.out.peaks, file(params.annotation_gtf), params.tss_bed ?: '')
        ANNOTATE_CONSENSUS_PEAKS(CONSENSUS_PEAKS.out.peaks, file(params.annotation_gtf), params.tss_bed ?: '')
        annotation_mqc_ch = ANNOTATE_SAMPLE_PEAKS.out.annotation
            .mix(ANNOTATE_GROUP_PEAKS.out.annotation)
            .mix(ANNOTATE_CONSENSUS_PEAKS.out.annotation)
    }

    FRIP(
        DEDUP_BAM.out.bam.map { meta, bam, bai -> bam }.collect(),
        DEDUP_BAM.out.bam.map { meta, bam, bai -> bai }.collect(),
        MERGE_GROUP_BAMS.out.bam.map { meta, bam, bai -> bam }.collect(),
        MERGE_GROUP_BAMS.out.bam.map { meta, bam, bai -> bai }.collect(),
        MACS3_SAMPLE_PEAKS.out.peaks,
        MACS3_GROUP_PEAKS.out.peaks,
        CONSENSUS_PEAKS.out.peaks,
        VALIDATE_ASSOCIATIONS.out.samplesheet,
        VALIDATE_ASSOCIATIONS.out.sample_peak_jobs,
        VALIDATE_ASSOCIATIONS.out.group_peak_jobs
    )

    tss_mqc_ch = Channel.empty()
    tss_deeptools_ch = Channel.empty()
    if (params.tss_bed) {
        DEEPTOOLS_TSS(BAMCOVERAGE_BW.out.bigwig.map { meta, bw -> bw }.collect(), file(params.tss_bed))
        TSS_ENRICHMENT(BAMCOVERAGE_BW.out.bigwig.map { meta, bw -> bw }.collect(), file(params.tss_bed))
        tss_mqc_ch = TSS_ENRICHMENT.out.files
        tss_deeptools_ch = DEEPTOOLS_TSS.out.files
    }

    DEEPTOOLS_PEAKS(
        BAMCOVERAGE_BW.out.bigwig.map { meta, bw -> bw }.collect(),
        CONSENSUS_PEAKS.out.peaks
    )

    MULTIQC_ALIGNMENT(
        alignment_stats_files
            .mix(filter_count_files)
            .mix(dedup_count_files)
            .mix(dedup_metric_files)
            .collect()
    )

    final_mqc_ch = raw_fastqc_files
        .mix(trimmed_fastqc_files)
        .mix(cutadapt_report_files)
        .mix(alignment_stats_files)
        .mix(filter_count_files)
        .mix(dedup_count_files)
        .mix(dedup_metric_files)
        .mix(READ_RETENTION.out.files.flatten())
        .mix(MACS3_SAMPLE_PEAKS.out.metrics.flatten())
        .mix(MACS3_GROUP_PEAKS.out.metrics.flatten())
        .mix(CONSENSUS_PEAKS.out.metrics.flatten())
        .mix(FRIP.out.files.flatten())
        .mix(TOOL_VERSIONS.out.versions)
        .mix(annotation_mqc_ch)
        .mix(tss_mqc_ch.flatten())
        .mix(tss_deeptools_ch.flatten())
        .mix(DEEPTOOLS_PEAKS.out.files.flatten())
        .collect()

    MULTIQC_FINAL(final_mqc_ch)
}

workflow {
    CUTNTAG_PE()
}
