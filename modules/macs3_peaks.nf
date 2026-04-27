process MACS3_SAMPLE_PEAKS {
    tag 'per_sample'
    label 'large'

    publishDir "${params.outdir}/08_peaks_per_sample", mode: 'copy'

    input:
    path bams
    path bais
    path sample_peak_jobs

    output:
    path 'per_sample_peaks', emit: peaks
    path '*peak_counts*', emit: metrics

    script:
    """
    python ${projectDir}/bin/macs3_batch_callpeak.py \\
      --level sample \\
      --jobs ${sample_peak_jobs} \\
      --bam-dir . \\
      --outdir per_sample_peaks \\
      --qvalue ${params.macs3_qvalue} \\
      --broad-cutoff ${params.macs3_broad_cutoff} \\
      --genome-size '${params.macs3_genome_size}' \\
      --extra '${params.macs3_extra}' \\
      --threads ${task.cpus} \\
      ${params.call_control_peaks ? '--call-control-peaks' : ''}
    cp per_sample_peaks/peak_counts.tsv per_sample_peak_counts.tsv
    cp per_sample_peaks/peak_counts_mqc.tsv per_sample_peak_counts_mqc.tsv
    """
}

process MACS3_GROUP_PEAKS {
    tag 'merged_groups'
    label 'large'

    publishDir "${params.outdir}/09_peaks_merged_groups", mode: 'copy'

    input:
    path bams
    path bais
    path group_peak_jobs

    output:
    path 'group_peaks', emit: peaks
    path '*peak_counts*', emit: metrics

    script:
    """
    python ${projectDir}/bin/macs3_batch_callpeak.py \\
      --level group \\
      --jobs ${group_peak_jobs} \\
      --bam-dir . \\
      --outdir group_peaks \\
      --qvalue ${params.macs3_qvalue} \\
      --broad-cutoff ${params.macs3_broad_cutoff} \\
      --genome-size '${params.macs3_genome_size}' \\
      --extra '${params.macs3_extra}' \\
      --threads ${task.cpus}
    cp group_peaks/peak_counts.tsv group_peak_counts.tsv
    cp group_peaks/peak_counts_mqc.tsv group_peak_counts_mqc.tsv
    """
}
