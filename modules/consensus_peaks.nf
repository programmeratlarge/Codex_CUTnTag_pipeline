process CONSENSUS_PEAKS {
    tag 'consensus'
    label 'small'

    publishDir "${params.outdir}/09_peaks_merged_groups/consensus", mode: 'copy'

    input:
    path group_peaks
    path consensus_jobs

    output:
    path 'consensus_peaks', emit: peaks
    path '*peak_counts*', emit: metrics

    script:
    """
    python ${projectDir}/bin/make_consensus_peaks.py \\
      --group-peaks ${group_peaks} \\
      --consensus-jobs ${consensus_jobs} \\
      --outdir consensus_peaks
    cp consensus_peaks/consensus_peak_counts.tsv consensus_peak_counts.tsv
    cp consensus_peaks/consensus_peak_counts_mqc.tsv consensus_peak_counts_mqc.tsv
    """
}
