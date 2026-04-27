process FRIP {
    tag 'frip'
    label 'large'

    publishDir "${params.outdir}/11_frip", mode: 'copy'

    input:
    path sample_bams
    path sample_bais
    path group_bams
    path group_bais
    path sample_peaks
    path group_peaks
    path consensus_peaks
    path samplesheet
    path sample_peak_jobs
    path group_peak_jobs

    output:
    path '*', emit: files

    script:
    """
    python ${projectDir}/bin/calculate_frip.py \\
      --samplesheet ${samplesheet} \\
      --sample-peak-jobs ${sample_peak_jobs} \\
      --group-peak-jobs ${group_peak_jobs} \\
      --sample-bam-dir . \\
      --group-bam-dir . \\
      --sample-peaks ${sample_peaks} \\
      --group-peaks ${group_peaks} \\
      --consensus-peaks ${consensus_peaks} \\
      --outdir .
    """
}
