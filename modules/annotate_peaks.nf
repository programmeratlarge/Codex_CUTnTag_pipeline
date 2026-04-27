process ANNOTATE_PEAKS {
    tag { peaks_dir.baseName }
    label 'medium'

    publishDir "${params.outdir}/10_peak_annotation", mode: 'copy'

    input:
    path peaks_dir
    path annotation_gtf
    val tss_bed

    output:
    path "${peaks_dir.baseName}_annotation", emit: annotation

    script:
    """
    python ${projectDir}/bin/annotate_peaks.py \\
      --peaks ${peaks_dir} \\
      --gtf ${annotation_gtf} \\
      --tss-bed '${tss_bed}' \\
      --outdir ${peaks_dir.baseName}_annotation \\
      --promoter-upstream ${params.promoter_upstream} \\
      --promoter-downstream ${params.promoter_downstream} \\
      --downstream-window ${params.downstream_window}
    """
}
