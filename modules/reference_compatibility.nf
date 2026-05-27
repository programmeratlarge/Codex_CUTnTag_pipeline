process REFERENCE_COMPATIBILITY {
    tag 'reference_compatibility'
    label 'small'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path chrom_sizes
    val blacklist_path
    val annotation_gtf_path
    val tss_bed_path

    output:
    path 'reference_contig_compatibility.tsv', emit: tsv
    path 'reference_contig_compatibility.txt', emit: report

    script:
    """
    python ${projectDir}/bin/check_reference_contigs.py \\
      --bam ${bam} \\
      --chrom-sizes ${chrom_sizes} \\
      --blacklist '${blacklist_path}' \\
      --annotation-gtf '${annotation_gtf_path}' \\
      --tss-bed '${tss_bed_path}' \\
      --tsv reference_contig_compatibility.tsv \\
      --report reference_contig_compatibility.txt
    """
}

