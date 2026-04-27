process READ_RETENTION {
    tag 'read_retention'
    label 'small'

    publishDir "${params.outdir}/05_filtering/read_retention", mode: 'copy'

    input:
    path raw_counts
    path trimmed_counts
    path filter_counts
    path dedup_counts

    output:
    path '*', emit: files

    script:
    """
    python ${projectDir}/bin/make_read_retention.py \\
      --output read_retention.tsv \\
      --multiqc read_retention_mqc.tsv \\
      --raw ${raw_counts.join(' ')} \\
      --trimmed ${trimmed_counts.join(' ')} \\
      --filter ${filter_counts.join(' ')} \\
      --dedup ${dedup_counts.join(' ')}
    """
}
