process COUNT_FASTQ {
    tag { "${meta.sample_id}:${stage}" }
    label 'small'

    publishDir { "${params.outdir}/logs/read_counts/${stage}" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    val stage

    output:
    tuple val(meta), path("${meta.sample_id}.${stage}.fastq_counts.tsv"), emit: counts

    script:
    """
    python ${projectDir}/bin/count_fastq_reads.py \\
      --sample-id '${meta.sample_id}' \\
      --r1 ${r1} \\
      --r2 ${r2} \\
      --stage '${stage}' \\
      --output '${meta.sample_id}.${stage}.fastq_counts.tsv'
    """
}
