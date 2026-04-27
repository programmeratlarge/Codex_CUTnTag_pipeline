process SCAN_FASTQ_PAIRS {
    tag 'scan_fastq_pairs'
    label 'small'

    publishDir "${params.outdir}/00_fastq_pairs", mode: 'copy'

    input:
    val input_dir
    val paired_pattern

    output:
    path 'fastq_pairs.tsv', emit: pairs

    script:
    """
    python ${projectDir}/bin/scan_fastq_pairs.py \\
      --input-dir '${input_dir}' \\
      --paired-pattern '${paired_pattern}' \\
      --output fastq_pairs.tsv
    """
}
