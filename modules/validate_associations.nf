process VALIDATE_ASSOCIATIONS {
    tag 'validate_associations'
    label 'small'

    publishDir "${params.outdir}/00_fastq_pairs", mode: 'copy'

    input:
    path fastq_pairs
    path association_csv
    val allow_extra_rows
    val allow_control_free
    val call_control_peaks

    output:
    path 'samplesheet.validated.tsv', emit: samplesheet
    path 'group_members.tsv', emit: group_members
    path 'sample_peak_jobs.tsv', emit: sample_peak_jobs
    path 'group_peak_jobs.tsv', emit: group_peak_jobs
    path 'consensus_jobs.tsv', emit: consensus_jobs
    path 'metadata.validation.json', emit: metadata_json
    path 'association_validation_report.txt', emit: report

    script:
    """
    python ${projectDir}/bin/validate_associations.py \\
      --fastq-pairs ${fastq_pairs} \\
      --association-csv ${association_csv} \\
      --outdir . \\
      ${allow_extra_rows ? '--allow-extra-association-rows' : ''} \\
      ${allow_control_free ? '--allow-control-free-peak-calling' : ''} \\
      ${call_control_peaks ? '--call-control-peaks' : ''}
    """
}
