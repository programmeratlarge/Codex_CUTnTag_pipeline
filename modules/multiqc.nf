process MULTIQC_RAW {
    tag 'raw'
    label 'small'
    publishDir "${params.outdir}/13_multiqc/raw", mode: 'copy'

    input:
    path inputs

    output:
    path 'multiqc_raw.html', emit: report
    path 'multiqc_raw_data', emit: data

    script:
    """
    multiqc . --force --filename multiqc_raw.html --outdir .
    """
}

process MULTIQC_TRIM {
    tag 'trimmed'
    label 'small'
    publishDir "${params.outdir}/13_multiqc/trimmed", mode: 'copy'

    input:
    path fastqc_inputs
    path cutadapt_reports

    output:
    path 'multiqc_trimmed.html', emit: report
    path 'multiqc_trimmed_data', emit: data

    script:
    """
    multiqc . --force --filename multiqc_trimmed.html --outdir .
    """
}

process MULTIQC_ALIGNMENT {
    tag 'alignment_filtering'
    label 'small'
    publishDir "${params.outdir}/13_multiqc/alignment_filtering", mode: 'copy'

    input:
    path inputs

    output:
    path 'multiqc_alignment_filtering.html', emit: report
    path 'multiqc_alignment_filtering_data', emit: data

    script:
    """
    multiqc . --force --filename multiqc_alignment_filtering.html --outdir .
    """
}

process MULTIQC_FINAL {
    tag 'final'
    label 'small'
    publishDir "${params.outdir}/13_multiqc/final", mode: 'copy'

    input:
    path inputs

    output:
    path 'multiqc_cutntag_final.html', emit: report
    path 'multiqc_cutntag_final_data', emit: data

    script:
    """
    multiqc . \\
      --force \\
      --config ${projectDir}/assets/multiqc_config.yaml \\
      --filename multiqc_cutntag_final.html \\
      --outdir .
    """
}

