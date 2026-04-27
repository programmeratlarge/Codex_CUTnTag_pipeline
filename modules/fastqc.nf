process FASTQC {
    tag { meta.sample_id }
    label 'small'

    publishDir { "${params.outdir}/${task.process.contains('FASTQC_RAW') ? '01_fastqc_raw' : '03_fastqc_trimmed'}" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path('*_fastqc.zip'), emit: zip
    tuple val(meta), path('*_fastqc.html'), emit: html

    script:
    """
    mkdir -p fastqc
    fastqc --threads ${task.cpus} --outdir fastqc ${r1} ${r2}
    mv fastqc/* .
    """
}
