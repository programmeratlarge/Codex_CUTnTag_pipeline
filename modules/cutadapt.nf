process CUTADAPT_PE {
    tag { meta.sample_id }
    label 'medium'

    publishDir "${params.outdir}/02_trimmed", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.sample_id}.trimmed.R1.fastq.gz"), path("${meta.sample_id}.trimmed.R2.fastq.gz"), emit: fastq
    tuple val(meta), path("${meta.sample_id}.cutadapt.log"), emit: report

    script:
    """
    cutadapt \\
      --cores ${task.cpus} \\
      -a '${params.adapter_fwd}' \\
      -A '${params.adapter_rev}' \\
      -m ${params.min_read_length} \\
      ${params.cutadapt_extra} \\
      -o ${meta.sample_id}.trimmed.R1.fastq.gz \\
      -p ${meta.sample_id}.trimmed.R2.fastq.gz \\
      ${r1} ${r2} \\
      > ${meta.sample_id}.cutadapt.log
    """
}

