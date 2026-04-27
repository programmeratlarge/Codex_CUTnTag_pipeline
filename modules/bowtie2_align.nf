process BOWTIE2_ALIGN {
    tag { meta.sample_id }
    label 'large'

    publishDir "${params.outdir}/04_alignment", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.sample_id}.sorted.bam"), path("${meta.sample_id}.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.sample_id}.bowtie2.log"), emit: log

    script:
    """
    bowtie2 \\
      -x '${params.bowtie2_index}' \\
      -1 ${r1} \\
      -2 ${r2} \\
      -p ${task.cpus} \\
      ${params.bowtie2_options} \\
      2> ${meta.sample_id}.bowtie2.log \\
    | samtools view -@ ${task.cpus} -bS - \\
    | samtools sort -@ ${task.cpus} -o ${meta.sample_id}.sorted.bam -

    samtools index -@ ${task.cpus} ${meta.sample_id}.sorted.bam
    """
}

