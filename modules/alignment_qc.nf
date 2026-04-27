process ALIGNMENT_QC {
    tag { meta.sample_id }
    label 'small'

    publishDir "${params.outdir}/04_alignment/qc", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "*.txt", emit: stats

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.sample_id}.flagstat.txt
    samtools idxstats ${bam} > ${meta.sample_id}.idxstats.txt
    samtools stats -@ ${task.cpus} ${bam} > ${meta.sample_id}.samtools_stats.txt
    """
}
