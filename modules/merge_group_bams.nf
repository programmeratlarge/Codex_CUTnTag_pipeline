process MERGE_GROUP_BAMS {
    tag { meta.group_id }
    label 'large'

    publishDir "${params.outdir}/09_peaks_merged_groups/merged_bams", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${meta.group_id}.merged.bam"), path("${meta.group_id}.merged.bam.bai"), emit: bam
    path "${meta.group_id}.members.tsv", emit: members

    script:
    def members = meta.sample_ids.split(';').collect { "${it}" }.join('\\n')
    """
    printf 'sample_id\\n${members}\\n' > ${meta.group_id}.members.tsv

    if [[ ${bams.size()} -eq 1 ]]; then
      cp ${bams[0]} ${meta.group_id}.merged.bam
      cp ${bais[0]} ${meta.group_id}.merged.bam.bai
    else
      samtools merge -@ ${task.cpus} -f ${meta.group_id}.unsorted.bam ${bams.join(' ')}
      samtools sort -@ ${task.cpus} -o ${meta.group_id}.merged.bam ${meta.group_id}.unsorted.bam
      samtools index -@ ${task.cpus} ${meta.group_id}.merged.bam
    fi
    """
}

