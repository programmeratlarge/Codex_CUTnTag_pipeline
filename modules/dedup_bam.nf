process DEDUP_BAM {
    tag { meta.sample_id }
    label 'medium'

    publishDir "${params.outdir}/05_filtering/deduplicated", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample_id}.final.bam"), path("${meta.sample_id}.final.bam.bai"), emit: bam
    tuple val(meta), path("${meta.sample_id}.dedup_counts.tsv"), emit: counts
    tuple val(meta), path("${meta.sample_id}.dedup_metrics.txt"), emit: metrics

    script:
    """
    sample='${meta.sample_id}'
    input_fragments=\$(samtools view -@ ${task.cpus} -c -f 64 ${bam})

    if [[ '${params.umi_enabled}' == 'true' || '${params.dedup_mode}' == 'umi_tools' ]]; then
      umi_tools dedup \\
        -I ${bam} \\
        -S \${sample}.final.bam \\
        --output-stats \${sample}.umi_tools \\
        ${params.umi_tools_options} \\
        > \${sample}.dedup_metrics.txt 2>&1
      samtools index -@ ${task.cpus} \${sample}.final.bam
    elif [[ '${params.dedup_mode}' == 'picard' ]]; then
      picard MarkDuplicates \\
        I=${bam} \\
        O=\${sample}.final.bam \\
        M=\${sample}.dedup_metrics.txt \\
        REMOVE_DUPLICATES=true \\
        ASSUME_SORT_ORDER=coordinate \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=SILENT
      [[ -s \${sample}.final.bai ]] && mv \${sample}.final.bai \${sample}.final.bam.bai
      [[ -s \${sample}.final.bam.bai ]] || samtools index -@ ${task.cpus} \${sample}.final.bam
    elif [[ '${params.dedup_mode}' == 'none' ]]; then
      cp ${bam} \${sample}.final.bam
      cp ${bai} \${sample}.final.bam.bai
      printf 'Deduplication disabled\\n' > \${sample}.dedup_metrics.txt
    else
      echo "Unsupported --dedup_mode '${params.dedup_mode}'" >&2
      exit 2
    fi

    final_fragments=\$(samtools view -@ ${task.cpus} -c -f 64 \${sample}.final.bam)
    removed=\$(( input_fragments - final_fragments ))
    if [[ "\${input_fragments}" -gt 0 ]]; then
      dup_rate=\$(awk -v r="\${removed}" -v t="\${input_fragments}" 'BEGIN{printf "%.6f", r/t}')
    else
      dup_rate=0
    fi

    cat > \${sample}.dedup_counts.tsv <<EOF
sample_id	pre_dedup_fragments	final_usable_fragments	duplicate_removed_fragments	duplicate_rate	dedup_mode	umi_enabled
${meta.sample_id}	\${input_fragments}	\${final_fragments}	\${removed}	\${dup_rate}	${params.dedup_mode}	${params.umi_enabled}
EOF
    """
}

