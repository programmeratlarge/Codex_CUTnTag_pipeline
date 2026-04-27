process FILTER_BAM {
    tag { meta.sample_id }
    label 'medium'

    publishDir "${params.outdir}/05_filtering/pre_dedup", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(bam), path(bai)
    val blacklist_path

    output:
    tuple val(meta), path("${meta.sample_id}.filtered.bam"), path("${meta.sample_id}.filtered.bam.bai"), emit: bam
    tuple val(meta), path("${meta.sample_id}.filter_counts.tsv"), emit: counts

    script:
    """
    sample='${meta.sample_id}'
    mapq='${params.mapq}'
    mito_names='${params.mito_chroms}'
    blacklist='${blacklist_path}'

    aligned_fragments=\$(samtools view -@ ${task.cpus} -c -f 64 -F 3852 ${bam})

    samtools view -@ ${task.cpus} -bh -f 2 -F 3852 -q \${mapq} ${bam} \\
      | samtools sort -@ ${task.cpus} -o \${sample}.mapq_proper.bam -
    samtools index -@ ${task.cpus} \${sample}.mapq_proper.bam
    properly_paired_fragments=\$(samtools view -@ ${task.cpus} -c -f 64 \${sample}.mapq_proper.bam)

    samtools view -@ ${task.cpus} -h \${sample}.mapq_proper.bam \\
      | awk -v mt="\${mito_names}" 'BEGIN{n=split(mt,a,","); for(i=1;i<=n;i++) mito[a[i]]=1} /^@/ {print; next} !(\$3 in mito)' \\
      | samtools view -@ ${task.cpus} -bh - \\
      | samtools sort -@ ${task.cpus} -o \${sample}.no_mito.bam -
    samtools index -@ ${task.cpus} \${sample}.no_mito.bam
    mito_removed_fragments=\$(samtools view -@ ${task.cpus} -c -f 64 \${sample}.no_mito.bam)

    if [[ -n "\${blacklist}" && -s "\${blacklist}" ]]; then
      bedtools intersect -abam \${sample}.no_mito.bam -b "\${blacklist}" -v \\
        | samtools sort -@ ${task.cpus} -o ${meta.sample_id}.filtered.bam -
    else
      cp \${sample}.no_mito.bam ${meta.sample_id}.filtered.bam
    fi
    samtools index -@ ${task.cpus} ${meta.sample_id}.filtered.bam
    blacklist_filtered_fragments=\$(samtools view -@ ${task.cpus} -c -f 64 ${meta.sample_id}.filtered.bam)

    cat > ${meta.sample_id}.filter_counts.tsv <<EOF
sample_id	aligned_fragments	properly_paired_fragments	mapq_proper_fragments	mitochondrial_removed_fragments	blacklist_filtered_fragments	mapq	mito_chroms	blacklist
${meta.sample_id}	\${aligned_fragments}	\${properly_paired_fragments}	\${properly_paired_fragments}	\${mito_removed_fragments}	\${blacklist_filtered_fragments}	\${mapq}	\${mito_names}	\${blacklist}
EOF
    """
}

