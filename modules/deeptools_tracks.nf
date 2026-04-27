process BAMCOVERAGE_BW {
    tag { meta.sample_id }
    label 'medium'

    publishDir "${params.outdir}/06_bigwig", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample_id}.bw"), emit: bigwig

    script:
    def eff = params.effective_genome_size ? "--effectiveGenomeSize ${params.effective_genome_size}" : ''
    """
    bamCoverage \\
      --bam ${bam} \\
      --outFileName ${meta.sample_id}.bw \\
      --outFileFormat bigwig \\
      --numberOfProcessors ${task.cpus} \\
      --binSize ${params.bigwig_bin_size} \\
      --normalizeUsing ${params.normalize_using} \\
      ${eff} \\
      ${params.bamcoverage_extra}
    """
}

process BAMCOVERAGE_BEDGRAPH {
    tag { meta.sample_id }
    label 'medium'

    publishDir "${params.outdir}/07_bedgraph", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample_id}.bedgraph.gz"), path("${meta.sample_id}.bedgraph.gz.tbi"), emit: bedgraph

    script:
    def eff = params.effective_genome_size ? "--effectiveGenomeSize ${params.effective_genome_size}" : ''
    """
    if [[ '${params.make_bedgraph}' == 'true' ]]; then
      bamCoverage \\
        --bam ${bam} \\
        --outFileName ${meta.sample_id}.bedgraph \\
        --outFileFormat bedgraph \\
        --numberOfProcessors ${task.cpus} \\
        --binSize ${params.bigwig_bin_size} \\
        --normalizeUsing ${params.normalize_using} \\
        ${eff} \\
        ${params.bamcoverage_extra}
      if [[ -s '${params.chrom_sizes}' ]]; then
        bedtools sort -faidx '${params.chrom_sizes}' -i ${meta.sample_id}.bedgraph | bgzip -c > ${meta.sample_id}.bedgraph.gz
      else
        sort -k1,1 -k2,2n ${meta.sample_id}.bedgraph | bgzip -c > ${meta.sample_id}.bedgraph.gz
      fi
      tabix -p bed ${meta.sample_id}.bedgraph.gz || touch ${meta.sample_id}.bedgraph.gz.tbi
    else
      touch ${meta.sample_id}.bedgraph.gz
      touch ${meta.sample_id}.bedgraph.gz.tbi
    fi
    """
}
