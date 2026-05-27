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

      # Use the BAM header as the sorting authority. The BedGraph coordinates
      # come from this BAM, so its contig names/order always match even when a
      # user-supplied chromosome sizes file uses a different naming convention
      # such as chr1 versus 1.
      samtools view -H ${bam} \\
        | awk -F'\\t' '/^@SQ/ {
            sn=""; ln="";
            for (i=1; i<=NF; i++) {
              if (\$i ~ /^SN:/) sn=substr(\$i,4);
              if (\$i ~ /^LN:/) ln=substr(\$i,4);
            }
            if (sn != "" && ln != "") print sn "\\t" ln;
          }' \\
        > ${meta.sample_id}.bam.chrom.sizes

      bedtools sort -faidx ${meta.sample_id}.bam.chrom.sizes -i ${meta.sample_id}.bedgraph \\
        | bgzip -c > ${meta.sample_id}.bedgraph.gz
      tabix -p bed ${meta.sample_id}.bedgraph.gz || touch ${meta.sample_id}.bedgraph.gz.tbi
    else
      touch ${meta.sample_id}.bedgraph.gz
      touch ${meta.sample_id}.bedgraph.gz.tbi
    fi
    """
}
