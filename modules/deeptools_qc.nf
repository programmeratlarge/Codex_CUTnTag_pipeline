process DEEPTOOLS_TSS {
    tag 'tss_matrix'
    label 'large'

    publishDir "${params.outdir}/12_deeptools", mode: 'copy'

    input:
    path bigwigs
    path tss_bed

    output:
    path '*', emit: files

    script:
    """
    mkdir -p matrices profiles heatmaps
    computeMatrix reference-point \\
      --referencePoint TSS \\
      -S ${bigwigs.join(' ')} \\
      -R ${tss_bed} \\
      --beforeRegionStartLength ${params.tss_window} \\
      --afterRegionStartLength ${params.tss_window} \\
      --numberOfProcessors ${task.cpus} \\
      --missingDataAsZero \\
      -out matrices/tss_matrix.gz
    plotProfile -m matrices/tss_matrix.gz -out profiles/tss_profile.pdf --perGroup
    plotHeatmap -m matrices/tss_matrix.gz -out heatmaps/tss_heatmap.pdf --colorMap RdBu
    """
}

process DEEPTOOLS_PEAKS {
    tag 'peak_matrix'
    label 'large'

    publishDir "${params.outdir}/12_deeptools", mode: 'copy'

    input:
    path bigwigs
    path consensus_peaks

    output:
    path '*', emit: files

    script:
    """
    mkdir -p matrices profiles heatmaps
    find ${consensus_peaks} -name '*.bed' -type f | sort > consensus_peak_files.txt
    if [[ -s consensus_peak_files.txt ]]; then
      computeMatrix reference-point \\
        --referencePoint center \\
        -S ${bigwigs.join(' ')} \\
        -R \$(cat consensus_peak_files.txt) \\
        --beforeRegionStartLength ${params.peak_center_window} \\
        --afterRegionStartLength ${params.peak_center_window} \\
        --numberOfProcessors ${task.cpus} \\
        --missingDataAsZero \\
        -out matrices/peak_center_matrix.gz
      plotProfile -m matrices/peak_center_matrix.gz -out profiles/peak_center_profile.pdf --perGroup
      plotHeatmap -m matrices/peak_center_matrix.gz -out heatmaps/peak_center_heatmap.pdf --colorMap RdBu
    else
      touch matrices/peak_center_matrix.gz profiles/peak_center_profile.pdf heatmaps/peak_center_heatmap.pdf
    fi
    """
}

process TSS_ENRICHMENT {
    tag 'tss_enrichment'
    label 'medium'

    publishDir "${params.outdir}/12_deeptools/tss_enrichment", mode: 'copy'

    input:
    path bigwigs
    path tss_bed

    output:
    path '*', emit: files

    script:
    """
    python ${projectDir}/bin/tss_enrichment.py \\
      --bigwigs ${bigwigs.join(' ')} \\
      --tss-bed ${tss_bed} \\
      --outdir . \\
      --window ${params.tss_window} \\
      --center-window ${params.tss_center_window}
    """
}
