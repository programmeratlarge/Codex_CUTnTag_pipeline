process TOOL_VERSIONS {
    tag 'versions'
    label 'small'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    path 'tool_versions.yml', emit: versions

    script:
    """
    set +e
    {
      echo 'tool_versions:'
      echo "  nextflow: \$(nextflow -version 2>/dev/null | head -n 1 | sed 's/: /:/' || true)"
      echo "  fastqc: \$(fastqc --version 2>&1 | head -n 1 || true)"
      echo "  multiqc: \$(multiqc --version 2>&1 | head -n 1 || true)"
      echo "  cutadapt: \$(cutadapt --version 2>&1 | head -n 1 || true)"
      echo "  bowtie2: \$(bowtie2 --version 2>&1 | head -n 1 || true)"
      echo "  samtools: \$(samtools --version 2>&1 | head -n 1 || true)"
      echo "  picard: \$(picard MarkDuplicates --version 2>&1 | head -n 1 || true)"
      echo "  umi_tools: \$(umi_tools --version 2>&1 | head -n 1 || true)"
      echo "  bedtools: \$(bedtools --version 2>&1 | head -n 1 || true)"
      echo "  deeptools: \$(bamCoverage --version 2>&1 | head -n 1 || true)"
      echo "  macs3: \$(macs3 --version 2>&1 | head -n 1 || true)"
    } > tool_versions.yml
    """
}

