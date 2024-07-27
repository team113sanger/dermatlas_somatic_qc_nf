
process CALCULATE_SAMPLE_TMB {
    publishDir "${params.outdir}/${params.release_version}/${meta.analysis_type}/plots_${file_id}", mode: params.publish_dir_mode
    input:
    tuple val(meta), path(maf_file)
    val(exome_size)

    output:
    path("mutations_per_Mb.tsv"), emit: tmb

    shell:
    file_id = maf_file.name.split("_caveman")[0]
    """
    for sample in \$(cut -f 11 !{maf_file} | grep PD | sort -u); do
        echo "Processing sample: \$sample"
        echo -ne "\${sample}\t" >> mutations_per_Mb.tsv
        muts="\$(grep "\${sample}" !{maf_file} | cut -f 4,5 | sort -u | wc -l )"
        echo "\${muts}/!{exome_size}" | bc -l >> mutations_per_Mb.tsv
    done 
    """

}

process MAF_TO_EXCEL {
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/maf:0.5.1"
    publishDir "${params.outdir}/${params.release_version}/${meta.analysis_type}", mode: params.publish_dir_mode
    
    input: 
    tuple val(meta), path(maf)

    output: 
    path("*.xlsx")

    script:
    """
    Rscript /opt/repo/maf2xlsx.R $maf
    """

}