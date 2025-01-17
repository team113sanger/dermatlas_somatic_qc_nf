process MAF_TO_EXCEL {
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/maf:0.6.1"
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