// Resolve a subcohort's publish-dir name: user-provided legacy map wins,
// otherwise fall back to the raw analysis_type key. Covariate mode adds a
// second layer (with_covariates / without_covariates) when present.
def dndscvPublishDir(meta) {
    def name = params.dndscv_subcohort_names?.get(meta.analysis_type) ?: meta.analysis_type
    def mode = meta.covariate_mode ? "/${meta.covariate_mode}" : ''
    "${params.dndscv_outdir}/${params.release_version}/${name}${mode}"
}


process MAF_TO_DNDSCV_INPUT {
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/maf:0.6.5"
    publishDir path: { dndscvPublishDir(meta) }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sig_maf)
    val(merge_by_patient)

    output:
    tuple val(meta), path("${meta.analysis_type}_dndscv.in"), emit: mut_table

    script:
    """
    Rscript /opt/repo/maf2dndscv.R \\
        --infile ${sig_maf} \\
        --outfile ${meta.analysis_type}_dndscv.in \\
        --by_patient ${merge_by_patient}
    """

    stub:
    """
    touch ${meta.analysis_type}_dndscv.in
    """
}


process DNDSCV_RUN {
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/dndscv:0.5.0"
    publishDir path: { dndscvPublishDir(meta) }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(mut_table), path(refdb), path(covariates)

    output:
    tuple val(meta), path("dndscv_genes_${meta.analysis_type}.tsv"), emit: genes
    tuple val(meta), path("dndscv.out"),                             emit: stdout_log

    script:
    def cov_arg = covariates.name == 'NO_COV' ? '' : "--cov ${covariates}"
    """
    Rscript /opt/repo/dndscv_grch38.R \\
        --mut_file ${mut_table} \\
        --suffix ${meta.analysis_type} \\
        --refdb ${refdb} \\
        ${cov_arg} > dndscv.out
    """

    stub:
    """
    touch dndscv_genes_${meta.analysis_type}.tsv dndscv.out
    """
}
