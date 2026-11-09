// Resolve a subcohort's OncodriveFML publish-dir name: user-provided legacy
// map wins, otherwise fall back to the raw analysis_type key.
def oncodrivefmlPublishDir(meta) {
    def name = params.oncodrivefml_subcohort_names?.get(meta.analysis_type) ?: meta.analysis_type
    "${params.oncodrivefml_outdir}/${params.release_version}/${name}"
}


process MAF_TO_ONCODRIVEFML_INPUT {
    container "docker.io/rocker/tidyverse:4.4.0"
    publishDir path: { oncodrivefmlPublishDir(meta) }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sig_maf)
    path(exclude_file)

    output:
    tuple val(meta), path("${meta.cohort_prefix}_oncodrivefml.tsv"), emit: variants

    script:
    def exclude_arg = exclude_file.name == 'NO_EXCLUDE' ? '' : "${exclude_file}"
    """
    maf2oncodrivefml.R ${sig_maf} ${meta.cohort_prefix}_oncodrivefml.tsv ${exclude_arg}
    """

    stub:
    """
    touch ${meta.cohort_prefix}_oncodrivefml.tsv
    """
}


process ONCODRIVEFML_RUN {
    container "gitlab-registry.internal.sanger.ac.uk/team113sanger/team113_software/oncodrivefml:2.4.0"
    publishDir path: { oncodrivefmlPublishDir(meta) }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(variants_tsv)
    path(elements_bed)
    path(cadd_tabix)
    path(cadd_tabix_tbi)
    path(analysis_conf_template)
    path(bgdata_conf_template)
    val(bgdata_repo)
    val(bgdata_cache)
    val(seed)

    output:
    tuple val(meta), path("*-oncodrivefml.tsv.gz"),               emit: results
    tuple val(meta), path("*-oncodrivefml.png"),  optional: true, emit: plot
    tuple val(meta), path("*-oncodrivefml.html"), optional: true, emit: report
    tuple val(meta), path("oncodrivefml.out"),                    emit: stdout_log

    script:
    """
    sed -e "s|__CADD_TABIX__|${cadd_tabix}|" ${analysis_conf_template} > oncodrive_v2.conf

    mkdir -p bbglab
    sed -e "s|__BGDATA_REPO__|${bgdata_repo}|" \\
        -e "s|__BGDATA_CACHE__|${bgdata_cache}|" \\
        ${bgdata_conf_template} > bbglab/bgdatav2.conf

    export XDG_CONFIG_HOME=\$PWD

    oncodrivefml \\
        -i ${variants_tsv} \\
        -e ${elements_bed} \\
        -c oncodrive_v2.conf \\
        -s wes \\
        -o . \\
        --seed ${seed} \\
        --debug > oncodrivefml.out 2>&1
    """

    stub:
    """
    touch ${meta.cohort_prefix}-oncodrivefml.tsv.gz
    touch ${meta.cohort_prefix}-oncodrivefml.png
    touch ${meta.cohort_prefix}-oncodrivefml.html
    touch oncodrivefml.out
    """
}
