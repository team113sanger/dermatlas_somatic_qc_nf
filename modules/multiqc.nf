
process MULTIQC {
    container "docker://quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"
    publishDir "${params.outdir}/${params.release_version}", mode: params.publish_dir_mode

    input:
    path(plot_dirs, stageAs: "plots/*")
    path(multiqc_config)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"),        emit: data

    script:
    """
    prepare_multiqc_inputs.py \
        --plot-dirs plots/* \
        --outdir mqc_inputs

    multiqc mqc_inputs \
        --config ${multiqc_config} \
        --title "Somatic Variant QC — Cohort Overview" \
        --force
    """
}
