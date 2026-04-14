
process MULTIQC {
    container "docker://quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"
    publishDir "${params.outdir}/${params.release_version}", mode: params.publish_dir_mode

    input:
    path(plot_dirs, stageAs: "plots/*")
    path(tmb_files, stageAs: "tmb/*")
    path(multiqc_config)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"),        emit: data

    script:
    """
    prepare_multiqc_inputs.py \
        --plot-dirs plots/* \
        --tmb-files tmb/* \
        --outdir mqc_inputs

    multiqc mqc_inputs \
        --config ${multiqc_config} \
        --config mqc_inputs/dynamic_mqc_config.yaml \
        --title "Somatic Variant QC — Cohort Overview" \
        --force
    """
}
