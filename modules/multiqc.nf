
process MULTIQC {
    container "docker://quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"
    publishDir "${params.outdir}/${params.release_version}/${meta.analysis_type}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(qc_tsvs), path(plot_dirs), path(tmb_tsvs)
    path(multiqc_config)

    output:
    tuple val(meta), path("multiqc_report.html"), emit: report
    tuple val(meta), path("multiqc_data"),        emit: data

    script:
    """
    prepare_multiqc_inputs.py \
        --qc-tsvs ${qc_tsvs} \
        --tmb-tsvs ${tmb_tsvs} \
        --plot-dirs ${plot_dirs} \
        --outdir mqc_inputs

    multiqc mqc_inputs \
        --config ${multiqc_config} \
        --title "Subcohort ${meta.analysis_type} - Somatic Variant QC" \
        --force
    """
}
