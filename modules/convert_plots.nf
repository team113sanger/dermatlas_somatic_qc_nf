
process CONVERT_PLOTS_TO_PNG {
    container "docker://quay.io/biocontainers/poppler:22.12.0--h091648b_0"

    input:
    tuple val(meta), path(plot_dirs)

    output:
    tuple val(meta), path("png_plots/*"), emit: plot_dirs

    script:
    """
    mkdir -p png_plots
    for d in ${plot_dirs}; do
        out="png_plots/\$(basename \$d)"
        mkdir -p "\$out"
        for pdf in "\$d"/*.pdf; do
            [ -e "\$pdf" ] || continue
            base=\$(basename "\$pdf" .pdf)
            pdftoppm -png -r 150 -singlefile "\$pdf" "\$out/\$base"
        done
    done
    """
}
