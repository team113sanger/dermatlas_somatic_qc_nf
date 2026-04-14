
process CONVERT_PLOTS_TO_PNG {
    container "docker://quay.io/biocontainers/poppler:25.07.0"

    input:
    tuple val(meta), path(plot_dirs)

    output:
    tuple val(meta), path("png_plots/*"), emit: plot_dirs

    script:
    def allowed = ["plots_keepPA_vaf_size_filt_matched", "plots_keep_vaf_size_filt_matched"]
    """
    mkdir -p png_plots
    for d in ${plot_dirs}; do
        name=\$(basename \$d)
        case "\$name" in
            ${allowed.join('|')}) ;;
            *) continue ;;
        esac
        out="png_plots/\$name"
        mkdir -p "\$out"
        for pdf in "\$d"/*.pdf; do
            [ -e "\$pdf" ] || continue
            base=\$(basename "\$pdf" .pdf)
            pdftoppm -png -r 150 -singlefile "\$pdf" "\$out/\$base"
        done
    done
    """
}
