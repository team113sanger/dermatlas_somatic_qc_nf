
process CONVERT_PLOTS_TO_PNG {
    container "docker://quay.io/biocontainers/poppler:25.07.0"

    input:
    tuple val(meta), path(plot_dirs)

    output:
    tuple val(meta), path("png_plots/*"), emit: plot_dirs

    script:
    def allowed = ["plots_keepPA_vaf_size_filt_matched", "plots_keep_vaf_size_filt_matched"]
    def sub = meta.analysis_type
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export FONTCONFIG_PATH="\$PWD/.cache/fontconfig"
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"

    mkdir -p png_plots
    for d in ${plot_dirs}; do
        name=\$(basename \$d)
        case "\$name" in
            ${allowed.join('|')}) ;;
            *) continue ;;
        esac
        out="png_plots/${sub}__\$name"
        mkdir -p "\$out"
        for pdf in "\$d"/*.pdf; do
            [ -e "\$pdf" ] || continue
            [ -s "\$pdf" ] || { echo "skip empty PDF: \$pdf" >&2; continue; }
            pages=\$(pdfinfo "\$pdf" 2>/dev/null | awk '/^Pages:/ {print \$2}')
            if [ -z "\$pages" ] || [ "\$pages" -lt 1 ]; then
                echo "skip 0-page PDF: \$pdf" >&2
                continue
            fi
            base=\$(basename "\$pdf" .pdf)
            pdftoppm -png -r 150 -singlefile "\$pdf" "\$out/\$base"
        done
    done
    """
}
