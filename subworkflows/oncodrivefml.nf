include { MAF_TO_ONCODRIVEFML_INPUT; ONCODRIVEFML_RUN } from "../modules/oncodrivefml.nf"

workflow ONCODRIVEFML {
    take:
    sig_maf_ch              // tuple(meta, sig_maf) — one per subcohort
    exclude_file            // path, or sentinel file named NO_EXCLUDE
    cadd_tabix              // path
    cadd_tabix_tbi          // path
    elements_bed            // path
    analysis_conf_template  // path
    bgdata_conf_template    // path
    bgdata_repo             // val (string; bind-mounted, not staged)
    bgdata_cache            // val (string; bind-mounted, not staged)
    seed                    // val

    main:

    MAF_TO_ONCODRIVEFML_INPUT(sig_maf_ch, exclude_file)

    ONCODRIVEFML_RUN(
        MAF_TO_ONCODRIVEFML_INPUT.out.variants,
        elements_bed,
        cadd_tabix,
        cadd_tabix_tbi,
        analysis_conf_template,
        bgdata_conf_template,
        bgdata_repo,
        bgdata_cache,
        seed
    )

    emit:
    results     = ONCODRIVEFML_RUN.out.results
    stdout_log  = ONCODRIVEFML_RUN.out.stdout_log
}
