include { MAF_TO_DNDSCV_INPUT; DNDSCV_RUN } from "../modules/dndscv.nf"

workflow DNDSCV {
    take:
    sig_maf_ch          // tuple(meta, sig_maf) — one per subcohort
    refdb               // path
    covariates          // path, or sentinel file named NO_COV when not provided
    merge_by_patient    // 'y' | 'n'

    main:

    MAF_TO_DNDSCV_INPUT(sig_maf_ch, merge_by_patient)

    // Fan out across covariate modes. When a real covariates file is supplied
    // we run dndscv twice per subcohort (with & without covariates) to match
    // the legacy manual analysis layout. When the sentinel NO_COV is used we
    // only run the without_covariates mode.
    def has_cov = covariates.name != 'NO_COV'

    MAF_TO_DNDSCV_INPUT.out.mut_table
    | flatMap { meta, mut ->
        def rows = [ tuple(meta + [covariate_mode: 'without_covariates'], mut, refdb, file("${projectDir}/assets/NO_COV")) ]
        if (has_cov) {
            rows << tuple(meta + [covariate_mode: 'with_covariates'], mut, refdb, covariates)
        }
        rows
    }
    | set { run_inputs }

    DNDSCV_RUN(run_inputs)

    emit:
    genes       = DNDSCV_RUN.out.genes
    stdout_log  = DNDSCV_RUN.out.stdout_log
}
