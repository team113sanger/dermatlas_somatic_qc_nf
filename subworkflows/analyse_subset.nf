include { QC_VARIANTS; CALCULATE_SAMPLE_TMB } from "../modules/filter_variants.nf"

workflow SUBSET_ANALYSIS {
    take: 
    file_list
    all_files
    sample_list
    genome_build
    filter_column
    filter_mode

    main:
     QC_VARIANTS(file_list,
                 all_files,
                 sample_list, 
                 genome_build, 
                 filter_column, 
                 filter_mode,
                 "test_maf")
    CALCULATE_SAMPLE_TMB(QC_VARIANTS.out.pass_maf)

    emit: 
    CALCULATE_SAMPLE_TMB.out

}