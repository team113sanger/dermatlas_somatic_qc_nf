include { FILTER_PASS_VARIANTS; ADD_COMMON_ANNOTATIONS } from "../modules/filter_variants.nf"

workflow PROCESS_VCFS {
    take:
    caveman_vcfs 
    pindel_vcfs
    baitset
    dbsnp_vars
    dbsnp_header

    main:
    caveman_vcfs
    | concat( pindel_vcfs )
    | set { raw_vcfs }

    FILTER_PASS_VARIANTS(raw_vcfs, baitset)
    ADD_COMMON_ANNOTATIONS(FILTER_PASS_VARIANTS.out, 
                           dbsnp_vars, 
                           dbsnp_header)
    ADD_COMMON_ANNOTATIONS.out.collect(flat: false)
    | set { all_files }
    
    emit:
    all_files
}