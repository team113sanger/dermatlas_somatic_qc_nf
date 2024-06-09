include { FILTER_PASS_VARIANTS; ADD_COMMON_ANNOTATIONS; QC_VARIANTS } from "../modules/filter_variants.nf"

workflow VCF_TO_ANNO_MAF {
    take:
    caveman_vcfs 
    pindel_vcfs

    
    main:
    caveman_vcfs
    | join(pindel_vcfs)
    | set { raw_vcfs }

    FILTER_PASS_VARIANTS(raw_vcfs)
    ADD_COMMON_ANNOTATIONS(FILTER_PASS_VARIANTS.out)
    QC_VARIANTS(ADD_COMMON_ANNOTATIONS.out)
    
    emit:
    QC_VARIANTS.out
    
}