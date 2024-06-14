include { FILTER_PASS_VARIANTS; ADD_COMMON_ANNOTATIONS; QC_VARIANTS; CALCULATE_SAMPLE_TMB } from "../modules/filter_variants.nf"

workflow VCF_TO_ANNO_MAF {
    take:
    caveman_vcfs 
    pindel_vcfs
    baitset
    dbsnp_vars
    dbsnp_header
    genome_build
    filter_column
    filter_mode

    main:
    caveman_vcfs
    | concat( pindel_vcfs )
    | set { raw_vcfs }

    raw_vcfs
    | map{meta, file -> meta}
    | collectFile(name: "sample_list.tsv"){
        meta -> 
        ["sample_list.tsv", "${meta["sample_id"]}\n"]} 
    | set {sample_list}

    FILTER_PASS_VARIANTS(raw_vcfs, baitset)
    ADD_COMMON_ANNOTATIONS(FILTER_PASS_VARIANTS.out, dbsnp_vars, dbsnp_header)
    ADD_COMMON_ANNOTATIONS.out.map{ files, indexes -> files}.collect()
    | set { annotated_files }

    ADD_COMMON_ANNOTATIONS.out.map{ files, indexes -> files.baseName}
    | collectFile(name: "vcf_list.tsv", storeDir: "results"){
        basenames ->
        ["vcfs_list.tsv", basenames + ".gz\n"]} 
    | set { file_list }

    ADD_COMMON_ANNOTATIONS.out.map{ files, indexes -> indexes}.collect()
    | set {indices}


    annotated_files.merge(indices) 
    | set {all_files}

    QC_VARIANTS(file_list,
                all_files,
                sample_list, 
                genome_build, 
                filter_column, 
                filter_mode,
                "test_maf")
    CALCULATE_SAMPLE_TMB(QC_VARIANTS.out.pass_maf)
    
    
    // emit:
    // QC_VARIANTS.out.pass_maf
    // QC_VARIANTS.out.voi_maf
    
}