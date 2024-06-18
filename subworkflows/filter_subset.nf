workflow FILTER_SAMPLES{
    take: 
    vcf_ch
    pair_identities

    main:

    pair_identities 
    | splitCsv(sep:"\t", header:['tumor', 'normal']) 
    | map{ meta -> 
        [sample_id: meta.tumor]}
    | set {sample_filter}

    
    vcf_ch
    | flatMap{n -> n}
    | map{ meta, file, index -> [meta.subMap('sample_id'), meta, file, index]}
    | groupTuple()
    | join(sample_filter)
    | transpose()
    | map{ sample_id, meta, file, index -> [meta, file, index] }
    | set { filtered_vcfs }

    // filtered_vcfs 
    // | map{ meta, file, index -> [file, index]}
    // | flatten()
    // | collect()
    // | view()


    filtered_vcfs
    | collectFile(name: "sample_list.tsv", storeDir: 'results'){
        meta, file, index -> 
        ["sample_list.tsv", "${meta["sample_id"]}\n"]} 
    | set {sample_list}

    filtered_vcfs
    | map{ meta, files, indexes -> files.baseName}
    | collectFile(name: "vcf_list.tsv", storeDir: "results"){
        basenames ->
        ["vcfs_list.tsv", basenames + ".gz\n"]} 
    | set { file_list }


    // QC_VARIANTS(file_list,
    //             filtered_vcfs,
    //             sample_list, 
    //             genome_build, 
    //             filter_column, 
    //             filter_mode,
    //             "test_maf")
    // CALCULATE_SAMPLE_TMB(QC_VARIANTS.out.pass_maf)




}