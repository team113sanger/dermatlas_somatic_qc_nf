include { QC_VARIANTS } from "../modules/filter_variants.nf"
include { CALCULATE_SAMPLE_TMB } from "../modules/calculate_tmb.nf"

workflow COHORT_ANALYSIS{
    take: 
    vcf_ch
    pair_identities
    genome_build
    filter_column
    filter_mode
    analysis_type
    exome_size

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
    | map { meta, file, index -> [meta + [analysis_type: analysis_type], file, index]}
    | set { filtered_vcfs }

    filtered_vcfs.map{ meta, file, index -> [file, index]}.flatten().collect()
    | map {file_list -> tuple([analysis_type: analysis_type], file_list)}
    | set { relevant_vcfs }
    
    filtered_vcfs
    | collectFile(storeDir: "${params.outdir}/${params.release_version}"){
        meta, file, index -> 
        new File("${params.outdir}/${params.release_version}/${meta.analysis_type}").mkdirs()
        def filename = "${meta.analysis_type}/sample_list.tsv"
        [filename, "${meta["sample_id"]}\n"]} 
    | set {sample_list}

    filtered_vcfs
    | map{ meta, files, indexes -> [meta, files.baseName]}
    | collectFile(storeDir: "${params.outdir}/${params.release_version}"){
        meta, basenames ->
        def filename = "${meta.analysis_type}/vcf_list.tsv"
        [filename, basenames + ".gz\n"]} 
    | set { file_list }


    QC_VARIANTS(relevant_vcfs,
                file_list,
                sample_list, 
                genome_build, 
                filter_column, 
                filter_mode)

    keep_ch = QC_VARIANTS.out.keep_maf.transpose()
    CALCULATE_SAMPLE_TMB(keep_ch, exome_size)




}