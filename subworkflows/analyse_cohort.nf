include { QC_VARIANTS } from "../modules/filter_variants.nf"
include { CALCULATE_SAMPLE_TMB } from "../modules/calculate_tmb.nf"
include { MAF_TO_EXCEL } from "../modules/maf_to_excel.nf"

workflow SUBCOHORT_ANALYSIS {
    take:
    subcohort_sample_sets  // channel of (subcohort_name, sample_list_file) tuples
    vcf_ch
    genome_build
    filter_column
    filter_mode
    exome_size
    alternative_transcripts

    main:

    // Build a set of (sample_id, subcohort_name) for filtering
    subcohort_sample_sets
    | flatMap { subcohort_name, sample_list ->
        sample_list.splitCsv(sep: "\t", header: ['tumor', 'normal'])
            .collect { row ->
                tuple(row.tumor, subcohort_name)
            }
    }
    | set { sample_to_subcohort }

    // Combine VCFs with subcohort membership, filter matches
    vcf_ch
    | flatMap { n -> n }
    | map { meta, file, index -> [meta.sample_id, meta, file, index] }
    | combine(sample_to_subcohort)
    | filter { sample_id, meta, file, index, tumor_id, subcohort_name ->
        sample_id == tumor_id
    }
    | map { sample_id, meta, file, index, _tumor_id, subcohort_name ->
        [meta + [analysis_type: subcohort_name], file, index]
    }
    | set { filtered_vcfs }
    filtered_vcfs.view()

    // Group VCFs by subcohort for QC_VARIANTS
    filtered_vcfs
    | map { meta, file, index -> [meta.analysis_type, file, index] }
    | groupTuple()
    | map { analysis_type, files, indexes ->
        tuple(analysis_type, [analysis_type: analysis_type], files.flatten() + indexes.flatten())
    }
    | set { relevant_vcfs }

    // Create sample list per subcohort
    filtered_vcfs
    | collectFile(storeDir: "${params.outdir}/${params.release_version}") {
        meta, file, index ->
        new File("${params.outdir}/${params.release_version}/${meta.analysis_type}").mkdirs()
        def filename = "${meta.analysis_type}/sample_list.tsv"
        [filename, "${meta["sample_id"]}\n"]
    }
    | map { sample_list_file ->
        def analysis_type = sample_list_file.parent.name
        tuple(analysis_type, sample_list_file)
    }
    | set { sample_list }

    // Create VCF list per subcohort
    filtered_vcfs
    | map { meta, files, indexes -> [meta, files.baseName] }
    | collectFile(storeDir: "${params.outdir}/${params.release_version}") {
        meta, basenames ->
        def filename = "${meta.analysis_type}/vcf_list.tsv"
        [filename, basenames + ".gz\n"]
    }
    | map { vcf_list_file ->
        def analysis_type = vcf_list_file.parent.name
        tuple(analysis_type, vcf_list_file)
    }
    | set { file_list }

    // Join all channels by analysis_type and restructure for QC_VARIANTS
    relevant_vcfs
    | join(file_list)
    | join(sample_list)
    | multiMap { analysis_type, meta, vcf_files, vcf_list_file, sample_list_file ->
        vcfs: tuple(meta, vcf_files)
        file_list: vcf_list_file
        sample_list: sample_list_file
    }
    | set { qc_channels }

    QC_VARIANTS(qc_channels.vcfs,
                qc_channels.file_list,
                qc_channels.sample_list,
                genome_build,
                filter_column,
                filter_mode,
                alternative_transcripts)

    keep_ch = QC_VARIANTS.out.keep_maf.transpose()
    CALCULATE_SAMPLE_TMB(keep_ch, exome_size)
    MAF_TO_EXCEL(keep_ch)

    emit:
    output_variants = QC_VARIANTS.out.keep_maf

}