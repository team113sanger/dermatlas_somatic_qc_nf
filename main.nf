#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { PROCESS_VCFS } from "./subworkflows/process_vcfs.nf"
include { SUBCOHORT_ANALYSIS } from "./subworkflows/analyse_cohort.nf"

workflow DERMATLAS_SOMATIC_VARIANT_QC {

    dbsnp_vars   = file(params.dbsnp_variants, checkIfExists: true)
    dbsnp_header = file(params.dbsnp_header, checkIfExists: true)
    baitset      = file(params.baitset, checkIfExists: true)
    metadata     = Channel.fromPath(params.metadata_manifest, checkIfExists: true) // Unused - for future extensions or deprecation


    caveman_vcf_ch = Channel.fromPath(params.caveman_vcfs)
    | map {file -> tuple([sample_id: file.simpleName,
                          caller: "caveman",
                          filename: file.name.split(".vcf")[0],
                          vcf_outdir: params.caveman_outdir], file)}

    pindel_vcf_ch = Channel.fromPath(params.pindel_vcfs)
    | map {file -> tuple([sample_id: file.simpleName,
                          caller: "pindel",
                          filename: file.name.split(".vcf")[0],
                          vcf_outdir: params.pindel_outdir], file)}

    PROCESS_VCFS(caveman_vcf_ch,
                 pindel_vcf_ch,
                 baitset,
                 dbsnp_vars,
                 dbsnp_header)

    if (params.subcohorts) {
        log.info("Running QC analysis for subcohorts: ${params.subcohorts.keySet().join(', ')}...")

        // Create channel of (subcohort_name, sample_list_file) tuples from params.subcohorts map
        subcohort_sample_sets = Channel.fromList(
            params.subcohorts.collect { subcohort, config ->
                tuple(subcohort, file(config.sample_list, checkIfExists: true))
            }
        )

        SUBCOHORT_ANALYSIS(
            subcohort_sample_sets,
            PROCESS_VCFS.out.all_files,
            params.genome_build,
            params.filtering_column,
            params.filter_option,
            params.exome_size,
            params.alternative_transcripts
        )
    }

    emit:
    output_ch = SUBCOHORT_ANALYSIS.out.output_variants
}

workflow {
    DERMATLAS_SOMATIC_VARIANT_QC()
}


