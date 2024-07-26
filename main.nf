#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { PROCESS_VCFS } from "./subworkflows/process_vcfs.nf"
include { COHORT_ANALYSIS as ANALYSE_OTPP} from "./subworkflows/analyse_cohort.nf"
include { COHORT_ANALYSIS as ANALYSE_INDEPENDENT} from "./subworkflows/analyse_cohort.nf"
include { COHORT_ANALYSIS as ANALYSE_ALL} from "./subworkflows/analyse_cohort.nf"

workflow {
    
    all_pairs          = Channel.fromPath(params.tumor_normal_pairs, checkIfExists: true)
    unique_pairs       = Channel.fromPath(params.one_per_patient, checkIfExists: true)
    independent_tumors = Channel.fromPath(params.independent, checkIfExists: true)
    patient_md         = Channel.fromPath(params.metadata_manifest, checkIfExists: true)
    dbsnp_vars         = file(params.dbsnp_variants, checkIfExists: true)
    dbsnp_header       = file(params.dbsnp_header, checkIfExists: true)
    baitset            = file(params.baitset, checkIfExists: true)
    metadata           = Channel.fromPath(params.metadata_manifest, checkIfExists: true)
    
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

    ANALYSE_ALL(PROCESS_VCFS.out.all_files, 
            all_pairs,
            params.genome_build,
            params.filtering_column,
            params.filter_option,
            "all")
    
    // if (!isEmpty(unique_pairs)){
    ANALYSE_OTPP(PROCESS_VCFS.out.all_files, 
                unique_pairs,
                params.genome_build,
                params.filtering_column,
                params.filter_option,
                "one_tumor_per_patient")
    // }

    // if (!isEmpty(independent_tumors)){
    ANALYSE_INDEPENDENT(PROCESS_VCFS.out.all_files, 
                        independent_tumors,
                        params.genome_build,
                        params.filtering_column,
                        params.filter_option,
                        "independent")
                        // }

}




