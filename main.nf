#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { PROCESS_VCFS } from "./subworkflows/process_vcfs.nf"
include { SUBSET_ANALYSIS } from "./subworkflows/analyse_subset.nf"
include { DERMATLAS_METADATA } from "./subworkflows/process_metadata.nf"

workflow {
    
    all_pairs          = Channel.fromPath(params.tumor_normal_pairs, checkIfExists: true)
    unique_pairs       = Channel.fromPath(params.one_per_patient, checkIfExists: true)
    independent_tumors = Channel.fromPath(params.independent, checkIfExists: true)
    patient_md         = Channel.fromPath(params.metadata_manifest, checkIfExists: true)
    dbsnp_vars         = file(params.dbsnp_variants, checkIfExists: true)
    dbsnp_header       = file(params.dbsnp_header, checkIfExists: true)
    baitset            = file(params.baitset, checkIfExists: true)
    metadata           = file(params.metadata_manifest, checkIfExists: true)
    all_pairs.view()
    caveman_vcf_ch = Channel.fromPath(params.caveman_vcfs)
    | map {file -> tuple([sample_id: file.simpleName, caller: "caveman"], file)}
    pindel_vcf_ch = Channel.fromPath(params.pindel_vcfs)
    | map {file -> tuple([sample_id: file.simpleName, caller: "pindel"], file)}

    // DERMATLAS_METADATA(independent_tumors, metadata)
    
    PROCESS_VCFS(caveman_vcf_ch, 
                 pindel_vcf_ch,
                 baitset, 
                 dbsnp_vars, 
                 dbsnp_header)
    // DIVIDE_COHORT(PROCESS_VCFS.out.file_list,
    //               PROCESS_VCFS.out.all_files,
    //               PROCESS_VCFS.out.sample_list,
    //               all_pairs
    //               unique_pairs
    //               independent_tumors
    //               )
   
   SUBSET_ANALYSIS(
    PROCESS_VCFS.out.file_list,
    PROCESS_VCFS.out.all_files,
    PROCESS_VCFS.out.sample_list,
    params.genome_build,
    params.filtering_column,
    params.filter_option
   )
    

}




