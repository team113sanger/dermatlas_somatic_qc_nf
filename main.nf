#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { PROCESS_VCFS } from "./subworkflows/process_vcfs.nf"
include { SUBSET_ANALYSIS } from "./subworkflows/analyse_subset.nf"

workflow {
    
    all_pairs          = Channel.fromPath(params.tumor_normal_pairs, checkIfExists: true)
    unique_pairs       = Channel.fromPath(params.one_per_patient, checkIfExists: true)
    independent_tumors = Channel.fromPath(params.independent, checkIfExists: true)
    patient_md         = Channel.fromPath(params.metadata_manifest, checkIfExists: true)
    dbsnp_vars         = file(params.dbsnp_variants, checkIfExists: true)
    dbsnp_header       = file(params.dbsnp_header, checkIfExists: true)
    baitset            = file(params.baitset, checkIfExists: true)
    
    caveman_vcf_ch = Channel.fromPath(params.caveman_vcfs)
    | map {file -> tuple([sample_id: file.simpleName, caller: "caveman"], file)}
    pindel_vcf_ch = Channel.fromPath(params.pindel_vcfs)
    | map {file -> tuple([sample_id: file.simpleName, caller: "pindel"], file)}

    
    PROCESS_VCFS(caveman_vcf_ch, 
                 pindel_vcf_ch,
                 baitset, 
                 dbsnp_vars, 
                 dbsnp_header)
   
   SUBSET_ANALYSIS(
    PROCESS_VCFS.out.file_list,
    PROCESS_VCFS.out.all_files,
    PROCESS_VCFS.out.sample_list,
    params.genome_build,
    params.filtering_column,
    params.filter_option
   )
    

}




