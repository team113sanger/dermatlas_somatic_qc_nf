#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { PROCESS_VCFS } from "./subworkflows/process_vcfs.nf"
include { COHORT_ANALYSIS as ONE_TUMOR_PER_PATIENT} from "./subworkflows/analyse_cohort.nf"
include { COHORT_ANALYSIS as INDEPENDENT_TUMORS} from "./subworkflows/analyse_cohort.nf"
include { COHORT_ANALYSIS as ALL_TUMORS} from "./subworkflows/analyse_cohort.nf"

workflow DERMATLAS_SOMATIC_VARIANT_QC {

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

    if (params.all_samples){
    all_pairs = Channel.fromPath(params.all_samples, checkIfExists: true)
    
    ALL_TUMORS(PROCESS_VCFS.out.all_files, 
            all_pairs,
            params.genome_build,
            params.filtering_column,
            params.filter_option,
            "all",
            params.exome_size,
            params.alternative_transcripts)
    }

    if( params.one_per_patient) {
    unique_pairs = Channel.fromPath(params.one_per_patient, checkIfExists: true)
    
    ONE_TUMOR_PER_PATIENT(PROCESS_VCFS.out.all_files, 
                unique_pairs,
                params.genome_build,
                params.filtering_column,
                params.filter_option,
                "onePerPatient",
                params.exome_size,
                params.alternative_transcripts)
    }

    if( params.independent) {
    independent_tumors = Channel.fromPath(params.independent, checkIfExists: true)
    
    INDEPENDENT_TUMORS(PROCESS_VCFS.out.all_files, 
                        independent_tumors,
                        params.genome_build,
                        params.filtering_column,
                        params.filter_option,
                        "independent",
                        params.exome_size,
                        params.alternative_transcripts)
    }
    emit: 
    all_ch = ALL_TUMORS.out.output_variants
}

workflow {
    DERMATLAS_SOMATIC_VARIANT_QC()
}


