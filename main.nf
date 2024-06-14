#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { VCF_TO_ANNO_MAF } from "./subworkflows/vcf_to_maf.nf"

workflow {
    caveman_vcf_ch = Channel.fromPath(params.caveman_vcfs)
    | map {file -> tuple([sample_id: file.simpleName, caller: "caveman"], file)}
    pindel_vcf_ch = Channel.fromPath(params.pindel_vcfs)
    | map {file -> tuple([sample_id: file.simpleName, caller: "pindel"], file)}

    dbsnp_vars = file(params.dbsnp_variants, checkIfExists: true)
    dbsnp_header = file(params.dbsnp_header, checkIfExists: true)
    baitset = file(params.baitset, checkIfExists: true)
    
    VCF_TO_ANNO_MAF(caveman_vcf_ch, 
                    pindel_vcf_ch,
                    baitset, 
                    dbsnp_vars, 
                    dbsnp_header,
                    params.genome_build,
                    params.filtering_column,
                    params.filter_option)

}




