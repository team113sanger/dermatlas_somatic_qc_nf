#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VCF_TO_ANNO_MAF } from "../subworkflows/vcf_to_maf.nf"
workflow {
    vcf_channel = Channel.fromPath(params.caveman_vcfs)
    VCF_TO_ANNO_MAF()

}


params.caveman_vcfs = 

