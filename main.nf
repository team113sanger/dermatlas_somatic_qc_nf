#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { VCF_TO_ANNO_MAF } from "../subworkflows/vcf_to_maf.nf"
workflow {
    caveman_vcf_ch = Channel.fromPath(params.caveman_vcfs)
    pindel_vcfs_vcf_ch = Channel.fromPath(params.pindel_vcfs)
    VCF_TO_ANNO_MAF(caveman_vcf_ch, 
                    pindel_vcfs_vcf_ch)

}




