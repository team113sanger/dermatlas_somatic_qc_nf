#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process FILTER_PASS_VARIANTS {
    container 
    input: 
    tuple val(meta), path(vcf)
    tuple val(meta), path(bedfile)

    def vcfout = "TBC"

    ouput: 
    tuple val(meta), path(vcfout), path("*.tbi")

    script:
    """
    bcftools view -f PASS -O z -o $vcfout -T $bedfile $vcf
    tabix -p vcf $vcfout
    """
    stub: 
    """
    echo stub > test.vcf
    echo stub > test.vcf.tbi
    """

}


process ADD_COMMON_ANNOTATIONS {
    container "gitlab//maf"
    input:
    tuple val(meta), path(vcf), path(index)

    output:
    tuple path("snpflagged.vcf.gz")

    script:
    def anno_vcf = "TBC"
    """
    add_commonSNPs2vcf.sh \
    -p $PROJECTDIR \
    -o $anno_vcf \
    -v $vcf
    """
    stub: 
    """
    echo stub > test.vcf
    echo stub > test.vcf.tbi
    """
}

process QC_VARIANTS {
    input:
    tuple val(meta), path(variants_list)

    output: 
    tuple val(meta), path("*.maf"), emit: maf_file

    script:
    """
    somatic_variants_qc.sh \
    -l variants_list \
    -m caveman_pindel_${f}.maf \
    -s /opt/repo \
    -b GRCh38 \
    -a gnomAD_AF \
    -f filter2
    """
    stub: 
    """
    echo stub > test.maf
    """
}

workflow VCF_TO_ANNO_MAF {
    take:
    caveman_vcfs 
    pindel_vcfs

    
    main:
    caveman_vcfs
    | join(pindel_vcfs)
    | set { raw_vcfs }

    FILTER_PASS_VARIANTS(raw_vcfs)
    ADD_COMMON_ANNOTATIONS(FILTER_PASS_VARIANTS.out)
    QC_VARIANTS(ADD_COMMON_ANNOTATIONS.out)
    
    emit:
    QC_VARIANTS.out
    
}