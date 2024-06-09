process FILTER_PASS_VARIANTS {
    container "biocontainers/bcftools:v1.9-1-deb_cv1"
    
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