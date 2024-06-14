process FILTER_PASS_VARIANTS {
    publishDir "results/${meta.caller}/${meta.sample_id}", mode: params.publish_dir_mode
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    
    input: 
    tuple val(meta), path(vcf)
    path(bedfile)

    output: 
    tuple val(meta), path("*.filt.vcf.gz"), path("*.tbi")
    
    script:
    def vcfout = "$meta.sample_id"
    """
    bcftools view -f PASS -O z -o ${vcfout}.filt.vcf.gz -T $bedfile $vcf
    tabix -p vcf ${vcfout}.filt.vcf.gz
    """
    
    stub: 
    """
    echo stub > test.vcf
    echo stub > test.vcf.tbi
    """

}


process ADD_COMMON_ANNOTATIONS {
    publishDir "results/${meta.caller}/${meta.sample_id}", mode: params.publish_dir_mode
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"    
    
    input:
    tuple val(meta), path(vcf), path(index)
    path(dbsnp_vars)
    path(header)

    output:
    tuple path("*snpflagged.vcf.gz"), path("*snpflagged.vcf.gz.tbi")

    script:
    def dbsnp_file = dbsnp_vars[0].name.split(".gz")[0]
    def vcfout = "$meta.sample_id" + "_" + "$meta.caller"
    """
    bcftools annotate \
    -a ${dbsnp_file}.gz \
    -h $header \
    -c CHROM,POS,-,REF,ALT,dbSNP \
    -O z -o ${vcfout}.snpflagged.vcf.gz $vcf
    tabix -p vcf ${vcfout}.snpflagged.vcf.gz
    """
    
    stub: 
    """
    echo stub > test.snpflagged.vcf.gz
    echo stub > test.snpflagged.vcf.gz.tbi
    """
}


process QC_VARIANTS {
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/maf/feature/dockerise:38bd755a"
    publishDir "results", mode: params.publish_dir_mode
    input:
    path(file_list)
    path(vcflist)
    path(sample_list)
    val(BUILD)
    val(AF_COL)
    val(filter)
    val(maf_file)

    output: 
    path("pass*.maf"), emit: pass_maf
    path("voi*.maf"), emit: voi_maf

    script:
    def f = "TBC"
    """
    /opt/repo/qc_somatic_variants.sh \
    -l $file_list \
    -m caveman_pindel_${f}.maf \
    -s /opt/repo \
    -b $BUILD \
    -a $AF_COL \
    -f $filter
    """
    stub: 
    """
    echo stub > test.maf
    echo stub > test.x
    """
}


process CALCULATE_SAMPLE_TMB {
    publishDir "results", mode: params.publish_dir_mode
    input:
    path(maf_file)

    output:
    path("mutations_per_Mb.tsv"), emit: tmb

    shell:
    f = "TBC"
    """
    sample_ids=`cut -f 11 !{maf_file} | grep PD |sort -u`
    for sample_id in "$sample_ids"; echo -ne "$sample_id\t"; muts=`grep "$sample_id" !{maf_file} | cut -f 4,5 | sort -u | wc -l`; echo "$muts"/48.225157 | bc -l; done > mutations_per_Mb.tsv
    """

}