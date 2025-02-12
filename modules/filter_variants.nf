process FILTER_PASS_VARIANTS {
    publishDir "${meta.vcf_outdir}/${meta.sample_id}", mode: params.publish_dir_mode
    container "quay.io/biocontainers/bcftools:1.9--ha228f0b_4"
    
    input: 
    tuple val(meta), path(vcf)
    path(bedfile)

    output: 
    tuple val(meta), path("*.filt.vcf.gz")
    
    script:
    def vcfout = "$meta.filename"
    """
    bcftools view -f PASS -O z -o ${vcfout}.filt.vcf.gz -T $bedfile $vcf
    """
    
    stub: 
    """
    echo stub > test.vcf
    """

}

process INDEX_PASS_VARIANTS {
    publishDir "${meta.vcf_outdir}/${meta.sample_id}", mode: params.publish_dir_mode
    container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"
    
    input: 
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path(vcf), path("*.tbi")

    script:
    """
    tabix -p vcf ${vcf}
    """
    stub:
    """
    echo stub > test_vcf.gz.tbi
    """

}


process ADD_COMMON_ANNOTATIONS {
    publishDir "${meta.vcf_outdir}/${meta.sample_id}", mode: params.publish_dir_mode
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    
    input:
    tuple val(meta), path(vcf), path(index)
    path(dbsnp_vars)
    path(header)

    output:
    tuple val(meta), path("*snpflagged.vcf.gz"), path("*snpflagged.vcf.gz.tbi")

    script:
    def dbsnp_file = dbsnp_vars[0].name.split(".gz")[0]
    def vcfout = "${meta.filename}" + ".filt"
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
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/qc:0.5.0"
    publishDir "${params.outdir}/${params.release_version}/${meta.analysis_type}", mode: params.publish_dir_mode
    
    input:
    tuple val(meta), path(vcf_files)
    path(file_list)
    path(sample_list)
    val(BUILD)
    val(AF_COL)
    val(filter)
    path(alternative_transcripts)

    output: 
    tuple val(meta), path("pass*.maf"), emit: pass_maf
    tuple val(meta), path("voi*.maf"), emit: voi_maf
    tuple val(meta), path("keep*.maf"), emit: keep_maf
    tuple val(meta), path("plots*"), emit: plot_dirs
    tuple val(meta), path(".tsv"), emit: qc_tsv


    script:
    def f = "${meta.analysis_type}"
    def use_alt = alternative_transcripts.name != "NO_FILE" ? "-t $alternative_transcripts": ''
    """
    /opt/repo/somatic_variants_qc.sh \
    -l $file_list \
    -m caveman_pindel_${f}.maf \
    -s /opt \
    -b $BUILD \
    -a $AF_COL \
    -f $filter \
    $use_alt
    """

    stub: 
    """
    echo stub > test.maf
    echo stub > test.x
    """
}

