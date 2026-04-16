process MAF_TO_TARGETS {
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/qc:0.5.1"
    publishDir "${params.outdir}/${params.release_version}/${meta.analysis_type}/sigprofiler", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(keep_maf)

    output:
    tuple val(meta), path("REGIONS/sbs-dbs_targets_*.tsv"), path("REGIONS/id_targets_*.tsv"), emit: targets

    script:
    """
    mkdir -p REGIONS
    cd REGIONS
    maf2targets.R ../${keep_maf}
    """

    stub:
    """
    mkdir -p REGIONS
    touch REGIONS/sbs-dbs_targets_STUB.tsv REGIONS/id_targets_STUB.tsv
    """
}


process BUILD_SAMPLE_VCF {
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"

    input:
    tuple val(meta),
          path(sbs_targets),
          path(id_targets),
          path(caveman_vcf), path(caveman_tbi),
          path(pindel_vcf),  path(pindel_tbi),
          path(keep_maf)

    output:
    tuple val(meta), path("${meta.sample_id}.vcf")

    script:
    """
    build_sample_vcf.sh \\
        -s ${meta.sample_id} \\
        -c ${caveman_vcf} \\
        -p ${pindel_vcf} \\
        -x ${sbs_targets} \\
        -i ${id_targets} \\
        -m ${keep_maf}
    """

    stub:
    """
    echo stub > ${meta.sample_id}.vcf
    """
}


// No container directive: this process is expected to run on the host with the
// `sigprofiler/1.1.21-virtual-environment` module loaded (see farm22 profile in
// nextflow.config). For non-farm22 environments, add a container or ensure
// SigProfilerExtractor + the GRCh38 reference bundle are on PATH.
process SIGPROFILER_EXTRACT {
    publishDir "${params.outdir}/${params.release_version}/${meta.analysis_type}/sigprofiler", mode: params.publish_dir_mode

    cpus 12
    memory { 12.GB * task.attempt }
    time { 48.h * task.attempt }

    input:
    tuple val(meta), path(vcfs, stageAs: "VCFS/*")
    val(genome_build)

    output:
    tuple val(meta), path("results"), emit: signatures

    script:
    def seed_arg = params.sigprofiler_seed ?: ''
    """
    run_sigprofilerextractor.py VCFS results ${genome_build} ${seed_arg}
    """

    stub:
    """
    mkdir -p results/SBS96/Suggested_Solution
    touch results/stub.txt
    """
}
