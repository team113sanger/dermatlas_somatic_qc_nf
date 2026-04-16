include { MAF_TO_TARGETS; BUILD_SAMPLE_VCF; SIGPROFILER_EXTRACT } from "../modules/sigprofiler.nf"

workflow SIGNATURES {
    take:
    keep_maf_ch        // tuple(meta, keep_maf)        — one per subcohort (meta.analysis_type)
    annotated_vcf_ch   // list of tuple(meta, vcf, tbi) — snpflagged VCFs (meta.sample_id + meta.caller)
    genome_build

    main:

    MAF_TO_TARGETS(keep_maf_ch)

    // Pair each subcohort's target files with its keep MAF (needed for warn-only count check).
    MAF_TO_TARGETS.out.targets
    | join(keep_maf_ch)
    | set { targets_with_maf }

    // Fan out per (subcohort, sample) — key on sample_id for downstream joins.
    targets_with_maf
    | flatMap { meta, sbs_list, id_list, keep_maf ->
        def sbs_files = [sbs_list].flatten()
        def id_files  = [id_list].flatten()
        sbs_files.collect { sbs ->
            def sid = sbs.name.replaceFirst(/^sbs-dbs_targets_/, '').replaceFirst(/\.tsv$/, '')
            def id_file = id_files.find { it.name == "id_targets_${sid}.tsv" }
            assert id_file : "MAF_TO_TARGETS did not emit id_targets file for sample ${sid}"
            tuple(sid, meta + [sample_id: sid], sbs, id_file, keep_maf)
        }
    }
    | set { targets_per_sample }

    // Split per-sample annotated VCFs by caller, rekey on sample_id.
    annotated_vcf_ch
    | flatMap { it }
    | branch { meta, vcf, tbi ->
        caveman: meta.caller == 'caveman'
            return tuple(meta.sample_id, vcf, tbi)
        pindel: meta.caller == 'pindel'
            return tuple(meta.sample_id, vcf, tbi)
    }
    | set { by_caller }

    // Combine-by-key: every (subcohort × sample) gets its matching caveman + pindel VCF.
    // `combine(by: 0)` is the right tool — `join` would drop duplicates when a sample appears
    // in multiple subcohorts (all, onePerPatient, independent).
    targets_per_sample
    | combine(by_caller.caveman, by: 0)
    | combine(by_caller.pindel,  by: 0)
    | map { _sid, meta, sbs, id_file, keep_maf, cv, cv_tbi, pv, pv_tbi ->
        tuple(meta, sbs, id_file, cv, cv_tbi, pv, pv_tbi, keep_maf)
    }
    | set { build_inputs }

    BUILD_SAMPLE_VCF(build_inputs)

    // Regroup per-sample VCFs back into a per-subcohort bundle for SigProfilerExtractor.
    // Keying on a fresh `[analysis_type: ...]` map drops sample_id so groupTuple collapses cleanly.
    BUILD_SAMPLE_VCF.out
    | map { meta, vcf -> tuple([analysis_type: meta.analysis_type], vcf) }
    | groupTuple()
    | set { subcohort_vcfs }

    SIGPROFILER_EXTRACT(subcohort_vcfs, genome_build)

    emit:
    signatures = SIGPROFILER_EXTRACT.out.signatures
}
