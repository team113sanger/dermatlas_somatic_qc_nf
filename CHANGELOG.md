# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Keywords

As of the release following 1.2.1 the following *keywords* are used at the start of each
changelog entry to indicate the impact of the change:

- **REPRODUCIBILITY** - a change to the pipeline's scientific processing that
  may cause the same input data to produce different scientific outputs or
  results, including changes to algorithms, tolerances, randomisation,
  scientific functionality, or output formats.
- **ROBUSTNESS** - a fix or improvement to the pipeline's scientific
  functionality that improves correctness, reliability, or the range of inputs
  that can be processed, without intentionally changing the scientific results
  of an equivalent successful analysis.
- **INTEGRATION** - a change to how the pipeline integrates with other systems
  or infrastructure, without changing its scientific processing or results.

## [Unreleased]
### Added
- **INTEGRATION** - run reporting. `lib/Utils.groovy` (shared verbatim with the other
  Dermatlas pipelines) is wired in by `workflow.onComplete { Utils.reportRun(workflow, params) }`
  and records each run in the Dermatlas website's analysis log (via `dermatlas-http`, >= 0.6.1)
  and/or posts a Slack message. Both are explicit opt-ins, gated by the
  `DERMATLAS_WEBSITE_LOGGING` / `DERMATLAS_SLACK_NOTIFICATIONS` environment toggles, and
  never fire on stub runs. `nextflow.config` gains `is_stub`,
  `analysis_pipeline_slug = 'somatic_pipe'` and `trace_file`.
- **INTEGRATION** - `nextflow.config` gains a `trace {}` block. The execution trace and
  report are named `execution_trace-<RUN_ID>.txt` / `execution_report-<RUN_ID>.html`
  under the launcher's `TRACE_DIR`, from the `RUN_ID` the launcher exports (a bare
  timestamp for a direct `nextflow run`).
- **INTEGRATION** - `assets/run_somatic_variants.sh` is rebuilt from the reference Dermatlas
  launcher (as in `dermatlas_copy_number_nf`): it sources the project `source_me.sh`
  (`SOURCE_ME`, `"none"` to skip), validates the environment before launch, reports a
  failed launch to stderr and (opt-in) Slack, holds an exclusive `flock` on
  `${PROJECT_DIR}/somatic_pipe/.lock` (a concurrent submission exits 75), writes a
  `.completed_successfully` / `.completed_with_error` sentinel, one log per nextflow
  command (`logs/nextflow-{pull,run}-<RUN_ID>.log`), per-revision `NXF_ASSETS` clones, a
  pinned `NXF_SINGULARITY_CACHEDIR`, and on success writes
  `stats/resource-stats-<RUN_ID>.txt`, reports the work-dir usage to the website
  (`dermatlas-http cohort analysis-workdir-stats`, >= 0.6.2, module-loaded via
  `DERMATLAS_HTTP_MODULE`) and deletes the work directory (`DERMATLAS_CLEANUP_WORK_DIR`).
  See "Without the website", "Toggles" and "Reclaiming disk space" in the README.
- **INTEGRATION** - `.update-version.sh` sets the release version in every file that
  records it; "Cutting a release" in the README now uses it.
- **REPRODUCIBILITY** - a fourth subcohort, `related`, is analysed from the tumours in
  `DNA_PAIR_LIST_RELATED_TUMOURS_ALL`. SigProfiler and dNdScv outputs for it are published
  under `related_tumours` (`sigprofiler_subcohort_names` / `dndscv_subcohort_names`). A
  cohort whose list is empty produces no `related` outputs and no error.
- **INTEGRATION** - `somatic_variants.config` feeds the `related` subcohort from
  `DNA_PAIR_LIST_RELATED_TUMOURS_ALL`, which `run_somatic_variants.sh` now checks is
  exported; the MANUAL ENVIRONMENT OVERRIDES block and the README's standalone contract
  list it too (ten required exports).

### Changed
- **INTEGRATION** - **Breaking:** the launcher now fails at launch, naming the variables,
  unless `./source_me.sh` (or the wrapper's overrides) exports `PROJECT_DIR COMMANDS_DIR
  ANALYSIS_DIR STUDY PROJECT COHORT DNA_PAIR_LIST_ANALYSED_ALL
  DNA_PAIR_LIST_INDEPENDENT_TUMOURS_ALL DNA_PAIR_LIST_ONE_TUMOUR_PER_PATIENT_ALL
  DNA_PAIR_LIST_RELATED_TUMOURS_ALL` (plus the website/Slack variables when those
  toggles are on). The generator that writes `source_me.sh` has to emit them before a
  project can run this release.
- **INTEGRATION** - `metadata_manifest` is optional (default `null`). No step reads it, so
  the pipeline now checks the file exists only when it is set, and the config takes it
  from `COHORT_METADATA_FILE` when that is exported.
- **INTEGRATION** - **Breaking:** `assets/somatic_variants.config` takes its sample lists from
  the variables dermanager exports for them (`DNA_PAIR_LIST_*_ALL`) instead of rebuilding
  their paths from a filename convention, `metadata_manifest` from `COHORT_METADATA_FILE`
  instead of `${PROJECT_DIR}/metadata/${STUDY}-biosample-manifest-completed.tsv`, and
  the VCF inputs and every output directory from `${ANALYSIS_DIR}` instead of
  `${PROJECT_DIR}/analysis`. Outputs land in the same place for a dermanager project.
- **INTEGRATION** - **Breaking:** the launcher's directory moves from
  `${PROJECT_DIR}/somatic_pipeline` to `${PROJECT_DIR}/somatic_pipe`, and its
  config from `commands/somatic_variants.config` to `commands/somatic_pipe/somatic_variants.config`,
  matching the slug dermanager unpacks the asset bundle under. A run started under the
  old directory cannot `-resume` in the new one.
- **INTEGRATION** - `.github/workflows/publish-assets.yml` no longer moves the rolling
  `main-latest` / `develop-latest` tags: each is created once and only its bundle is
  replaced, so `git hf release finish` no longer fails on a moved tag.
- **INTEGRATION** - the repository is GitHub-primary: `manifest.homePage`, the README and
  the workflow header no longer point at or defer to GitLab. The GitLab container
  registry is unchanged.
- **INTEGRATION** - `docs/source/conf.py` records the pipeline version (was stale at 1.0.0).

### Removed
- **INTEGRATION** - the timestamp-only report name in `nextflow.config`.
- **INTEGRATION** - the stale launcher and config copies embedded in the README and the
  docs, which now document and link to `assets/`.

## [1.2.1] - 2026-08-27
### Added
- `.github/workflows/publish-assets.yml` publishes `assets/` to GitHub Releases as
  `projectify_asset_bundle.tar.gz` (and a `.sha256` of it) on every push to `main` and
  `develop` - as the rolling `main-latest` and `develop-latest` pre-releases - and on
  every `X.Y.Z` tag. `dermanager projectify` fetches assets from those release URLs
  instead of the GitHub API, which needs no token and is not rate limited. See
  "Asset release bundles" in the README.

## [1.2.0] - 2026-04-22
### Added
- `DNDSCV` subworkflow: per-subcohort significantly-mutated-gene discovery with dNdScv.
  - `MAF_TO_DNDSCV_INPUT` — converts the subcohort `sig_maf` into the 5-column dNdScv mutation table, with optional patient-level merging of sibling tumours sharing a PDXXXX prefix.
  - `DNDSCV_RUN` — runs dNdScv per subcohort; when a covariates file is supplied, runs twice (with/without covariates) and publishes into `with_covariates` / `without_covariates` sub-directories to match the manual-analysis layout.
- New params: `run_dndscv` (default `true`), `dndscv_outdir`, `dndscv_refdb`, `dndscv_covariates`, `dndscv_merge_by_patient`, `dndscv_subcohort_names`.
- `farm22` profile: defaults for `dndscv_refdb` / `dndscv_covariates` and an 8.GB memory boost for `DNDSCV_RUN`.
- `dndscv_outdir` set to `${PROJECT_DIR}/analysis/dndscv` in `assets/somatic_variants.config`.
- `run_dndscv: false` added to `tests/testdata/test_params.json` so the stub/test suite keeps skipping the new subworkflow.

## [1.1.0] - 2026-04-20
### Added
- `SIGNATURES` subworkflow: per-subcohort mutational-signature extraction with SigProfilerExtractor.
  - `MAF_TO_TARGETS` — derives SBS/DBS and ID target regions from the subcohort keep MAF.
  - `BUILD_SAMPLE_VCF` — pulls PASS variants from caveman (SBS/DBS) and pindel (ID) per sample, concats into a SigProfiler-ready VCF, warn-only MAF-vs-VCF count check.
  - `GROUP_SUBCOHORT_VCFS` — bgzip/tabix per-sample VCFs and `bcftools concat` into `VCFS_GROUPED/all.vcf` for the cohort.
  - `SIGPROFILER_EXTRACT` — runs SigProfilerExtractor on the subcohort VCFs (host module `sigprofiler/1.1.21-virtual-environment`); emits `results/`, and the matrix-generator `VCFS/input` / `VCFS/output` artefacts.
- `QC_VARIANTS` now emits a `sig_maf` channel (`keep_vaf_size_filt_matched_caveman_pindel_*.maf`) used as the signature-calling input.
- New params: `run_signatures` (default `true`), `sigprofiler_outdir`, `sigprofiler_seed`, `sigprofiler_subcohort_names`.
- `sigprofiler_subcohort_names` maps subcohort keys → legacy publish-dir names (e.g. `onePerPatient` → `one_tumour_per_patient`, `independent` → `independent_tumours`) to match the manual-analysis layout; unmapped keys fall back to the raw key.
- `sigprofiler_outdir` set to `${PROJECT_DIR}/analysis/sigprofiler` in `assets/somatic_variants.config`.
- `farm22` profile: HPC modules for `BUILD_SAMPLE_VCF` / `GROUP_SUBCOHORT_VCFS` and `long` queue for `SIGPROFILER_EXTRACT`.


## [1.0.0] - 2026-01-05
### Added
- Removed opinionated cohort grouping for running on multiple arbitary sample subsets 
- Updated SOP documentation to build within repository




## [0.6.4] - 2025-10-06
### Added
- Asset file for initialisation with dermanager and multi-pipeline running


## [0.6.3] - 2025-07-11
### Changed
- Updated cacheDir
- Updated metadata and default parameters
- Updated README to use new github repo
- Fixes to tests

## [0.6.2] - 2025-03-11
### Changed
- Typo in OPPT to correct publication with MAF release scripts

## [0.6.1] - 2025-02-14
### Changed
- Updated QC to fix AF filtering for alternative transcripts

## [0.6.0] - 2025-02-13
### Fixed
- Omission of sample_matched and other kinds of file lists from QC step. 
- Documentation updates to ensure things work smoothly around existing Dermatlas manual QC.
- Fixed publication directory for log 


## [0.5.2] - 2025-01-17
### Changed
- Fixes for pulling from a pipeline address and integration testing.

## [0.5.1] - 2024-12-13
### Changed
- Updated pipeline name and docs to better reflect what it does

## [0.5.0] - 2024-09-18
### Added
- A CI running on secure lustre for nf-test
### Changed
- Tests updated for the new CI.

## [0.4.0] - 2024-09-12
### Added
- Fixes to the container conflicts caused by QC/MAF cross calling 
- Added support for alternative canonical transcript labels 
- Corrected to use same module call whilst on farm22 as previously

## [0.3.0] - 2024-07-29
### Added
- Support for tsv -> .xlsx conversion
- Optional execution of independent/one tumor per patient/all sample analyses

## [0.2.0] - 2024-07-26
### Added
- Calculate and output tumour mutation burden
- Harmonised script name inside and outside of container (MAF) added here

## [0.1.1] - 2024-07-16
### Changed
- Update readme for reliable Farm singularity

## [0.1.0] - 2024-07-15
### Changed
- BCFtools container changed to `bcftools:1.20` for dockerising steps.

### Added 
Initial version of this pipeline. Support for 
- Basic variant filtering and adding annotations for dbSNP155 common variants (Caveman & Pindel)
- Make MAF files and run variant call QC

