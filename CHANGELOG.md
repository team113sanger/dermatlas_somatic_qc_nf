# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# 1.1.0 [April 20th 2026]
### Added
- `DNDSCV` subworkflow: per-subcohort significantly-mutated-gene discovery with dNdScv.
  - `MAF_TO_DNDSCV_INPUT` — converts the subcohort `sig_maf` into the 5-column dNdScv mutation table, with optional patient-level merging of sibling tumours sharing a PDXXXX prefix.
  - `DNDSCV_RUN` — runs dNdScv per subcohort; when a covariates file is supplied, runs twice (with/without covariates) and publishes into `with_covariates` / `without_covariates` sub-directories to match the manual-analysis layout.
- New params: `run_dndscv` (default `true`), `dndscv_outdir`, `dndscv_refdb`, `dndscv_covariates`, `dndscv_merge_by_patient`, `dndscv_subcohort_names`.
- `farm22` profile: defaults for `dndscv_refdb` / `dndscv_covariates` and an 8.GB memory boost for `DNDSCV_RUN`.
- `dndscv_outdir` set to `${PROJECT_DIR}/analysis/dndscv` in `assets/somatic_variants.config`.

# 1.1.0 [April 20th 2026]
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


# 1.0.0 [January 5th 2026]
### Added
- Removed opinionated cohort grouping for running on multiple arbitary sample subsets 
- Updated SOP documentation to build within repository




# 0.6.4 [October 6th 2025]
### Added
- Asset file for initialisation with dermanager and multi-pipeline running


# 0.6.3 [July 11th 2025]
### Changed
- Updated cacheDir
- Updated metadata and default parameters
- Updated README to use new github repo
- Fixes to tests

# 0.6.2 [March 11th 2025]
### Changed
- Typo in OPPT to correct publication with MAF release scripts

# 0.6.1 [February 14th 2025]
### Changed
- Updated QC to fix AF filtering for alternative transcripts

# 0.6.0 [February 12th 2025]
### Fixed
- Omission of sample_matched and other kinds of file lists from QC step. 
- Documentation updates to ensure things work smoothly around existing Dermatlas manual QC.
- Fixed publication directory for log 


# 0.5.2 [December 13th 2024]
### Changed
- Fixes for pulling from a pipeline address and integration testing.

# 0.5.1 [December 13th 2024]
### Changed
- Updated pipeline name and docs to better reflect what it does

# 0.5.0 [September 18th 2024]
### Added
- A CI running on secure lustre for nf-test
### Changed
- Tests updated for the new CI.

# 0.4.0 [September 11th 2024]
### Added
- Fixes to the container conflicts caused by QC/MAF cross calling 
- Added support for alternative canonical transcript labels 
- Corrected to use same module call whilst on farm22 as previously

# 0.3.0 [July 2024]
### Added
- Support for tsv -> .xlsx conversion
- Optional execution of independent/one tumor per patient/all sample analyses

# 0.2.0 [July 2024]
### Added
- Calculate and output tumour mutation burden
- Harmonised script name inside and outside of container (MAF) added here

# 0.1.1 [July 2024]
### Changed
- Update readme for reliable Farm singularity

# 0.1.0 [June 2024]
### Changed
- BCFtools container changed to `bcftools:1.20` for dockerising steps.

### Added 
Initial version of this pipeline. Support for 
- Basic variant filtering and adding annotations for dbSNP155 common variants (Caveman & Pindel)
- Make MAF files and run variant call QC

