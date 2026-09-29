# Nextflow: Somatic variant calling pipeline

Somatic variant calling and post-processing for DERMATLAS can be run mostly with a single nextflow pipeline in a largely "set-and-forget" manner to reproduce the manual steps detailed in [DERMATLAS - Post-processing CaVEMan and Pindel calls]()

This document contains an overview of how to configure and run the pipeline; for a more detailed explanation of the pipeline inputs and requirements for running can be found within the project [README](https://github.com/team113sanger/dermatlas_somatic_qc_nf)

## Workflow Overview

1. **Generating pipeline inputs for Dermatlas studies** 
2. **Generating the cohort config file**
3. **Running the pipeline**
4. **Make a release folder**


### 1. Generating pipeline inputs for Dermatlas studies

The study-specific inputs required to run the somatic QC pipeline are: the sample vcfs for a study's matched tumour-normal pairs; sets of the sample groupings; and the study metadata sheet. See [dematlas analysis setup](https://dermatlas-docs-4861e2.pages.internal.sanger.ac.uk/user_docs/dermatlas_analysis_setup.html) for instructions on how to prepare these.

Somatic QC is run on several cohort sample lists: all analysed tumours, one tumour per patient, all independent tumours and related tumours (which is empty for most cohorts, and then produces no outputs).

### 2. Generating the cohort config file

The nextflow pipeline's config file encodes all the options and inputs we might want to pass to the pipeline. When a project is provisioned from the [Dermatlas cohorts page](https://team113.sanger.ac.uk/dermatlas/cohorts/) this is written for you to `commands/somatic_pipe/somatic_variants.config`, alongside the launcher `commands/somatic_pipe/run_somatic_variants.sh` and the project `source_me.sh` both of them read their locations from.

The defaults provided should be suitable for running out of the box but for some pipeline runs there are a handful of parameters that you might consider changing:


- The `cohort prefix` (which will be used in the labelling of output files)
- The path to the `caveman_vcfs`  and the path to the `pindel_vcfs`
- The paths to write the outputs of caveman and pindel filtering to (`caveman_outdir` and `pindel_outdir`; in Dermatlas this is normally the same directory that `caveman_vcfs` and `pindel_vcfs` reside in)
- The output directory to publish somatic variant QC results into
- The `sigprofiler_outdir` for SigProfilerExtractor signature-calling outputs (set `run_signatures = false` to skip that subworkflow)
- The `dndscv_outdir` for dNdScv significantly-mutated-gene outputs (set `run_dndscv = false` to skip that subworkflow)
- The cohort metadata manifest (`COHORT_METADATA_FILE`; optional and not yet read by any step)
- The all-tumour, one-tumour-per-patient, independent-tumour and related-tumour sample lists generated in stage 1 (`DNA_PAIR_LIST_*_ALL`)

The reference config is [`assets/somatic_variants.config`](https://github.com/team113sanger/dermatlas_somatic_qc_nf/blob/main/assets/somatic_variants.config);
every cohort-specific value in it is read from a variable exported by `source_me.sh`.

**Launching the nextflow pipeline**

From the project directory, submit the launcher to the Sanger farm22 LSF:

```bash
cd <project_dir> && bsub -o somatic.out -e somatic.err < commands/somatic_pipe/run_somatic_variants.sh
```

The bsub magic at the start of the wrapper script will send a nextflow "master job", which looks after all other jobs to the oversubscribed queue (where it can live in peace running for a long period without fear of termination). Nextflow will shortly start submitting jobs on your behalf to the relevant queues.

If you haven't initialised your project with dermanager, see "Without the website" in the
[README](https://github.com/team113sanger/dermatlas_somatic_qc_nf#without-the-website): the reference launcher is
[`assets/run_somatic_variants.sh`](https://github.com/team113sanger/dermatlas_somatic_qc_nf/blob/main/assets/run_somatic_variants.sh),
and its `REVISION` selects which release of the pipeline runs.

### Troubleshooting problem nextflow runs:

There are several reasons the somatic QC pipeline might fail including bugs in the pipeline; issues with LSF; or misconfiguration.  In most cases (especially when you suspect a farm/ LSF failure), simply re-submitting the pipeline with the same command will trigger the nextflow `-resume` directive and the pipeline will pick up where it left off. A resubmission while the previous run is still alive exits immediately with status 75.

It is often worth taking a glance at the pipeline logs (`<YOUR_PROJECT_DIR>/somatic_pipe/logs/nextflow-run-<RUN_ID>.log`) to follow and see what's going on, especially if things have failed. A successful run deletes its work directory, but a failed one keeps it.

When jobs fail, nextflow will provide the path to the directory a failed job was run in. I'd recommend inspecting the files in here with `ls -la` and printing some of the log files for the job with

```bash
cat .command.err
cat .command.out
cat .command.sh
```

**Make a release folder**

"Creating a release" allows you to copy out the key results from Somatic variant QC into a seperate folder. This can be performed in exactly the same way as detailed in [DERMATLAS - Post-processing CaVEMan and Pindel calls](/spaces/CAS/pages/116131145/DERMATLAS+-+Post-processing+CaVEMan+and+Pindel+calls)  

```bash
export PROJECT_DIR="/lustre/scratch125/casm/team113da/projects/dermatlas_pu11_project_dir/7651_3442_DERMATLAS_Superficial_acral_fibromyxoma_WES"
export STUDY=7651
export PROJECT=3442
export COHORT="Superficial_acral_fibromyxoma"
PREFIX="${STUDY}_${PROJECT}"
  
cd $PROJECT_DIR/analysis/variants_combined
i=1
mkdir release_v${i}
source ${PROJECT_DIR}/scripts/maf/source_me.sh
bash ${PROJECT_DIR}/scripts/maf/make_variant_release.sh \
$PROJECT_DIR ${COHORT} version${i} release_v${i} > release_v${i}/files.log
```
