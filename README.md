# dermatlas_somatic_qc_nf

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.04.5-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

dermatlas_somatic_qc_nf is a bioinformatics pipeline written in [Nextflow](http://www.nextflow.io) for performing processing and QC on somatic variant calls from cohorts of FFPE tumors within the Dermatlas project. 

## Pipeline summary

In brief, the pipeline takes the Caveman and Pindel VCF files for a set samples – which have been pre-processed by the Dermatlas ingestion pipeline – and then:
- Links each sample vcf to it's associated metadata.
- Filters `PASS` variants from the file.
- Adds an annotation field to variants present in dbSNP 
- Performs Dermatlas variant-QC and generates Dermatlas diagnostic plots 
- Calculates the TMB of Dermatlas `keep` samples produced by Dermatlas variant-QC
- Creates `.xlsx` file outputs from mafs for releasing to project scientists

## Inputs 

### Cohort-dependent variables

- `caveman_vcfs`: path to a set of Caveman vcf files (using **.vcf expansion)
- `pindel_vcfs`: path to a set of Pindel vcf files (using **.vcf expansion)
- `metadata_manifest`: path to a tab-delimited manifest containing sample PD IDs and information about sample phenotype/preparation.
- `cohort_prefix`: Prefix to add to output file names
- `exome_size`: Size in Mb of the baitset (for Dermatlas this is `48.225157`)
- `outdir`: Directory to publish results 
- `caveman_outdir`: Directory to publish the results of variant processing to match Dermatlas conventions (typically the `analysis` dir for a publishable unit).
- `pindel_outdir`: Directory to publish the results of variant processing to match Dermatlas conventions (typically the `analysis` dir for a publishable unit).
- `release_version`: Directory to release results into within an output directory (e.g.`version1`)

**Optional**
- `all_samples`: path to a file containing a tab-delimited list of all matched tumour-normal pairs in a cohort.
- `one_per_patient`: path to a file containing a tab-delimited list of matched tumour-normal pairs with one tmor selected per-patient.
- `independent`: path to a file containing a tab-delimited list of matched tumour-normal pairs with all independent comparisons to perform.
- `alternative_transcripts`: path to a file containing a tab-delimited list of HUGO gene symbol - transcript ID pairs for correcting the transcript considered canonical.


### Cohort-independent variables
Reference files that are reused across pipeline executions have been placed within the pipeline's default `nextflow.config` file to simplify configuration and can be ommited from setup. Behind the scences, the following reference files are required for a run: 
- `dbsnp_variants`: path to DBSNP vcf file and it's `.tbi` index file (`dbSNP155_common.tsv.gz{,.tbi}`)
- `dbsnp_header`: Path to a file detailing dbsnp header info
- `genome_build`: Genome build string (`GRCh38`) to use in somatic QC steps
- `filtering_column`: Column within VCF to use in filtering likely germline variants (default: `gnomAD_AF`)
- `filter_option`: String to determine the mode of filter applied to variants (One of `filter1` or `filter2`)

Default reference file values supplied within the `nextflow.config` file can be overided by adding them to the params `.json` file. An example complete params file `example_params.json` is supplied within this repo for demonstation.

## Usage 

The recommended way to launch this pipeline is using a wrapper script (e.g. `bsub < my_wrapper.sh`) that submits nextflow as a job and records the version (**e.g.** `-r 0.6.0`)  and the `.json` parameter file supplied for a run.

An example wrapper script:
```
#!/bin/bash
#BSUB -q oversubscribed
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000
#BSUB -oo logs/somatic_variants_pipeline%J.o
#BSUB -eo logs/somatic_variants_pipeline%J.e

PARAMS_FILE="/lustre/scratch125/casm/team113da/users/jb63/home/ubuntu/projects/dermatlas_somatic_qc_nf/example_params.json"

# Load module dependencies
module load nextflow-23.10.0
module load /software/modules/ISG/singularity/3.11.4

# Create a nextflow job that will spawn other jobs

nextflow run 'https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/dermatlas_mafqc_nf' \
-r 0.6.0 \
-params-file $PARAMS_FILE \
-profile farm22 
```


When running the pipeline for the first time on the farm you will need to provide credentials to pull singularity containers from the team113 sanger gitlab. You should be able to do this by running
```
module load singularity/3.11.4 
singularity remote login --username $(whoami) docker://gitlab-registry.internal.sanger.ac.uk
```

The pipeline can configured to run on either Sanger OpenStack secure-lustre instances or farm22 by changing the profile speicified:
`-profile secure_lustre` or `-profile farm22`. 

## Pipeline visualisation 
Created using nextflow's in-built visualitation features.
```nextflow run main.nf -preview -with-dag -params-file tests/testdata/test_params.json flowchart.mmd```

```mermaid
flowchart TB
subgraph " "
    v0["Channel.fromPath"]
    v2["Channel.fromPath"]
    v4["Channel.fromPath"]
    v6["Channel.fromPath"]
    v9["bedfile"]
    v11["dbsnp_vars"]
    v12["header"]
    v15["Channel.fromPath"]
    v32["BUILD"]
    v33["AF_COL"]
    v34["filter"]
    v40["exome_size"]
    v45["Channel.fromPath"]
    v62["BUILD"]
    v63["AF_COL"]
    v64["filter"]
    v70["exome_size"]
    v75["Channel.fromPath"]
    v92["BUILD"]
    v93["AF_COL"]
    v94["filter"]
    v100["exome_size"]
    end
    subgraph " "
    v1["patient_md"]
    v3["metadata"]
    v36[" "]
    v37[" "]
    v38[" "]
    v42[" "]
    v44[" "]
    v66[" "]
    v67[" "]
    v68[" "]
    v72[" "]
    v74[" "]
    v96[" "]
    v97[" "]
    v98[" "]
    v102[" "]
    v104[" "]
    end
    subgraph PROCESS_VCFS
    v10([FILTER_PASS_VARIANTS])
    v13([ADD_COMMON_ANNOTATIONS])
    v5(( ))
    v14(( ))
    end
    subgraph ALL_TUMORS
    v35([QC_VARIANTS])
    v41([CALCULATE_SAMPLE_TMB])
    v43([MAF_TO_EXCEL])
    v39(( ))
    end
    subgraph ONE_TUMOR_PER_PATIENT
    v65([QC_VARIANTS])
    v71([CALCULATE_SAMPLE_TMB])
    v73([MAF_TO_EXCEL])
    v69(( ))
    end
    subgraph INDEPENDENT_TUMORS
    v95([QC_VARIANTS])
    v101([CALCULATE_SAMPLE_TMB])
    v103([MAF_TO_EXCEL])
    v99(( ))
    end
    v0 --> v1
    v2 --> v3
    v4 --> v5
    v6 --> v5
    v9 --> v10
    v5 --> v10
    v10 --> v13
    v11 --> v13
    v12 --> v13
    v13 --> v14
    v15 --> v14
    v32 --> v35
    v33 --> v35
    v34 --> v35
    v14 --> v35
    v35 --> v38
    v35 --> v37
    v35 --> v36
    v35 --> v39
    v40 --> v41
    v39 --> v41
    v41 --> v42
    v39 --> v43
    v43 --> v44
    v45 --> v14
    v62 --> v65
    v63 --> v65
    v64 --> v65
    v14 --> v65
    v65 --> v68
    v65 --> v67
    v65 --> v66
    v65 --> v69
    v70 --> v71
    v69 --> v71
    v71 --> v72
    v69 --> v73
    v73 --> v74
    v75 --> v14
    v92 --> v95
    v93 --> v95
    v94 --> v95
    v14 --> v95
    v95 --> v98
    v95 --> v97
    v95 --> v96
    v95 --> v99
    v100 --> v101
    v99 --> v101
    v101 --> v102
    v99 --> v103
    v103 --> v104
```

## Testing

This pipeline has been developed with the [nf-test](http://nf-test.com) testing framework. Unit tests and small test data are provided within the pipeline `test` subdirectory. A snapshot has been taken of the outputs of most steps in the pipeline to help detect regressions when editing. You can run all tests on openstack with:

```
nf-test test 
```
and individual tests with:
```
nf-test test tests/modules/x.nf.test
```

For faster testing of the flow of data through the pipeline **without running any of the tools involved**, stubs have been provided to mock the results of each succesful step.
```
nextflow run main.nf \
-params-file params.json \
-c tests/nextflow.config \
--stub-run
```


