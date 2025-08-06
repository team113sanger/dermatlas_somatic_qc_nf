#!/bin/bash
#BSUB -q oversubscribed
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000

set -euo pipefail
 
source source_me.sh 
 
export REVISION="0.6.3"
 
# Load module dependencies
module load nextflow-23.10.0
module load /software/modules/ISG/singularity/3.11.4
  
# Create a nextflow job that will spawn other jobs
nextflow pull 'https://github.com/team113sanger/dermatlas_somatic_qc_nf'

nextflow run 'https://github.com/team113sanger/dermatlas_somatic_qc_nf' \
-resume \
-r ${REVISION} \
-c ${PROJECT_DIR}/commands/somatic_variants.config \
-profile farm22