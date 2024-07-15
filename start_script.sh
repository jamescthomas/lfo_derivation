#!/bin/bash

# starts the pipeline for Biobank analysis. updated by JCT on 2023/04/05

#SBATCH -J full_stream
#SBATCH -p short
#SBATCH -a 1-41415%450
#SBATCH -o /well/webb/users/vox025/logs/full_stream_%j.out 

module load FSL/5.0.11-foss-2018b-Python-3.6.6
source $FSLDIR/etc/fslconf/fsl.sh

module load AFNI/18.1.09-foss-2018b-Python-3.6.6

module load matlab/2019a

SEEDFILE=/well/webb/users/vox025/subj_all.txt
eid=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" $SEEDFILE)

echo `date`: Executing task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID} on `hostname` as user ${USER}

bash /well/webb/users/vox025/full_stream_.sh $eid
