#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 72:00:00
#SBATCH -p single
#SBATCH -A hpc_ceds3d
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J cox_potential_model

date

module purge
module load python/3.9.7-anaconda

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .

python cox_alefem.py

exit 0

