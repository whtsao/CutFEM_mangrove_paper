#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 00:10:00
#SBATCH -p workq
#SBATCH -A hpc_ceds3d
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values

date

module purge
module load python/3.9.7-anaconda

mkdir -p $WORK/beji_alefem.$SLURM_JOBID
cd $WORK/beji_alefem.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .

python cox_alefem.py

exit 0

