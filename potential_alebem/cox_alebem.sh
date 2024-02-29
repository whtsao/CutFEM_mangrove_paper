#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 1# specify 6 threads per process
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH -A hpc_ceds3d
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J cox_potential_alebem

date

module purge
module load intel/2021.5.0 # run on mike
#module load intel/19.0.5 # run on qbc

export HOME_DIR=/home/$USER
export WORK_DIR=/work/$USER

##SLURM_JOBID: Job ID number given to this job
##SLURM_JOB_NODELIST: List of nodes allocated to the job
##SLURM_SUBMIT_DIR: Directory where the sbatch command was executed

mkdir -p $WORK_DIR/wavebar_alebem.$SLURM_JOBID
cd $WORK_DIR/wavebar_alebem.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.sh .
cp $SLURM_SUBMIT_DIR/*.ipt .
cp $SLURM_SUBMIT_DIR/*.f90 .
cp $SLURM_SUBMIT_DIR/*.out .
cp $SLURM_SUBMIT_DIR/*.err .
ifort -qopenmp cox_alebem.f90 -mkl

./a.out

# Mark the time it finishes.
date
# exit the job
exit 0
