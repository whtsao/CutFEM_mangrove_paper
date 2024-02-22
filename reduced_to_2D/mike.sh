#!/bin/bash
#SBATCH -N 16
#SBATCH -n 1024
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH -A hpc_ceds3d
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J cox_reduced_to_2D

date

module purge
module load proteus/fct
module load intel/2021.5.0
module load mvapich2/2.3.7/intel-2021.5.0
module load gcc/11.2.0
export LD_LIBRARY_PATH=/home/packages/compilers/intel/compiler/2022.0.2/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
export MV2_HOMOGENEOUS_CLUSTER=1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.csv .
#cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
#cp $SLURM_SUBMIT_DIR/petsc.options.asm .
cp $SLURM_SUBMIT_DIR/*.sh .

srun parun -l5 -v -p --TwoPhaseFlow cox_flume2DV.py -C "he=0.04 mangrove_porous=True filename='inp_HD_TR1.csv'"

date

exit 0

