#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 00:10:00
#SBATCH -p workq
#SBATCH -A hpc_ceds3d
#SBATCH -J waveforce_tree
#SBATCH -o o.out
#SBATCH -e e.err
#load proteus module and ensure proteus's python is in path

date

module purge
module load intel/2021.5.0 
module load mvapich2/2.3.7/intel-2021.5.0   
module load gcc/11.2.0
module load proteus/fct
export LD_LIBRARY_PATH=/home/packages/compilers/intel/compiler/2022.0.2/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
export MV2_HOMOGENEOUS_CLUSTER=1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cd "$SLURM_SUBMIT_DIR"
cp *.py $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.pyx $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.sh $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID

python setup.py build_ext -i
srun parun --TwoPhaseFlow waveforce_tree.py -F -l 5 -C "he=0.2"
exit 0
