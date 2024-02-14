#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 00:10:00
#SBATCH -p workq
#SBATCH -A hpc_ceds3d
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J cox_2D_HD_020824_CS

date

module purge
module load intel/2021.5.0
module load mvapich2/2.3.7/intel-2021.5.0
module load gcc/11.2.0
module load proteus/fct
export LD_LIBRARY_PATH=/home/packages/compilers/intel/compiler/2022.0.2/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
export MV2_HOMOGENEOUS_CLUSTER=1


mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.csv .
#cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
#cp $SLURM_SUBMIT_DIR/petsc.options.asm .
cp $SLURM_SUBMIT_DIR/*.sh .

#srun parun -l5 -v -p --TwoPhaseFlow cox_flume2DV.py -C "he=0.5 wave_type='Time'" -O petsc.options.superlu_dist
#srun parun -l5 -v -p --TwoPhaseFlow cox_flume2DV.py -C "he=0.5 wave_type='Time'" -O petsc.options.asm

#srun parun -l5 -v -p --TwoPhaseFlow cox_flume2DV.py -C "he=0.02 mangrove_porous=True filename='inp_HD_TR1.csv'"
#srun parun -l5 -v -p --hotStart --TwoPhaseFlow cox_flume2DV.py -C "he=0.04 mangrove_porous=True filename='inp_HD_TR1.csv'"
srun parun -l5 -v -p --TwoPhaseFlow cox_flume2DV.py -C "he=0.2 mangrove_porous=True filename='inp_HD_TR1.csv'"


date

exit 0

