#!/bin/bash
#SBATCH -N 4
#SBATCH -n 192
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH -A loni_ceds3d
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J cox_2D_BL_tr2

date

module purge
module load proteus/1.8.1
#which python
#unset PYTHONPATH
#which python
#export PROTEUS_ARCH=linux
#export PYTHONPATH=/project/cekees/cekees/proteus/1.7.5/lib/python3.7/site-packages:$PYTHONPATH
#export PATH=${PROJECT}/bin:${PROJECT}/proteus/1.7.5/bin:${PATH}
#export LD_LIBRARY_PATH=${PROJECT}/proteus/1.7.5/lib:${PROJECT}/proteus/1.7.5/lib64:${LD_LIBRARY_PATH}

#. ~/.bashrc
mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.csv .
#cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
#cp $SLURM_SUBMIT_DIR/petsc.options.asm .
cp $SLURM_SUBMIT_DIR/*.sh .

#conda activate proteus-dev
#unset PYTHONPATH

#mpiexec -n 16 srun -N1 -n16 -c1 beji_battjes_so.py -l 5 -C "he=0.01 T=50.0" -O petsc.options.superlu_dist
#srun parun beji_battjes_so.py -l 5 -C "he=0.006 T=60.0 free_slip=True" -O petsc.options.superlu_dist
#srun parun cox_flume3D.py -F -l 5 -C "he=0.01" -O petsc.options.asm
#srun parun -l5 -v -p --TwoPhaseFlow cox_flume2DV.py -C "he=0.5 wave_type='Time'" -O petsc.options.superlu_dist
#srun parun -l5 -v -p --TwoPhaseFlow cox_flume2DV.py -C "he=0.5 wave_type='Time'" -O petsc.options.asm

srun parun -l5 -v -p --TwoPhaseFlow cox_flume2DV.py -C "he=0.02 mangrove_porous=False filename='inp_BL_TR2.csv' final_Time=49.99"

date

exit 0

