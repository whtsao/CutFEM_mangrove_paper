#!/bin/bash
#SBATCH -N 8
#SBATCH -n 384
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 48:00:00
#SBATCH -p workq
#SBATCH -A loni_ceds3d
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J CS_adj_p_gauges

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
srun parun -l5 -v -p --TwoPhaseFlow cox_flume2DV_qbc.py -C "he=0.1 mangrove_porous=True filename='inp_HD_TR1.csv'"


date

exit 0

