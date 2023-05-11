#!/bin/bash -l


#SBATCH -p debug
####SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:15:00  
#SBATCH -J tmatch
#SBATCH -e tmatch-reorg-jun4-1540-2tmp_%j.err
#SBATCH -o tmatch-reorg-jun4-1540-2tmp_%j.out
#SBATCH -D /project/projectdirs/m1248/dassa-new
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use Haswell nodes


#export DATA_DIR=/global/cscratch1/sd/dbin
export HDF5_USE_FILE_LOCKING=FALSE

#rm tmatch-verify-data/output/tmach-result-2-file-after-reorg.h5

srun -n 1 ./tmatch-reoganized -c  tmatch-reorg-demo-jun4-1540-2tmp.config

#h5diff -v -d 0.0000000001 tmatch-verify-data/output/tmach-result-2-file-after-reorg.h5 tmatch-verify-data/output/tmach-result-2-file-before-reorg.h5

