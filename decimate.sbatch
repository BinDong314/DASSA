#!/bin/bash
# Job name:
#SBATCH --job-name=decimate
# QoS:
#SBATCH --qos=lr_debug
#
# Partition:
#SBATCH --partition=lr6
# 
#
# Processors:
#SBATCH --nodes=3
#SBATCH -n 3
#SBATCH --ntasks-per-node=1
#
# Wall clock limit:
#SBATCH --time=00:30:00
#

## Run command
##module load openmpi
module load gcc/7.4.0 hdf5/1.10.5-gcc-p fftw boost
mpirun -n 3 ./decimate -c decimate-vero-10f.config

