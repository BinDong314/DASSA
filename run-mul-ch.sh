#!/bin/bash

#
# ##Before are job script 
# 
# Job name:
#SBATCH --job-name=xcorr
#
# Partition:
#SBATCH --partition=lr5
#
# QoS:
#SBATCH --qos=lr_normal
#
# Processors:
#SBATCH --ntasks=64
#
# Wall clock limit:
#SBATCH --time=00:30:00


#SBATCH -e xcorrelation_%j.err
#SBATCH -o xcorrelation_%j.out



 
config_file=xcorrelation-multi-channel.config
output_fn_header=/global/scratch/dbin/westSac_180103000031_xcorr


output_fn_tail=.h5

##More pattern to run loop
MastChannelSet=(1 3 5 7)
for MC in ${MastChannelSet[*]}
do
    change-ch.sh xcorrelation-test.config $MC ./westSac_180103000031_Xcorr_Master$MC.h5
    sed -i -e "s%\(master_index\).*%\1=$MC%"   $config_file
    sed -i -e "s%\(output_file\).*%\1=$output_fn_header$MC$output_fn_tail%"  $config_file
    
    ##Replace the below to more process
    ## mpirun -n 8 ./xcorrelation -c xcorrelation-test.config
    ## mpirun -n 64 ./xcorrelation -c xcorrelation-test.config

    ./xcorrelation -c xcorrelation-test.config
done

