#!/bin/bash

echo "Clear test data ..."
##	#./create_test_data -f  ./test-data/test-das.h5  -d /dat -n 2 -s 12,16 -t 2 -r -m
##Clean data
rm ./test-data/dir-output/*
rm ./test-data/dir-decimate/*
###Create test data ends

echo "Start to run test  ..."

# To test of running a command
# $1 command name
# $2 output data of the command (HDF5 format)
# $3 correct output of the command (HDF5 format)
function run_command(){
    $1 > /dev/null 2>&1
    if [ $? -eq 0 ]
    then
        echo "Test run $1  ... [PASSED]"
    else
        echo "Test run $1  ... [FAILED]" >&2
    fi

    if [ "$#" -gt 1 ]; then
        h5diff  -v ./test-data/$5/$2 ./test-data-good/$4 > /dev/null 2>&1
        if [ $? -eq 0 ]
        then
            echo "Checked output of $2 [PASSED]"
        else
            echo "Checked output of $2 [FAILED]" >&2
        fi   

        h5diff ./test-data/$5/$3 ./test-data-good/$4 > /dev/null 2>&1
 	if [ $? -eq 0 ]
        then
            echo "Checked output of $3 [PASSED]"
        else
            echo "Checked output of $3 [FAILED]" >&2
        fi   
    fi
}


## Run the test code
run_command ./xcorrelation test-das-1.h5  test-das-2.h5 test-das-xcorr.h5 dir-output
run_command ./decimate     test-das-1.h5  test-das-2.h5 test-das-dec.h5 dir-decimate

