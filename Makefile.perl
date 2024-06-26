#Bin's Mac
#AU_DIR=/Users/dbin/work/test-fasttensor/fasttensor/build
#HDF5_DIR=/Users/dbin/work/test-fasttensor/hdf5-1.12.0/build
#FFTW_DIR=/opt/local
#DASH_DIR=/Users/dbin/opt/dash-0.4.0
#EIGEN3_DIR=/Users/dbin/work/test-fasttensor/fasttensor/examples/dassa/Tools3rd

#Cori
#FT_DIR=/project/projectdirs/m1248/fasttensor-new/build
#DASH_DIR=/project/projectdirs/m1248/fasttensor-new/tools3rd/dash/build/install
#EIGEN3_DIR=${PWD}/Tools3rd/
#FFTW_DIR=/project/projectdirs/m1248/fasttensor-new/tools3rd/fftw-3.3.10/build
#SHELL:=/bin/bash

#Cori
FT_DIR=/global/cfs/projectdirs/m1248/fasttensor-new/build
#DASH_DIR=/global/cfs/projectdirs/m1248/fasttensor-new/tools3rd/dash/build/install
EIGEN3_DIR=${PWD}/Tools3rd/
FFTW_DIR=/global/cfs/projectdirs/m1248/fasttensor-new/tools3rd/fftw-3.3.10/build
SHELL:=/bin/bash



CCC=CC
OTHER_FLAGS=-O3  -fopenmp  -std=c++17 -I./DasLib/

FT_FLAG=-I$(FT_DIR)/include -L$(FT_DIR)/lib -lFastensor
HDF5_FLAG=-I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5  -lhdf5_hl -lz
FFTW_FLAG=-I$(FFTW_DIR)/include -L$(FFTW_DIR)/lib -lfftw3 
#DASH_FLAG=-I$(DASH_DIR)/include -L$(DASH_DIR)/lib -ldash-mpi -ldart-mpi -ldart-base -lpthread -DDASH_ENABLE_HDF5 -DDASH_MPI_IMPL_ID='mpich'
EIGEN3_FLAG=-I$(EIGEN3_DIR)
ALL_FLAGS= $(OTHER_FLAGS) $(FT_FLAG) $(HDF5_FLAG) $(FFTW_FLAG)  $(DASH_FLAG) $(EIGEN3_FLAG) 



.PHONY:all
all:stack tdms2h5 xcorrelation decimate similarity  create_test_data tmatch tmatch-reoganized

stack:stack.cpp
	$(CCC) -o stack stack.cpp $(ALL_FLAGS)

tdms2h5:tdms2h5.cpp
	$(CCC) -o  tdms2h5 tdms2h5.cpp $(ALL_FLAGS)

xcorrelation:xcorrelation.cpp
	$(CCC) -o xcorrelation xcorrelation.cpp  $(ALL_FLAGS)

decimate:decimate.cpp
	$(CCC) -o decimate decimate.cpp  $(ALL_FLAGS)

similarity:similarity.cpp
	$(CCC) -o similarity similarity.cpp  $(ALL_FLAGS)

tmatch:tmatch.cpp
	$(CCC) -o tmatch tmatch.cpp $(ALL_FLAGS)

tmatch-reoganized:tmatch-reoganized.cpp
	$(CCC) -o tmatch-reoganized tmatch-reoganized.cpp $(ALL_FLAGS)


create_test_data: create_test_data.cpp
	$(CCC) -o create_test_data create_test_data.cpp $(ALL_FLAGS)

clean:
	rm stack tdms2h5 xcorrelation decimate similarity create_test_data tmatch-reoganized

test:
	srun -n 1  ./test.sh
