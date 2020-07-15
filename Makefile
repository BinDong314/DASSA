#Bin's Mac
#AU_DIR=/Users/dbin/work/test-fasttensor/fasttensor/build
#HDF5_DIR=/Users/dbin/work/test-fasttensor/hdf5-1.12.0/build
#FFTW_DIR=/opt/local
#DASH_DIR=/Users/dbin/opt/dash-0.4.0
#EIGEN3_DIR=/Users/dbin/work/test-fasttensor/fasttensor/examples/dassa/Tools3rd

#Lawrencium
AU_DIR=/clusterfs/bear/BinDong_DAS_Data/fasttensor/build
DASH_DIR=/clusterfs/bear/BinDong_DAS_Data/fasttensor/tools3rd/dash/build/opt/dash-0.4.0
EIGEN3_DIR=${PWD}/Tools3rd/


SHELL:=/bin/bash

CCC=mpicxx
OTHER_FLAGS=-O3  -std=c++17 -I./DasLib/

AU_FLAG=-I$(AU_DIR)/include -L$(AU_DIR)/lib -lArrayUDF
HDF5_FLAG=-I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5  -lhdf5_hl -lz
FFTW_FLAG=-I$(FFTW_DIR)/include -L$(FFTW_DIR)/lib -lfftw3 
DASH_FLAG=-I$(DASH_DIR)/include -L$(DASH_DIR)/lib -ldash-mpi -ldart-mpi -ldart-base -lpthread -DDASH_ENABLE_HDF5 -DDASH_MPI_IMPL_ID='mpich'
EIGEN3_FLAG=-I$(EIGEN3_DIR)
ALL_FLAGS= $(OTHER_FLAGS) $(AU_FLAG) $(HDF5_FLAG) $(FFTW_FLAG)  $(DASH_FLAG) $(EIGEN3_FLAG) 

.PHONY:all
all:stack tdms2h5 xcorrelation

stack:stack.cpp
	module load gcc/7.4.0 hdf5/1.10.5-gcc-p fftw; \
	$(CCC) -o stack stack.cpp $(ALL_FLAGS)

tdms2h5:tdms2h5.c
	module load gcc/7.4.0 hdf5/1.10.5-gcc-p fftw; \
	$(CCC) -o  tdms2h5 tdms2h5.c $(ALL_FLAGS) 

xcorrelation:xcorrelation.cpp
	module load gcc/7.4.0 hdf5/1.10.5-gcc-p fftw boost; \
	$(CCC) -o xcorrelation xcorrelation.cpp  $(ALL_FLAGS) 
	
clean:
	rm stack tdms2h5 xcorrelation

