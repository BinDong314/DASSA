#Bin's Mac
#FT_DIR=/Users/dbin/work/test-fasttensor/fasttensor/build
#HDF5_DIR=/Users/dbin/work/test-fasttensor/hdf5-1.12.0/build
#FFTW_DIR=/opt/local
#DASH_DIR=/Users/dbin/opt/dash-0.4.0
#EIGEN3_DIR=/Users/dbin/work/test-fasttensor/fasttensor/examples/dassa/Tools3rd

#Lawrencium
FT_DIR=/clusterfs/bear/BinDong_DAS_Data/fasttensor/build
DASH_DIR=/clusterfs/bear/BinDong_DAS_Data/fasttensor/tools3rd/dash/build/opt/dash-0.4.0
EIGEN3_DIR=${PWD}/Tools3rd/


SHELL:=/bin/bash

CCC=mpicxx
OTHER_FLAGS=-O3  -std=c++17 -I./DasLib/

AU_FLAG=-I$(FT_DIR)/include -L$(FT_DIR)/lib -lArrayUDF
HDF5_FLAG=-I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5  -lhdf5_hl -lz
FFTW_FLAG=-I$(FFTW_DIR)/include -L$(FFTW_DIR)/lib -lfftw3 
DASH_FLAG=-I$(DASH_DIR)/include -L$(DASH_DIR)/lib -ldash-mpi -ldart-mpi -ldart-base -lpthread -DDASH_ENABLE_HDF5 -DDASH_MPI_IMPL_ID='mpich'
EIGEN3_FLAG=-I$(EIGEN3_DIR)
ALL_FLAGS= $(OTHER_FLAGS) $(AU_FLAG) $(HDF5_FLAG) $(FFTW_FLAG)  $(DASH_FLAG) $(EIGEN3_FLAG) 

.PHONY:all
all:stack tdms2h5 xcorrelation decimate similarity  create_test_data

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

create_test_data: create_test_data.cpp
	$(CCC) -o create_test_data create_test_data.cpp $(ALL_FLAGS)

clean:
	rm stack tdms2h5 xcorrelation decimate similarity create_test_data

test:
	./test.sh

