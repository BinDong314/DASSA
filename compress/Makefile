
#Bin's Mac

#AU_DIR=/Users/dbin/work/fastensor-dassa/build
AU_DIR=/Users/dbin/work/fastensor/build
HDF5_DIR=/Users/dbin/work/soft/hdf5-1.10.7/build
FFTW_DIR=/Users/dbin/work/soft/fftw-3.3.9/build
DASH_DIR=/Users/dbin/work/soft/dash/build/install
EIGEN3_DIR=/Users/dbin/work/dassa/Tools3rd
SZ_DIR=/Users/dbin/work/soft/SZ/build


#Lawrencium
#AU_DIR=/clusterfs/bear/BinDong_DAS_Data/fasttensor/build
#DASH_DIR=/clusterfs/bear/BinDong_DAS_Data/fasttensor/tools3rd/dash/build/opt/dash-0.4.0
#EIGEN3_DIR=${PWD}/Tools3rd/


SHELL:=/bin/bash

CCC=mpicxx
OTHER_FLAGS=-O3  -std=c++17 -I./DasLib/

SZ_FLAG=
#-I$(SZ_DIR)/include -L$(SZ_DIR)/lib -Wl,-rpath,$(SZ_DIR)/lib  -Wl,-rpath,$(SZ_DIR) -lSZ   -lzlib -lhdf5sz -ldl
AU_FLAG=-I$(AU_DIR)/include -L$(AU_DIR)/lib -lFastensor
HDF5_FLAG=-I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5  -lhdf5_hl -lz
FFTW_FLAG=-I$(FFTW_DIR)/include -L$(FFTW_DIR)/lib -lfftw3 
DASH_FLAG=-I$(DASH_DIR)/include -L$(DASH_DIR)/lib -ldash-mpi -ldart-mpi -ldart-base -lpthread -DDASH_ENABLE_HDF5 -DDASH_MPI_IMPL_ID='mpich' -DHAS_DASH_ENDPOINT
EIGEN3_FLAG=-I$(EIGEN3_DIR)
ALL_FLAGS= $(OTHER_FLAGS) $(AU_FLAG) $(HDF5_FLAG) $(FFTW_FLAG)  $(DASH_FLAG) $(EIGEN3_FLAG)  $(SZ_FLAG)

.PHONY:all
all:  compress compress-ushort

compress:compress.cpp
	$(CCC)   -o compress compress.cpp $(ALL_FLAGS) 

compress-ushort:compress-ushort.cpp
	$(CCC)   -o compress-ushort compress-ushort.cpp $(ALL_FLAGS) 
decompress:decompress.cpp
	$(CCC)  -o decompress decompress.cpp $(ALL_FLAGS)
 	
clean:
	rm -f compress decompress compress-ushort

