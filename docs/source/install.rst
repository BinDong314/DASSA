.. _install2:


Installation
============


1, Install Dependency
-----------------------

For the below dependent software,  a quick way is to use existing ones  the systems where you run the code. For example, on the Cori system @ NERSC, you may use 'module load fftw3' to use the pre-installed on there. Here are general steps to install these dependent software from scratch. 

1.1, FFTW 
^^^^^^^^^
For example, fftw-3.3.8 can be installed via ::

   > tar zxvf fftw-3.3.8.tar.gz
   > cd fftw-3.3.8
   > ./configure  --prefix=$PWD/build  #$PWD capture current working directory
   > make install   #the fftw will be installed with $PWD/build

2.2 FasTensor
^^^^^^^^^^^^^
Please see the Readme file from FasTensor https://bitbucket.org/dbin_sdm/fasttensor/src/master/ to install it and its Dependencies, e.g. HDF5 and DASH. Both HDF5 and DASH are compulsory to use DASSA.  
  
2, Compile DASSA
----------------

Edit the Makefile within DASSA source code to have correct setting for below items::

   FT_DIR=
   HDF5_DIR=
   FFTW_DIR=
   DASH_DIR=

'FT_DIR' points to the installation directory of FasTensor. 'HDF5_DIR' points to the  installation directory of HDF5. 'FFTW_DIR' points to the  installation directory of FFTW. 'DASH_DIR' points to the  installation directory of DASH.
Note that these setting may not be necessary when the system exports these variables. 

After setting proper variables in Makefile, it uses the 'make' command to compile the code ::

 > make
 
All executable code will be produced under the current directory.


