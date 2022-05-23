.. _install:


Installation Guide
==================


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



3, Notes to install DASSA on LAWRENCIUM cluster @ LBNL
-----------------------------------------------------

If you are the users of the `LAWRENCIUM cluster @ LBNL <https://sites.google.com/a/lbl.gov/high-performance-computing-services-group/lbnl-supercluster/lawrencium>`_ , below steps and softwares may be used to compile DASSA (as well as FasTensor) ::

   > module load gcc/7.4.0 hdf5/1.10.5-gcc-p fftw
   > edit Makefile to have correct setting FT_DIR/HDF5_DIR/FFTW_DIR/DASH_DIR
   > make


4, Notes to install DASSA on Cori Supercomputer @ NERSC
-------------------------------------------------------
If you are the users of the `Cori @ NERSC <https://www.nersc.gov/systems/cori/>`_ , below items are needed to be considered to compile DASSA (as well as FasTensor).


**Do NOT use "module load fftw" as the fftw on Cori is broken.**

Build out own fftw ::
   > wget http://www.fftw.org/fftw-3.3.10.tar.gz
   > tar zxvf fftw-3.3.10.tar.gz
   > cd fftw-3.3.10
   > ./configure --prefix=$PWD/build --enable-threads --disable-fortran --enable-shared --enable-sse2
   > make & make install

Then, we can compile DASSA: 
   > module load cray-hdf5-parallel/1.10.5.2
   > export HDF5_USE_FILE_LOCKING=FALSE
   > cd dassa
   > make -f Makefile.cori

Notes: "cray-hdf5-parallel/1.10.5.2" and "cray-hdf5-parallel/1.12.1.1" are available on Cori. But the  "cray-hdf5-parallel/1.12.1.1" will report below error in compiling. So please use "cray-hdf5-parallel/1.10.5.2" ::


   /usr/lib64/gcc/x86_64-suse-linux/7/../../../../x86_64-suse-linux/bin/ld:
   /global/project/projectdirs/m1248/fasttensor-new/src/ft_endpoint_hdf5.cpp:1007: 
   undefined reference to `H5Ovisit'

See more about H5Ovisit here
https://forum.hdfgroup.org/t/second-hdf5-1-12-0-alpha-release-is-available-for-testing/6480
