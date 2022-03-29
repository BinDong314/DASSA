
This is simple README file for DASSA, which now has the same Copyright as ArrayUDF now.
Please see the Coryright at the end.


**For the tips on Lawrencium, see the section 6 at the end** 


**1, Install FastTensor (ArrayUDF)** 

  See Readme file at the below link
  https://bitbucket.org/dbin_sdm/fasttensor/src/master/

**2, Install Dependency** 

   2.1, FFTW etc.

   You can either use existing modules (e.g., ``apt get insall fftw3``, ``module load fftw3``, etc.) or install your own (see DIY tips below in Section 5)
 

**3, Compile DASSA**

Edit the Makefile to to have correct setting for below items
```bash
   AU_DIR=/Users/dbin/work/FastTensor/build
   HDF5_DIR=/Users/dbin/work/soft/hdf5-1.10.5/build
   FFTW_DIR=/Users/dbin/work/test-fasttensor/fasttensor/examples/dassa/Tools3rd/fftw-3.3.8/build
   DASH_DIR=/Users/dbin/opt/dash-0.4.
```

Then
```bash
> make
``` 

**4, Simple test**

```bash
./stack -i /clusterfs/bear/BinDong_DAS_Data/xcorr_examples_h5
```

-i: the directory for input files 
  
Output files will be in current directoy

```
    xcorr_examples_h5_stack_data_in_sum.h5
    xcorr_examples_h5_stack_phaseWeight.h5
    xcorr_examples_h5_stack_semblanceWeight.h5
    xcorr_examples_h5_stack_semblance_denom_sum.h5
    xcorr_examples_h5_stack_final_pwstack.h5
```

Please see the section 7 for more about control for input, output control, and runtime parameters.

**5, How to install FFTW:**

```bash
> cd Tools3rd
> tar zxvf fftw-3.3.8.tar.gz
> cd fftw-3.3.8
> ./configure  --prefix=$PWD/buildproperties
> make install
```


**6, How to install on Lawrencium**
   
6.1 Install FastTensor: 

   You need to first install FastTensor following the method in Section 1.  

   Or, if you have access, you can use the FastTensor installed by Bin at /clusterfs/bear/BinDong_DAS_Data/fasttensor/build
   
   Then, check out the Dassa code and try the make.

```bash
   > git clone https://dbin_sdm@bitbucket.org/dbin_sdm/dassa.git
   > cd dassa
   > make
```
   
   If it works, go back to Section 5 above to run the test code

   If it does not work, following below steps to check environment or edit your Makefile

```bash 
   > module load gcc/7.4.0 hdf5/1.10.5-gcc-p fftw
   > edit Makefile to have correct setting AU_DIR/HDF5_DIR/FFTW_DIR/DASH_DIR
   > make
```

**7, Input/Output/Runtime Parameters**

The stack.config file provides an example to specify the Input/Output/Runtime parameters.
To use the config file:

```bash
  > stack -c stack.config
```

The stack.config is the example configure files.
I hope that their names are self-identified. 
There are some ways I encode the parameters for sub functions, e.g.

```bash
CausalityFlagging_ButterLow_order = 3
```

This line is to set the order parameter for Butterworth filter (low pass) and it is called by  CausalityFlagging.

Another assumption in the input parameter are dataset (variable) name (i.e.,  below two parameters):

```bash
xcorr_input_dataset_name = /xcoor ; dataset name within all XCORR files
stack_output_file_dataset_name = /data ;dataset name for all output files
```

To simplify the case, the two dataset names apply to all input files and output files. 
Their file names can  be specified for different variables. 


** 8 Notes on NERSC
Do NOT use "module load fftw" as the fftw on Cori is broken.
We need to build them by ourself

```bash
> wget http://www.fftw.org/fftw-3.3.10.tar.gz
> tar zxvf fftw-3.3.10.tar.gz
> cd fftw-3.3.10
> ./configure --prefix=$PWD/build --enable-threads --disable-fortran --enable-shared --enable-sse2
> make & make install
```

Then, compile DASSA
```dash
> module load cray-hdf5-parallel/1.10.5.2
> export HDF5_USE_FILE_LOCKING=FALSE
> cd dassa
> make -f Makefile.cori
> ./xcorrelation
```
Notes: "cray-hdf5-parallel/1.10.5.2" and "cray-hdf5-parallel/1.12.1.1"
But the  "cray-hdf5-parallel/1.12.1.1" will report below error in compiling
```dash
/usr/lib64/gcc/x86_64-suse-linux/7/../../../../x86_64-suse-linux/bin/ld:
/global/project/projectdirs/m1248/fasttensor-new/src/ft_endpoint_hdf5.cpp:1007: 
undefined reference to `H5Ovisit'
```
See more about H5Ovisit here
https://forum.hdfgroup.org/t/second-hdf5-1-12-0-alpha-release-is-available-for-testing/6480


Example on Cori
```dash
> sbatch xcorrelation-cori.bs
  which uses "xcorrelation-cori.config" as config file
  The "xcorrelation-cori.config" run corr on two files in my $SCRATCH /global/cscratch1/sd/dbin/DarkFiber_ambientTest_201110/
  check progress from file xcorr_56717603.out and xcorr_56717603.err
> h5dump  /global/cscratch1/sd/dbin/DarkFiber_ambientTest_201110_sample_output/DF__UTC_20201112_001132.602.h5
  check result at /global/cscratch1/sd/dbin/DarkFiber_ambientTest_201110_sample_output/
```


****************************


DASSA Copyright (c) 2021, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.


****************************


*** License Agreement ***

DASSA Copyright (c) 2021, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches,
or upgrades to the features, functionality or performance of the source
code ("Enhancements") to anyone; however, if you choose to make your
Enhancements available either publicly, or directly to Lawrence Berkeley
National Laboratory, without imposing a separate written license agreement
for such Enhancements, then you hereby grant the following license: a
non-exclusive, royalty-free perpetual license to install, use, modify,
prepare derivative works, incorporate into other computer software,
distribute, and sublicense such enhancements or derivative works thereof,
in binary and source code form.