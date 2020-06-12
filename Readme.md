
This is simple README file for DASSA, which now has the same Copyright as ArrayUDF now.
Please see the Coryright at the end.


**For the tips on Lawrencium, see the section 6 at the end** 


1, **Install FastTensor (ArrayUDF)** 

  See Readme file at the below link
  https://bitbucket.org/dbin_sdm/fasttensor/src/master/

2, **Install Dependency** 

   2.1, FFTW etc.

   You can either use existing modules (e.g., ``apt get insall fftw3``, ``module load fftw3``, etc.) or install your own (see DIY tips below in Section 5)
 

3, **Compile DASSA**

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

4, **Simple test**

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

5, **How to install FFTW:**

```bash
> cd Tools3rd
> tar zxvf fftw-3.3.8.tar.gz
> cd fftw-3.3.8
> ./configure  --prefix=$PWD/buildproperties
> make install
```


6, **How to install on Lawrencium**
   
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

7, **Input/Output/Runtime Parameters**


**Copyright Notice ArrayUDF Copyright (c) 2017**

The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved. If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov. NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
