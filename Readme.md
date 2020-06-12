
This is simple README file for DASSA, which now has the same Copyright as ArrayUDF now.
Please see the Coryright at the end.


For the tips on Lawrencium, see the section 6 at the end 


1, **Install FastTensor (ArrayUDF)** 

  See Readme file at https://bitbucket.org/dbin_sdm/arrayudf-test/src/master/

2, **Install Dependency** 

   2.1, FFTW etc.

   You can either use existing modules (e.g., apt get insall fftw3, module load fftw3, etc.)
  or install your own (see DIY tips below in Section 7)
 

3, *Compile DASSA*

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
```

5, **DIY Tips on FFTW:**

```bash
> cd Tools3rd
> tar zxvf fftw-3.3.8.tar.gz
> cd fftw-3.3.8
> ./configure  --prefix=$PWD/buildproperties
> make install
```


6, **Tips on Lawrencium**
   
6.1 A quick start: 
   It uses the FastTensor installed by Bin at /clusterfs/bear/BinDong_DAS_Data/fasttensor/build

```bash
   > git clone https://dbin_sdm@bitbucket.org/dbin_sdm/dassa.git
   > cd dassa
   > make
```
   
   Then, go back to 5. above to run the test code

6.2 Some compile details for DIY

```bash 
   > git clone https://dbin_sdm@bitbucket.org/dbin_sdm/dassa.git
   > module load gcc/7.4.0 hdf5/1.10.5-gcc-p fftw
   > edit Makefile to have correct setting AU_DIR/HDF5_DIR/FFTW_DIR/DASH_DIR
   > make
```




**Copyright Notice ArrayUDF Copyright (c) 2017**

The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved. If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov. NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
