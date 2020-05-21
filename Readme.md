
This is simple README file for DASSA, which now have the same copyright as ArrayUDF (https://bitbucket.org/dbin_sdm/arrayudf-test/src/master/)


For the tips on Lawrencium, see the section 8 at the end 


1, Install FastTensor (ArrayUDF) 

  See Readme file at https://bitbucket.org/dbin_sdm/arrayudf-test/src/master/

2, Install Dependency 

   2.1, FFTW etc.

   You can either use existing modules (e.g., apt get insall fftw3, module load fftw3, etc.)
  or install your own (see DIY tips below in Section 7)
 

3, Compile DASSA

Edit the Makefile to to have correct setting for below items
```properties
   AU_DIR=/Users/dbin/work/FastTensor/build
   HDF5_DIR=/Users/dbin/work/soft/hdf5-1.10.5/build
   FFTW_DIR=/Users/dbin/work/test-fasttensor/fasttensor/examples/dassa/Tools3rd/fftw-3.3.8/build
   DASH_DIR=/Users/dbin/opt/dash-0.4.
```

Then
```properties
> make
``` 

4, Simple test

```properties
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

5, DIY Tips on FFTW:

```properties
> cd Tools3rd
> tar zxvf fftw-3.3.8.tar.gz
> cd fftw-3.3.8
> ./configure  --prefix=$PWD/buildproperties
> make install
```


6, Tips on Lawrencium
   
   6.1 A quick start:
   
   It uses the FastTensor installed by Bin at /clusterfs/bear/BinDong_DAS_Data/fasttensor/build

   ```properties
   > git clone https://dbin_sdm@bitbucket.org/dbin_sdm/dassa.git
   > cd dassa
   > make
   ```
   
   Then, go back to 5. above to run the test code

   6.2 Some compile details for DIY

   ```properties 
   > git clone https://dbin_sdm@bitbucket.org/dbin_sdm/dassa.git
   > module load gcc/7.4.0 hdf5/1.10.5-gcc-p fftw
   > edit Makefile to have correct setting AU_DIR/HDF5_DIR/FFTW_DIR/DASH_DIR
   > make
   ```


