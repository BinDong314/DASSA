
dsi_file='/Users/dbin/work/dassa/vincent-stack-debug/test-dsi/1minXcorr_180104000031_ch4650_4850.mat';
h5_file_name='/Users/dbin/work/dassa/vincent-stack-debug/test-h5/1minXcorr_180104000031_ch4650_4850.h5';
h5_dset_name='dat'
h5_data_type='double' % int32 int16, float/single, double
transpose_flag=1 % from col-major (matlab) to row-major (hdf5)
from_mem_flag = 0 %

%importdata(test_dsi_file);
dsi2h5(dsi_file, h5_file_name , h5_dset_name, h5_data_type, transpose_flag, from_mem_flag);

