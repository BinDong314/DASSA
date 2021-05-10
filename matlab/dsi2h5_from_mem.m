
%%
% Convert dsi_variable to HDF5 file
%  -- The "dsi_variable" is the matlabe variable, which has a cell namely
%    'dat' for the data.
%  -- The data in "dsi_file.dat" is stored as "double"
%  -- The "h5_type" specifies the type to store the data in HDF5.
%     Thus, there is a conversion from "double" to the "h5_type".
%     The "h5_type" accepts 'single', 'int16', 'int32' 'float', and 'double'
%  -- The "h5_file" and "h5_dataset" contains the HDF5 file name and dataset
%     name for HDF5 file.
%  
%  ! Above four parameter are necessay to call dsi2h5
%
%  The "transpose_flag" is optional parameter, by default it is set to be 0
%  "transpose_flag=0" it stores the data in dsi as it is in dsi_file.   
%
%  About metadata:
%     The dsi2h5 will attach below metadata to record the semantic of the
%     data within the HDF5.
%     - SamplingFrequency[Hz]
%     - MeasureLength[m] and SpatialResolution[m], 
%     - nTrace is the number of channel, size(dsi_file.dat)[1] 
%     - nPoint is the number of points in time series, , size(dsi_file.dat)[1] 
%     - fh (as string)
%     - dh (as dataset)

function dsi2h5_from_mem(dsi_variable, h5_file, h5_dataset, h5_type, transpose_flag)
    if nargin < 4
      error('At least four parameters are needed to call dsi2h5_from_mem, dsi_file, h5_file, h5_dataset and h5_type.');
    end
    
    if nargin < 5
        transpose_flag =  0;
    end
   dsi2h5_write(dsi_variable, h5_file, h5_dataset, h5_type, transpose_flag, 1);
end
