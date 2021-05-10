%% Script to convert dsi format DAS data into hdf5 format to use as input for cross-correlation code 
%% Written by Bin Dong - Sent to V. Rodríguez Tribaldos on February 23rd, 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Input
% -- Dsi file after application of pre-processing for cross-correlation preparation (i.e. remove mean and trend, decimation, temporal normalization and spectral whitening)

%%% Output
% -- Cross-correlation on hdf5 format

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%%% History

% - Original version: 
%   Bin Dong (Computational Research Division, Data Science & Technology Department, LBNL) -- Sent to V. Rodríguez Tribaldos 

% - Modifications
% V. Rodríguez Tribaldos -- add some minor comments, paths to data files.
% 
% Add from_mem_flag  Apr 5  Bin Dong
%     



%%
% Convert dsi_file to HDF5 file
%  -- The "dsi_file" is the file name to a DSI file, which has a cell namely
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
%     The dsi2h5 will attach two metadata to record the semantic of the
%     data within the HDF5.
%     - nTrace is the number of channel, size(dsi_file.dat)[1] 
%     - nPoint is the number of points in time series, , size(dsi_file.dat)[1] 
function dsi2h5(dsi_file, h5_file, h5_dataset, h5_type, transpose_flag)
    if nargin < 4
      error('At least four parameters are needed to call dsi2h5, dsi_file, h5_file, h5_dataset and h5_type.');
    end
    
    if nargin < 5
        transpose_flag =  0;
    end
    dsi2h5_write(dsi_file, h5_file, h5_dataset, h5_type, transpose_flag, 0);
end











