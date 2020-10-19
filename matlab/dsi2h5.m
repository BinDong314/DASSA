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


%%
% Convert dsi_file to HDF5 file
% h5_type_str: only allow 'single', 'int16', 'int32' 'float', and 'double'

function dsi2h5(dsi_file, h5_file, h5_dataset, h5_type_str, transpose_flag)
dsi_data_raw = importdata(dsi_file);

fcpl = H5P.create('H5P_FILE_CREATE');
fapl = H5P.create('H5P_FILE_ACCESS');
fid  = H5F.create(h5_file, 'H5F_ACC_TRUNC', fcpl, fapl);

%Get the data and its size, type, etc.
%Force to use the  single/float/double/int16 type 
if(strcmp(h5_type_str, 'single'))
    dsi_data =single(cell2mat(dsi_data_raw.dat));
elseif(strcmp(h5_type_str, 'float'))
    dsi_data =float(cell2mat(dsi_data_raw.dat));
elseif(strcmp(h5_type_str, 'double'))
    dsi_data = double(cell2mat(dsi_data_raw.dat));
elseif(strcmp(h5_type_str, 'int16'))
    dsi_data =int16(cell2mat(dsi_data_raw.dat));
elseif(strcmp(h5_type_str, 'int32'))
    dsi_data =int32(cell2mat(dsi_data_raw.dat));
else
    disp('Not known h5_type_str, I only understand single, int16, int32, float, double.')
end

if(transpose_flag == 1)
    dsi_data = transpose(dsi_data);
end

h5_dims = size(dsi_data);
 
%Create the space size on disk /w size of float_data_set
space_id = H5S.create_simple(2, h5_dims, h5_dims);


%Create the actual dataset to store the 2D data
%Get the data and its size, type, etc.
if(strcmp(h5_type_str, 'single'))
    type_id = H5T.copy('H5T_NATIVE_FLOAT');
    order = H5ML.get_constant_value('H5T_ORDER_BE');
    H5T.set_order(type_id,order);
elseif(strcmp(h5_type_str, 'int16'))
    type_id = H5T.copy('H5T_NATIVE_SHORT');
    order = H5ML.get_constant_value('H5T_ORDER_BE');
    H5T.set_order(type_id,order);
elseif(strcmp(h5_type_str, 'int32'))
    type_id = H5T.copy('H5T_NATIVE_INT');
    order = H5ML.get_constant_value('H5T_ORDER_BE');
    H5T.set_order(type_id,order);
elseif(strcmp(h5_type_str, 'float'))
    type_id = H5T.copy('H5T_NATIVE_FLOAT');
    order = H5ML.get_constant_value('H5T_ORDER_BE');
    H5T.set_order(type_id,order);
elseif(strcmp(h5_type_str, 'double'))
    type_id = H5T.copy('H5T_NATIVE_DOUBLE');
    order = H5ML.get_constant_value('H5T_ORDER_BE');
    H5T.set_order(type_id,order);
else
    disp('Not known type_str, I only understand single(float), double or int16.')
end

dset_id = H5D.create(fid, h5_dataset, type_id, space_id, 'H5P_DEFAULT');

%Write data
%Here we transpose the column-major to row-major
H5D.write(dset_id, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',transpose(dsi_data));
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);
H5F.close(fid);
end










