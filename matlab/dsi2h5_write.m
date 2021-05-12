%
% A common function used by dsi2h5.m and dsi2h5_from_mem.m
%
function dsi2h5_write(dsi_file, h5_file, h5_dataset, h5_type_str, transpose_flag, from_mem_flag)

if(from_mem_flag == 0)
    dsi_data_raw = importdata(dsi_file);
else
    dsi_data_raw =   dsi_file;  
end

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
    error('Not known h5_type_str, I only understand single, int16, int32, float, double.')
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
    disp('Not known type_str, I only understand single(float), double, int16 and int32.')
end

dset_id = H5D.create(fid, h5_dataset, type_id, space_id, 'H5P_DEFAULT');
H5D.write(dset_id, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',transpose(dsi_data));


acpl = H5P.create('H5P_ATTRIBUTE_CREATE');
%attr_type_id = H5T.copy('H5T_NATIVE_DOUBLE');
attr_space_id = H5S.create('H5S_SCALAR');

%nTrace
% is the number of channels
%It is stored as string in HDF5
%  e.g.,
%    ATTRIBUTE "nTrace" {
%          DATATYPE  H5T_STRING {
%             STRSIZE 3;
%             STRPAD H5T_STR_NULLTERM;
%             CSET H5T_CSET_ASCII;
%             CTYPE H5T_C_S1;
%          }
%          DATASPACE  SCALAR
%          DATA {
%          (0): "375"
%          }
%       }
attr_type_id = H5T.copy('H5T_C_S1');
nTrace_str=int2str(h5_dims(2));
H5T.set_size(attr_type_id,strlength(nTrace_str));
encoding = H5ML.get_constant_value('H5T_CSET_ASCII');
H5T.set_cset(attr_type_id, encoding);
attr_id = H5A.create(dset_id,'nTrace', attr_type_id, attr_space_id, acpl);
H5A.write(attr_id,'H5ML_DEFAULT',nTrace_str)
H5A.close(attr_id);
H5T.close(attr_type_id);

% nPoint
%  is the number of point in time series 
% e.g.
% ATTRIBUTE "nPoint" {
%          DATATYPE  H5T_STRING {
%             STRSIZE 6;
%             STRPAD H5T_STR_NULLTERM;
%             CSET H5T_CSET_ASCII;
%             CTYPE H5T_C_S1;
%          }
%          DATASPACE  SCALAR
%          DATA {
%          (0): "149999"
%          }
%       }
nPoint_str=int2str(h5_dims(1));
attr_type_id_2 = H5T.copy('H5T_C_S1');
H5T.set_size(attr_type_id_2,strlength(nPoint_str));
H5T.set_cset(attr_type_id_2, encoding);
attr_id_2 = H5A.create(dset_id,'nPoint', attr_type_id_2, attr_space_id, acpl);
H5A.write(attr_id_2,'H5ML_DEFAULT',nPoint_str)
H5A.close(attr_id_2);
H5T.close(attr_type_id_2);


% SamplingFrequency[Hz]
%  is the SamplingFrequency[Hz] of the orginal data
% Saved as string in HDF5
% e.g.
%   ATTRIBUTE "SamplingFrequency[Hz]" {
%          DATATYPE  H5T_STRING {
%             STRSIZE 3;
%             STRPAD H5T_STR_NULLTERM;
%             CSET H5T_CSET_ASCII;
%             CTYPE H5T_C_S1;
%          }
%          DATASPACE  SCALAR
%          DATA {
%          (0): "125"
%          }
%       }
%if (isfield(dsi_data_raw, 'fh')) 
%    SamplingFrequency_int = 1/dsi_data_raw.fh{8};
    %
    %Fixme: we conver HZ to be normalized ones (based on 1 minute's data)
    %       This convert only happens when the data spans multiple
    %       munites
    %
    % We do this to let DASSA understand data layout
    %
%    SamplingFrequency_int =  SamplingFrequency_int * h5_dims(1) / (500 * 60)
%    SamplingFrequency_str=int2str(SamplingFrequency_int)
%     attr_type_id_22 = H5T.copy('H5T_C_S1');
%     H5T.set_size(attr_type_id_22,strlength(SamplingFrequency_str));
%     H5T.set_cset(attr_type_id_22, encoding);
%     attr_id_22 = H5A.create(dset_id,'SamplingFrequency[Hz]', attr_type_id_22, attr_space_id, acpl);
%     H5A.write(attr_id_22,'H5ML_DEFAULT',SamplingFrequency_str)
%     H5A.close(attr_id_22);
%     H5T.close(attr_type_id_22);
% end

%"MeasureLength[m]"
% In this case, we set the MeasureLength[m] = fh{1}
% The MeasureLength[m] is not the same as the MeasureLength[m] in orginal HDF5 or
% TDMS file from iDAS devices
% We set MeasureLength[m] = fh{1} to let
% MeasureLength[m]/SpatialResolution[m] = nTrace
% Both MeasureLength[m] and SpatialResolution[m] are used by DASSA's
% XCORR'code to figure out column-vector or so.
% if (isfield(dsi_data_raw, 'fh'))
%     %
%     %Fixme: we treate the nTrace as the MeasureLength to let DASSA understand data layout
%     %       
%     %
%     MeasureLength_str=int2str(dsi_data_raw.fh{1});
%     attr_type_id_23 = H5T.copy('H5T_C_S1');
%     H5T.set_size(attr_type_id_23,strlength(MeasureLength_str));
%     H5T.set_cset(attr_type_id_23, encoding);
%     attr_id_23 = H5A.create(dset_id,'MeasureLength[m]', attr_type_id_23, attr_space_id, acpl);
%     H5A.write(attr_id_23,'H5ML_DEFAULT',MeasureLength_str)
%     H5A.close(attr_id_23);
%     H5T.close(attr_type_id_23);
% end

%SpatialResolution[m]
% In this case, we set the SpatialResolution[m] = 1
% This is not the same number as the SpatialResolution[m] in orginal HDF5 or
% TDMS file from iDAS devices
% We set SpatialResolution[m] = 1 to let
% MeasureLength[m]/SpatialResolution[m] = nTrace
% Both MeasureLength[m] and SpatialResolution[m] are used by DASSA's
% XCORR'code to figure out column-vector or so.
% if (isfield(dsi_data_raw, 'fh'))  
%     %
%     %Fixme: we only put '1' here to let DASSA understand data layout
%     %       
%     %
%     SpatialResolution_str='1';
%     attr_type_id_24 = H5T.copy('H5T_C_S1');
%     H5T.set_size(attr_type_id_24,strlength(SpatialResolution_str));
%     H5T.set_cset(attr_type_id_24, encoding);
%     attr_id_24 = H5A.create(dset_id,'SpatialResolution[m]', attr_type_id_24, attr_space_id, acpl);
%     H5A.write(attr_id_24,'H5ML_DEFAULT',SpatialResolution_str)
%     H5A.close(attr_id_24);
%     H5T.close(attr_type_id_24);
% end

%
%fh header, is saved as attribute 
% we convert fh to be string, which ignore the empy ones, fh{2, 3, 4, 11}
%  fh{5} and fh{6} are flatted
% e.g.,
%  fh = {375,[],[],[],[2020,11,11,10,9,0.580000000000000],[2020,11,11,10,10,0.580000000000000],300000,0.00200000000000000,0,59.9980000000000,[],1,1}
%  It is stored as
%    ATTRIBUTE "fh" {
%          DATATYPE  H5T_STRING {
%             STRSIZE 80;
%             STRPAD H5T_STR_NULLTERM;
%             CSET H5T_CSET_ASCII;
%             CTYPE H5T_C_S1;
%          }
%          DATASPACE  SCALAR
%          DATA {
%          (0): "375 2020 11 11 10 9 0.58 2020 11 11 10 10 0.58 149999 0.008 -599.992 599.992 1 1"
%          }
%       }
if (isfield(dsi_data_raw, 'fh'))
   fh_str = convertStringsToChars(join(string(cell2mat(dsi_data_raw.fh))));
   attr_type_id_3 = H5T.copy('H5T_C_S1');
   H5T.set_size(attr_type_id_3, strlength(fh_str));
   H5T.set_cset(attr_type_id_3, encoding);
   attr_id_3 = H5A.create(dset_id,'fh', attr_type_id_3, attr_space_id, acpl);
   H5A.write(attr_id_3,'H5ML_DEFAULT',fh_str)
   H5A.close(attr_id_3);
   H5T.close(attr_type_id_3);
end

% th header is stored as a 2D data set at it is too large for attribute
%  see blow code for error information
if (isfield(dsi_data_raw, 'th'))
    th_mat = cell2mat(dsi_data_raw.th);
    th_dims = size(th_mat);
    th_space_id = H5S.create_simple(2, th_dims, th_dims);
    th_type_id = H5T.copy('H5T_NATIVE_DOUBLE');
    th_order = H5ML.get_constant_value('H5T_ORDER_BE');
    H5T.set_order(th_type_id,th_order);
    th_dset_id = H5D.create(fid, 'th', th_type_id, th_space_id, 'H5P_DEFAULT');
    H5D.write(th_dset_id, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',transpose(th_mat));
    H5S.close(th_space_id);
    H5T.close(th_type_id);
    H5D.close(th_dset_id);
end
%th
% if (isfield(dsi_data_raw, 'th'))
%    th_mat = cell2mat(dsi_data_raw.th);
%    th_vec = th_mat(:);
%    th_str = convertStringsToChars(join(string(th_vec)));
%    attr_type_id_4 = H5T.copy('H5T_C_S1');   
%    H5T.set_size(attr_type_id_4,strlength(th_str));
%    H5T.set_cset(attr_type_id_4, encoding);
%    attr_id_4 = H5A.create(dset_id,'th', attr_type_id_4, attr_space_id, acpl);
%    H5A.write(attr_id_4,'H5ML_DEFAULT',th_str);
%    H5A.close(attr_id_4);
%    H5T.close(attr_type_id_4);
% end
%Error using hdf5lib2
%The HDF5 library encountered an error and produced the following stack trace information:
%    H5O_alloc              object header message is too large
%    H5O_msg_alloc          unable to allocate space for message
%    H5O_msg_append_real    unable to create new message
%    H5O_attr_create        unable to create new attribute in header
%    H5A_create             unable to create attribute in object header
%    H5Acreate1             unable to create attribute
    

H5S.close(attr_space_id);

%Write data
%Here we transpose the column-major to row-major
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);
H5F.close(fid);
end

