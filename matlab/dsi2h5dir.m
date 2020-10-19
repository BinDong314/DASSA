%
% Users intput
%

dsi_dir = "/Users/dbin/work/dassa/vincent-stack-debug/test-dsi/"
dsi_dir_filter_pattern="*.mat"

h5_dir = "/Users/dbin/work/dassa/vincent-stack-debug/test-h5/"
h5_dset_name='dat'
h5_data_type='double' % int32 int16, float/single, double
h5_file_append='.h5'  % h5_file_name is not needed, the append '.h5' to file name

transpose_flag=1 % from col-major (matlab) to row-major (hdf5)


%
% Some code to work
%
dsi_files = [dir(dsi_dir+dsi_dir_filter_pattern)];
dsi_file_list_path = struct2cell(dsi_files);
cols = size(dsi_file_list_path,2);

for i = 1:cols
    if (i<2)
        filelist = dsi_file_list_path{1,i};
    else
        filelist = [filelist; dsi_file_list_path{1,i}];
    end
end

filelist_cell = cellstr(filelist);
dsiPath = 'xcorr_examples/';
nStack = length(filelist_cell);


for k = 1 : nStack    
    dsi_file = char(filelist_cell{k});
    dsi_fullpath = fullfile(dsi_dir, dsi_file);
    h5_fullpath=fullfile(h5_dir, dsi_file) + h5_file_append;
    dsi2h5(dsi_fullpath, h5_fullpath, h5_dset_name, h5_data_type, transpose_flag)
end    






