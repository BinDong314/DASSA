%
% A test function to read the data saved with dsi2h5_write
% Author: Bin Dong May 2021

% Simple test code 1
% h5_file = '/Users/dbin/work/dassa/veronica-10-files/data/vero-10-file-merge-issue/mat/IV_10min_20201111100058_20201111100958_ch1700_3200.mat.h5';
% h5_dataset = '/dat';
% dsi = dsi2h5_read(h5_file, h5_dataset);

% Simple test code 2
h5_file = '/Users/dbin/work/dassa/veronica-10-files/data/vero-10-file-merge-issue/mat-dsi2h5-xcorr/IV_10min_20201111100058_20201111100958_ch1700_3200.mat.h5';
h5_dataset = '/dat';
dsi_from_h5 = dsi2h5_read(h5_file, h5_dataset);

%
% The function dsi2h5_read accepts two parameters,  h5_file and h5_dataset
% 
function [Dsi_placeholder]=dsi2h5_read(h5_file, h5_dataset)
    if nargin < 2
      error('At least two parameters are needed to call h5_file and h5_dataset.');
    end
    
    %Get the data
    dsi_dat = h5read(h5_file,h5_dataset)';
    Dsi_placeholder.dat{1} = dsi_dat;
%     size(dsi_dat)
%     nTrace = size(dsi_dat, 2);
%     nPoint = size(dsi_dat, 1);
    
    %Get the fh
%     dsi_fh_str = h5readatt(h5_file,h5_dataset, 'fh')
%     dsi_fh_vec = str2double(split(dsi_fh_str));
%  
%     assert(dsi_fh_vec(1) == nTrace);
%     assert(dsi_fh_vec(14) == nPoint);
% 
%     Dsi_placeholder.fh = cell(1, 13);
%     Dsi_placeholder.fh{1} = nTrace;
%     Dsi_placeholder.fh{7} = nPoint;
%     Dsi_placeholder.fh{12} = 1;
%     Dsi_placeholder.fh{13} = 1;
%     
%     tb_datevec = dsi_fh_vec(2:7)';
%     te_datevec = dsi_fh_vec(8:13)';
%     Dsi_placeholder.fh{5} = tb_datevec;
%     Dsi_placeholder.fh{6} = te_datevec;
%     
%     Dsi_placeholder.fh{8} = dsi_fh_vec(15);
% 
%     Dsi_placeholder.fh{10} = dsi_fh_vec(16);     % fh Word #10 = Trace end time (s)
%     Dsi_placeholder.fh{11} = dsi_fh_vec(17); % fh Word #11 = P value in tdms header

    %Get the th
%     dsi_th_dat = h5read(h5_file, '/th');
%     Dsi_placeholder.th{1} = dsi_th_dat;

 
end

