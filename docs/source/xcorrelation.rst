.. _xcorrelation:


X-correlation
=============

X-correlation contains a series of ambient noise related data analysis operations to generate empirical Green’s functions used for imaging the shallow subsurface. The processing approach involves a long pipeline to convert the raw DAS data into shear-wave velocity profiles. It converts raw data into a noise- correlation in the frequency domain, after applying a series of pre-processing and filtering operations.

Tables of contents of this page

- :ref:`1, Quick Start Example`
- :ref:`2, Configuration File`
- :ref:`3, Documentation for Functions within the Xcorrelation`


1, Quick Start Example 
----------------------

Below is a quick start example to run `xcorrelation` with the existing configuration file ( xcorrelation.config ) in the repo. 
The xcorrelation accepts one parameter `-c` which comes with the configuration file name.

.. code-block:: bash

    > ./xcorrelation -c xcorrelation.config

    Parameters to run the Xcorrelation:
    Input parameters:
            input_dir_file = test-data/test-das.h5
            input_file_type = EP_HDF5
            input_data_type = short
            is_column_vector = true
            is_channel_range = false
            is_channel_stride = false
            n_files_to_concatenate = 1

    Runtime parameters:
            is_space_decimate = 0
            space_decimate_chs = 32
            space_decimate_operation = ave
            DT = 0.002
            DT_NEW = 0.008

    Output parameters:
            is_output_single_file = 1
            output_type = EP_HDF5
            output_file_dir = test-data/test-das-xcorr.h5
            output_dataset = /dat

    Rank 0, Xcorr: chunk_size = : 12,16
    In init_xcorr, nPoint = 3, lts_per_file =12, DT_NEW = 0.008, DT =  0.002 , nfft = 8, round_dt_dew_dt = 4
    Rank 0, BUTTER_A = : 0.0316893,0.095068,0.095068,0.0316893
    Rank 0, BUTTER_B = : 1,-1.45903,0.910369,-0.197825

    Rank 0,        data size: 12,16
    Rank 0,       chunk size: 12,16
    Rank 0,     overlap size: 0,0
    Rank 0, current chunk id:: : 0
    Rank 0,     total chunks: : 1
    Rank 0, max_offset_upper = : 11,15
    After time-domain decimate, with 16 channels,  each with 3 points
    Rank 0, Output vector_shape: : 5,16
    Timing Results of All
    Read      time (s) : max = 0.001005, min = 0.001005, ave = 0.001005
    UDF       time (s) : max = 0.003153, min = 0.003153, ave = 0.003153
    Write     time (s) : max = 0.002369, min = 0.002369, ave = 0.002369

2, Configuration File
---------------------

.. code-block:: bash


    #
    # Author Bin Dong Apr 7 2021
    # This file demonstrates the usage of config file
    #  It perform the xcorr on a single file "test-data/test-das.h5"
    #  It stores the result into another file "test-data/test-das-xcorr.h5"
    #
    # Note: 
    #  Each config file starts with the "[parameter]"
    #  The comment lines starting with "#" ";" and "%" will be skipped 
    #  The ";" can also be used to comment the line at the end
    #
    #  In the blow text,
    #  -- I will use the "#" to comment a group of lines.
    #     We may have three groups, input parameter, output parameter and runtime parameter.
    #      
    #  -- I will use the ";" to comment each line
    #
    # Also, all parameters have default values (even for the input file)
    #       So, some missing will not report the error. 
    #       Todo: we may classify the parameter as optional ones or essential ones
    #
    #

    [parameter]

    #############################
    #   Input data's parameters #
    #############################

    is_input_single_file = true   ; true / false (default),
    input_file_type = EP_HDF5     ; only takes  EP_HDF5(default) / EP_TMDS
    input_dir_file = test-data/test-das.h5 ; when is_input_single_file = true,  it points to a file
                                        ; when is_input_single_file = false, it points to a directory
    input_dataset=/dat    ; dataset name in HDF5, use (h5dump -A) to get it

    input_data_type = short ;short (default)/double/float the data element type in dataset,  

    is_column_vector = true ; true (default)/false
                            ; column vector: each column is a time series (most cases)
                            ; row vector: each row is a time series 

    n_files_to_concatenate = 1 ; the number of files to concatenate
                            ; only works when input is a directory
                            ; 1 : not to concatenate
                            ; 2 : concatenate every two files ...
                            ; Note: better not to have leftover 

    is_input_search_rgx = false  ; true / false(default) 
    input_search_rgx = (.*?)[1](\.h5) ; filter the input file names as input
                                    ; See : https://www.cplusplus.com/reference/regex/ECMAScript/


    is_channel_range = false   ;  true / false(default) 
    channel_range_start = 0    ;  Select a few channels to run xcorr
    channel_range_end = 2      ;  channel_range_start is "0" based. 

    is_channel_stride = false     ; true / false (default)
    channel_stride_size = 1       ; Only used when is_ch_stride = true
                        ; Pick every [ch_stride_size] channel from the first (zero based)
                        ; E.g.,  ch_stride_size = 99,
                        ; It picks channels 0, 99, 198, ....



    is_file_range = false      ; false or 0,  true/1 (by default 0) 
                            ; only works when  is_input_single_file = flase, i.e., a directory
                            ; pick the [file_range_start_index]th file to  [file_range_end_index]th file
                            ; All files are sorted by the filenames (kind of time order)
    file_range_start_index = 0 ; Note: zero based and inclusive
    file_range_end_index = 4

    #####################################
    #        Output data's parameters   #
    #####################################

    is_output_single_file = true                     ; true / false(default)
    output_type = EP_HDF5                            ; only takes EP_HDF5 now
    output_file_dir = test-data/test-das-xcorr.h5    ; when is_output_single_file = true,  it points to a file
                                                    ; when is_output_single_file = false,  it points to a directory
    output_dataset = /dat                            ; dataset for output file

    is_dir_output_match_replace_rgx = false          ; true / false(default), only works in directory mode
                                                    ; whether it has a way to auto generate the output file 
                                                    ; name from the input file name
    output_file_regex_match = ^(.*)\.h5$            ; regex pattern to match original file name
    output_file_regex_replace = $1-xcorr.h5          ; regex pattern to replace original file name



    ################################
    #      Runtime parameters      #
    ################################

    master_index   = 0     ; master channel index ("0" based)
                        ; if is_channel_range = true, "master_index = 0" means the "channel_range_start"
                        ; so, it is relative channel index

    butter_order = 3       ; The order for Butterworth filter
    dt = 0.002             ; original sampling interval (1/frequency) 
    dt_new = 0.008         ; new sampling interval (1/frequency) , for resample

    winLen_sec = 0.5              ;window length (in second) for the move mean operation
    z = 0, 0.5, 1.0, 1.0, 0.5, 0  ;for interp1
    F1 = 0
    F2 = 0.002
    F3 = 0.006
    F4 = 14.5
    F5 = 15
    eCoeff = 1                    ;for whitening 


    is_space_decimate = false       ; it may have space_decimate after (resample)
    space_decimate_chs = 32         ; the number of channels to decimate 
    space_decimate_operation = ave  ; ave(default)/median/min/max


    ##########################
    # Other Parameters       #
    ##########################

    #
    # These parameters are only needed when you want to specify the attribute names used for auto-layout detection
    # Another option is to set "is_column_vector = true/false", which will ignore the auto-layout detection
    #
    attribute_name_measure_length = MeasureLength[m]          ; the length of the fiber  
    attribute_name_spatial_resolution = SpatialResolution[m]  ; the resolution of the fiber 
    attribute_name_sampling_frequency = SamplingFrequency[Hz] ; the sampling frequency


- is_input_single_file
- input_file_type

3, Documentation for Functions within the Xcorrelation
------------------------------------------------------

:download:`Documentation <xcorrdoc.pdf>`