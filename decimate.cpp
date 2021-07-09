/**
 *
 * Author: Bin Dong
 * Email questions to dbin@lbl.gov
 * Scientific Data Management Research Group
 * Lawrence Berkeley National Laboratory
 *
 */

//
// This  example code is for decimating DAS data
//     detread
//     filtfilt
//     resample
//

#include <iostream>
#include <stdarg.h>
#include <vector>
#include <stdlib.h>

#include "ft.h"
#include "DasLib.h"

using namespace std;
using namespace AU;
using namespace DasLib;

int read_config_file(std::string file_name, int mpi_rank);
std::string config_file = "./decimate.config";

bool is_input_single_file = false;
std::string input_dir_file = "/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir";
std::string input_file_type = "EP_TDMS";
bool is_input_search_rgx = false;
std::string input_search_rgx = "^(.*)[1234]\\.tdms$";

int chs_per_file = 11648;
int lts_per_file = 30000;
int n_files_to_concatenate = 2;

bool is_output_single_file = false;
std::string output_type = "EP_HDF5";
std::string output_file_dir = "./tdms-dir-dec/";
std::string output_dataset = "/dat";

bool is_dir_output_match_replace_rgx = true;
std::string dir_output_match_rgx = "^(.*)\\.tdms$";
std::string dir_output_replace_rgx = "$1.h5";

bool is_space_decimate = false;
int space_decimate_rows = 32;
/**
 * 
 * Name convertion for the space_decimate_operation
*     ave, average, mean
*     med, median
*     max, maximum
*     min, minimum
*     sum, summary
 */
std::string space_decimate_operation = "ave"; //"median", "min", "max"

double DT = 0.002;
double DT_NEW = 0.008;

void printf_help(char *cmd);

//
vector<double> BUTTER_A;
vector<double> BUTTER_B;

//int nfft;
//FIND_M_POWER2(nPoint, nfft);

int butter_order = 3;
double cut_frequency_low = 0.25;

bool is_channel_range = false;
int channel_range_start = 0;
int channel_range_end = 1;
bool is_many_files = false;
int many_files_split_n = 10;

bool is_column_major = true;
bool is_column_major_from_config = false;

std::string MeasureLengthName = "MeasureLength[m]";
std::string SpatialResolutionName = "SpatialResolution[m]";
std::string SamplingFrequencyName = "SamplingFrequency[Hz]";

bool is_channel_stride = false;
int channel_stride_size = 0;

bool is_stencil_tag_once = false;

void InitDecimate()
{
    //int nPoint = ceil(lts_per_file * n_files_to_concatenate / (DT_NEW / DT));
    int nPoint = ceil(lts_per_file / (DT_NEW / DT));
    cut_frequency_low = (0.5 / DT_NEW) / (0.5 / DT);
    ButterLow(butter_order, cut_frequency_low, BUTTER_A, BUTTER_B);
    if (!ft_rank)
        std::cout << "After decimate, nPoint = " << nPoint << "\n";
}

inline Stencil<std::vector<double>> udf_decimate(const Stencil<short> &iStencil)
{
    std::vector<int> max_offset_upper;
    iStencil.GetOffsetUpper(max_offset_upper);
    PrintVector("max_offset_upper = ", max_offset_upper);

    //int chs_per_file_udf = max_offset_upper[0] + 1, lts_per_file_udf = max_offset_upper[1] + 1;
    int chs_per_file_udf, lts_per_file_udf;
    std::vector<int> start_offset = {0, 0};
    std::vector<int> end_offset = {max_offset_upper[0], max_offset_upper[1]};

    if (is_column_major)
    {
        chs_per_file_udf = max_offset_upper[1] + 1;
        lts_per_file_udf = max_offset_upper[0] + 1;
    }
    else
    {
        chs_per_file_udf = max_offset_upper[0] + 1;
        lts_per_file_udf = max_offset_upper[1] + 1;
    }

    if (is_channel_range)
    {
        assert(channel_range_end - channel_range_start + 1 <= chs_per_file_udf);
        chs_per_file_udf = channel_range_end - channel_range_start + 1;
        if (is_column_major)
        {
            start_offset[1] = channel_range_start;
            end_offset[1] = channel_range_end;
        }
        else
        {
            start_offset[0] = channel_range_start;
            end_offset[0] = channel_range_end;
        }
        PrintVector("start_offset = ", start_offset);
        PrintVector("end_offset = ", end_offset);
    }

    // //Todo: not done
    // if (is_many_files)
    // {
    //     if (is_space_decimate)
    //     {
    //         many_files_split_n = space_decimate_rows;
    //     }
    //     else
    //     {
    //         many_files_split_n = chs_per_file_udf / n_files_to_concatenate;
    //     }

    //     if (!ft_rank)
    //         std::cout << "Using the is_many_files, many_files_split_n = " << many_files_split_n << " \n";
    // }

    std::vector<short> ts_short;
    iStencil.ReadNeighbors(start_offset, end_offset, ts_short);

    //Convert to row-vector here if it is column-vector
    //Because all the following code are built as row-vector (2D vector)
    //Each row is a time series
    if (is_column_major)
    {
        std::vector<short> ts_short_temp;
        ts_short_temp.resize(ts_short.size());
        transpose(ts_short.data(), ts_short_temp.data(), lts_per_file_udf, chs_per_file_udf);
        ts_short = ts_short_temp;
    }

    std::vector<std::vector<double>> ts2d = DasLib::Vector1D2D<short, double>(lts_per_file_udf, ts_short);

    //We may update the data based on is_channel_stride/channel_stride_size
    //e.g.,  channel_range_start = 0, channel_range_end = 99
    //       channel_stride_size = 99
    if (is_channel_stride)
    {
        std::vector<std::vector<double>> ts2d_temp;
        for (int iiii = 0; iiii < chs_per_file_udf; iiii++)
        {
            if (iiii % channel_stride_size == 0)
            {
                std::cout << iiii << "\n";
                //ts2d.erase(ts2d.begin() + iiii);
                ts2d_temp.push_back(ts2d[iiii]);
            }
        }
        ts2d = ts2d_temp;
        chs_per_file_udf = ts2d.size();
        std::cout << "chs_per_file_udf (after skip ) = " << ts2d.size() << ", " << ts2d[0].size() << "\n";
    }

    std::vector<std::vector<double>> ts2d_ma;

    std::cout << "ts2d.size() = " << ts2d.size() << ",ts2d[0].size() = " << ts2d[0].size() << ", lts_per_file_udf =" << lts_per_file_udf << ", ts_short.size() = " << ts_short.size() << "\n";

    //std::cout << "Got data ! at rank " << ft_rank << " \n";

    std::vector<double> ts_temp2;
    //Resample in time-domain
    for (int i = 0; i < chs_per_file_udf; i++)
    {
        detrend(ts2d[i].data(), lts_per_file_udf);       //Detread
        filtfilt(BUTTER_A, BUTTER_B, ts2d[i], ts_temp2); //filtfilt
        resample(1, DT_NEW / DT, ts_temp2, ts2d[i]);     //resample
    }
    DasLib::clear_vector(ts_temp2);
    if (!ft_rank)
        std::cout << "Finish time-domain decimate ! \n";

    if (is_space_decimate)
    {
        //decimate in space-domain
        int ma_batches = (chs_per_file_udf % space_decimate_rows) ? (chs_per_file_udf / space_decimate_rows + 1) : (chs_per_file_udf / space_decimate_rows);
        int start_row, end_row;
        for (int i = 0; i < ma_batches; i++)
        {
            start_row = i * space_decimate_rows;
            end_row = (((i + 1) * space_decimate_rows - 1) < chs_per_file_udf) ? ((i + 1) * space_decimate_rows - 1) : chs_per_file_udf - 1;
            //std::cout << "ma_batches =" << ma_batches << ", start_row = " << start_row << ", end_row =  " << end_row << "\n";
            ts2d_ma.push_back(spacedecimate(ts2d, start_row, end_row, space_decimate_operation));
        }
    }
    else
    {
        ts2d_ma = ts2d;
    }
    DasLib::clear_vector(ts2d);
    if (!ft_rank)
        std::cout << "Finish space-domain decimate ! \n";

    std::vector<double> ts_temp = Convert2DVTo1DV(ts2d_ma);
    Stencil<std::vector<double>> oStencil;
    std::vector<size_t> vector_shape(2);

    if (is_column_major)
    {
        std::vector<double> ts_temp_column;
        ts_temp_column.resize(ts_temp.size());
        transpose(ts_temp.data(), ts_temp_column.data(), ts2d_ma.size(), ts2d_ma[0].size());
        ts_temp = ts_temp_column;
        vector_shape[1] = ts2d_ma.size();
        vector_shape[0] = ts2d_ma[0].size();
    }
    else
    {
        vector_shape[0] = ts2d_ma.size();
        vector_shape[1] = ts2d_ma[0].size();
    }

    //
    // Deal with tag
    //
    if (iStencil.HasTagMap() && !is_stencil_tag_once)
    {
        std::map<std::string, std::string> tag_map;
        iStencil.GetTagMap(tag_map);
        for (std::map<std::string, std::string>::iterator it = tag_map.begin(); it != tag_map.end(); ++it)
        {
            //std::cout << " key : " << it->first << ", value:" << it->second << " \n";
            if (it->first == "SamplingFrequency[Hz]")
            {
                int temp_sf = 1 / DT_NEW;
                it->second = std::to_string(temp_sf); // "125";
            }
        }
        oStencil.SetTagMap(tag_map);
        if (is_output_single_file) //We only deal with meta
            is_stencil_tag_once = true;
    }
    //vector_shape[0] = ts2d_ma.size();
    //vector_shape[1] = ts2d_ma[0].size();
    PrintVector("vector_shape: ", vector_shape);
    std::cout << "vector_shape[0] = " << vector_shape[0] << ",vector_shape[1] = " << vector_shape[1] << "\n";
    DasLib::clear_vector(ts2d_ma);
    oStencil.SetShape(vector_shape);
    oStencil = ts_temp;
    return oStencil;
}

int main(int argc, char *argv[])
{
    int copt;
    bool has_config_file_flag = false;
    while ((copt = getopt(argc, argv, "c:h")) != -1)
        switch (copt)
        {
        case 'c':
            config_file.assign(optarg);
            has_config_file_flag = true;
            break;
        case 'h':
            printf_help(argv[0]);
            exit(0);
        default:
            printf("Wrong option [%c] for %s \n", copt, argv[0]);
            printf_help(argv[0]);
            exit(-1);
            break;
        }

    //Init the MPICH, etc.
    AU_Init(argc, argv);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int numprocs, rank, namelen;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);
    printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

    if (has_config_file_flag)
        read_config_file(config_file, ft_rank);

    // set up the chunk size and the overlap size
    // 11648, 30000 for each dataset
    std::vector<int> chunk_size(2);
    std::vector<int> overlap_size = {0, 0};

    //Input data
    AU::Array<short> *A;
    if (!is_input_single_file)
    {
        A = new AU::Array<short>("EP_DIR:" + input_file_type + ":" + input_dir_file);
        std::cout << "EP_DIR:" + input_file_type + ":" + input_dir_file << "\n";

        std::vector<std::string> file_size_str;
        A->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str);
        String2Vector(file_size_str[0], chunk_size);
    }
    else
    {
        A = new AU::Array<short>(input_file_type + ":" + input_dir_file);
        std::cout << input_file_type + ":" + input_dir_file << "\n";
        std::vector<unsigned long long> array_size;
        A->GetArraySize(array_size);
        chunk_size[0] = array_size[0];
        chunk_size[1] = array_size[1];
    }
    PrintVector("chunk_size = ", chunk_size);

    //By default, the TDMS file is column major
    if (input_file_type == "EP_TDMS")
    {
        is_column_major_from_config = true; //Set to skip the check from file
        is_column_major = true;
    }

    if (!is_column_major_from_config)
    {
        /*
         * Here we try to read  MeasureLengthName/SpatialResolutionName/SamplingFrequencyName or nTrace/nPoint
         *    to detect the layout of the file.
         * 
         */
        int meta_chs = -1, meta_time_series_points = -1;
        std::string MeasureLength, SpatialResolution, SamplingFrequency;
        A->GetTag(MeasureLengthName, MeasureLength);
        A->GetTag(SpatialResolutionName, SpatialResolution);
        if (MeasureLength.empty() || SpatialResolution.empty())
        {
            A->GetTag("nTrace", MeasureLength);
            if (MeasureLength.empty())
            {
                std::cout << "Error: Can not find out [" << MeasureLengthName << "/" << SpatialResolutionName << "] or [nTrace] in the file to detect organization. You can also set [is_column_vector = true/false] in config file to disable this detection ! \n";
                exit(-1);
            }
            meta_chs = std::stoi(MeasureLength);
        }
        else
        {
            meta_chs = std::stoi(MeasureLength) / std::stoi(SpatialResolution);
        }

        A->GetTag(SamplingFrequencyName, SamplingFrequency);
        if (SamplingFrequency.empty())
        {
            A->GetTag("nPoint", SamplingFrequency);
            if (SamplingFrequency.empty())
            {
                std::cout << "Error: Can not find out [" << SamplingFrequencyName << "] or [nPoint] in the file to detect organization. You can also set [is_column_vector = true/false] in config file to disable this detection ! \n";
                exit(-1);
            }
            meta_time_series_points = std::stoi(SamplingFrequency);
        }
        else
        {
            meta_time_series_points = 60 * std::stoi(SamplingFrequency);
        }

        //std::cout << "meta_time_series_points = " << meta_time_series_points << " , meta_chs =  " << meta_chs << " \n";
        if (chunk_size[0] == meta_time_series_points && chunk_size[1] == meta_chs)
        {
            is_column_major = true;
            if (!ft_rank)
            {
                std::cout << termcolor::reset << "\n";
                std::cout << termcolor::magenta << "Found data organization = " << termcolor::green << " column vector (time by channel) \n";
                std::cout << termcolor::reset << "\n";
            }
        }
        else if (chunk_size[0] == meta_chs && chunk_size[1] == meta_time_series_points)
        {
            is_column_major = false;
            if (!ft_rank)
            {
                std::cout << termcolor::reset << "\n";
                std::cout << termcolor::magenta << "Found data organization = " << termcolor::green << " row vector (channel by time)\n";
                std::cout << termcolor::reset << "\n";
            }
        }
        else
        {
            std::cout << "Metadata and data are inconsistent in the size ! " << std::endl;
            exit(-1);
        }
    }

    if (is_column_major)
    {
        chunk_size[0] = chunk_size[0] * n_files_to_concatenate;
        lts_per_file = chunk_size[0];
        chs_per_file = chunk_size[1];
    }
    else
    {
        chunk_size[1] = chunk_size[1] * n_files_to_concatenate;
        lts_per_file = chunk_size[1];
        chs_per_file = chunk_size[0];
    }
    if (!ft_rank)
        std::cout << "lts_per_file = " << lts_per_file << ",chs_per_file = " << chs_per_file << "\n";

    ///Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir", chunk_size, overlap_size);
    A->SetChunkSize(chunk_size);
    A->SetOverlapSize(overlap_size);

    //Extract tag to be along with Input Stencil
    A->GetStencilTag();

    std::vector<std::string> aug_merge_index, aug_dir_sub_cmd, aug_input_search_rgx;

    //Set the index to merge file
    if (is_column_major)
    {
        aug_merge_index.push_back("0");
        A->EndpointControl(DIR_MERGE_INDEX, aug_merge_index);
    }
    else
    {
        aug_merge_index.push_back("1");
        A->EndpointControl(DIR_MERGE_INDEX, aug_merge_index);
    }

    //Set fhe search reges on file
    aug_input_search_rgx.push_back(input_search_rgx);
    if (is_input_search_rgx)
        A->EndpointControl(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx);

    if (!is_column_major && input_file_type == "EP_TDMS")
    {
        aug_dir_sub_cmd.push_back("BINARY_ENABLE_TRANSPOSE_ON_READ");
        A->EndpointControl(DIR_SUB_CMD_ARG, aug_dir_sub_cmd);
    }

    InitDecimate();

    //Result data
    AU::Array<double> *B;
    if (is_output_single_file)
    {
        //Store into a single file
        B = new AU::Array<double>(output_type + ":" + output_file_dir + ":" + output_dataset);
    }
    else
    {
        //Store into multiple file
        B = new AU::Array<double>("EP_DIR:" + output_type + ":" + output_file_dir + ":" + output_dataset);
        //Use the below rgx pattern to name the file
        std::vector<std::string> aug_output_replace_arg;
        aug_output_replace_arg.push_back(dir_output_match_rgx);
        aug_output_replace_arg.push_back(dir_output_replace_rgx);
        B->EndpointControl(DIR_MERGE_INDEX, aug_merge_index);

        if (is_dir_output_match_replace_rgx)
            B->EndpointControl(DIR_OUPUT_REPLACE_RGX, aug_output_replace_arg);
    }

    //Stride on execution
    //Each chunk only runs the udf_decimate once
    A->EnableApplyStride(chunk_size);

    A->SetVectorDirection(AU_FLAT_OUTPUT_ROW);

    //Run
    A->Transform(udf_decimate, B);

    A->ReportTime();
    //Clear
    delete A;
    delete B;

    AU_Finalize();

    return 0;
}

void printf_help(char *cmd)
{
    char *msg = (char *)"Usage: %s [OPTION]\n\
      	  -h help (--help)\n\
         -c config file for parameters  \n\
          Example: mpirun -n 1 %s -c  decimate.config\n";

    fprintf(stdout, msg, cmd, cmd);
}

int read_config_file(std::string file_name, int mpi_rank)
{
    INIReader reader(file_name);

    if (reader.ParseError() < 0)
    {
        std::cout << "Can't load [" << file_name << " ]\n";
        return 1;
    }

    input_dir_file = reader.Get("parameter", "input_dir_file", "/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir");

    input_file_type = reader.Get("parameter", "input_file_type", "EP_TDMS");

    std::string temp_str = reader.Get("parameter", "is_input_search_rgx", "false");

    is_input_search_rgx = (temp_str == "false") ? false : true;

    if (is_input_search_rgx)
    {
        input_search_rgx = reader.Get("parameter", "input_search_rgx", "^(.*)[1234]\\.tdms$");
    }

    chs_per_file = reader.GetInteger("parameter", "chs_per_file", 11648);
    lts_per_file = reader.GetInteger("parameter", "lts_per_file", 30000);

    temp_str = reader.Get("parameter", "is_channel_range", "false");
    is_channel_range = (temp_str == "false" || temp_str == "0") ? false : true;
    if (is_channel_range)
    {
        channel_range_start = reader.GetInteger("parameter", "channel_range_start", 0);
        channel_range_end = reader.GetInteger("parameter", "channel_range_end", 1);
    }

    temp_str = reader.Get("parameter", "is_channel_stride", "false");
    is_channel_stride = (temp_str == "false" || temp_str == "0") ? false : true;
    if (is_channel_stride)
    {
        channel_stride_size = reader.GetInteger("parameter", "channel_stride_size", 1);
    }

    temp_str = reader.Get("parameter", "is_column_vector", "NULL");
    if (temp_str == "NULL")
    {
        is_column_major_from_config = false;
    }
    else
    {
        is_column_major = (temp_str == "false" || temp_str == "0") ? false : true;
        is_column_major_from_config = true;
    }

    n_files_to_concatenate = reader.GetInteger("parameter", "n_files_to_concatenate", 2);

    MeasureLengthName = reader.Get("parameter", "attribute_name_measure_length", "MeasureLength[m]");
    SpatialResolutionName = reader.Get("parameter", "attribute_name_spatial_resolution", "SpatialResolution[m]");
    SamplingFrequencyName = reader.Get("parameter", "attribute_name_sampling_frequency", "SamplingFrequency[Hz]");

    temp_str = reader.Get("parameter", "is_output_single_file", "false");

    is_output_single_file = (temp_str == "false") ? false : true;

    output_type = reader.Get("parameter", "output_type", "EP_HDF5");

    output_file_dir = reader.Get("parameter", "output_file_dir", "./tdms-dir-dec/test.h5");

    output_dataset = reader.Get("parameter", "output_dataset", "/data");

    temp_str = reader.Get("parameter", "is_dir_output_match_replace_rgx", "false");

    is_dir_output_match_replace_rgx = (temp_str == "false") ? false : true;

    if (is_dir_output_match_replace_rgx)
    {
        dir_output_match_rgx = reader.Get("parameter", "output_file_regex_match", "^(.*)\\.tdms$");

        dir_output_replace_rgx = reader.Get("parameter", "output_file_regex_replace", "$1.h5");
    }

    temp_str = reader.Get("parameter", "is_space_decimate", "false");

    is_space_decimate = (temp_str == "false") ? false : true;

    if (is_space_decimate)
    {

        space_decimate_rows = reader.GetInteger("parameter", "space_decimate_rows", 32);

        space_decimate_operation = reader.Get("parameter", "space_decimate_operation", "ave");
    }

    DT = reader.GetReal("parameter", "dt", 0.002);

    DT_NEW = reader.GetReal("parameter", "dt_new", 0.008);

    if (!mpi_rank)
    {
        std::cout << "\n\n";
        std::cout << termcolor::red << "Parameters to run the Decimate: ";

        std::cout << termcolor::blue << "\n\n Input parameters: ";
        std::cout << termcolor::magenta << "\n        input_dir_file = " << termcolor::green << input_dir_file;
        std::cout << termcolor::magenta << "\n        input_file_type = " << termcolor::green << input_file_type;

        if (is_input_search_rgx)
        {
            std::cout << termcolor::magenta << "\n        input_search_rgx = " << termcolor::green << input_search_rgx;
        }
        std::cout << termcolor::blue << "\n\n Runtime parameters: ";
        std::cout << termcolor::magenta << "\n\n        lts_per_file = " << termcolor::green << lts_per_file;
        std::cout << termcolor::magenta << "\n        chs_per_file = " << termcolor::green << chs_per_file;
        std::cout << termcolor::magenta << "\n        n_files_to_concatenate = " << termcolor::green << n_files_to_concatenate;

        std::cout << termcolor::magenta << "\n        is_space_decimate = " << termcolor::green << is_space_decimate;

        std::cout << termcolor::magenta << "\n        space_decimate_rows = " << termcolor::green << space_decimate_rows;

        std::cout << termcolor::magenta << "\n        space_decimate_operation = " << termcolor::green << space_decimate_operation;

        std::cout << termcolor::magenta << "\n        DT = " << termcolor::green << DT;

        std::cout << termcolor::magenta << "\n        DT_NEW = " << termcolor::green << DT_NEW;

        std::cout << termcolor::blue << "\n\n Output parameters: ";

        std::cout << termcolor::magenta << "\n        is_output_single_file = " << termcolor::green << is_output_single_file;
        std::cout << termcolor::magenta << "\n        output_type = " << termcolor::green << output_type;
        std::cout << termcolor::magenta << "\n        output_file_dir = " << termcolor::green << output_file_dir;
        std::cout << termcolor::magenta << "\n        output_dataset = " << termcolor::green << output_dataset;

        if (is_dir_output_match_replace_rgx)
        {
            std::cout << termcolor::magenta << "\n        dir_output_match_rgx = " << termcolor::green << dir_output_match_rgx;
            std::cout << termcolor::magenta << "\n        dir_output_replace_rgx = " << termcolor::green << dir_output_replace_rgx;
        }

        std::cout << termcolor::reset << "\n\n";
    }
    fflush(stdout);

    return 0;
}
