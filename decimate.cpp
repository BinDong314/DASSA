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

#include "au.h"
#include "DasLib.h"

using namespace std;
using namespace AU;
using namespace DasLib;
// help functions

std::string input_dir = "/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir";
std::string input_file_type = "EP_TDMS";
bool is_input_search_rgx = false;
std::string input_search_rgx = "^(.*)[1234]\\.tdms$";

int chs_per_file = 11648;
int lts_per_file = 30000;
int time_decimate_files = 2;

bool is_output_single_file = true;
std::string output_type = "EP_HDF5";
std::string output_file_dir = "./tdms-dir-dec/test.h5";
std::string output_dataset = "/DataCT";

bool is_dir_output_match_replace_rgx = false;

std::string dir_output_match_rgx = "^(.*)\\.tdms$";

std::string dir_output_replace_rgx = "$1.h5";

bool is_space_decimate = false;
int space_decimate_rows = 32;
std::string space_decimate_operation = "mean"; //"ave", "min", "max"

double DT = 0.002;
double DT_NEW = 0.008;

void printf_help(char *cmd);

//
int nPoint = ceil(chunk_size_col / (DT_NEW / DT));

vector<double> BUTTER_A;
vector<double> BUTTER_B;

//int nfft;
//FIND_M_POWER2(nPoint, nfft);

int butter_order = 3;
double cut_frequency_low = 0.25;

void InitDecimate()
{
    cut_frequency_low = (0.5 / DT_NEW) / (0.5 / DT);
    ButterLow(butter_order, cut_frequency_low, BUTTER_A, BUTTER_B);
    if (!au_rank)
        std::cout << " nPoint = " << nPoint << "\n";
}

inline Stencil<std::vector<double>> udf_decimate(const Stencil<short> &iStencil)
{
    std::vector<int> start_offset{0, 0}, end_offset{chunk_size_row - 1, chunk_size_col - 1};
    std::vector<short> ts_short = iStencil.Read(start_offset, end_offset);
    std::vector<double> ts(ts_short.begin(), ts_short.end());
    std::vector<std::vector<double>> ts2d = DasLib::Vector1D2D(chunk_size_col, ts), ts2d_ma;
    //Resample in time-domain
    for (int i = 0; i < chunk_size_row; i++)
    {
        detrend(ts2d[i].data(), chunk_size_col);        //Detread
        filtfilt(BUTTER_A, BUTTER_B, ts2d[i], ts_temp); //filtfilt
        resample(1, DT_NEW / DT, ts_temp, ts2d[i]);     //resample
    }
    if (is_space_decimate)
    {
        //Moving mean in space-domain
        int ma_batches = chunk_size_row / rows_to_mean, ma_batches_remainder = chunk_size_row % rows_to_mean;
        for (int i = 0; i < ma_batches; i++)
        {
            ts2d_ma.push_back(SpaceMoveMean(ts2d, i * rows_to_mean, (i + 1) * rows_to_mean - 1));
        }
        if (ma_batches_remainder != 0)
        {
            ts2d_ma.push_back(SpaceMoveMean(ts2d, chunk_size_row - ma_batches_remainder, chunk_size_row - 1));
        }
    }
    else
    {
        ts2d_ma = ts2d;
    }
    std::vector<double> ts_temp = Convert2DVTo1DV(ts2d_ma);
    Stencil<std::vector<double>> oStencil;
    std::vector<size_t> vector_shape(2);
    vector_shape[0] = ts2d_ma.size();
    vector_shape[1] = ts2d_ma[0].size();
    oStencil.SetOutputVectorShape(vector_shape);
    oStencil = ts_temp;
    return oStencil;
}

int main(int argc, char *argv[])
{
    //Init the MPICH, etc.
    AU_Init(argc, argv);

    // set up the chunk size and the overlap size
    // 11648, 30000 for each dataset
    std::vector<int> chunk_size(2);
    chunk_size[0] = chs_per_file;
    chunk_size[1] = lts_per_file * time_decimate_files;
    std::vector<int> overlap_size = {0, 0};

    //Input data
    AU::Array<short> *A = new AU::Array<short>("EP_DIR:" + input_file_type + ":" + input_dir, chunk_size, overlap_size);

    ///Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir", chunk_size, overlap_size);

    std::vector<std::string> aug_merge_index, aug_dir_sub_cmd, aug_input_search_rgx;

    aug_merge_index.push_back("1");
    aug_dir_sub_cmd.push_back("BINARY_ENABLE_TRANSPOSE_ON_READ");

    aug_input_search_rgx.push_back(input_search_rgx);

    A->EndpointControl(DIR_MERGE_INDEX, aug_merge_index);
    //A->EndpointControl(DIR_SUB_CMD_ARG, aug_dir_sub_cmd1); //Not needed
    A->EndpointControl(DIR_SUB_CMD_ARG, aug_dir_sub_cmd);
    if (is_input_search_rgx)
        A->EndpointControl(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx);

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
    A->Apply(udf_decimate, B);

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
          -i input dir\n\
          -o output file name \n\
	      -t dataset name for intput time series \n\
          -d dataset name for output decimation \n\
          -r chunk size row \n\
          -l chunk size col \n\
          -c config file for parameters (has high priority than commands if existing) \n\
          Example: mpirun -n 1 %s -i  -o fft-test.arrayudf.h5  -g / -t /DataCT -x /Xcorr\n";

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

    input_dir = reader.Get("parameter", "input_dir", "/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir");

    input_file_type = reader.Get("parameter", "input_file_type", "EP_TDMS");

    std::string temp_str = reader.Get("parameter", "is_input_search_rgx", "false");

    is_input_search_rgx = (temp_str == "false") ? false : true;

    if (is_input_search_rgx)
    {
        input_search_rgx = reader.Get("parameter", "input_search_rgx", "^(.*)[1234]\\.tdms$");
    }

    chs_per_file = reader.GetInteger("parameter", "chs_per_file", 11648);
    lts_per_file = reader.GetInteger("parameter", "lts_per_file", 30000);

    time_decimate_files = reader.GetInteger("parameter", "time_decimate_files", 2);

    temp_str = reader.Get("parameter", "is_output_single_file", "false");

    is_output_single_file = (temp_str == "false") ? false : true;

    output_type = reader.Get("parameter", "output_type", "EP_HDF5");

    output_file_dir = reader.Get("parameter", "output_file_dir", "./tdms-dir-dec/test.h5");

    output_dataset = reader.Get("parameter", "output_dataset", "/data");

    temp_str = reader.Get("parameter", "is_dir_output_match_replace_rgx", "false");

    is_dir_output_match_replace_rgx = (temp_str == "false") ? false : true;

    if (is_dir_output_match_replace_rgx)
    {
        output_file_regex_match = reader.Get("parameter", "output_file_regex_match", "^(.*)\\.tdms$");

        output_file_regex_replace = reader.Get("parameter", "output_file_regex_replace", "$1.h5");
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
    return 0;
}