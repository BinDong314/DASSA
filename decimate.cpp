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
std::string output_dataset = "/dat";

bool is_dir_output_match_replace_rgx = false;

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

void InitDecimate()
{
    int nPoint = ceil(lts_per_file * time_decimate_files / (DT_NEW / DT));
    cut_frequency_low = (0.5 / DT_NEW) / (0.5 / DT);
    ButterLow(butter_order, cut_frequency_low, BUTTER_A, BUTTER_B);
    if (!au_rank)
        std::cout << "After decimate, nPoint = " << nPoint << "\n";
}

inline Stencil<std::vector<double>> udf_decimate(const Stencil<short> &iStencil)
{
    std::vector<int> max_offset_upper;
    iStencil.GetOffsetUpper(max_offset_upper);
    PrintVector("max_offset_upper = ", max_offset_upper);
    int chs_per_file_udf = max_offset_upper[0] + 1, lts_per_file_udf = max_offset_upper[1] + 1;
    std::vector<int> start_offset{0, 0}, end_offset{chs_per_file_udf - 1, lts_per_file_udf - 1};
    std::vector<short> ts_short = iStencil.ReadNeighbors(start_offset, end_offset);
    //std::vector<double> ts(ts_short.begin(), ts_short.end());
    std::vector<std::vector<double>> ts2d = DasLib::Vector1D2D<short, double>(lts_per_file_udf, ts_short), ts2d_ma;

    std::vector<double> ts_temp2;
    //Resample in time-domain
    for (int i = 0; i < chs_per_file_udf; i++)
    {
        detrend(ts2d[i].data(), lts_per_file_udf);       //Detread
        filtfilt(BUTTER_A, BUTTER_B, ts2d[i], ts_temp2); //filtfilt
        resample(1, DT_NEW / DT, ts_temp2, ts2d[i]);     //resample
    }
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
    std::vector<double> ts_temp = Convert2DVTo1DV(ts2d_ma);
    Stencil<std::vector<double>> oStencil;
    std::vector<size_t> vector_shape(2);
    vector_shape[0] = ts2d_ma.size();
    vector_shape[1] = ts2d_ma[0].size();
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

    if (has_config_file_flag)
        read_config_file(config_file, au_rank);

    // set up the chunk size and the overlap size
    // 11648, 30000 for each dataset
    std::vector<int> chunk_size(2);
    chunk_size[0] = chs_per_file;
    chunk_size[1] = lts_per_file * time_decimate_files;
    std::vector<int> overlap_size = {0, 0};

    std::cout << "EP_DIR:" + input_file_type + ":" + input_dir << "\n";
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
        std::cout << termcolor::magenta << "\n        input_dir = " << termcolor::green << input_dir;
        std::cout << termcolor::magenta << "\n        input_file_type = " << termcolor::green << input_file_type;

        if (is_input_search_rgx)
        {
            std::cout << termcolor::magenta << "\n        input_search_rgx = " << termcolor::green << input_search_rgx;
        }
        std::cout << termcolor::blue << "\n\n Runtime parameters: ";
        std::cout << termcolor::magenta << "\n\n        lts_per_file = " << termcolor::green << lts_per_file;
        std::cout << termcolor::magenta << "\n        chs_per_file = " << termcolor::green << chs_per_file;
        std::cout << termcolor::magenta << "\n        time_decimate_files = " << termcolor::green << time_decimate_files;

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