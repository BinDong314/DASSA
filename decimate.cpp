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

using namespace std;
using namespace AU;

// help functions
void printf_help(char *cmd);

bool is_result_multiple_files = false;
int chunk_size_row = 11648;
int chunk_size_col = 60000;

double DT = 0.002;
double DT_NEW = 0.008;
int nPoint = ceil(chunk_size_col / (DT_NEW / DT));
int nfft;
FIND_M_POWER2(nPoint, nfft);

//UDF One: duplicate the original data
inline Stencil<std::vector<double>> udf_decimate(const Stencil<short> &iStencil)
{
    std::vector<int> start_offset{0, 0}, end_offset{chunk_size_row - 1, chunk_size_col - 1};
    std::vector<double> ts = iStencil.Read(start_offset, end_offset);
    std::vector<std::vector<double>> ts2d = DasLib::Vector1D2D(chunk_size_col, ts);
    std::vector<double> ts_temp;

    //Resample in time-domain
    for (int i = 0; i < chunk_size_row; i++)
    {
        //Detread
        detrend(ts2d[i].data(), chunk_size_col);

        //filtfilt
        filtfilt(BUTTER_A, BUTTER_B, ts2d[i], ts_temp);

        //resample
        resample(1, DT_NEW / DT, ts_temp, ts2d[i]);
    }

    ts_temp = Vector2D1D(ts2d);
    Stencil<std::vector<double>> oStencil;

    std::vector<size_t> vector_shape(2);
    vector_shape[0] = ts2d.size();
    vector_shape[1] = ts2d[0].size();
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
    chunk_size[0] = chunk_size_row;
    chunk_size[1] = chunk_size_col;
    std::vector<int> overlap_size = {0, 0};

    //Input data
    Array<short> *A = new Array<short>("EP_DIR:EP_TDMS:/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir", chunk_size, overlap_size);
    std::vector<std::string> aug_merge_index, aug_dir_sub_cmd, aug_input_search_rgx;

    aug_merge_index.push_back("1");
    aug_dir_sub_cmd.push_back("BINARY_ENABLE_TRANSPOSE_ON_READ");
    aug_input_search_rgx.push_back("^(.*)[135]\\.tdms$");

    A->EndpointControl(DIR_MERGE_INDEX, aug_merge_index);
    //A->EndpointControl(DIR_SUB_CMD_ARG, aug_dir_sub_cmd1); //Not needed
    A->EndpointControl(DIR_SUB_CMD_ARG, aug_dir_sub_cmd);
    A->EndpointControl(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx);

    //Result data

    Array<double> *B;
    if (!is_result_multiple_files)
    {
        //Store into a single file
        B = new Array<double>("EP_HDF5:./tdms-dir-dec/test.h5:/DataCT");
    }
    else
    {
        //Store into multiple file
        B = new Array<double>("EP_DIR:EP_HDF5:./tdms-dir-dec/:/DataCT");
        //Use the below rgx pattern to name the file
        std::vector<std::string> aug_output_replace_arg;
        aug_output_replace_arg.push_back("^(.*)\\.tdms$");
        aug_output_replace_arg.push_back("$1.h5");
        B->EndpointControl(DIR_MERGE_INDEX, aug_merge_index);
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
