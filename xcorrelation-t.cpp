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
std::string input_dir_file = "./test-data/dir";
std::string input_h5_dataset = "/dat";
std::string input_file_type = "EP_HDF5";
bool is_input_search_rgx = false;
std::string input_search_rgx = "^(.*)[1234]\\.tdms$";

int chs_per_file = 11648;
int lts_per_file = 30000;
int n_files_time_decimate = 1;

bool is_output_single_file = false;
std::string output_type = "EP_HDF5";
std::string output_file_dir = "./test-data/dir-output/";
std::string output_dataset = "/dat";

bool is_dir_output_match_replace_rgx = true;

std::string dir_output_match_rgx = "^(.*)\\.h5$";

std::string dir_output_replace_rgx = "$1-xcorr.h5";

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

//For moveing mean
double WINLEN_SEC = 0.5;
int nPoint_hal_win;

//For interp1
std::vector<double> shapingFilt;
double fNyquist = 250; //250
std::vector<double> INTERP_Z{0, 0.5, 1, 1, 0.5, 0};
std::vector<double> INTERP_ZF{0, 0.002, 0.006, 14.5, 15, fNyquist};
double df;
double eCoeff = 1.0;

unsigned long long MASTER_INDEX = 0;

bool is_test_flag = false;
std::vector<int> output_chunk_size(2);

bool is_column_major = true;

/**
 * @brief input data type
 *    0 : short (by default)
 *    1 : double 
 *    2 : float 
 */
int input_data_type = 0;

void init_xcorr()
{
    int nPoint = ceil(lts_per_file * n_files_time_decimate / (DT_NEW / DT));
    cut_frequency_low = (0.5 / DT_NEW) / (0.5 / DT);
    ButterLow(butter_order, cut_frequency_low, BUTTER_A, BUTTER_B);

    nPoint_hal_win = floor((2 * floor(WINLEN_SEC / DT_NEW / 2) + 1) / 2);

    size_t nfft;
    minpower2(nPoint, nfft);
    shapingFilt.resize(nfft);
    fNyquist = 0.5 / DT_NEW;
    INTERP_ZF[5] = fNyquist;

    if (!ft_rank)
        std::cout << "After decimate, nPoint = " << nPoint << ", lts_per_file =" << lts_per_file << ", DT_NEW = " << DT_NEW << ", DT =  " << DT << " , nfft = " << nfft << "\n";

    df = 2.0 * fNyquist / (double)nfft;

    std::vector<double> FF_LHS, LHS;
    FF_LHS.resize(nfft / 2);
    LHS.resize(nfft / 2);
    for (int i = 0; i < nfft / 2; i++)
    {
        FF_LHS[i] = df * (i + 1);
    }
    interp1(INTERP_ZF, INTERP_Z, FF_LHS, LHS);
    int nfft_half = nfft / 2;
    for (int i = 0; i < nfft_half; i++)
    {
        shapingFilt[i] = LHS[i];
        shapingFilt[i + nfft_half] = LHS[nfft_half - i - 1];
    }
    FF_LHS.clear();
    LHS.clear();
    std::vector<double>().swap(FF_LHS);
    std::vector<double>().swap(LHS);
}

template <class TT>
inline Stencil<std::vector<double>> udf_xcorr(const Stencil<TT> &iStencil)
{
    std::vector<int> max_offset_upper;
    iStencil.GetOffsetUpper(max_offset_upper);
    PrintVector("max_offset_upper = ", max_offset_upper);

    int chs_per_file_udf = max_offset_upper[0] + 1, lts_per_file_udf = max_offset_upper[1] + 1;
    std::vector<int> start_offset = {0, 0};
    std::vector<int> end_offset = {chs_per_file_udf - 1, lts_per_file_udf - 1};

    if (is_channel_range)
    {
        assert(channel_range_end - channel_range_start + 1 <= chs_per_file_udf);
        chs_per_file_udf = channel_range_end - channel_range_start + 1;
        start_offset[0] = channel_range_start;
        end_offset[0] = channel_range_end;
        PrintVector("start_offset = ", start_offset);
        PrintVector("end_offset = ", end_offset);
    }

    if (is_many_files)
    {
        if (is_space_decimate)
        {
            many_files_split_n = space_decimate_rows;
        }
        else
        {
            many_files_split_n = chs_per_file_udf / n_files_time_decimate;
        }

        if (!ft_rank)
            std::cout << "Using the is_many_files, many_files_split_n = " << many_files_split_n << " \n";
    }

    std::vector<TT> ts_short;
    iStencil.ReadNeighbors(start_offset, end_offset, ts_short);

    if (is_column_major)
    {
        std::vector<TT> ts_short_temp;
        ts_short_temp.resize(ts_short.size());
        transpose(ts_short.data(), ts_short_temp.data(), chs_per_file_udf, lts_per_file_udf);
        ts_short = ts_short_temp;
        int temp = chs_per_file_udf;
        chs_per_file_udf = lts_per_file_udf;
        lts_per_file_udf = temp;
    }

    std::vector<std::vector<double>> ts2d = DasLib::Vector1D2D(lts_per_file_udf, ts_short);
    //PrintVV("ts2d :", ts2d);

    std::vector<std::vector<double>> ts2d_ma;

    std::cout << "ts2d.size() = " << ts2d.size() << ",ts2d[0].size() = " << ts2d[0].size() << ", lts_per_file_udf =" << lts_per_file_udf << ", ts_short.size() = " << ts_short.size() << "\n";
    std::cout << "Got data ! at rank " << ft_rank << " \n";

    std::vector<double> ts_temp2;
    //Resample in time-domain
    for (int i = 0; i < chs_per_file_udf; i++)
    {
        detrend(ts2d[i].data(), lts_per_file_udf); //Detread
        //PrintVector("ts2d[i] = ", ts2d[i]);
        filtfilt(BUTTER_A, BUTTER_B, ts2d[i], ts_temp2); //filtfilt
        resample(1, DT_NEW / DT, ts_temp2, ts2d[i]);     //resample
    }

    //PrintVV("ts2d :", ts2d);

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

        if (!ft_rank)
            std::cout << "Finish space-domain decimate ! \n";
    }
    else
    {
        ts2d_ma = ts2d;
    }
    DasLib::clear_vector(ts2d);

    //***************
    //Move Aveage
    //MOVING_MEAN(X_l, C_l, nPoint_hal_win);

    ts2d.resize(ts2d_ma.size());
    std::vector<std::complex<double>> fft_out, fft_in, master_fft;
    size_t nPoint_before_fft;
    double temp_real, temp_imag, temp_f;
    for (int i_row = 0; i_row < ts2d_ma.size(); i_row++)
    {
        moving_mean(ts2d_ma[i_row], ts2d[i_row], nPoint_hal_win);
        nPoint_before_fft = ts2d[i_row].size();
        fftv_forward_p2(ts2d[i_row], fft_out);
        fft_in.resize(fft_out.size());
        for (int ii = 0; ii < fft_out.size(); ii++)
        {
            temp_f = pow(sqrt(fft_out[ii].real() * fft_out[ii].real() + fft_out[ii].imag() * fft_out[ii].imag()), eCoeff) + 0.001;
            temp_real = (fft_out[ii].real() + 0.001) / temp_f * shapingFilt[ii];
            fft_in[ii].real(temp_real);
            temp_imag = (fft_out[ii].imag()) / temp_f * shapingFilt[ii];
            fft_in[ii].imag(temp_imag);
        }
        fftv_backward(fft_in, fft_out);

        ts2d[i_row].resize(nPoint_before_fft);
        for (int i = 0; i < nPoint_before_fft; i++)
        {
            ts2d[i_row][i] = fft_out[i].real();
        }

        fftv_forward_p2(ts2d[i_row], fft_out);

        if (i_row == MASTER_INDEX)
        {
            master_fft = fft_out;
            /*
            for (int ij = 0; ij < master_fft.size(); ij++)
            {
                std::cout << " master_fft [" << ij << "] =  " << master_fft[ij] << " \n";
            }*/
        }

        fft_in.clear();
        fft_in.resize(fft_out.size());
        for (int j = 0; j < fft_out.size(); j++)
        {
            fft_in[j].real(master_fft[j].real() * fft_out[j].real() + master_fft[j].imag() * fft_out[j].imag());
            fft_in[j].imag(master_fft[j].imag() * fft_out[j].real() - master_fft[j].real() * fft_out[j].imag());

            //std::cout << " fft_in [" << j << "] =  " << fft_in[j] << " \n";
        }

        fftv_backward(fft_in, fft_out);

        ts2d_ma[i_row].resize(2 * nPoint_before_fft - 1);
        int ts2d_ma_index = 0;
        for (int k = fft_out.size() - nPoint_before_fft + 1; k < fft_out.size(); k++)
        {
            ts2d_ma[i_row][ts2d_ma_index] = fft_out[k].real();
            ts2d_ma_index++;
        }

        for (int l = 0; l < nPoint_before_fft; l++)
        {
            ts2d_ma[i_row][ts2d_ma_index] = fft_out[l].real();
            ts2d_ma_index++;
        }
    }

    //****************
    std::vector<double> ts_temp = Convert2DVTo1DV(ts2d_ma);
    Stencil<std::vector<double>> oStencil;
    std::vector<size_t> vector_shape(2);
    vector_shape[0] = ts2d_ma.size();
    vector_shape[1] = ts2d_ma[0].size();
    //PrintVector("vector_shape: ", vector_shape);
    //std::cout << "vector_shape[0] = " << vector_shape[0] << ",vector_shape[1] = " << vector_shape[1] << "\n";
    DasLib::clear_vector(ts2d_ma);
    oStencil.SetShape(vector_shape);
    oStencil = ts_temp;

    //
    // Deal with tag
    //
    if (iStencil.HasTagMap())
    {
        std::map<std::string, std::string> tag_map;
        iStencil.GetTagMap(tag_map);
        for (std::map<std::string, std::string>::iterator it = tag_map.begin(); it != tag_map.end(); ++it)
        {
            std::cout << " key : " << it->first << ", value:" << it->second << " \n";
            if (it->first == "SamplingFrequency[Hz]")
            {
                it->second = "125";
            }
        }
        oStencil.SetTagMap(tag_map);
    }

    return oStencil;
}

int main(int argc, char *argv[])
{
    int copt;
    bool has_config_file_flag = false;
    while ((copt = getopt(argc, argv, "c:ht")) != -1)
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

    std::cout << "EP_DIR:" + input_file_type + ":" + input_dir_file << ":" << input_h5_dataset << "\n";
    std::string A_endpoint_id;

    if (!is_input_single_file)
    {
        A_endpoint_id = "EP_DIR:" + input_file_type + ":" + input_dir_file + ":" + input_h5_dataset;
    }
    else
    {
        A_endpoint_id = input_file_type + ":" + input_dir_file + ":" + input_h5_dataset;
    }

    //Input data
    //AU::Array<short> *

    FT::ArrayBase *A;
    AuEndpointDataType t;
    if (input_data_type == 0)
    {
        t = AuEndpointDataType::AU_SHORT;
        A = new FT::Array<short>(A_endpoint_id);
    }
    else if (input_data_type == 1)
    {
        t = AuEndpointDataType::AU_DOUBLE;
        A = new FT::Array<double>(A_endpoint_id);
    }
    else if (input_data_type == 2)
    {
        t = AuEndpointDataType::AU_FLOAT;
        A = new FT::Array<float>(A_endpoint_id);
    }
    else
    {
        std::cout << "Not supported input_data_type \n";
        exit(-1);
    }

    A->GetStencilTag();

    if (!is_input_single_file)
    {
        std::vector<std::string> file_size_str;
        A->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str);
        String2Vector(file_size_str[0], chunk_size);
    }
    else
    {
        std::vector<unsigned long long> array_size;
        A->GetArraySize(array_size);
        chunk_size[0] = array_size[0];
        chunk_size[1] = array_size[1];
    }

    std::cout << "A_endpoint_id = " << A_endpoint_id << "\n";

    if (!is_column_major)
    {
        chunk_size[1] = chunk_size[1] * n_files_time_decimate;
        chs_per_file = chunk_size[0];
        lts_per_file = chunk_size[1];
    }
    else
    {
        chunk_size[0] = chunk_size[0] * n_files_time_decimate;
        chs_per_file = chunk_size[1];
        lts_per_file = chunk_size[0];
    }

    A->SetChunkSize(chunk_size);
    A->SetOverlapSize(overlap_size);

    std::cout << "chunk_size = " << chunk_size[0] << " , " << chunk_size[1] << " \n";

    std::vector<std::string> aug_merge_index;
    aug_merge_index.push_back("1");
    A->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);

    if (is_input_search_rgx)
    {
        std::vector<std::string> aug_input_search_rgx;
        aug_input_search_rgx.push_back(input_search_rgx);
        A->ControlEndpoint(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx);
    }

    init_xcorr();
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
        B->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);

        if (is_dir_output_match_replace_rgx)
            B->ControlEndpoint(DIR_OUPUT_REPLACE_RGX, aug_output_replace_arg);
    }

    //Stride on execution
    //Each chunk only runs the udf_decimate once
    A->EnableApplyStride(chunk_size);

    A->SetVectorDirection(AU_FLAT_OUTPUT_ROW);

    //Run
    //A->Transform<std::vector<double>>(udf_xcorr, B);
    TRANSFORM(A, udf_xcorr, B, t, std::vector<double>);
    A->ReportCost();
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
          -c config file for parameters (has high priority than commands if existing) \n\
          -t do a quick test with a short input\n\
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

    input_dir_file = reader.Get("parameter", "input_dir_file", "/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir");

    input_h5_dataset = reader.Get("parameter", "input_dataset", "/dat");

    input_file_type = reader.Get("parameter", "input_file_type", "EP_TDMS");

    input_data_type = reader.GetInteger("parameter", "input_data_type", 0);

    std::string temp_str = reader.Get("parameter", "is_input_single_file", "false");
    is_input_single_file = (temp_str == "false") ? false : true;

    temp_str = reader.Get("parameter", "is_input_search_rgx", "false");

    is_input_search_rgx = (temp_str == "false") ? false : true;

    if (is_input_search_rgx)
    {
        input_search_rgx = reader.Get("parameter", "input_search_rgx", "^(.*)[1234]\\.tdms$");
    }

    //chs_per_file = reader.GetInteger("parameter", "chs_per_file", 11648);
    //lts_per_file = reader.GetInteger("parameter", "lts_per_file", 30000);

    temp_str = reader.Get("parameter", "is_channel_range", "false");

    is_channel_range = (temp_str == "false" || temp_str == "0") ? false : true;
    if (is_channel_range)
    {
        channel_range_start = reader.GetInteger("parameter", "channel_range_start", 0);
        channel_range_end = reader.GetInteger("parameter", "channel_range_end", 1);
    }

    temp_str = reader.Get("parameter", "is_column_major", "true");
    is_column_major = (temp_str == "false" || temp_str == "0") ? false : true;

    n_files_time_decimate = reader.GetInteger("parameter", "n_files_time_decimate", 1);

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

    std::string temp_str2 = reader.Get("parameter", "z", "0, 0.5, 1, 1, 0.5, 0");
    std::stringstream iss(temp_str2);
    double number;
    std::vector<double> Z;
    while (iss >> number)
    {
        Z.push_back(number);
        if (iss.peek() == ',')
            iss.ignore();
    }
    INTERP_Z = Z;
    INTERP_ZF[0] = reader.GetReal("parameter", "F1", 0);
    INTERP_ZF[1] = reader.GetReal("parameter", "F2", 0.002);
    INTERP_ZF[2] = reader.GetReal("parameter", "F3", 0.006);
    INTERP_ZF[3] = reader.GetReal("parameter", "F4", 14.5);
    INTERP_ZF[4] = reader.GetReal("parameter", "F5", 15);

    WINLEN_SEC = reader.GetReal("parameter", "winLen_sec", 0.5);

    eCoeff = reader.GetReal("parameter", "eCoeff", 1.0);
    MASTER_INDEX = reader.GetInteger("parameter", "master_index", 0);

    butter_order = reader.GetInteger("parameter", "butter_order", 3);

    if (!mpi_rank)
    {
        std::cout << "\n\n";
        std::cout << termcolor::red << "Parameters to run the Decimate: ";

        std::cout << termcolor::blue << "\n\n Input parameters: ";
        std::cout << termcolor::magenta << "\n        input_dir_file = " << termcolor::green << input_dir_file;
        std::cout << termcolor::magenta << "\n        input_file_type = " << termcolor::green << input_file_type;
        std::cout << termcolor::magenta << "\n        input_data_type = " << termcolor::green << input_data_type;
        std::cout << termcolor::magenta << "\n        is_column_major = " << termcolor::green << is_column_major;

        if (is_input_search_rgx)
        {
            std::cout << termcolor::magenta << "\n        input_search_rgx = " << termcolor::green << input_search_rgx;
        }
        std::cout << termcolor::blue << "\n\n Runtime parameters: ";
        // std::cout << termcolor::magenta << "\n\n        lts_per_file = " << termcolor::green << lts_per_file;
        // std::cout << termcolor::magenta << "\n        chs_per_file = " << termcolor::green << chs_per_file;
        std::cout << termcolor::magenta << "\n        n_files_time_decimate = " << termcolor::green << n_files_time_decimate;

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
