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
std::string config_file = "./xcorrelation.config";

bool is_input_single_file = false;
std::string input_dir_file = "./test-data/dir";
std::string input_h5_dataset = "/dat";
std::string input_file_type = "EP_HDF5";

bool is_input_search_rgx = false;
std::string input_search_rgx = "^(.*)[1234]\\.tdms$";

int chs_per_file = 11648;
int lts_per_file = 30000;
int n_files_to_concatenate = 1;

bool is_output_single_file = false;
std::string output_type = "EP_HDF5";
std::string output_file_dir = "./test-data/dir-output/";
std::string output_dataset = "/dat";

bool is_dir_output_match_replace_rgx = false;

std::string dir_output_match_rgx = "^(.*)\\.h5$";

std::string dir_output_replace_rgx = "$1-xcorr.h5";

bool is_space_decimate = false;
int space_decimate_chs = 32;
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
int round_dt_dew_dt;
void printf_help(char *cmd);

//
vector<double> BUTTER_A;
vector<double> BUTTER_B;

// int nfft;
// FIND_M_POWER2(nPoint, nfft);

int butter_order = 3;
double cut_frequency_low = 0.25;

bool is_channel_range = false;
int channel_range_start = 0;
int channel_range_end = 1;
bool is_many_files = false;
int many_files_split_n = 10;

// For moveing mean
double WINLEN_SEC = 0.5;
int nPoint_hal_win;

// For interp1
std::vector<double> shapingFilt;
double fNyquist = 250; // 250
std::vector<double> INTERP_Z{0, 0.5, 1, 1, 0.5, 0};
std::vector<double> INTERP_ZF{0, 0.002, 0.006, 14.5, 15, fNyquist};
double df;
double eCoeff = 1.0;

unsigned long long MASTER_INDEX = 0;

bool is_test_flag = false;
std::vector<int> output_chunk_size(2);

bool is_column_major = true;
bool is_column_major_from_config = false;

/**
 * @brief input data type
 *    0 : short (by default)
 *    1 : double
 *    2 : float
 */
int input_data_type = 0;

std::string MeasureLengthName = "MeasureLength[m]";
std::string SpatialResolutionName = "SpatialResolution[m]";
std::string SamplingFrequencyName = "SamplingFrequency[Hz]";

bool is_channel_stride = false;
int channel_stride_size = 0;

bool is_file_range = false;
int file_range_start_index = 0;
int file_range_end_index = 1;
std::string file_range_indexes_str;
std::vector<std::string> aug_merge_index;

double time_begin, time_end; // the time for the begin and the end of result series

bool is_setview_in_main = true;
void init_xcorr()
{
    round_dt_dew_dt = round(DT_NEW / DT);
    // int nPoint = ceil(lts_per_file / (DT_NEW / DT));
    int nPoint = ceil(lts_per_file / round_dt_dew_dt);
    cut_frequency_low = (0.5 / DT_NEW) / (0.5 / DT);
    ButterLow(butter_order, cut_frequency_low, BUTTER_A, BUTTER_B);

    nPoint_hal_win = floor((2 * floor(WINLEN_SEC / DT_NEW / 2) + 1) / 2);

    size_t nfft;
    minpower2(nPoint, nfft);
    shapingFilt.resize(nfft);
    fNyquist = 0.5 / DT_NEW;
    INTERP_ZF[5] = fNyquist;

    if (!ft_rank)
    {
        std::cout << "In init_xcorr, nPoint = " << nPoint << ", lts_per_file =" << lts_per_file << ", DT_NEW = " << DT_NEW << ", DT =  " << DT << " , nfft = " << nfft << ", round_dt_dew_dt = " << round_dt_dew_dt << "\n";
        PrintVector("BUTTER_A = ", BUTTER_A);
        PrintVector("BUTTER_B = ", BUTTER_B);
    }

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
    if (!ft_rank)
        PrintVector("max_offset_upper = ", max_offset_upper);

    // Here, we assume data is row-vector
    // int chs_per_file_udf = max_offset_upper[0] + 1, lts_per_file_udf = max_offset_upper[1] + 1;
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

    if (!is_setview_in_main && is_channel_range)
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
        PrintVector("start_offset after channel_range  = ", start_offset);
        PrintVector("end_offset after channel_range = ", end_offset);
    }

    std::vector<TT> ts_short;
    iStencil.ReadNeighbors(start_offset, end_offset, ts_short);

    // Convert to row-vector here if it is column-vector
    // Because all the following code are built as row-vector (2D vector)
    // Each row is a time series
    if (is_column_major)
    {
        std::vector<TT> ts_short_temp;
        ts_short_temp.resize(ts_short.size());
        transpose(ts_short.data(), ts_short_temp.data(), lts_per_file_udf, chs_per_file_udf);
        ts_short = ts_short_temp;
        //    int temp = chs_per_file_udf;
        //    chs_per_file_udf = lts_per_file_udf;
        //    lts_per_file_udf = temp;
    }

    // std::cout << "Before Vector1D2D ), chs_per_file_udf = " << chs_per_file_udf << ", lts_per_file_udf = " << lts_per_file_udf << "\n";
    std::vector<std::vector<double>> ts2d = DasLib::Vector1D2D(lts_per_file_udf, ts_short);

    // PrintVV("ts2d before space decimate:", ts2d);

    if (is_space_decimate)
    {
        std::vector<std::vector<double>> ts2d_temp;
        // decimate in space-domain
        int ma_batches = (chs_per_file_udf % space_decimate_chs) ? (chs_per_file_udf / space_decimate_chs + 1) : (chs_per_file_udf / space_decimate_chs);
        int start_row, end_row;
        for (int i = 0; i < ma_batches; i++)
        {
            start_row = i * space_decimate_chs;
            end_row = (((i + 1) * space_decimate_chs - 1) < chs_per_file_udf) ? ((i + 1) * space_decimate_chs - 1) : chs_per_file_udf - 1;
            // std::cout << "ma_batches =" << ma_batches << ", start_row = " << start_row << ", end_row =  " << end_row << "\n";
            ts2d_temp.push_back(spacedecimate(ts2d, start_row, end_row, space_decimate_operation));
        }

        if (!ft_rank)
            std::cout << "Finish space-domain decimate with [" << ts2d_temp.size() << "] channels  with [ " << ts2d_temp[0].size() << " ] points\n";

        ts2d = ts2d_temp;
        chs_per_file_udf = ts2d_temp.size();
    }

    //
    // Starts from here: ts_short is row-major order,
    // e.g. ts = {time series 1 , time series 2}
    //

    // std::cout << "chs_per_file_udf (before skip ) = " << ts2d.size() << ", " << ts2d[0].size() << ", chs_per_file_udf = " << chs_per_file_udf << ", lts_per_file_udf = " << lts_per_file_udf << "\n";

    // PrintVV("ts2d before stride after space-dec: \n ", ts2d);

    // We may update the data based on is_channel_stride/channel_stride_size
    // e.g.,  channel_range_start = 0, channel_range_end = 99
    //        channel_stride_size = 99
    if (is_channel_stride)
    {
        std::vector<std::vector<double>> ts2d_temp;
        for (int iiii = 0; iiii < chs_per_file_udf; iiii++)
        {
            if (iiii % channel_stride_size == 0)
            {
                // std::cout << iiii << "\n";
                // ts2d.erase(ts2d.begin() + iiii);
                ts2d_temp.push_back(ts2d[iiii]);
            }
        }
        ts2d = ts2d_temp;
        chs_per_file_udf = ts2d.size();
        if (!ft_rank)
            std::cout << "chs_per_file_udf (after skip ) = " << ts2d.size() << ", each with " << ts2d[0].size() << " points\n";
    }

    // std::cout << "ts2d.size() = " << ts2d.size() << ",ts2d[0].size() = " << ts2d[0].size() << ", lts_per_file_udf =" << lts_per_file_udf << ", ts_short.size() = " << ts_short.size() << "\n";
    // std::cout << "Got data ! at rank " << ft_rank << " \n";
    // PrintVV("ts2d before detrend: ", ts2d);

    std::vector<double> ts_temp2;
    // Resample in time-domain
    for (int i = 0; i < chs_per_file_udf; i++)
    {
        detrend(ts2d[i].data(), lts_per_file_udf); // Detread
        // PrintVector("ts2d[i] = ", ts2d[i]);
        filtfilt(BUTTER_A, BUTTER_B, ts2d[i], ts_temp2); // filtfilt
        resample(1, round_dt_dew_dt, ts_temp2, ts2d[i]); // resample
    }
    DasLib::clear_vector(ts_temp2);
    if (!ft_rank)
        std::cout << "After time-domain decimate, with " << ts2d.size() << " channels,  each with " << ts2d[0].size() << " points \n";
    // PrintVV("ts2d :", ts2d);

    // Find the time_begin, time_end for the result series
    time_end = (ts2d[0].size() - 1) * DT_NEW;
    time_begin = -time_end;
    // if (!ft_rank)
    //     std::cout << "Finish time-domain decimate ! \n";

    std::vector<std::vector<double>> ts2d_ma = ts2d;

    // Moved the blow code to be before
    // if (is_space_decimate)
    // {
    //     //decimate in space-domain
    //     int ma_batches = (chs_per_file_udf % space_decimate_chs) ? (chs_per_file_udf / space_decimate_chs + 1) : (chs_per_file_udf / space_decimate_chs);
    //     int start_row, end_row;
    //     for (int i = 0; i < ma_batches; i++)
    //     {
    //         start_row = i * space_decimate_chs;
    //         end_row = (((i + 1) * space_decimate_chs - 1) < chs_per_file_udf) ? ((i + 1) * space_decimate_chs - 1) : chs_per_file_udf - 1;
    //         //std::cout << "ma_batches =" << ma_batches << ", start_row = " << start_row << ", end_row =  " << end_row << "\n";
    //         ts2d_ma.push_back(spacedecimate(ts2d, start_row, end_row, space_decimate_operation));
    //     }

    //     if (!ft_rank)
    //         std::cout << "Finish space-domain decimate ! \n";
    // }
    // else
    // {
    //     ts2d_ma = ts2d;
    // }
    DasLib::clear_vector(ts2d);

    //***************
    // Move Aveage
    // MOVING_MEAN(X_l, C_l, nPoint_hal_win);

    ts2d.resize(ts2d_ma.size());
    std::vector<std::complex<double>> fft_out, fft_in, master_fft;
    size_t nPoint_before_fft;
    double temp_real, temp_imag, temp_f;

    if (MASTER_INDEX != 0)
    {
        moving_mean(ts2d_ma[MASTER_INDEX], ts2d[MASTER_INDEX], nPoint_hal_win);
        nPoint_before_fft = ts2d[MASTER_INDEX].size();
        fftv_forward_p2(ts2d[MASTER_INDEX], fft_out);
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

        ts2d[MASTER_INDEX].resize(nPoint_before_fft);
        for (int i = 0; i < nPoint_before_fft; i++)
        {
            ts2d[MASTER_INDEX][i] = fft_out[i].real();
        }

        fftv_forward_p2(ts2d[MASTER_INDEX], fft_out);
        master_fft = fft_out;
    }

    for (int i_row = 0; i_row < ts2d_ma.size(); i_row++)
    {
        moving_mean(ts2d_ma[i_row], ts2d[i_row], nPoint_hal_win);
        nPoint_before_fft = ts2d[i_row].size();

        fftv_forward_p2(ts2d[i_row], fft_out);
        fft_in.resize(fft_out.size());
        // std::cout << "nPoint_before_fft =" << nPoint_before_fft << ", ts2d_ma[i_row].size() = " << ts2d_ma[i_row].size() << ", ts2d[i_row].size() = " << ts2d[i_row].size() << ", fft_out.size = " << fft_out.size() << ", shapingFilt.size() =" << shapingFilt.size() << " \n";
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

        if (i_row == 0 && i_row == MASTER_INDEX)
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

            // std::cout << " fft_in [" << j << "] =  " << fft_in[j] << " \n";
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

    if (!ft_rank)
        PrintVector("Output vector_shape: ", vector_shape);
    // std::cout << "vector_shape[0] = " << vector_shape[0] << ",vector_shape[1] = " << vector_shape[1] << "\n";
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
            // std::cout << " key : " << it->first << ", value:" << it->second << " \n";
            if (it->first == "SamplingFrequency[Hz]")
            {
                int temp_sf = 1 / DT_NEW;
                it->second = std::to_string(temp_sf); // "125";
            }
        }
        tag_map.insert({"time_begin", std::to_string(time_begin)});
        tag_map.insert({"time_end", std::to_string(time_end)});
        oStencil.SetTagMap(tag_map);
        // To add Dsi_out.fh{9} = tb; Dsi_out.fh{10} = te;
        // nPoint = size(seis_resamp, 1);
        // t_resamp = dt_new .* (0 : (nPoint - 1)); tb = t_resamp(1); te = t_resamp(end);
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

    // Init the MPICH, etc.
    AU_Init(argc, argv);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int numprocs, rank, namelen;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // MPI_Get_processor_name(processor_name, &namelen);
    // printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

    if (has_config_file_flag)
        read_config_file(config_file, ft_rank);

    // set up the chunk size and the overlap size
    // 11648, 30000 for each dataset
    std::vector<int> chunk_size(2);
    std::vector<int> overlap_size = {0, 0};

    // std::cout << "EP_DIR:" + input_file_type + ":" + input_dir_file << ":" << input_h5_dataset << "\n";
    std::string A_endpoint_id;

    if (!is_input_single_file)
    {
        A_endpoint_id = "EP_DIR:" + input_file_type + ":" + input_dir_file + ":" + input_h5_dataset;
    }
    else
    {
        A_endpoint_id = input_file_type + ":" + input_dir_file + ":" + input_h5_dataset;
    }

    // Input data
    // AU::Array<short> *

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

    if (is_channel_range)
    {
        if (is_column_major)
        {
            A->SetView(channel_range_start, channel_range_end - channel_range_start + 1, 1);
        }
        else
        {
            A->SetView(channel_range_start, channel_range_end - channel_range_start + 1, 0);
        }
    }

    A->GetStencilTag();

    if (is_input_search_rgx && is_input_single_file == false)
    {
        std::vector<std::string> aug_input_search_rgx;
        aug_input_search_rgx.push_back(input_search_rgx);
        A->ControlEndpoint(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx);
    }

    if (is_file_range && is_input_single_file == false)
    {
        std::vector<size_t> file_range_index;
        for (size_t i = file_range_start_index; i <= file_range_end_index; i++)
        {
            file_range_index.push_back(i);
        }
        file_range_indexes_str = Vector2String(file_range_index);
        // if (!ft_rank)
        //     std::cout << " file_range_indexes_str =" << file_range_indexes_str << "\n";

        std::vector<std::string> index_param;
        index_param.push_back(file_range_indexes_str);
        A->ControlEndpoint(DIR_FILE_SORT_INDEXES, index_param);

        // index_param.clear();
        // A->ControlEndpoint(DIR_SKIP_SIZE_CHECK, index_param);
    }

    if (!is_input_single_file)
    {
        if (is_column_major)
        {
            aug_merge_index.push_back("0");
        }
        else
        {
            aug_merge_index.push_back("1");
        }
        A->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);
    }

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

    if (!ft_rank)
        PrintVector("Xcorr: chunk_size = ", chunk_size);

    // By default, the TDMS file is column major
    if (input_file_type == "EP_TDMS")
    {
        is_column_major_from_config = true; // Set to skip the check from file
        is_column_major = true;
    }

    // std::cout << "A_endpoint_id = " << A_endpoint_id << "\n";
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

        // std::cout << "MeasureLength= " << MeasureLength << " , SpatialResolution =" << SpatialResolution << ", SamplingFrequency = " << SamplingFrequency << std::endl;
        // if (!MeasureLength.empty() && !SpatialResolution.empty() && !SamplingFrequency.empty())
        //{
        //     meta_chs = std::stoi(MeasureLength) / std::stoi(SpatialResolution);
        //     meta_time_series_points = 60 * std::stoi(SamplingFrequency);
        // }
        // else
        //{
        //     std::cout << "Metadata can not be found for, " << MeasureLengthName << ", " << SpatialResolutionName << ", and " << SamplingFrequencyName << ", please specify the is_column_vector in config file, check attribute_name_*, ... " << std::endl;
        //     exit(-1);
        // }

        // std::cout << "meta_time_series_points = " << meta_time_series_points << " , meta_chs =  " << meta_chs << " \n";
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
        if (is_channel_range)
        {
            chunk_size[1] = channel_range_end - channel_range_start + 1;
            // A->SetView(channel_range_start, chunk_size[1], 1);
        }
        chunk_size[0] = chunk_size[0] * n_files_to_concatenate;
        chs_per_file = chunk_size[1];
        lts_per_file = chunk_size[0];
    }
    else
    {
        if (is_channel_range)
        {
            chunk_size[0] = channel_range_end - channel_range_start + 1;
            // A->SetView(channel_range_start, chunk_size[0], 0);
        }
        chunk_size[1] = chunk_size[1] * n_files_to_concatenate;
        chs_per_file = chunk_size[0];
        lts_per_file = chunk_size[1];
    }

    A->SetChunkSize(chunk_size);
    A->SetOverlapSize(overlap_size);

    // std::cout << "chunk_size = " << chunk_size[0] << " , " << chunk_size[1] << " \n";

    init_xcorr();
    // Result data

    AU::Array<double> *B;
    if (is_output_single_file)
    {
        // Store into a single file
        B = new AU::Array<double>(output_type + ":" + output_file_dir + ":" + output_dataset);
    }
    else
    {
        // Store into multiple file
        B = new AU::Array<double>("EP_DIR:" + output_type + ":" + output_file_dir + ":" + output_dataset);
        // Use the below rgx pattern to name the file
        std::vector<std::string> aug_output_replace_arg;
        aug_output_replace_arg.push_back(dir_output_match_rgx);
        aug_output_replace_arg.push_back(dir_output_replace_rgx);
        B->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);

        if (is_dir_output_match_replace_rgx)
            B->ControlEndpoint(DIR_OUPUT_REPLACE_RGX, aug_output_replace_arg);
    }

    // Stride on execution
    // Each chunk only runs the udf_decimate once
    A->EnableApplyStride(chunk_size);

    A->SetVectorDirection(AU_FLAT_OUTPUT_ROW);

    // Run
    // A->Transform<std::vector<double>>(udf_xcorr, B);
    TRANSFORM(A, udf_xcorr, B, t, std::vector<double>);
    A->ReportCost();

    // Clear
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
          Example: mpirun -n 1 %s -c %s.config \n";

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

    std::string temp_str_str;
    temp_str_str = reader.Get("parameter", "input_data_type", "short");
    if (temp_str_str == "short" || temp_str_str == "0")
    {
        input_data_type = 0;
    }
    else if (temp_str_str == "double" || temp_str_str == "1")
    {
        input_data_type = 1;
    }
    else if (temp_str_str == "float" || temp_str_str == "2")
    {
        input_data_type = 2;
    }
    else
    {
        AU_EXIT("Don't understand the [input_data_type] in config file. It can take short/double/float or 0/1/2\n");
    }

    std::string temp_str = reader.Get("parameter", "is_input_single_file", "false");
    is_input_single_file = (temp_str == "false") ? false : true;

    temp_str = reader.Get("parameter", "is_input_search_rgx", "false");

    is_input_search_rgx = (temp_str == "false") ? false : true;

    if (is_input_search_rgx)
    {
        input_search_rgx = reader.Get("parameter", "input_search_rgx", "^(.*)[1234]\\.tdms$");
    }

    // chs_per_file = reader.GetInteger("parameter", "chs_per_file", 11648);
    // lts_per_file = reader.GetInteger("parameter", "lts_per_file", 30000);

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

    std::string is_file_range_str = reader.Get("parameter", "is_file_range", "false");
    if (is_file_range_str == "false" || is_file_range_str == "0")
    {
        is_file_range = false;
    }
    else if (is_file_range_str == "true" || is_file_range_str == "1")
    {
        is_file_range = true;
    }
    else
    {
        AU_EXIT("Don't read the is_file_range's value " + is_file_range_str);
    }

    if (is_file_range)
    {
        file_range_start_index = reader.GetInteger("parameter", "file_range_start_index", 0);
        file_range_end_index = reader.GetInteger("parameter", "file_range_end_index", 1);
    }

    MeasureLengthName = reader.Get("parameter", "attribute_name_measure_length", "MeasureLength[m]");
    SpatialResolutionName = reader.Get("parameter", "attribute_name_spatial_resolution", "SpatialResolution[m]");
    SamplingFrequencyName = reader.Get("parameter", "attribute_name_sampling_frequency", "SamplingFrequency[Hz]");

    n_files_to_concatenate = reader.GetInteger("parameter", "n_files_to_concatenate", 1);

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

        space_decimate_chs = reader.GetInteger("parameter", "space_decimate_chs", 32);

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
        std::cout << termcolor::red << "Parameters to run the Xcorrelate : ";

        std::cout << termcolor::blue << "\n\n Input parameters: ";
        std::cout << termcolor::magenta << "\n        input_dir_file = " << termcolor::green << input_dir_file;
        std::cout << termcolor::magenta << "\n        input_file_type = " << termcolor::green << input_file_type;
        if (input_data_type == 0)
        {
            std::cout << termcolor::magenta << "\n        input_data_type = " << termcolor::green << "short";
        }
        else if (input_data_type == 1)
        {
            std::cout << termcolor::magenta << "\n        input_data_type = " << termcolor::green << "double";
        }
        else if (input_data_type == 2)
        {
            std::cout << termcolor::magenta << "\n        input_data_type = " << termcolor::green << "float";
        }

        if (is_column_major_from_config)
        {
            if (is_column_major == true)
            {
                std::cout << termcolor::magenta << "\n        is_column_vector = " << termcolor::green << "true";
            }
            else
            {
                std::cout << termcolor::magenta << "\n        is_column_vector = " << termcolor::green << "false";
            }
        }

        if (is_channel_range)
        {
            std::cout << termcolor::magenta << "\n        is_channel_range = " << termcolor::green << "true";
            std::cout << termcolor::magenta << "\n        channel_range_start = " << termcolor::green << channel_range_start;
            std::cout << termcolor::magenta << "\n        channel_range_end = " << termcolor::green << channel_range_end;
        }
        else
        {
            std::cout << termcolor::magenta << "\n        is_channel_range = " << termcolor::green << "false";
        }

        if (is_channel_stride)
        {
            std::cout << termcolor::magenta << "\n        is_channel_stride = " << termcolor::green << "true";
            std::cout << termcolor::magenta << "\n        channel_stride_size = " << termcolor::green << channel_stride_size;
        }
        else
        {
            std::cout << termcolor::magenta << "\n        is_channel_stride = " << termcolor::green << "false";
        }

        if (is_input_search_rgx)
        {
            std::cout << termcolor::magenta << "\n        input_search_rgx = " << termcolor::green << input_search_rgx;
        }
        std::cout << termcolor::magenta << "\n        n_files_to_concatenate = " << termcolor::green << n_files_to_concatenate;

        std::cout << termcolor::blue << "\n\n Runtime parameters: ";
        // std::cout << termcolor::magenta << "\n\n        lts_per_file = " << termcolor::green << lts_per_file;
        // std::cout << termcolor::magenta << "\n        chs_per_file = " << termcolor::green << chs_per_file;

        std::cout << termcolor::magenta << "\n        is_space_decimate = " << termcolor::green << is_space_decimate;

        std::cout << termcolor::magenta << "\n        space_decimate_chs = " << termcolor::green << space_decimate_chs;

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
