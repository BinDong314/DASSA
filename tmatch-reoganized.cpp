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
#include <math.h>
#include <sys/time.h>
#include <deque>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "ft.h"
#include "DasLib.h"

using namespace std;
using namespace AU;
using namespace DasLib;

struct timeval begin_time, end_time;

void printf_help(char *cmd);

int read_config_file(std::string file_name, int mpi_rank);
std::string config_file = "./tmatch.config";

// Input Parameter
bool is_input_single_file = false;
std::string das_dir = "/Users/dbin/work/dassa/template-match/template-match-data";
std::string das_h5_dataset = "/Acoustic";
std::string das_file_type = "EP_HDF5";

std::string template_dir = "/Users/dbin/work/dassa/template-match/template_dir/";

bool is_das_input_search_rgx = false;
std::string das_input_search_rgx = "^(.*)[1234]\\.tdms$";

bool is_template_input_search_rgx = false;
std::string template_input_search_rgx = "(ci37327652|ci37329164)";

int n_files_to_concatenate = 20;

bool is_file_range = true;
int file_range_start_index = 0;
int file_range_end_index = 20;
std::string file_range_indexes_str;
std::vector<std::string> aug_merge_index;

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

// Output parameter
bool is_output_single_file = false;
std::string output_type = "EP_HDF5";
std::string output_file_dir = "./test-data/dir-output/";
std::string output_dataset = "/dat";

bool is_dir_output_match_replace_rgx = false;
std::string dir_output_match_rgx = "^(.*)\\.h5$";

std::string dir_output_replace_rgx = "$1-xcorr.h5";

std::string template_list_output_file = "template_final_list.txt";

bool is_space_decimate = false;
int space_decimate_chs = 32;

// Run-time paramter
double DT = 0.002;        // dt0
int chs_per_file = 6912;  // nchan0
int lts_per_file = 30000; // npts0
int &npts0 = lts_per_file;
int &nchan0 = chs_per_file;
double &dt0 = DT;
int npts1;

double taperwidth = 5;
int round_dt_dew_dt;
int decifac = 10;
double DT_NEW; // = decifac * DT
double &dt1 = DT_NEW;

//
vector<double> BUTTER_A;
vector<double> BUTTER_B;
int butter_order = 2;
std::vector<double> fbands = {0.5, 16};

bool is_channel_range = false;
int channel_range_start = 0;
int channel_range_end = 1;

bool is_channel_stride = true;
int channel_stride_start = 2;
int channel_stride_size = 3;

std::vector<double> bfsos1 = {1, 2, 1, 1, 0.502206458992975, 0.225842674724860, 1, -2, 1, 1, -1.91118566554443, 0.915114848946725};

std::vector<double> bfg1 = {0.412719392485785};

int nchan1;

double tsegment;
int nseg;
int nof1;

int nlen_template = 2; // UI
int npts1_template;
std::vector<double> ctap0, ctap_template1, ctap_template2;

int ntemplates;
std::vector<std::vector<std::vector<double>>> template_data; //[template index][channel][points]
std::vector<double> template_winlen;
std::vector<std::vector<double>> template_tstart;  //[template index][channel]
std::vector<std::vector<double>> template_weights; //[template index][channel]

std::vector<double> cheby1_b = {3.58632426674156e-09, 2.86905941339325e-08, 1.00417079468764e-07, 2.00834158937527e-07, 2.51042698671909e-07, 2.00834158937527e-07, 1.00417079468764e-07, 2.86905941339325e-08, 3.58632426674156e-09};
std::vector<double> cheby1_a = {1, -7.39772047094363, 24.0727277609670, -44.9989146659833, 52.8434667669224, -39.9159143524684, 18.9376658097605, -5.15917942637567, 0.617869501520363};

bool is_template_file_range = false;
int template_file_range_start_index = 0;
int template_file_range_end_index = 1;
std::string template_file_indexes_str;

int omp_num_threads_p = 32;

void all_gather_vector(const std::vector<double> &v_to_send, std::vector<double> &v_to_receive);

// bool is_correlation_method = false;
int correlation_method = 0; // correlation_method = 0 (dot-product), 1 (xcorr-max), 2 (fft-max)
#define CORR_DOT_PRODUCT 0
#define CORR_XCORR_MAX 1
#define CORR_FFT_MAX 2

void init_xcorr()
{
    std::vector<std::string> index_param;
    index_param.push_back(template_file_indexes_str);

    DT_NEW = decifac * DT;
    if (!ft_rank)
        PrintScalar("DT_NEW (dt1) = ", DT_NEW);

    // fbands=[0.5 16];
    //  fs1=1/DT_NEW; % NEW SAMPLING FREQUENCY 50 (from 500 to 50)
    //[bfz1,bfp1,bfk1]=butter(2,fbands./(fs1/2),'bandpass'); % bgz1: 1:1;-1-1]
    for (int i = 0; i < butter_order; i++)
    {
        fbands[i] = fbands[i] / ((1 / DT_NEW) / 2);
    }
    ButterPass(butter_order, fbands, BUTTER_A, BUTTER_B);

    if (!ft_rank)
    {
        PrintVector("fbands = ", fbands);
        PrintVector("BUTTER_A = ", BUTTER_A);
        PrintVector("BUTTER_B = ", BUTTER_B);
    }

    tsegment = 1200;                            // UI
    nseg = round(86400 / tsegment);             // IF EACH SEGMENT IS 1200 S, WE WILL PROCESS 72 SEGMENTS IN A DAY
                                                // 86400 = (24 hours *60 minutes *60 seconds)
    nof1 = round(tsegment / (dt0 * npts0)) + 1; // 21, NUMBER OF FILES IT WILL PROCESS AT ONE TIME, ADDING 1

    // ctap0=tukeywin((nof1*npts0),5/(nof1*npts0*dt0));
    tukeywin(ctap0, nof1 * npts0, 5 / (nof1 * npts0 * dt0));

    if (!ft_rank)
    {
        PrintScalar("taperwidth / (nof1 * npts0 * dt0) = ", 5 / (nof1 * npts0 * dt0));
        PrintScalar("nof1 * npts0 = ", nof1 * npts0);
        PrintVector("ctap0 = ", ctap0);
    }

    // npts1=round((nof1*npts0)/decifac);
    npts1 = round((nof1 * npts0) / decifac);

    // ctap_template1=tukeywin((npts0*nlen_template),5/(npts0*nlen_template*dt0));
    tukeywin(ctap_template1, (npts0 * nlen_template), 5 / (npts0 * nlen_template * dt0));
    if (!ft_rank)
        PrintVector("ctap_template1 = ", ctap_template1);

    npts1_template = round(npts0 * nlen_template / decifac);

    // tstart_ci39534271.txt
    AU::Array<double> *T_tsstart;
    T_tsstart = new AU::Array<double>("EP_DIR:EP_CSV:" + template_dir);

    std::vector<double> T_tstart_weight;
    std::vector<double> T_tstart;
    std::vector<double> T_weight;

    std::vector<std::string> control_para_ve, aug_merge_index, aug_input_search_rgx, file_size_str, null_str;
    std::vector<int> chunk_size_tsstart, overlap_size_tsstart = {0, 0};

    if (is_template_file_range)
    {
        T_tsstart->ControlEndpoint(DIR_FILE_SORT_INDEXES, index_param);
    }

    T_tsstart->ControlEndpoint(DIR_SKIP_SIZE_CHECK, null_str);

    control_para_ve.push_back(std::to_string(CSV_SET_DELIMITER));
    control_para_ve.push_back(" ");
    T_tsstart->ControlEndpoint(DIR_SUB_CMD_ARG, control_para_ve);

    if (!is_template_input_search_rgx)
    {
        aug_input_search_rgx.push_back("(.*)tstart(.*)txt$");
    }
    else
    {
        aug_input_search_rgx.push_back("(.*)tstart(.*)" + template_input_search_rgx + "(.*)txt$");
    }
    T_tsstart->EndpointControl(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx);

    T_tsstart->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str);
    String2Vector(file_size_str[0], chunk_size_tsstart);
    PrintVector("chunk_size_tsstart = ", chunk_size_tsstart);
    T_tsstart->SetChunkSize(chunk_size_tsstart);
    T_tsstart->SetOverlapSize(overlap_size_tsstart);
    T_tsstart->SetChunkSchedulingMethod(CHUNK_SCHEDULING_CR);

    // winlen_ci39534271.txt
    AU::Array<double> *T_winlen;
    T_winlen = new AU::Array<double>("EP_DIR:EP_CSV:" + template_dir);

    T_winlen->ControlEndpoint(DIR_SKIP_SIZE_CHECK, null_str);

    std::vector<double> T_winlen_data;
    std::vector<int> chunk_size_winlen, overlap_size_winlen = {0, 0};
    std::vector<std::string> aug_input_search_rgx_winlen, file_size_str_winlen;

    if (is_template_file_range)
    {
        T_winlen->ControlEndpoint(DIR_FILE_SORT_INDEXES, index_param);
    }

    if (!is_template_input_search_rgx)
    {
        aug_input_search_rgx_winlen.push_back("(.*)winlen(.*)txt$"); // tstart_ci39534271.txt
    }
    else
    {
        aug_input_search_rgx_winlen.push_back("(.*)winlen(.*)" + template_input_search_rgx + "(.*)txt$"); // tstart_ci39534271.txt
    }
    T_winlen->EndpointControl(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx_winlen);

    T_winlen->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str_winlen);
    String2Vector(file_size_str_winlen[0], chunk_size_winlen);
    PrintVector("chunk_size_winlen = ", chunk_size_winlen);
    T_winlen->SetChunkSize(chunk_size_winlen);
    T_winlen->SetOverlapSize(overlap_size_winlen);
    T_winlen->SetChunkSchedulingMethod(CHUNK_SCHEDULING_CR);

    // ci39534271.h5
    // winlen_ci39534271.txt
    AU::Array<short> *T_h5;
    T_h5 = new AU::Array<short>("EP_DIR:EP_HDF5:" + template_dir + ":/Acoustic");

    T_h5->ControlEndpoint(DIR_SKIP_SIZE_CHECK, null_str);

    int T_pts, T_chs;
    std::vector<short> T_h5_data;
    std::vector<int> chunk_size_h5, overlap_size_h5 = {0, 0};
    std::vector<std::string> aug_input_search_rgx_h5, file_size_str_h5;

    if (is_template_file_range)
    {
        T_h5->ControlEndpoint(DIR_FILE_SORT_INDEXES, index_param);
    }

    if (!is_template_input_search_rgx)
    {
        aug_input_search_rgx_h5.push_back("(.*)h5$"); // tstart_ci39534271.txt
    }
    else
    {
        aug_input_search_rgx_h5.push_back("(.*)" + template_input_search_rgx + "(.*)h5$"); // tstart_ci39534271.txt
    }
    T_h5->EndpointControl(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx_h5);

    T_h5->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str_h5);
    String2Vector(file_size_str_h5[0], chunk_size_h5);
    if (!ft_rank)
        PrintVector("chunk_size_h5 = ", chunk_size_h5);
    T_pts = chunk_size_h5[0];
    T_chs = chunk_size_h5[1];
    T_h5->SetChunkSize(chunk_size_h5);
    T_h5->SetOverlapSize(overlap_size_h5);

    int ntemplates_on_my_rank = 0;
    if (ft_size > 1)
    {
        T_h5->SetChunkSchedulingMethod(CHUNK_SCHEDULING_CR);
        unsigned long long my_chunk_start, my_chunk_end;
        T_h5->GetMyChunkStartEnd(my_chunk_start, my_chunk_end);
        if (!ft_rank)
            std::cout << "rank [" << ft_rank << "]: my_chunk_start = " << my_chunk_start << ", my_chunk_end = " << my_chunk_end << "\n";
        ntemplates_on_my_rank = my_chunk_end - my_chunk_start;
    }
    T_h5->ControlEndpoint(DIR_N_FILES, file_size_str_h5);
    ntemplates = std::stoi(file_size_str_h5[0]);

    int ntemplates_to_go;
    if (ft_size > 1)
    {
        ntemplates_to_go = ntemplates_on_my_rank;
    }
    else
    {
        ntemplates_to_go = ntemplates;
    }

    T_h5->ControlEndpoint(DIR_SAVE_FINAL_FILE_LIST, template_list_output_file);

    if (!ft_rank)
        std::cout << "Rank " << ft_rank << "ntemplates = " << ntemplates << ", ntemplates_to_go = " << ntemplates_to_go << std::endl;

    template_data.resize(ntemplates_to_go);
    template_winlen.resize(ntemplates_to_go);
    template_tstart.resize(ntemplates_to_go);
    template_weights.resize(ntemplates_to_go);

    template_winlen.clear();
    std::vector<std::vector<double>> T_ts2d; //[Channels][Points]
    double T_weight_sum = 0;
    size_t T_tstart_weight_size;
    for (int rc2 = 0; rc2 < ntemplates_to_go; rc2++)
    {
        // Read winlen data
        T_winlen->ReadNextChunk(T_winlen_data);
        // PrintVector("T_winlen_data = ", T_winlen_data);
        template_winlen.push_back(std::round(T_winlen_data[0] / dt1));
        // PrintScalar("template_winlen[0] = ", template_winlen[0]);

        // Read tstart_weight data and split it into tstart and weight
        T_tsstart->ReadNextChunk(T_tstart_weight);
        // PrintVector("T_tstart_weight = ", T_tstart_weight);
        T_tstart.clear();
        T_weight.clear();
        T_weight_sum = 0;
        T_tstart_weight_size = T_tstart_weight.size();
        for (int i = 0; i < T_tstart_weight_size; i = i + 2)
        {
            // template_tstart(:,rc2)=round((tstart1(:,1)+30)./dt1)+1;
            T_tstart.push_back(std::round((T_tstart_weight[i] + 30) / dt1) + 1);
            T_weight.push_back(T_tstart_weight[i + 1]);
            T_weight_sum = T_tstart_weight[i + 1] + T_weight_sum;
        }
        T_tstart_weight_size = T_tstart_weight_size / 2;
        for (int i = 0; i < T_tstart_weight_size; i++)
        {
            T_weight[i] = T_weight[i] / T_weight_sum;
        }
        // PrintVector("T_tstart = ", T_tstart);
        // PrintVector("T_weight = ", T_weight);
        template_tstart[rc2] = T_tstart;
        template_weights[rc2] = T_weight;

        // Read template data and pro-process it
        T_h5->ReadNextChunk(T_h5_data);
        // PrintVector("T_h5_data = ", T_h5_data);

        tukeywin(ctap_template2, template_winlen[rc2], 0.04); // Put outside if template_winlen[rc2] is constent across tempaltes
                                                              // PrintVector("ctap_template2 = ", ctap_template2);

        // std::cout << " template_tstart[" << rc2 << "].size() = " << template_tstart[rc2].size() << " \n";
        // std::cout << " template_weights[" << rc2 << "].size() = " << template_weights[rc2].size() << " \n";
        // std::cout << "T_chs = " << T_chs << " \n";

        T_ts2d = DasLib::Vector1D2DByColStride(T_chs, T_h5_data, 2, 3); // filter the data starting at ch (2-1) and every 3 chs
                                                                        // std::cout << "T_ts2d.size() = " << T_ts2d.size() << "\n";

        // PrintVV("T_ts2d of template = ", T_ts2d);

#if defined(_OPENMP)
#pragma omp parallel for
#endif
        for (int i = 0; i < T_ts2d.size(); i++)
        {
#if defined(_OPENMP)
            if ((!ft_rank && !omp_get_thread_num() && (!i)))
            {
                printf("Init Inside the OpenMP parallel region thread 0, we have %d threads, at template %d of MPI rank %d .\n", omp_get_num_threads(), i, ft_rank);
            }
#endif
            //  detrend(T_ts2d[i].data(), T_pts); // Detread
            // for (int j = 0; j < ctap_template1.size(); j++)
            //{
            //     T_ts2d[i][j] = T_ts2d[i][j] * ctap_template1[j];
            // }
            // T_ts2d[i].push_back(0); // add one zero
            //  decimate(T_ts2d[i].data());
            // decimate(T_ts2d[i], 10);

            // FILTER [bfsos1,bfg1]=zp2sos(bfz1,bfp1,bfk1);
            // [bfsos1,bfg1]=zp2sos(bfz1,bfp1,bfk1);
            // bfsos1 = [1,2,1,1,0.502206458992975,0.225842674724860;
            //        1,-2,1,1,-1.91118566554443,0.915114848946725]
            // bfg1 = [0.412719392485785]
            // atemp2=filtfilt(bfsos1,bfg1,atemp2);
            // filtfilt(bfsos1, bfg1, T_ts2d[i], atemp2);
            //  PrintVector("Before ddff T_ts2d[0] = ", T_ts2d[i]);
            //  Apply detrend, decimate, filtfilt on the data, 10 is decimate factor
            std::vector<double> atemp2;
            atemp2 = ddff(T_ts2d[i], ctap_template1, 10, BUTTER_A, BUTTER_B, cheby1_b, cheby1_a);
            //    PrintVector("After ddff atemp2 = ", atemp2);
            //   % SELECTING TEMPLATE-DEPENDENT WINDOW STARTING FROM CHANNEL
            //% DEPENDENT BEGIN TIME, DETRENDING, MULTPLYING WITH TAPER
            //  atemp3=detrend(atemp2(template_tstart(rc1,rc2):template_tstart(rc1,rc2)+template_winlen(rc2)-1)).*ctap_template2;
            slice(atemp2, template_tstart[rc2][i] - 1, template_tstart[rc2][i] - 1 + template_winlen[rc2] - 1);
            // std::cout << " atemp2.size() = " << atemp2.size() << ", template_tstart[rc2][i] = " << template_tstart[rc2][i] << ", template_winlen[rc2] = " << template_winlen[rc2] << ",  ctap_template2.size= " << ctap_template2.size() << "\n";
            detrend(atemp2.data(), atemp2.size());
            // VectorElementMulti(atemp2, ctap_template2);
            // double atemp2_norm = norm_matlab(atemp2);
            // VectorDivideByScalar(atemp2, atemp2_norm);
            VectorElementMultiNormal(atemp2, ctap_template2);
            T_ts2d[i] = atemp2;
            // std::cout << " end for  channel #" << i << "\n";
        } // end for all channels of each template
        if (!ft_rank)
            std::cout << " end for all channels of each template \n";
        // PrintVV("T_ts2d cha x points = ", T_ts2d);
        template_data[rc2] = T_ts2d;

        double template_tstart_min = *(std::min_element(template_tstart[rc2].begin(), template_tstart[rc2].end()));
        VectorMinusByScalar(template_tstart[rc2], template_tstart_min);
        // if (!ft_rank)
        //     PrintVector("template_tstart after norm = ", template_tstart[rc2]);
    } // end for each template
    if (!ft_rank)
        std::cout << " end for each template\n";

    if (ft_size > 1)
    {
        // we need to merge across nodes
        // std::vector<std::vector<std::vector<double>>> template_data; //[template index][channel][points]
        // std::vector<double> template_winlen;
        // std::vector<std::vector<double>> template_tstart;  //[template index][channel]
        // std::vector<std::vector<double>> template_weights; //[template index][channel]

        // std::cout << "Rank " << ft_rank << " template_data.size() = " << template_data.size() << "\n";
        // if (template_data.size() > 0)
        // {
        //     std::cout << "Rank " << ft_rank << ": template_data[0].size() = " << template_data[0].size() << "\n";
        // }

        int channels, max_channels, timepoints, max_timepoints;
        if ((!template_data.empty()) && (!template_data[0].empty()))
        {
            channels = template_data[0].size();
            timepoints = template_data[0][0].size();
        }

        MPI_Allreduce(&channels, &max_channels, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&timepoints, &max_timepoints, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        std::vector<double> template_data_1D = Convert3DVTo1DV(template_data);
        std::vector<double> template_data_1D_merged;
        all_gather_vector(template_data_1D, template_data_1D_merged);
        template_data = Convert1DVTo3DV(template_data_1D_merged, ntemplates, max_channels, max_timepoints);

        std::vector<double> template_tstart_1D = Convert2DVTo1DV(template_tstart);
        std::vector<double> template_tstart_1D_merged;
        all_gather_vector(template_tstart_1D, template_tstart_1D_merged);
        template_tstart = ConvertVector1DTo2D(template_tstart_1D_merged, ntemplates, max_channels, true);

        std::vector<double> template_weights_1D = Convert2DVTo1DV(template_weights);
        std::vector<double> template_weights_1D_merged;
        all_gather_vector(template_weights_1D, template_weights_1D_merged);
        template_weights = ConvertVector1DTo2D(template_weights_1D_merged, ntemplates, max_channels, true);

        std::vector<double> template_winlen_merged;
        all_gather_vector(template_winlen, template_winlen_merged);
        template_winlen = template_winlen_merged;
    }

    // std::cout << " template_data.size =" << template_data.size() << " , template_data[0].size = " << template_data[0].size() << " , template_data[0][0].size = " << template_data[0][0].size();

    // PrintVV("template_data[0] = ", template_data[0]);
}

template <class TT>
inline Stencil<std::vector<double>> udf_template_match(const Stencil<TT> &iStencil)
{
    double init_xcorr_t_start = AU_WTIME;

    // Input pramters
    std::vector<int> max_offset_upper; // Size of input data for the whole chunk, zero based
    iStencil.GetOffsetUpper(max_offset_upper);
    if (ft_rank == 0 || ft_rank == (ft_size - 1))
        PrintVector("max_offset_upper = ", max_offset_upper);
    int chs_per_file_udf, lts_per_file_udf; // lts_per_file_udf is the points per channel in the chunk
    std::vector<int> start_offset = {0, 0};
    std::vector<int> end_offset = {max_offset_upper[0], max_offset_upper[1]};
    std::vector<TT> ts_short;              // Buffer to read data in
    std::vector<std::vector<double>> ts2d; // Buffer to stored the converted from ts_short to 2D
                                           // ts2d[0] is a channel ....

    // double template_tstart_max;           // temporary value for the tstart
    // size_t dx1;                           // npts2 is the size of vector for cross-correlation
    // dx1 is the start index for the points
    std::vector<std::vector<double>> xc0; // [template index][correlation]

    std::vector<std::vector<double>> amat1; // filter data , [channel][time points]
    std::vector<double> ts_temp2;           // temporary value for each ch during filter, try to remove it

    // *************
    // Output Parameters
    // ****************
    std::vector<double> ts_temp;
    Stencil<std::vector<double>> oStencil;
    std::vector<size_t> vector_shape(2);

    // ***********************************
    // Start to read data in
    // Here, we assume data is row-vector read from HDF5
    // We can transter it to be used column-major
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
        PrintVector("start_offset after channel_range  = ", start_offset);
        PrintVector("end_offset after channel_range = ", end_offset);
    }

    if (lts_per_file_udf != (nof1 * npts0))
    {

        std::cout << "Skip this chunk, it has = " << lts_per_file_udf << " points, but we initialize for [" << nof1 * npts0 << "]points \n";
        // goto ReturnEmpty;
        oStencil.SetEmpty();
        if (oStencil.IsEmpty())
        {
            std::cout << "I am IsEmpty \n";
        }
        return oStencil;
    }

    iStencil.ReadNeighbors(start_offset, end_offset, ts_short);

    if (!ft_rank)
        std::cout << "ReadNeighbors (s) = " << AU_WTIME - init_xcorr_t_start << std::endl;
    init_xcorr_t_start = AU_WTIME;

    ts2d = DasLib::Vector1D2DByColStride(chs_per_file_udf, ts_short, 2, 3);
    if (!ft_rank || ft_rank == (ft_size - 1))
        std::cout << "chs = " << ts2d.size() << ", each with " << ts2d[0].size() << " points\n";

    // PrintVV("ts2d = ", ts2d);
    if (!ft_rank)
        std::cout << "Vector1D2DByColStride (s) = " << AU_WTIME - init_xcorr_t_start << std::endl;
    init_xcorr_t_start = AU_WTIME;

    chs_per_file_udf = ts2d.size();
    lts_per_file_udf = ts2d[0].size();

    amat1.resize(chs_per_file_udf);

    int npts1_new = round(((nof1 - 1) * npts0) / decifac) + round(taperwidth / dt1) + MaxVector(template_winlen) + MaxVectorVector(template_tstart) + 1;
    if (npts1_new < npts1)
        npts1 = npts1_new;
// npts1 = min([npts1 npts1_new]) ;

// if (!ft_rank)
// {
//     std::cout << "round(((nof1 - 1) * npts0) / decifac) = " << round(((nof1 - 1) * npts0) / decifac)
//               << ", round(taperwidth / dt1) = " << round(taperwidth / dt1)
//               << ", MaxVector(template_winlen) = " << MaxVector(template_winlen)
//               << ", MaxVectorVector(template_tstart) = " << MaxVectorVector(template_tstart) << "\n";
//     PrintVV("template_tstart = ", template_tstart);
// }

// PrintVV("ts2d  of das data ", ts2d);

// Resample in time-domain
#if defined(_OPENMP)
#pragma omp parallel for firstprivate(ts_temp2)
#endif
    for (int ii = 0; ii < chs_per_file_udf; ii++)
    {
        // if (ii % 1000 == 0)
        //     std::cout << "ts2d[" << ii << " ], chs_per_file_udf = " << chs_per_file_udf << "\n";
        ts_temp2 = ddff(ts2d[ii], ctap0, 10, BUTTER_A, BUTTER_B, cheby1_b, cheby1_a);
        // ts_temp2.pop_back();
        ts_temp2.resize(npts1); // get rid of the last one
        amat1[ii] = ts_temp2;
    }

    // std::cout << " amat1.size =" << amat1.size() << " , amat1[0].size = " << amat1[0].size() << " \n";
    // PrintVV("amat1  ", amat1);

    // PrintVV("template_data[0] = ", template_data[0]);

    if (!ft_rank)
        std::cout << "ddff (s) = " << AU_WTIME - init_xcorr_t_start << std::endl;
    init_xcorr_t_start = AU_WTIME;

    // PrintVV("amat1 =", amat1);
    //  ************************
    //  Cross with Template    *
    //  ************************
    xc0.resize(ntemplates);
    nchan1 = chs_per_file_udf;
    // double micro_init_xcorr_t_start = AU_WTIME;

    size_t npts2_max = 0, temp_value;
    std::vector<size_t> npts2_vector;
    npts2_vector.resize(ntemplates, 0);
    double template_tstart_max;
    for (int rc2 = 0; rc2 < ntemplates; rc2++)
    {
        temp_value = MaxVector(template_tstart[rc2]); //*(std::max_element(std::begin(template_tstart[rc2]), std::end(template_tstart[rc2])));
        npts2_vector[rc2] = npts1 - template_winlen[rc2] - temp_value + 1;
        if (npts2_vector[rc2] > npts2_max)
        {
            npts2_max = npts2_vector[rc2];
        }
    }

    for (int rc2 = 0; rc2 < ntemplates; rc2++)
    {
        xc0[rc2].resize(npts2_max, 0);
    }

    if (!ft_rank)
        std::cout << "npts1 = " << npts1 << ", npts2_max = " << npts2_max << ", chs_per_file_udf =" << chs_per_file_udf << ", nchan1 = " << nchan1 << ", mpi_rank =" << ft_rank << ", ntemplates = " << ntemplates << ", nchan1 = " << nchan1 << ", npts2_vector[0] = " << npts2_vector[0] << ", npts1_new = " << npts1_new << "\n";

    // std::deque<double> ch_window_buffer(template_winlen[0], 0); firstprivate(ch_window_buffer)
    //  #if defined(_OPENMP)
    //  #endif
    ////#pragma omp parallel for
    for (int rc2 = 0; rc2 < ntemplates; rc2++)
    {
        double micro_init_xcorr_t_start = AU_WTIME;
        // std::vector<std::vector<double>> xc_channel_time;
        // xc_channel_time.resize(nchan1);
        // https://stackoverflow.com/questions/15349695/pre-allocated-private-stdvector-in-openmp-parallelized-for-loop-in-c
        // #if defined(_OPENMP)
        // #pragma omp parallel for schedule(static, 1)
        // #endif
        for (int rc1 = 0; rc1 < nchan1; rc1++)
        {
            if (template_weights[rc2][rc1] > 0)
            {

                // cross correlation per channel Channels rc1=1:nchan1
                // std::vector<double> xc1(npts2_vector[rc2], 0);
                double xmean, Sxx = 0, xsum;
                xmean = (template_winlen[rc2] - 1) / 2;
                xsum = (template_winlen[rc2] - 1) * template_winlen[rc2] / 2;
                for (int iiii = 0; iiii < template_winlen[rc2]; iiii++)
                {
                    Sxx += (iiii - xmean) * (iiii - xmean);
                }
                // #if defined(_OPENMP)
                // #pragma omp parallel {
                // #endif
                std::vector<double> sdcn_v(template_winlen[rc2], 0);
                size_t dx1;
                // https://stackoverflow.com/questions/15349695/pre-allocated-private-stdvector-in-openmp-parallelized-for-loop-in-c
#if defined(_OPENMP)
#pragma omp parallel for firstprivate(sdcn_v, dx1)
#endif
                for (int rc3 = 0; rc3 < npts2_vector[rc2]; rc3++)
                {
                    dx1 = rc3 + template_tstart[rc2][rc1];
                    // Replace below line with the following to find difference of two version
                    detrend_range_one_pass_std(amat1[rc1], dx1, template_winlen[rc2], ctap_template2, xmean, xsum, Sxx, sdcn_v);
                    // sdcn(amat1[rc1], sdcn_v, dx1, template_winlen[rc2], ctap_template2);
                    // double temp_xcorr;
                    switch (correlation_method)
                    {
                    case CORR_DOT_PRODUCT:
                        xc0[rc2][rc3] = xc0[rc2][rc3] + template_weights[rc2][rc1] * dot_product(sdcn_v, template_data[rc2][rc1]);
                        break;
                    case CORR_XCORR_MAX:
                        xc0[rc2][rc3] = xc0[rc2][rc3] + template_weights[rc2][rc1] * xcross_max(sdcn_v, template_data[rc2][rc1]);
                        break;
                    case CORR_FFT_MAX:
                        if (!rc3 && !rc1)
                        {
                            PrintVector("sdcn_v: ", sdcn_v);
                            PrintVector("template_data[rc2][rc1]: ", template_data[rc2][rc1]);
                            std::cout << "xcross_fft = " << xcross_fft(sdcn_v, template_data[rc2][rc1]) << "\n";
                        }

                        xc0[rc2][rc3] = xc0[rc2][rc3] + template_weights[rc2][rc1] * xcross_fft(sdcn_v, template_data[rc2][rc1]);
                        break;
                    default:
                        std::cout << "Unsupported Correlation Method code " << correlation_method << "\n";
                        std::cout << "Please set correlation_method = 0 (dot-product), 1 (xcorr-max), 2 (fft-max) " << correlation_method << "\n";

                        exit(-1);
                        break;
                    }
                    // xc0[rc2][rc3] = xc0[rc2][rc3] + template_weights[rc2][rc1] * temp_xcorr;
                }
                // #if defined(_OPENMP)
                //             }
                // #endif
            }
        }
        // Stack of all channels at time rc3 [template index][time] for template rc2
        // sum_weight_by_time(xc_channel_time, template_weights[rc2], xc0[rc2]);
        // xc_channel_time[rc1] = xc1;
    }

    if (!ft_rank)
        std::cout << "sdcn (for loop of all templates ) (s) = " << AU_WTIME - init_xcorr_t_start << std::endl;
    init_xcorr_t_start = AU_WTIME;

    ts_temp = Convert2DVTo1DV(xc0);
    if (is_column_major)
    {
        std::vector<double> ts_temp_column;
        ts_temp_column.resize(ts_temp.size());
        transpose(ts_temp.data(), ts_temp_column.data(), xc0.size(), xc0[0].size());
        ts_temp = ts_temp_column;
        vector_shape[1] = xc0.size();
        vector_shape[0] = xc0[0].size();
    }
    else
    {
        vector_shape[0] = xc0.size();
        vector_shape[1] = xc0[0].size();
    }

    if (ft_rank == 0 || ft_rank == (ft_size - 1))
        PrintVector("Output vector_shape: ", vector_shape);

    DasLib::clear_vector(xc0);
    oStencil.SetShape(vector_shape);
    oStencil = ts_temp;
    // exit(0);
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

    // MPI_Get_processor_name(processor_name, &namelen);
    // printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

    if (has_config_file_flag)
        read_config_file(config_file, ft_rank);

        // char processor_name[MPI_MAX_PROCESSOR_NAME];
        // int numprocs, rank, namelen;
        // MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#if defined(_OPENMP)
    // unsigned int thread_qty = atoi(std::getenv("OMP_NUM_THREADS"));
    omp_set_num_threads(omp_num_threads_p);
    if (!ft_rank)
        std::cout << "Set [omp_set_num_threads] to " << omp_num_threads_p << "\n";
#endif

    /*if (correlation_method == 2)
    {
        fftw_make_planner_thread_safe();
    }*/

    gettimeofday(&begin_time, 0);

    // set up the chunk size and the overlap size
    // 11648, 30000 for each dataset
    std::vector<int> chunk_size(2);
    std::vector<int> overlap_size = {0, 0};

    // std::cout << "EP_DIR:" + das_file_type + ":" + das_dir << ":" << das_h5_dataset << "\n";
    std::string A_endpoint_id;

    if (!is_input_single_file)
    {
        A_endpoint_id = "EP_DIR:" + das_file_type + ":" + das_dir + ":" + das_h5_dataset;
    }
    else
    {
        A_endpoint_id = das_file_type + ":" + das_dir + ":" + das_h5_dataset;
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

    // A->GetStencilTag();

    std::vector<std::string> null_str;
    if (!is_input_single_file)
    {
        A->SkipFileTail();
        A->ExecuteUDFOnce();
        A->ControlEndpoint(DIR_SKIP_SIZE_CHECK, null_str);
    }
    if (is_das_input_search_rgx && is_input_single_file == false)
    {
        std::vector<std::string> aug_input_search_rgx;
        aug_input_search_rgx.push_back(das_input_search_rgx);
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
        if (!ft_rank)
            std::cout << " file_range_indexes_str =" << file_range_indexes_str << "\n";

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
    if (das_file_type == "EP_TDMS")
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

    if (!is_input_single_file)
    {
        if (is_column_major)
        {
            if (is_channel_range)
            {
                chunk_size[1] = channel_range_end - channel_range_start + 1;
                // A->SetView(channel_range_start, chunk_size[1], 1);
            }
            overlap_size[1] = 0;
            overlap_size[0] = chunk_size[0];
            chunk_size[0] = chunk_size[0] * n_files_to_concatenate;
            chs_per_file = chunk_size[1];
            // lts_per_file = chunk_size[0];
        }
        else
        {
            if (is_channel_range)
            {
                chunk_size[0] = channel_range_end - channel_range_start + 1;
                // A->SetView(channel_range_start, chunk_size[0], 0);
            }
            overlap_size[1] = chunk_size[0];
            overlap_size[0] = 0;
            chunk_size[1] = chunk_size[1] * n_files_to_concatenate;
            chs_per_file = chunk_size[0];
            // lts_per_file = chunk_size[1];
        }
    }
    else
    {
        overlap_size[1] = 0;
        overlap_size[0] = 0;
    }
    A->SetChunkSize(chunk_size);
    A->SetOverlapSize(overlap_size);

    if (!is_input_single_file)
    {
        A->DisableOverlapLower(); // Only have one extra data in upper side
    }
    if (!ft_rank)
    {
        PrintVector("chunk_size = ", chunk_size);
        PrintVector("overlap_size = ", overlap_size);
    }

    std::vector<int> tempalte_file_indexes;
    if (is_template_file_range)
    {
        for (int i = template_file_range_start_index; i <= template_file_range_end_index; i++)
        {
            tempalte_file_indexes.push_back(i);
        }
    }
    template_file_indexes_str = Vector2String(tempalte_file_indexes);

    double init_xcorr_t_start = AU_WTIME;
    init_xcorr();
    if (!ft_rank)
        std::cout << "init_xcorr time (s) = " << AU_WTIME - init_xcorr_t_start << std::endl;
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
    TRANSFORM_NO_MP(A, udf_template_match, B, t, std::vector<double>);
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

    template_dir = reader.Get("parameter", "template_dir", "/Users/dbin/work/dassa/template-match/template_dir/");

    das_dir = reader.Get("parameter", "das_dir", "/Users/dbin/work/dassa/template-match/template-match-data");

    das_h5_dataset = reader.Get("parameter", "input_dataset", "/Acoustic");

    das_file_type = reader.Get("parameter", "das_file_type", "EP_HDF5");

    std::string temp_str_str;
    temp_str_str = reader.Get("parameter", "das_data_type", "short");
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

    temp_str = reader.Get("parameter", "is_das_input_search_rgx", "false");
    is_das_input_search_rgx = (temp_str == "false") ? false : true;

    if (is_das_input_search_rgx)
    {
        das_input_search_rgx = reader.Get("parameter", "das_input_search_rgx", "^(.*)[1234]\\.tdms$");
    }

    temp_str = reader.Get("parameter", "is_template_input_search_rgx", "false");
    is_template_input_search_rgx = (temp_str == "false") ? false : true;

    if (is_template_input_search_rgx)
    {
        template_input_search_rgx = reader.Get("parameter", "template_input_search_rgx", "(ci37327652|ci37329164)");
    }

    std::string is_file_range_str = reader.Get("parameter", "is_das_file_range", "false");
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
        file_range_start_index = reader.GetInteger("parameter", "das_file_range_start_index", 0);
        file_range_end_index = reader.GetInteger("parameter", "das_file_range_end_index", 1);
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

    correlation_method = reader.GetInteger("parameter", "correlation_method", 0);
    if (correlation_method < 0 || correlation_method > 2)
    {
        AU_EXIT("Don't understand the correlation_method's value " + correlation_method);
    }

    std::string is_template_file_range_str = reader.Get("parameter", "is_template_file_range", "false");
    if (is_template_file_range_str == "false" || is_template_file_range_str == "0")
    {
        is_template_file_range = false;
    }
    else if (is_template_file_range_str == "true" || is_template_file_range_str == "1")
    {
        is_template_file_range = true;
    }
    else
    {
        AU_EXIT("Don't read the is_template_file_range's value " + is_template_file_range);
    }

    if (is_template_file_range)
    {
        template_file_range_start_index = reader.GetInteger("parameter", "template_file_range_start_index", 0);
        template_file_range_end_index = reader.GetInteger("parameter", "template_file_range_end_index", 1);
    }

    MeasureLengthName = reader.Get("parameter", "attribute_name_measure_length", "MeasureLength[m]");
    SpatialResolutionName = reader.Get("parameter", "attribute_name_spatial_resolution", "SpatialResolution[m]");
    SamplingFrequencyName = reader.Get("parameter", "attribute_name_sampling_frequency", "SamplingFrequency[Hz]");

    taperwidth = reader.GetInteger("parameter", "taperwidth", 5);

    n_files_to_concatenate = reader.GetInteger("parameter", "n_files_to_concatenate", 1);

    temp_str = reader.Get("parameter", "is_output_single_file", "false");

    is_output_single_file = (temp_str == "false") ? false : true;

    output_type = reader.Get("parameter", "output_type", "EP_HDF5");

    template_list_output_file = reader.Get("parameter", "template_list_output_file", "template_final_list.txt");

    output_file_dir = reader.Get("parameter", "output_file_dir", "./tdms-dir-dec/test.h5");

    output_dataset = reader.Get("parameter", "output_dataset", "/data");

    temp_str = reader.Get("parameter", "is_dir_output_match_replace_rgx", "false");

    is_dir_output_match_replace_rgx = (temp_str == "false") ? false : true;

    if (is_dir_output_match_replace_rgx)
    {
        dir_output_match_rgx = reader.Get("parameter", "output_file_regex_match", "^(.*)\\.tdms$");

        dir_output_replace_rgx = reader.Get("parameter", "output_file_regex_replace", "$1.h5");
    }

    DT = reader.GetReal("parameter", "dt", 0.002);

    decifac = reader.GetReal("parameter", "decifac", 10);

    butter_order = reader.GetInteger("parameter", "butter_order", 2);

    omp_num_threads_p = reader.GetInteger("parameter", "omp_num_threads", 32);
    if (omp_num_threads_p < 0)
    {
        AU_EXIT("omp_num_threads must be positive integer : " + std::to_string(omp_num_threads_p));
    }

    // fbands

    temp_str = reader.Get("parameter", "is_output_single_file", "false");
    if (temp_str == "true")
    {
        fbands.resize(1);
        fbands[0] = reader.GetReal("parameter", "fbands_low", 0.5);
        fbands[1] = reader.GetReal("parameter", "fbands_high", 16);
    }

    if (!mpi_rank)
    {
        std::cout << "\n\n";
        std::cout << termcolor::red << "Parameters to run the Decimate: ";

        // template_dir

        std::cout << termcolor::blue << "\n\n Input parameters: ";
        std::cout << termcolor::magenta << "\n        template_dir = " << termcolor::green << template_dir;

        std::cout << termcolor::magenta << "\n        das_dir = " << termcolor::green << das_dir;
        std::cout << termcolor::magenta << "\n        das_file_type = " << termcolor::green << das_file_type;
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

        std::cout << termcolor::magenta << "\n        is_template_file_range = " << termcolor::green << is_template_file_range;
        if (is_template_file_range)
        {
            std::cout << termcolor::magenta << "\n        template_file_range_start_index = " << termcolor::green << template_file_range_start_index;
            std::cout << termcolor::magenta << "\n        template_file_range_end_index = " << termcolor::green << template_file_range_end_index;
        }

        std::cout << termcolor::magenta << "\n        is_das_file_range = " << termcolor::green << is_file_range;
        if (is_file_range)
        {
            std::cout << termcolor::magenta << "\n        das_file_range_start_index = " << termcolor::green << file_range_start_index;
            std::cout << termcolor::magenta << "\n        das_file_range_end_index = " << termcolor::green << file_range_end_index;
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

        if (is_das_input_search_rgx)
        {
            std::cout << termcolor::magenta << "\n        das_input_search_rgx = " << termcolor::green << das_input_search_rgx;
        }
        std::cout << termcolor::magenta << "\n        n_files_to_concatenate = " << termcolor::green << n_files_to_concatenate;

        std::cout << termcolor::blue << "\n\n Runtime parameters: ";
        // std::cout << termcolor::magenta << "\n\n        lts_per_file = " << termcolor::green << lts_per_file;
        // std::cout << termcolor::magenta << "\n        chs_per_file = " << termcolor::green << chs_per_file;

        std::cout << termcolor::magenta << "\n        DT = " << termcolor::green << DT;

        std::cout << termcolor::magenta << "\n        decifac = " << termcolor::green << decifac;
        std::cout << termcolor::magenta << "\n        OpenMP_num_threads = " << termcolor::green << omp_num_threads_p;
        switch (correlation_method)
        {
        case CORR_DOT_PRODUCT:
            std::cout << termcolor::magenta << "\n        correlation_method = " << termcolor::green << CORR_DOT_PRODUCT << " dot_product";
            break;
        case CORR_XCORR_MAX:
            std::cout << termcolor::magenta << "\n        correlation_method = " << termcolor::green << CORR_XCORR_MAX << " xcorr_max";
            break;
        case CORR_FFT_MAX:
            std::cout << termcolor::magenta << "\n        correlation_method = " << termcolor::green << CORR_FFT_MAX << " fft_max";
            break;
        default:
            break;
        }

        std::cout << termcolor::blue << "\n\n Output parameters: ";

        std::cout << termcolor::magenta << "\n    template_list_output_file = " << termcolor::green << template_list_output_file;
        std::cout << termcolor::magenta << "\n        is_output_single_file = " << termcolor::green << is_output_single_file;
        std::cout << termcolor::magenta << "\n        output_type = " << termcolor::green << output_type;
        std::cout << termcolor::magenta << "\n        output_file_dir = " << termcolor::green << output_file_dir;
        std::cout << termcolor::magenta << "\n        output_dataset = " << termcolor::green << output_dataset;
        std::cout << termcolor::magenta << "\n        taperwidth     = " << termcolor::green << taperwidth;

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

void all_gather_vector(const std::vector<double> &v_to_send, std::vector<double> &v_to_receive)
{
    // https://rookiehpc.github.io/mpi/docs/mpi_allgatherv/index.html
    int local_sum = v_to_send.size();

    std::vector<int> local_sum_vector;
    local_sum_vector.resize(ft_size);
    MPI_Allgather(&local_sum, 1, MPI_INT, local_sum_vector.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int global_sum = 0;
    std::vector<int> displacements;
    displacements.resize(ft_size);
    displacements[0] = 0;
    for (int i = 0; i < local_sum_vector.size(); i++)
    {
        global_sum = global_sum + local_sum_vector[i];
        if (i > 0)
            displacements[i] = displacements[i - 1] + local_sum_vector[i - 1];
    }
    v_to_receive.resize(global_sum);
    MPI_Allgatherv(v_to_send.data(), local_sum, MPI_DOUBLE, v_to_receive.data(), local_sum_vector.data(), displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
}
