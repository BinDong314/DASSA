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

#include "ft.h"
#include "DasLib.h"

using namespace std;
using namespace AU;
using namespace DasLib;

struct timeval begin_time, end_time;

void printf_help(char *cmd);

int read_config_file(std::string file_name, int mpi_rank);
std::string config_file = "./template-match.config";

// Input Parameter
bool is_input_single_file = false;
std::string input_dir_file = "/Users/dbin/work/dassa/template-match/template-match-data";
std::string input_h5_dataset = "/Acoustic";
std::string input_file_type = "EP_HDF5";

std::string template_file_tsstart = "/Users/dbin/work/dassa/template-match/template_dir/";

bool is_input_search_rgx = false;
std::string input_search_rgx = "^(.*)[1234]\\.tdms$";

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
std::vector<std::vector<double>> template_tstart;
std::vector<std::vector<double>> template_weights; //[template index][channel]

std::vector<double> cheby1_b = {3.58632426674156e-09, 2.86905941339325e-08, 1.00417079468764e-07, 2.00834158937527e-07, 2.51042698671909e-07, 2.00834158937527e-07, 1.00417079468764e-07, 2.86905941339325e-08, 3.58632426674156e-09};
std::vector<double> cheby1_a = {1, -7.39772047094363, 24.0727277609670, -44.9989146659833, 52.8434667669224, -39.9159143524684, 18.9376658097605, -5.15917942637567, 0.617869501520363};

void init_xcorr()
{
    DT_NEW = decifac * DT;
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
        PrintScalar("5 / (nof1 * npts0 * dt0) = ", 5 / (nof1 * npts0 * dt0));
        PrintScalar("nof1 * npts0 = ", nof1 * npts0);
        PrintVector("ctap0 = ", ctap0);
    }

    // npts1=round((nof1*npts0)/decifac);
    npts1 = round((nof1 * npts0) / decifac);

    // ctap_template1=tukeywin((npts0*nlen_template),5/(npts0*nlen_template*dt0));
    tukeywin(ctap_template1, (npts0 * nlen_template), 5 / (npts0 * nlen_template * dt0));
    PrintVector("ctap_template1 = ", ctap_template1);

    npts1_template = round(npts0 * nlen_template / decifac);

    ntemplates = 1;
    template_data.resize(ntemplates);
    template_winlen.resize(ntemplates);
    template_tstart.resize(ntemplates);
    template_weights.resize(ntemplates);

    // tstart_ci39534271.txt
    AU::Array<double> *T_tsstart;
    T_tsstart = new AU::Array<double>("EP_DIR:EP_CSV:" + template_file_tsstart);

    std::vector<double> T_tstart_weight;
    std::vector<double> T_tstart;
    std::vector<double> T_weight;

    std::vector<std::string> control_para_ve, aug_merge_index, aug_input_search_rgx, file_size_str;
    std::vector<int> chunk_size_tsstart, overlap_size_tsstart = {0, 0};

    control_para_ve.push_back(std::to_string(CSV_SET_DELIMITER));
    control_para_ve.push_back(" ");
    T_tsstart->ControlEndpoint(DIR_SUB_CMD_ARG, control_para_ve);

    aug_input_search_rgx.push_back("(.*)tstart(.*)txt$");
    T_tsstart->EndpointControl(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx);

    T_tsstart->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str);
    String2Vector(file_size_str[0], chunk_size_tsstart);
    PrintVector("chunk_size_tsstart = ", chunk_size_tsstart);
    T_tsstart->SetChunkSize(chunk_size_tsstart);
    T_tsstart->SetOverlapSize(overlap_size_tsstart);

    // winlen_ci39534271.txt
    AU::Array<double> *T_winlen;
    T_winlen = new AU::Array<double>("EP_DIR:EP_CSV:" + template_file_tsstart);

    std::vector<double> T_winlen_data;
    std::vector<int> chunk_size_winlen, overlap_size_winlen = {0, 0};
    std::vector<std::string> aug_input_search_rgx_winlen, file_size_str_winlen;

    aug_input_search_rgx_winlen.push_back("(.*)winlen(.*)txt$"); // tstart_ci39534271.txt
    T_winlen->EndpointControl(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx_winlen);

    T_winlen->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str_winlen);
    String2Vector(file_size_str_winlen[0], chunk_size_winlen);
    PrintVector("chunk_size_winlen = ", chunk_size_winlen);
    T_winlen->SetChunkSize(chunk_size_winlen);
    T_winlen->SetOverlapSize(overlap_size_winlen);

    // ci39534271.h5
    // winlen_ci39534271.txt
    AU::Array<short> *T_h5;
    T_h5 = new AU::Array<short>("EP_DIR:EP_HDF5:" + template_file_tsstart + ":/Acoustic");

    int T_pts, T_chs;
    std::vector<short> T_h5_data;
    std::vector<int> chunk_size_h5, overlap_size_h5 = {0, 0};
    std::vector<std::string> aug_input_search_rgx_h5, file_size_str_h5;

    aug_input_search_rgx_h5.push_back("(.*)h5$"); // tstart_ci39534271.txt
    T_h5->EndpointControl(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx_h5);

    T_h5->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str_h5);
    String2Vector(file_size_str_h5[0], chunk_size_h5);
    PrintVector("chunk_size_h5 = ", chunk_size_h5);
    T_pts = chunk_size_h5[0];
    T_chs = chunk_size_h5[1];
    T_h5->SetChunkSize(chunk_size_h5);
    T_h5->SetOverlapSize(overlap_size_h5);

    template_winlen.clear();
    std::vector<std::vector<double>> T_ts2d; //[Channels][Points]
    double T_weight_sum = 0;
    size_t T_tstart_weight_size;
    for (int rc2 = 0; rc2 < ntemplates; rc2++)
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
        PrintVector("T_tstart = ", T_tstart);
        PrintVector("T_weight = ", T_weight);
        template_tstart[rc2] = T_tstart;
        template_weights[rc2] = T_weight;

        // Read template data and pro-process it
        T_h5->ReadNextChunk(T_h5_data);
        // PrintVector("T_h5_data = ", T_h5_data);

        tukeywin(ctap_template2, template_winlen[rc2], 0.04); // Put outside if template_winlen[rc2] is constent across tempaltes
        // PrintVector("ctap_template2 = ", ctap_template2);

        std::vector<double> atemp2;

        // std::cout << " template_tstart[" << rc2 << "].size() = " << template_tstart[rc2].size() << " \n";
        // std::cout << " template_weights[" << rc2 << "].size() = " << template_weights[rc2].size() << " \n";
        // std::cout << "T_chs = " << T_chs << " \n";

        T_ts2d = DasLib::Vector1D2DByColStride(T_chs, T_h5_data, 2, 3); // filter the data starting at ch (2-1) and every 3 chs
        // std::cout << "T_ts2d.size() = " << T_ts2d.size() << "\n";
        for (int i = 0; i < T_ts2d.size(); i++)
        {
            // detrend(T_ts2d[i].data(), T_pts); // Detread

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
            atemp2 = ddff(T_ts2d[i], ctap_template1, 10, BUTTER_A, BUTTER_B, cheby1_b, cheby1_a);
            //    PrintVector("After ddff atemp2 = ", atemp2);
            //   % SELECTING TEMPLATE-DEPENDENT WINDOW STARTING FROM CHANNEL
            //% DEPENDENT BEGIN TIME, DETRENDING, MULTPLYING WITH TAPER
            //  atemp3=detrend(atemp2(template_tstart(rc1,rc2):template_tstart(rc1,rc2)+template_winlen(rc2)-1)).*ctap_template2;
            slice(atemp2, template_tstart[rc2][i] - 1, template_tstart[rc2][i] - 1 + template_winlen[rc2] - 1);
            // std::cout << " atemp2.size() = " << atemp2.size() << ", template_tstart[rc2][i] = " << template_tstart[rc2][i] << ", template_winlen[rc2] = " << template_winlen[rc2] << ",  ctap_template2.size= " << ctap_template2.size() << "\n";
            detrend(atemp2.data(), atemp2.size());
            VectorElementMulti(atemp2, ctap_template2);
            double atemp2_norm = norm_matlab(atemp2);
            VectorDivideByScalar(atemp2, atemp2_norm);
            T_ts2d[i] = atemp2;
            // std::cout << " end for  channel #" << i << "\n";
        } // end for all channels of each template
        std::cout << " end for all channels of each template \n";
        // PrintVV("T_ts2d cha x points = ", T_ts2d);
        template_data[rc2] = T_ts2d;

        double template_tstart_min = *(std::min_element(template_tstart[rc2].begin(), template_tstart[rc2].end()));
        VectorMinusByScalar(template_tstart[rc2], template_tstart_min);
        PrintVector("template_tstart after norm = ", template_tstart[rc2]);
    } // end for each template
    std::cout << " end for each template\n";

    // std::cout << " template_data.size =" << template_data.size() << " , template_data[0].size = " << template_data[0].size() << " , template_data[0][0].size = " << template_data[0][0].size();
}

template <class TT>
inline Stencil<std::vector<double>> udf_template_match(const Stencil<TT> &iStencil)
{
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

    double template_tstart_max;           // temporary value for the tstart
    size_t npts2, dx1;                    // npts2 is the size of vector for cross-correlation
                                          // dx1 is the start index for the points
    std::vector<double> xc1;              // cross correlation per channel
    std::vector<std::vector<double>> xc0; // [template index][correlation]

    std::vector<std::vector<double>> amat1; // filter data , [channel][time points]
    std::vector<double> ts_temp2;           // temporary value for each ch during filter, try to remove it
    std::vector<double> sdcn_v;

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

    ts2d = DasLib::Vector1D2DByColStride(chs_per_file_udf, ts_short, 2, 3);
    if (!ft_rank)
        std::cout << "chs = " << ts2d.size() << ", each with " << ts2d[0].size() << " points\n";

    chs_per_file_udf = ts2d.size();
    lts_per_file_udf = ts2d[0].size();

    amat1.resize(chs_per_file_udf);
    // Resample in time-domain
    for (int ii = 0; ii < chs_per_file_udf; ii++)
    {
        // if (ii % 1000 == 0)
        //     std::cout << "ts2d[" << ii << " ], chs_per_file_udf = " << chs_per_file_udf << "\n";
        ts_temp2 = ddff(ts2d[ii], ctap0, 10, BUTTER_A, BUTTER_B, cheby1_b, cheby1_a);
        ts_temp2.pop_back();
        amat1[ii] = ts_temp2;
    }

    // PrintVV("amat1 =", amat1);
    //  ************************
    //  Cross with Template    *
    //  ************************
    xc1.resize(chs_per_file_udf);
    xc0.resize(ntemplates);
    nchan1 = chs_per_file_udf;
    for (int rc2 = 0; rc2 < ntemplates; rc2++)
    {
        template_tstart_max = *(std::max_element(std::begin(template_tstart[rc2]), std::end(template_tstart[rc2])));
        npts2 = npts1 - template_winlen[rc2] - template_tstart_max;
        // npts2=npts1-template_winlen(rc2)-max(template_tstart(:,rc2))+1; % [62182]
        // % VECTOR WITH CROSS-CORRELATION RESULTS
        //     xc0=zeros(1,npts2); %
        xc0[rc2].resize(npts2);
        // Points rc3=1:npts2
        for (int rc3 = 0; rc3 < npts2; rc3++)
        {
            // Channels rc1=1:nchan1
            for (int rc1 = 0; rc1 < nchan1; rc1++)
            {

                if (template_weights[rc2][rc1] > 0)
                {
                    // dx1=rc3+template_tstart(rc1,rc2);
                    dx1 = rc3 + template_tstart[rc2][rc1] + 1;
                    // atemp3=(detrend(amat1(idx1:(idx1+template_winlen(rc2)-1),rc1))).*ctap_template2;
                    // atemp3=atemp3./norm(atemp3);
                    // if (rc1 < 3)
                    // {
                    //     std::cout << "++dx1 =" << dx1 << " \n";
                    //     PrintVector("Before sdcn amat1[rc1] =", amat1[rc1]);
                    //     std::vector<double> debug_v(amat1[rc1].begin() + dx1, amat1[rc1].begin() + dx1 + template_winlen[rc2]);
                    //     PrintVector("amat1[rc1] debug_v =", debug_v);
                    // }
                    sdcn(amat1[rc1], sdcn_v, dx1, template_winlen[rc2], ctap_template2);

                    xc1[rc1] = dot_product(sdcn_v, template_data[rc2][rc1]);
                    // if (rc1 < 3)
                    // {
                    //     PrintVector("sdcn_v =", sdcn_v);
                    //     PrintVector("template_data[rc2][rc1]", template_data[rc2][rc1]);
                    //     std::cout << "xc1[rc1] = " << xc1[rc1] << "\n";
                    // }
                }
            }
            xc0[rc2][rc3] = sum_weight(xc1, template_weights[rc2]);
            // if (rc3 < 3)
            // {
            //     PrintVector(" * xc1 =", xc1);
            //     PrintVector(" * template_weights[rc2] =", template_weights[rc2]);
            //     std::cout << " * xc0[rc2][rc3] = " << xc0[rc2][rc3] << "\n\n\n";
            // }
            // else
            // {
            //     exit(0);
            // }
        }
    }

    // Set output
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

    gettimeofday(&begin_time, 0);

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

    // A->GetStencilTag();

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

    A->SetChunkSize(chunk_size);
    A->SetOverlapSize(overlap_size);
    A->DisableOverlapLower(); // Only have one extra data in upper side

    if (!ft_rank)
    {
        PrintVector("chunk_size = ", chunk_size);
        PrintVector("overlap_size = ", overlap_size);
    }

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
    TRANSFORM(A, udf_template_match, B, t, std::vector<double>);
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

    DT = reader.GetReal("parameter", "dt", 0.002);

    DT_NEW = reader.GetReal("parameter", "dt_new", 0.008);

    butter_order = reader.GetInteger("parameter", "butter_order", 3);

    if (!mpi_rank)
    {
        std::cout << "\n\n";
        std::cout << termcolor::red << "Parameters to run the Decimate: ";

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