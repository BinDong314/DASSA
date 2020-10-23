/*
 *ArrayUDF Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

 *If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

* NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so. 
 */
/**
 *
 * Email questions to dbin@lbl.gov
 * Scientific Data Management Research Group
 * Lawrence Berkeley National Laboratory
 *
 */

#include <iostream>
#include <stdarg.h>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include "ft.h"
#include "DasLib.h"
#include <stdlib.h>

using namespace std;
using namespace FT;

/**
 * @brief some help function at the end
 * 
 * @param cmd 
 */
void printf_help(char *cmd);
int stack_config_reader(std::string file_name, int mpi_rank);

extern int au_size;
extern int au_rank;

int lts_per_file = 14999;
int chs_per_file = 201;

std::string config_file = "./stack.config";

//Output file name
std::string xcorr_input_dir = "/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/stacking_files/xcorr_examples_h5";
//"/clusterfs/bear/BinDong_DAS_Data/xcorr_examples_h5/";
std::string xcorr_input_dataset_name = "/xcoor";

std::string stack_output_dir = "./";

//Output file name
std::string stack_output_file_data_in_sum_name = "xcorr_examples_h5_stack_data_in_sum.h5";
std::string stack_output_file_final_pwstack = "xcorr_examples_h5_stack_final_pwstack.h5";
std::string stack_output_file_phaseWeight = "xcorr_examples_h5_stack_phaseWeight.h5";
std::string stack_output_file_semblanceWeight = "xcorr_examples_h5_stack_semblanceWeight.h5";
std::string stack_output_file_semblance_denom_sum = "xcorr_examples_h5_stack_semblance_denom_sum.h5";
std::string stack_output_file_dataset_name = "/data";

double t_start = -59.9920000000000, t_end = 59.9920000000000, sample_rate = 0.00800000000000000;
double sub_start_t = -59, sub_end_t = 59;
double pow_u = 0.3;
double CausalityFlagging_tmin = 0.05;
double CausalityFlagging_tmax = 3.0;
double CausalityFlagging_fmax = 10;
double CausalityFlagging_ButterLow_order = 3;
double CausalityFlagging_ButterLow_fcf = 0.16; //filter cutoff frequency

bool is_flipud_flag = true;

std::vector<unsigned long long> sc_size(2);
FT::Array<double> *semblance_denom_sum;
FT::Array<std::complex<double>> *coherency_sum;
FT::Array<double> *data_in_sum;
FT::Array<double> *phaseWeight;
FT::Array<double> *semblanceWeight;

FT::Array<double> *final_pwstack;

double nStack = 0;

double IO_Read_Time = 0, IO_Reduce_Time = 0, CPU_Time = 0, IO_Read_Time_Con = 0, DetMean_tim = 0, Subset_tim = 0, CausalityFlagging_tim = 0, instanPhaseEstimator_Time = 0, instanPhaseEstimator_Time2 = 0, sum_Time = 0, sum_micro = 0, sum_micro_sub = 0;
double temp_time, temp_time_large, micro_time, micro_time_sub;

bool is_ml_weight = false;
bool is_ml_weight_ordered = true;
std::string ml_weight_file = "ml_wight.txt";
int n_weighted_to_stack = -1; // -1 mean to stack all
std::vector<double> ml_weight;
double ml_weight_sum;
std::vector<size_t> sorted_indexes; //sorted index
std::string sorted_indexes_str;     //string of sort_indexes after cut to n_weighted_to_stack

void read_ml_weight(const std::string ml_weight_file_p, std::vector<double> &ml_weight_p, std::vector<size_t> &sort_indexes_p);

inline Stencil<double>
stack_udf(const Stencil<double> &iStencil)
{
    nStack++;

    std::vector<int> start_offset{0, 0}, end_offset{chs_per_file - 1, lts_per_file - 1};
    std::vector<double> ts;
    iStencil.ReadNeighbors(start_offset, end_offset, ts);
    std::vector<std::vector<double>> ts2d = DasLib::Vector1D2D(lts_per_file, ts);

    //PrintVV("ts2d ", ts2d);
    /*
    std::vector<std::vector<double>> ts2d;
    ts2d.resize(chs_per_file);
    for (int i = 0; i < chs_per_file; i++)
    {
        ts2d[i].resize(lts_per_file);
        for (int j = 0; j < lts_per_file; j++)
        {
            ts2d[i][j] = iStencil(i, j);
        }
    }
    */

    DasLib::DeleteMedian(ts2d);

    DetMean_tim = DetMean_tim + (AU_WTIME - temp_time);
    temp_time = AU_WTIME;

    //Remove the media
    for (int i = 0; i < chs_per_file; i++)
    {
        //Subset
        ts2d[i] = DasLib::TimeSubset(ts2d[i], t_start, t_end, sub_start_t, sub_end_t, sample_rate);
    }

    if (is_flipud_flag)
    {
        bool flag = DasLib::CausalityFlagging(ts2d, CausalityFlagging_tmin, CausalityFlagging_tmax, CausalityFlagging_fmax, sub_start_t, sub_end_t, sample_rate, CausalityFlagging_ButterLow_order, CausalityFlagging_ButterLow_fcf);
        if (flag == false)
        {
            std::cout << "flipud  for the " << nStack << "th file  at process rank " << au_rank << "\n";
            for (int i = 0; i < chs_per_file; i++)
            {
                std::reverse(ts2d[i].begin(), ts2d[i].end());
            }
        }
    }

    //PrintVV("ts2d after flipud ", ts2d);

    size_t chunk_chunk_id = iStencil.GetChunkID();

    //cout << "Current Chunk ID: =" << chunk_chunk_id << ", ml weight: " << ml_weight[chunk_chunk_id] << "\n";

    size_t LTS_new = ts2d[0].size();
    std::vector<std::vector<double>> semblance_denom;
    std::vector<std::vector<std::complex<double>>> coherency;
    semblance_denom.resize(chs_per_file);
    //coherency.resize(chs_per_file);
    for (int i = 0; i < chs_per_file; i++)
    {
        semblance_denom[i].resize(LTS_new);
        //coherency[i].resize(LTS_new);
        for (int j = 0; j < LTS_new; j++)
        {
            if (is_ml_weight)
            {
                ts2d[i][j] = ts2d[i][j] * ml_weight[chunk_chunk_id];
            }

            semblance_denom[i][j] = ts2d[i][j] * ts2d[i][j];
        }
        //coherency[i] = DasLib::instanPhaseEstimator(ts2d[i]);
    }
    coherency = DasLib::instanPhaseEstimatorVector(ts2d);

    //PrintVV("ts2d after add ml weight ", ts2d);

    CPU_Time = CPU_Time + (AU_WTIME - temp_time_large);
    temp_time_large = AU_WTIME;

    std::vector<unsigned long long> H_start{0, 0}, H_end{static_cast<unsigned long long>(chs_per_file) - 1, static_cast<unsigned long long>(LTS_new) - 1};

    std::vector<double> semblance_denom_sum_v;
    semblance_denom_sum->ReadArray(H_start, H_end, semblance_denom_sum_v);

    std::vector<std::complex<double>> coherency_sum_v;
    coherency_sum->ReadArray(H_start, H_end, coherency_sum_v);

    std::vector<double> data_in_sum_v;
    data_in_sum->ReadArray(H_start, H_end, data_in_sum_v);

    //PrintVector("coherency_sum_v (before) " + std::to_string(au_rank), data_in_sum_v);
    int offset;
    for (int i = 0; i < chs_per_file; i++)
    {
        for (int j = 0; j < LTS_new; j++)
        {
            offset = i * LTS_new + j;
            coherency_sum_v[offset] = coherency_sum_v[offset] + coherency[i][j];
            semblance_denom_sum_v[offset] = semblance_denom_sum_v[offset] + semblance_denom[i][j];
            data_in_sum_v[offset] = data_in_sum_v[offset] + ts2d[i][j];
        }
    }
    semblance_denom_sum->WriteArray(H_start, H_end, semblance_denom_sum_v);
    coherency_sum->WriteArray(H_start, H_end, coherency_sum_v);
    data_in_sum->WriteArray(H_start, H_end, data_in_sum_v);

    //std::cout << "finish one file, temp_index " << std::endl;
    return 0;
}

int main(int argc, char *argv[])
{
    int copt, mpi_rank, mpi_size;
    bool has_config_file_flag = false;
    while ((copt = getopt(argc, argv, "i:o:c:h")) != -1)
        switch (copt)
        {
        case 'i':
            xcorr_input_dir.assign(optarg);
            std::cout << xcorr_input_dir << "\n";
            break;
        case 'o':
            stack_output_dir.assign(optarg);
            break;
        case 'c':
            config_file.assign(optarg);
            has_config_file_flag = true;
            break;
        default:
            printf("Wrong option [%c] for %s \n", copt, argv[0]);
            printf_help(argv[0]);
            exit(-1);
            break;
        }

    //Init the MPICH, etc.
    FT_Init(argc, argv);

    if (has_config_file_flag)
        stack_config_reader(config_file, au_rank);

    if (is_ml_weight)
    {
        cout.precision(17);
        std::vector<double> ml_weight_temp;
        read_ml_weight(ml_weight_file, ml_weight_temp, sorted_indexes);
        std::vector<size_t> sorted_indexes_cut;
        ml_weight_sum = 0;

        if (n_weighted_to_stack > 0 && n_weighted_to_stack <= ml_weight_temp.size())
        {
            for (int i = 0; i < n_weighted_to_stack; i++)
            {
                ml_weight.push_back(ml_weight_temp[sorted_indexes[i]]);
                ml_weight_sum = ml_weight_sum + ml_weight_temp[sorted_indexes[i]];
                sorted_indexes_cut.push_back(sorted_indexes[i]);
            }
        }
        else
        {
            sorted_indexes_cut = sorted_indexes;
            for (int i = 0; i < ml_weight_temp.size(); i++)
            {
                ml_weight.push_back(ml_weight_temp[sorted_indexes[i]]);
                ml_weight_sum = ml_weight_sum + ml_weight_temp[sorted_indexes[i]];
            }
        }
        sorted_indexes_str = Vector2String(sorted_indexes_cut);

        std::cout << " sorted_indexes_str =" << sorted_indexes_str << "\n";
        std::cout << " sum of weight =" << ml_weight_sum << "\n";
        PrintVector("  ml_weight(ordered): ", ml_weight);
    }

    // set up the chunk size and the overlap size
    std::vector<int> chunk_size = {chs_per_file, lts_per_file};
    std::vector<int> overlap_size = {0, 0};

    size_t size_after_subset = DasLib::InferTimeSubsetSize(t_start, t_end, sub_start_t, sub_end_t, sample_rate);
    sc_size[0] = chs_per_file;
    sc_size[1] = size_after_subset;

    std::cout << "size_after_subset = " << size_after_subset << "\n";

    semblance_denom_sum = new FT::Array<double>("EP_MEMORY", sc_size);
    coherency_sum = new FT::Array<std::complex<double>>("EP_MEMORY", sc_size);
    data_in_sum = new FT::Array<double>("EP_MEMORY", sc_size);

    if (!au_rank)
        std::cout << "EP_HDF5:" + stack_output_dir + "/" + stack_output_file_final_pwstack + ":" + stack_output_file_dataset_name << "\n";

    final_pwstack = new AU::Array<double>("EP_HDF5:" + stack_output_dir + "/" + stack_output_file_final_pwstack + ":" + stack_output_file_dataset_name, sc_size);
    if (!au_rank)
        std::cout << "init  final_pwstack all"
                  << "\n";

    semblanceWeight = new AU::Array<double>("EP_HDF5:" + stack_output_dir + "/" + stack_output_file_semblanceWeight + ":" + stack_output_file_dataset_name, sc_size);
    if (!au_rank)
        std::cout << "init  semblanceWeight "
                  << "\n";

    phaseWeight = new AU::Array<double>("EP_HDF5:" + stack_output_dir + "/" + stack_output_file_phaseWeight + ":" + stack_output_file_dataset_name, sc_size);
    if (!au_rank)
        std::cout << "init  phaseWeight "
                  << "\n";

    if (!au_rank)
        std::cout << "init   all"
                  << "\n";

    final_pwstack->EndpointControl(OP_DISABLE_MPI_IO, std::vector<std::string>());
    final_pwstack->EndpointControl(OP_DISABLE_COLLECTIVE_IO, std::vector<std::string>());

    semblanceWeight->EndpointControl(OP_DISABLE_COLLECTIVE_IO, std::vector<std::string>());
    semblanceWeight->EndpointControl(OP_DISABLE_MPI_IO, std::vector<std::string>());

    phaseWeight->EndpointControl(OP_DISABLE_COLLECTIVE_IO, std::vector<std::string>());
    phaseWeight->EndpointControl(OP_DISABLE_MPI_IO, std::vector<std::string>());

    if (!au_rank)
        std::cout << "disable collective IO"
                  << "\n";

    //semblanceWeight->Nonvolatile("EP_HDF5:./xcorr_examples_h5_stack_semblanceWeight.h5:/data");
    //phaseWeight->Nonvolatile("EP_HDF5:./xcorr_examples_h5_stack_phaseWeight.h5:/data");
    //Input data,

    FT::Array<double> *A = new FT::Array<double>("EP_DIR:EP_HDF5:" + xcorr_input_dir + ":" + xcorr_input_dataset_name, chunk_size, overlap_size);

    std::vector<int> skip_size = {chs_per_file, lts_per_file};
    A->EnableApplyStride(skip_size);
    if (is_ml_weight)
    {
        std::vector<std::string> index_param;
        index_param.push_back(sorted_indexes_str);

        A->EndpointControl(DIR_FILE_SORT_INDEXES, index_param);
    }
    //std::cout << "Pre clone \n";
    //Clone to create local copy
    //std::complex<double> complex_zero(0, 0);
    //coherency_sum->Fill(complex_zero);
    //std::cout << "Fill \n";

    semblance_denom_sum->Clone();
    //std::cout << "Clone semblance_denom_sum  \n";
    coherency_sum->Clone();
    //std::cout << "Clone coherency_sum  \n";
    data_in_sum->Clone();
    //std::cout << "Pre apply \n";

    if (!au_rank)
        std::cout << "Run apply"
                  << "\n";
    //Run
    A->Apply(stack_udf);

    //std::vector<unsigned long long> H_start_test{0, 0}, H_end_test{static_cast<unsigned long long>(chs_per_file) - 1, static_cast<unsigned long long>(size_after_subset) - 1};
    //std::vector<double> data_in_sum_v_test;
    //data_in_sum->ReadArray(H_start_test, H_end_test, data_in_sum_v_test);
    //PrintVector("data_in_sum_v_test", data_in_sum_v_test);

    semblance_denom_sum->Merge(AU_SUM);
    coherency_sum->Merge(AU_SUM);
    data_in_sum->Merge(AU_SUM);

    double TotalStack;
    AU_Reduce(&nStack, &TotalStack, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (!au_rank)
        std::cout << "Total nStack = " << TotalStack << "\n";
    //semblance_denom_sum->Nonvolatile("EP_HDF5:/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/stacking_files/xcorr_examples_h5_stack_semblance_denom_sum.h5:/semblance_denom_sum");
    //coherency_sum->Nonvolatile("EP_HDF5:/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/stacking_files/xcorr_examples_h5_stack_coherency_sum.h5:/coherency_sum");

    //Clear
    delete A;

    if (!au_rank)
    {
        std::vector<unsigned long long> H_start{0, 0}, H_end{static_cast<unsigned long long>(chs_per_file) - 1, static_cast<unsigned long long>(size_after_subset) - 1};

        std::vector<double> semblance_denom_sum_v;
        semblance_denom_sum->ReadArray(H_start, H_end, semblance_denom_sum_v);
        std::vector<std::complex<double>> coherency_sum_v;
        coherency_sum->ReadArray(H_start, H_end, coherency_sum_v);
        std::vector<double> data_in_sum_v;
        data_in_sum->ReadArray(H_start, H_end, data_in_sum_v);

        std::vector<double> phaseWeight_v(coherency_sum_v.size());
        std::vector<double> semblanceWeight_v(coherency_sum_v.size());

        //PrintVector("data_in_sum_v", data_in_sum_v);

        //PrintVector("coherency_sum_v", coherency_sum_v);

        //PrintVector("semblance_denom_sum_v", semblance_denom_sum_v);

        //std::cout << "Write ... EP_HDF5:" + stack_output_dir + "/" + stack_output_file_semblance_denom_sum + ":" + stack_output_file_dataset_name << "\n"
        //          << std::flush;
        semblance_denom_sum->Nonvolatile("EP_HDF5:" + stack_output_dir + "/" + stack_output_file_semblance_denom_sum + ":" + stack_output_file_dataset_name);

        for (int i = 0; i < chs_per_file * size_after_subset; i++)
        {
            semblanceWeight_v[i] = data_in_sum_v[i] * data_in_sum_v[i] / semblance_denom_sum_v[i];
            phaseWeight_v[i] = std::pow(std::abs(coherency_sum_v[i] / TotalStack), pow_u);
        }

        //PrintVector("semblanceWeight_v", semblanceWeight_v);
        //PrintVector("phaseWeight_v", phaseWeight_v);
        std::vector<double> final_pwstack_v;
        final_pwstack_v.resize(data_in_sum_v.size());
        for (int i = 0; i < data_in_sum_v.size(); i++)
        {
            if (!is_ml_weight)
            {
                data_in_sum_v[i] = data_in_sum_v[i] / TotalStack;
            }
            else
            {
                data_in_sum_v[i] = data_in_sum_v[i] / ml_weight_sum;
            }
            final_pwstack_v[i] = data_in_sum_v[i] * phaseWeight_v[i];
        }

        std::cout << "Store data_in_sum... \n ";

        data_in_sum->WriteArray(H_start, H_end, data_in_sum_v);
        data_in_sum->Nonvolatile("EP_HDF5:" + stack_output_dir + "/" + stack_output_file_data_in_sum_name + ":" + stack_output_file_dataset_name);

        std::cout << "Store semblanceWeight... \n ";

        semblanceWeight->WriteArray(H_start, H_end, semblanceWeight_v);
        phaseWeight->WriteArray(H_start, H_end, phaseWeight_v);
        // std::cout << "Store final_pwstack... \n ";
        final_pwstack->WriteArray(H_start, H_end, final_pwstack_v);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete semblance_denom_sum;
    delete coherency_sum;
    delete data_in_sum;
    delete phaseWeight;

    /*
    PrintScalar("IO_Read_Time: ", IO_Read_Time);
    PrintScalar("CPU_Time: ", CPU_Time);
    PrintScalar("IO_Reduce_Time: ", IO_Reduce_Time);
    PrintScalar("IO_Read/sum/Write_Time_Con: ", IO_Read_Time_Con);
    // DetMean_tim = 0, Subset_tim = 0, CausalityFlagging_tim = 0, instanPhaseEstimator_Time = 0
    PrintScalar("CUP_DetMean_tim :", DetMean_tim);
    PrintScalar("CPU_Subset_tim :", Subset_tim);
    PrintScalar("CPU_CausalityFlagging_tim :", CausalityFlagging_tim);
    PrintScalar("CPU_instanPhaseEstimator_Time :", instanPhaseEstimator_Time);
    PrintScalar("CPU_instanPhaseEstimator_Time2 :", instanPhaseEstimator_Time2);

    PrintScalar("sum_Time :", sum_Time);

    PrintScalar("sum_micro :", sum_micro);
    PrintScalar("sum_micro_sub :", sum_micro_sub);
    */

    FT_Finalize();

    return 0;
}

void read_ml_weight(const std::string ml_weight_file_p, std::vector<double> &ml_weight_p, std::vector<size_t> &sort_indexes_p)
{
    std::string filename_str;
    ifstream inf;
    inf.open(ml_weight_file_p, std::ifstream::in);
    if (!inf)
    {
        AU_EXIT("failed to open weight file : " + ml_weight_file_p);
    }
    string line;
    double val;
    while (getline(inf, line))
    {
        if (line.length() > 1)
        {
            istringstream iss(line);
            iss >> filename_str; //drop the file name
            iss >> val;
            ml_weight_p.push_back(val);
        }
    }
    if (is_ml_weight_ordered)
    {
        sort_indexes_p = DasLib::sort_indexes(ml_weight_p);
    }
    else
    {
        //fill 0, ... ml_weight_p.size()
        sort_indexes_p.resize(ml_weight_p.size());
        std::iota(sort_indexes_p.begin(), sort_indexes_p.end(), 0);
    }
}

void printf_help(char *cmd)
{
    char *msg = (char *)"Usage: %s [OPTION]\n\
      	  -h help (--help)\n\
          -i input directory (xcorr files) \n\
          -o output directory (./ by default) \n\
          Example: mpirun -n 1 %s -i /clusterfs/bear/BinDong_DAS_Data/xcorr_examples_h5/\n";
    fprintf(stdout, msg, cmd, cmd);
}

int stack_config_reader(std::string file_name, int mpi_rank)
{
    INIReader reader(file_name);

    if (reader.ParseError() < 0)
    {
        std::cout << "Can't load [" << file_name << " ]\n";
        return 1;
    }

    std::string xcorr_input_dir_temp = reader.Get("parameter", "xcorr_input_dir", "/clusterfs/bear/BinDong_DAS_Data/xcorr_examples_h5/");

    xcorr_input_dir = realpathEx(xcorr_input_dir_temp);
    //std::string xcorr_input_dir_temp = real_path;
    //free(real_path);
    //xcorr_input_dir = xcorr_input_dir_temp;

    xcorr_input_dataset_name = reader.Get("parameter", "xcorr_input_dataset_name", "/xcoor");

    std::string is_ml_weight_str = reader.Get("parameter", "is_ml_weight", "false");
    if (is_ml_weight_str == "false" || is_ml_weight_str == "0")
    {
        is_ml_weight = false;
    }
    else
    {
        is_ml_weight = true;
    }

    if (is_ml_weight)
    {
        ml_weight_file = reader.Get("parameter", "ml_weight_file", "false");
        n_weighted_to_stack = reader.GetInteger("parameter", "n_weighted_to_stack", -1);
        std::string is_ml_weight_ordered_str = reader.Get("parameter", "is_ml_weight_ordered", "true");
        //std::cout << "is_ml_weight_ordered_str =" << is_ml_weight_ordered_str << "\n";
        if (is_ml_weight_ordered_str == "false" || is_ml_weight_ordered_str == "0")
        {
            is_ml_weight_ordered = false;
        }
        else
        {
            is_ml_weight_ordered = true;
        }
    }

    stack_output_dir = reader.Get("parameter", "stack_output_dir", "./");

    std::string is_flipud_flag_str = reader.Get("parameter", "is_flipud", "true");
    if (is_flipud_flag_str == "false")
    {
        is_flipud_flag = false;
    }
    else
    {
        is_flipud_flag = true;
    }

    chs_per_file = reader.GetInteger("parameter", "chs_per_file", 201);
    lts_per_file = reader.GetInteger("parameter", "lts_per_file", 14999);

    t_start = reader.GetReal("parameter", "t_start", -59.9920000000000);
    t_end = reader.GetReal("parameter", "t_end", 59.9920000000000);
    sub_start_t = reader.GetReal("parameter", "sub_start_t", -59);
    sub_end_t = reader.GetReal("parameter", "sub_end_t", 59);
    sample_rate = reader.GetReal("parameter", "sample_rate", 0.00800000000000000);
    CausalityFlagging_tmin = reader.GetReal("parameter", "CausalityFlagging_tmin", 0.05);
    CausalityFlagging_tmax = reader.GetReal("parameter", "CausalityFlagging_tmax", 3.0);
    CausalityFlagging_fmax = reader.GetReal("parameter", "CausalityFlagging_fmax", 10);
    CausalityFlagging_ButterLow_order = reader.GetReal("parameter", "CausalityFlagging_ButterLow_order", 3);
    CausalityFlagging_ButterLow_fcf = reader.GetReal("parameter", "CausalityFlagging_ButterLow_fcf", 0.16);
    pow_u = reader.GetReal("parameter", "dsiStackFileSet_versatile_pow_u", 0.3);

    stack_output_file_data_in_sum_name = reader.Get("parameter", "stack_output_file_data_in_sum", "xcorr_examples_h5_stack_data_in_sum.h5");
    stack_output_file_final_pwstack = reader.Get("parameter", "stack_output_file_final_pwstack", "xcorr_examples_h5_stack_final_pwstack.h5");
    stack_output_file_phaseWeight = reader.Get("parameter", "stack_output_file_phaseWeight", "xcorr_examples_h5_stack_phaseWeight.h5");
    stack_output_file_semblanceWeight = reader.Get("parameter", "stack_output_file_semblanceWeight", "xcorr_examples_h5_stack_semblanceWeight.h5");
    stack_output_file_semblance_denom_sum = reader.Get("parameter", "stack_output_file_semblance_denom_sum", "xcorr_examples_h5_stack_semblance_denom_sum.h5");
    stack_output_file_dataset_name = reader.Get("parameter", "stack_output_file_dataset_name", "/data");
    if (!mpi_rank)
    {
        std::cout << "\n\n";
        std::cout << termcolor::red << "Parameters to run the Stack: ";

        std::cout << termcolor::blue << "\n\n Input parameters: ";
        std::cout << termcolor::magenta << "\n        xcorr_input_dir = " << termcolor::green << xcorr_input_dir;
        std::cout << termcolor::magenta << "\n        xcorr_input_dataset_name = " << termcolor::green << xcorr_input_dataset_name;
        std::cout << termcolor::magenta << "\n        is_ml_weight = " << termcolor::green << is_ml_weight;
        if (is_ml_weight)
        {
            std::cout << termcolor::magenta << "\n        ml_weight_file = " << termcolor::green << ml_weight_file;
            std::cout << termcolor::magenta << "\n        n_weighted_to_stack = " << termcolor::green << n_weighted_to_stack;
            std::cout << termcolor::magenta << "\n        is_ml_weight_ordered = " << termcolor::green << is_ml_weight_ordered;
        }
        std::cout << termcolor::blue << "\n\n Runtime parameters: ";
        std::cout << termcolor::magenta << "\n\n        lts_per_file = " << termcolor::green << lts_per_file;
        std::cout << termcolor::magenta << "\n\n        chs_per_file = " << termcolor::green << chs_per_file;
        std::cout << termcolor::magenta << "\n        t_start = " << termcolor::green << t_start;
        std::cout << termcolor::magenta << "\n        t_end = " << termcolor::green << t_end;
        std::cout << termcolor::magenta << "\n        sample_rate = " << termcolor::green << sample_rate;
        std::cout << termcolor::magenta << "\n        sub_start_t = " << termcolor::green << sub_start_t;
        std::cout << termcolor::magenta << "\n        sub_end_t = " << termcolor::green << sub_end_t;
        std::cout << termcolor::magenta << "\n        CausalityFlagging_tmin = " << termcolor::green << CausalityFlagging_tmin;
        std::cout << termcolor::magenta << "\n        CausalityFlagging_tmax = " << termcolor::green << CausalityFlagging_tmax;
        std::cout << termcolor::magenta << "\n        CausalityFlagging_fmax = " << termcolor::green << CausalityFlagging_fmax;
        std::cout << termcolor::magenta << "\n        CausalityFlagging_ButterLow_order = " << termcolor::green << CausalityFlagging_ButterLow_order;
        std::cout << termcolor::magenta << "\n        CausalityFlagging_ButterLow_fcf = " << termcolor::green << CausalityFlagging_ButterLow_fcf;
        std::cout << termcolor::magenta << "\n        dsiStackFileSet_versatile_pow_u = " << termcolor::green << pow_u;
        std::cout << termcolor::magenta << "\n        is_flipud_flag = " << termcolor::green << is_flipud_flag;

        std::cout << termcolor::blue << "\n\n Output parameters: ";
        std::cout << termcolor::magenta << "\n\n        stack_output_dir = " << termcolor::green << stack_output_dir;
        std::cout << termcolor::magenta << "\n        output_file_data_in_sum_name = " << termcolor::green << stack_output_file_data_in_sum_name;
        std::cout << termcolor::magenta << "\n        output_file_final_pwstack = " << termcolor::green << stack_output_file_final_pwstack;
        std::cout << termcolor::magenta << "\n        output_file_phaseWeight = " << termcolor::green << stack_output_file_phaseWeight;
        std::cout << termcolor::magenta << "\n        output_file_semblanceWeight = " << termcolor::green << stack_output_file_semblanceWeight;
        std::cout << termcolor::magenta << "\n        output_file_semblance_denom_sum = " << termcolor::green << stack_output_file_semblance_denom_sum;
        std::cout << termcolor::magenta << "\n        output_file_dataset_name = " << termcolor::green << stack_output_file_dataset_name;

        std::cout << termcolor::reset << "\n\n";
    }
    fflush(stdout);

    return 0;
}