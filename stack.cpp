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
#include "au.h"
#include "DasLib.h"
#include <stdlib.h>

using namespace std;
using namespace AU;

/**
 * @brief some help function at the end
 * 
 * @param cmd 
 */
void printf_help(char *cmd);
int stack_config_reader(std::string file_name, int mpi_rank);

extern int au_size;
extern int au_rank;

//#define LTS 14999 //length of time series
//#define CHS 201   //channels

int lts_per_file = 14999;
int chs_per_file = 201;

std::string config_file = "./stack.config";

std::string xcorr_input_dir = "/clusterfs/bear/BinDong_DAS_Data/xcorr_examples_h5/";
std::string xcorr_input_dataset_name = "/xcoor";
std::string stack_output_dir = "./";

double t_start = -59.9920000000000, t_end = 59.9920000000000, sample_rate = 0.00800000000000000;
double sub_start_t = -59, sub_end_t = 59;
double pow_u = 0.3;
double CausalityFlagging_tmin = 0.05;
double CausalityFlagging_tmax = 3.0;
double CausalityFlagging_fmax = 10;
double CausalityFlagging_ButterLow_order = 3;
double CausalityFlagging_ButterLow_fcf = 0.16; //filter cutoff frequency

std::vector<unsigned long long> sc_size(2);
AU::Array<double> *semblance_denom_sum;
AU::Array<std::complex<double>> *coherency_sum;
AU::Array<double> *data_in_sum;
AU::Array<double> *phaseWeight;
AU::Array<double> *semblanceWeight;

AU::Array<double> *final_pwstack;

double nStack = 0;

inline Stencil<double>
stack_udf(const Stencil<double> &iStencil)
{
    nStack++;
    //std::cout << "nStack: " << nStack++ << " at " << au_mpi_rank_global << "\n";

    std::vector<int> start_offset{0, 0}, end_offset{chs_per_file - 1, lts_per_file - 1};
    //std::vector<double> ts = iStencil.Read(start_offset, end_offset);
    //std::vector<double> ts(CHS * LTS);
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
    //PrintVV("ts2d at rank " + std::to_string(au_rank), ts2d);

    DasLib::DeleteMedian(ts2d);

    //Remove the media
    for (int i = 0; i < chs_per_file; i++)
    {
        //Subset
        ts2d[i] = DasLib::TimeSubset(ts2d[i], t_start, t_end, sub_start_t, sub_end_t, sample_rate);
    }

    bool flag = DasLib::CausalityFlagging(ts2d, CausalityFlagging_tmin, CausalityFlagging_tmax, CausalityFlagging_fmax, sub_start_t, sub_end_t, sample_rate, CausalityFlagging_ButterLow_order, CausalityFlagging_ButterLow_fcf);
    if (flag == false)
    {
        std::cout << "flipud  at channel " << nStack << "\n";
        for (int i = 0; i < chs_per_file; i++)
        {
            std::reverse(ts2d[i].begin(), ts2d[i].end());
        }
        //PrintVector("Before reverse ts2d at " + std::to_string(nStack), ts2d[0]);
        //PrintVector("After reverse ts2d at " + std::to_string(nStack), ts2d[0]);
    }

    size_t LTS_new = ts2d[0].size();
    std::vector<std::vector<double>> semblance_denom;
    std::vector<std::vector<std::complex<double>>> coherency;
    semblance_denom.resize(chs_per_file);
    coherency.resize(chs_per_file);
    for (int i = 0; i < chs_per_file; i++)
    {
        semblance_denom[i].resize(LTS_new);
        //coherency[i].resize(LTS_new);
        for (int j = 0; j < LTS_new; j++)
        {
            semblance_denom[i][j] = ts2d[i][j] * ts2d[i][j];
        }

        coherency[i] = DasLib::instanPhaseEstimator(ts2d[i]);
    }

    std::vector<unsigned long long> H_start{0, 0}, H_end{static_cast<unsigned long long>(chs_per_file) - 1, static_cast<unsigned long long>(LTS_new) - 1};
    std::vector<double> semblance_denom_sum_v = semblance_denom_sum->ReadArray(H_start, H_end);
    std::vector<std::complex<double>> coherency_sum_v = coherency_sum->ReadArray(H_start, H_end);
    std::vector<double> data_in_sum_v = data_in_sum->ReadArray(H_start, H_end);

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

    //PrintVector("coherency_sum_v (after)" + std::to_string(au_rank), coherency_sum_v);

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
    AU_Init(argc, argv);

    if (has_config_file_flag)
        stack_config_reader(config_file, au_rank);

    // set up the chunk size and the overlap size
    std::vector<int> chunk_size = {chs_per_file, lts_per_file};
    std::vector<int> overlap_size = {0, 0};

    size_t size_after_subset = DasLib::InferTimeSubsetSize(t_start, t_end, sub_start_t, sub_end_t, sample_rate);
    sc_size[0] = chs_per_file;
    sc_size[1] = size_after_subset;

    std::cout << "size_after_subset = " << size_after_subset << "\n";

    semblance_denom_sum = new AU::Array<double>("EP_MEMORY", sc_size);
    coherency_sum = new AU::Array<std::complex<double>>("EP_MEMORY", sc_size);
    data_in_sum = new AU::Array<double>("EP_MEMORY", sc_size);
    final_pwstack = new AU::Array<double>("EP_HDF5:./xcorr_examples_h5_stack_final_pwstack.h5:/data", sc_size);
    semblanceWeight = new AU::Array<double>("EP_HDF5:./xcorr_examples_h5_stack_semblanceWeight.h5:/data", sc_size);
    phaseWeight = new AU::Array<double>("EP_HDF5:./xcorr_examples_h5_stack_phaseWeight.h5:/data", sc_size);

    //semblanceWeight->Nonvolatile("EP_HDF5:./xcorr_examples_h5_stack_semblanceWeight.h5:/data");
    //phaseWeight->Nonvolatile("EP_HDF5:./xcorr_examples_h5_stack_phaseWeight.h5:/data");
    //Input data,

    AU::Array<double> *A = new AU::Array<double>("EP_DIR:EP_HDF5:" + xcorr_input_dir + ":" + xcorr_input_dataset_name, chunk_size, overlap_size);

    std::vector<int> skip_size = {chs_per_file, lts_per_file};
    A->EnableApplyStride(skip_size);

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

    //Run
    A->Apply(stack_udf);

    std::vector<unsigned long long> H_start_test{0, 0}, H_end_test{static_cast<unsigned long long>(chs_per_file) - 1, static_cast<unsigned long long>(size_after_subset) - 1};
    std::vector<double> data_in_sum_v_test = data_in_sum->ReadArray(H_start_test, H_end_test);
    PrintVector("data_in_sum_v_test", data_in_sum_v_test);

    semblance_denom_sum->Merge(AU_SUM);
    coherency_sum->Merge(AU_SUM);
    data_in_sum->Merge(AU_SUM);

    double TotalStack;
    AU_Reduce(&nStack, &TotalStack, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    std::cout << "Total nStack = " << TotalStack << "\n";
    //semblance_denom_sum->Nonvolatile("EP_HDF5:/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/stacking_files/xcorr_examples_h5_stack_semblance_denom_sum.h5:/semblance_denom_sum");
    //coherency_sum->Nonvolatile("EP_HDF5:/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/stacking_files/xcorr_examples_h5_stack_coherency_sum.h5:/coherency_sum");

    //Clear
    delete A;

    if (!au_rank)
    {
        std::vector<unsigned long long> H_start{0, 0}, H_end{static_cast<unsigned long long>(chs_per_file) - 1, static_cast<unsigned long long>(size_after_subset) - 1};

        std::vector<double> semblance_denom_sum_v = semblance_denom_sum->ReadArray(H_start, H_end);
        std::vector<std::complex<double>> coherency_sum_v = coherency_sum->ReadArray(H_start, H_end);
        std::vector<double> data_in_sum_v = data_in_sum->ReadArray(H_start, H_end);

        std::vector<double> phaseWeight_v(coherency_sum_v.size());
        std::vector<double> semblanceWeight_v(coherency_sum_v.size());

        //PrintVector("data_in_sum_v", data_in_sum_v);

        //PrintVector("coherency_sum_v", coherency_sum_v);

        //PrintVector("semblance_denom_sum_v", semblance_denom_sum_v);

        semblance_denom_sum->Nonvolatile("EP_HDF5:./xcorr_examples_h5_stack_semblance_denom_sum.h5:/data");

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
            data_in_sum_v[i] = data_in_sum_v[i] / TotalStack;
            final_pwstack_v[i] = data_in_sum_v[i] * phaseWeight_v[i];
        }
        data_in_sum->WriteArray(H_start, H_end, data_in_sum_v);

        data_in_sum->Nonvolatile("EP_HDF5:./xcorr_examples_h5_stack_data_in_sum.h5:/data");

        semblanceWeight->WriteArray(H_start, H_end, semblanceWeight_v);
        phaseWeight->WriteArray(H_start, H_end, phaseWeight_v);

        final_pwstack->WriteArray(H_start, H_end, final_pwstack_v);
    }
    delete semblance_denom_sum;
    delete coherency_sum;
    delete data_in_sum;
    delete phaseWeight;

    AU_Finalize();

    return 0;
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

    stack_output_dir = reader.Get("parameter", "stack_output_dir", "./");

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
    pow_u = reader.GetReal("parameter", "pow_u", 0.3);

    if (!mpi_rank)
    {
        std::cout << "\n\n";
        std::cout << "Configurations to run the Stack: ";
        std::cout << "\n        xcorr_input_dir = " << xcorr_input_dir;
        std::cout << "\n        xcorr_input_dataset_name = " << xcorr_input_dataset_name;
        std::cout << "\n        stack_output_dir = " << stack_output_dir;
        std::cout << "\n        chs_per_file = " << chs_per_file;
        std::cout << "\n        t_start = " << t_start;
        std::cout << "\n        t_end = " << t_end;
        std::cout << "\n        sample_rate = " << sample_rate;
        std::cout << "\n        sub_start_t = " << sub_start_t;
        std::cout << "\n        sub_end_t = " << sub_end_t;
        std::cout << "\n        CausalityFlagging_tmin = " << CausalityFlagging_tmin;
        std::cout << "\n        CausalityFlagging_tmax = " << CausalityFlagging_tmax;
        std::cout << "\n        CausalityFlagging_fmax = " << CausalityFlagging_fmax;
        std::cout << "\n        CausalityFlagging_ButterLow_order = " << CausalityFlagging_ButterLow_order;
        std::cout << "\n        CausalityFlagging_ButterLow_fcf = " << CausalityFlagging_ButterLow_fcf;
        std::cout << "\n        pow_u = " << pow_u;
        std::cout << "\n\n";
    }
    fflush(stdout);

    return 0;
}