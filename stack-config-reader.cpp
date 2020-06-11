// Example that shows simple usage of the INIReader class

#include <iostream>
#include <sstream>
#include <sstream>
#include <vector>
#include <string>
#include "INIReader.h"
#include <iomanip> // std::setprecision

extern std::string xcorr_input_dir;
extern std::string stack_output_dir;
extern int chs_per_file; //channels (per file)
extern int lts_per_file; //length of time series (per file)
extern double t_start, t_end, sample_rate;
extern double sub_start_t, sub_end_t;
extern double pow_u;
extern double CausalityFlagging_tmin;
extern double CausalityFlagging_tmax;
extern double CausalityFlagging_fmax;
extern double CausalityFlagging_ButterLow_order;
extern double CausalityFlagging_ButterLow_fcf;

int stack_config_reader(std::string file_name, int mpi_rank)
{
    INIReader reader(file_name);

    if (reader.ParseError() < 0)
    {
        std::cout << "Can't load [" << file_name << " ]\n";
        return 1;
    }

    std::string temp_str;
    temp_str = reader.Get("parameter", "xcorr_input_dir", "UNKNOWN");
    if (temp_str != "UNKNOWN")
    {
        xcorr_input_dir = temp_str;
    }

    temp_str = reader.Get("parameter", "stack_output_dir", "UNKNOWN");
    if (temp_str != "UNKNOWN")
    {
        stack_output_dir = temp_str;
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
    pow_u = reader.GetReal("parameter", "pow_u", 0.3);

    if (!mpi_rank)
    {
        std::cout << "Configurations to run the Stack: ";
        std::cout << "\n        xcorr_input_dir = " << xcorr_input_dir;
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
    }
    fflush(stdout);

    return 0;
}