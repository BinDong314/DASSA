#ifndef __XCORRELATION_H__
#define __XCORRELATION_H__

#include <vector>
#include <exception>
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>
#include <stdarg.h>
#include <vector>
#include <stdlib.h>
#include <math.h> /* ceil  and floor*/
#include <cstring>
#include <numeric>
#include <iomanip> // std::setprecision
#define NAME_LENGTH 1024

//#include "DasLib/das-fft-help.cc"
//#include "DasLib/liir.c"
//#include "DasLib/bwlp.cc"
//#include "DasLib/config-reader.cc"

#include "DasLib.h"

// n = order of the filter
// fc = filter cutoff frequency as a fraction of Pi [0,1]
int butter_low(int n, double fcf, std::vector<double> &A, std::vector<double> &B)
{
    int sff = 1;  // scale flag: 1 to scale, 0 to not scale ccof
    int i;        // loop variables
    double sf;    // scaling factor
    double *dcof; // d coefficients
    int *ccof;    // c coefficients

    /* calculate the d coefficients */
    dcof = dcof_bwlp(n, fcf);
    if (dcof == NULL)
    {
        perror("Unable to calculate d coefficients");
        return (-1);
    }

    /* calculate the c coefficients */
    ccof = ccof_bwlp(n);
    if (ccof == NULL)
    {
        perror("Unable to calculate c coefficients");
        return (-1);
    }

    sf = sf_bwlp(n, fcf); /* scaling factor for the c coefficients */

    A.resize(n + 1);
    B.resize(n + 1);
    /* Output the c (A) coefficients */
    if (sff == 0)
    {
        for (i = 0; i <= n; ++i)
        {
            A[i] = ccof[i];
        }
    }
    else
    {
        for (i = 0; i <= n; ++i)
        {
            A[i] = (double)ccof[i] * sf;
        }
    }

    /* Output the d(B) coefficients */
    for (i = 0; i <= n; ++i)
    {
        B[i] = dcof[i];
    }

    free(dcof);
    free(ccof);

    return 0;
}

//#define DEBUG_VERBOSE 0
using namespace std;
using namespace DasLib;

//N0: the size of orignal vector
//nPoint: size after resample = N0 * 1/R (R = 4)
//nfft: size of vector for FFT, power of 2 and > 2*nPoint - 1
//nXCORR: size of result gatherXcorr vector,  = 2*nPoint - 1
int n0;
int nPoint;
int nfft;
int nXCORR;

//Paramters fom butter(3, 2*0.25, 'low'), for filtfilt
int butter_order = 3;
double cut_frequency_low = 0.25;
vector<double> BUTTER_A;
vector<double> BUTTER_B;

//vector<double> BUTTER_A{0.031689343849711, 0.095068031549133, 0.095068031549133, 0.031689343849711};
//vector<double> BUTTER_B{1.000000000000000, -1.459029062228061, 0.910369000290069, -0.197825187264319};

//For resample, R = dt_new/dt
double DT = 0.002;
double DT_NEW = 0.008;

//For moveing mean
double WINLEN_SEC = 0.5;
int nPoint_hal_win;

//For interp1
double fNyquist; //250
std::vector<double> INTERP_Z{0, 0.5, 1, 1, 0.5, 0};
std::vector<double> INTERP_ZF{0, 0.002, 0.006, 14.5, 15, fNyquist};
double df;
double eCoeff = 1.0;

//Maste channel
unsigned long long MASTER_INDEX = 0;

//Parameters for ArrayUDF
int auto_chunk_dims_index = 0;
std::vector<int> strip_size(2);
std::vector<int> chunk_size(2);
int chunked_batch_factor = 32; //To enable a large chunk

//Parameters to enable "window" operation
int user_window_size = 1;
int set_window_size_flag = 0;
int window_batch = 1;

//Flag to have FFT in row-direction (first dimension)
//By default, it work in column-direction (second dimension)
int row_major_flag = 1;

//View's parameter
bool enable_view_flag = 0;
std::vector<unsigned long long> view_start{0, 0};
std::vector<unsigned long long> view_count{30000, 11648};
std::vector<int> view_os_size{0, 0};

int decimation_flag = 0;
int decimation_size;

//Using n0 as intial intput
//MR: mpi rank
//MS: mpi size
//DATADIMS: size of 2D data
#define INIT_PARS(MR, MS, DATADIMS)                                                     \
    {                                                                                   \
        if (row_major_flag == 0)                                                        \
        {                                                                               \
            auto_chunk_dims_index = 0;                                                  \
            strip_size[0] = DATADIMS[0];                                                \
            strip_size[1] = 1;                                                          \
            if (set_window_size_flag == 0)                                              \
            {                                                                           \
                n0 = DATADIMS[0];                                                       \
            }                                                                           \
            else                                                                        \
            {                                                                           \
                assert(user_window_size > 0 && user_window_size < DATADIMS[0]);         \
                n0 = user_window_size;                                                  \
                if (DATADIMS[0] % user_window_size == 0)                                \
                {                                                                       \
                    window_batch = DATADIMS[0] / user_window_size;                      \
                }                                                                       \
                else                                                                    \
                {                                                                       \
                    window_batch = DATADIMS[0] / user_window_size + 1;                  \
                }                                                                       \
            }                                                                           \
        }                                                                               \
        else                                                                            \
        {                                                                               \
            auto_chunk_dims_index = 1;                                                  \
            strip_size[0] = 1;                                                          \
            strip_size[1] = DATADIMS[1];                                                \
            if (set_window_size_flag == 0)                                              \
            {                                                                           \
                n0 = DATADIMS[1];                                                       \
            }                                                                           \
            else                                                                        \
            {                                                                           \
                assert(user_window_size > 0 && user_window_size < DATADIMS[1]);         \
                n0 = user_window_size;                                                  \
                if (DATADIMS[1] % user_window_size == 0)                                \
                {                                                                       \
                    window_batch = DATADIMS[1] / user_window_size;                      \
                }                                                                       \
                else                                                                    \
                {                                                                       \
                    window_batch = DATADIMS[1] / user_window_size + 1;                  \
                }                                                                       \
            }                                                                           \
        }                                                                               \
        nPoint = ceil(n0 / (DT_NEW / DT));                                              \
        FIND_M_POWER2(nPoint, nfft);                                                    \
        nXCORR = 2 * nPoint - 1;                                                        \
        fNyquist = 0.5 / DT_NEW;                                                        \
        nPoint_hal_win = floor((2 * floor(WINLEN_SEC / DT_NEW / 2) + 1) / 2);           \
        INTERP_ZF[5] = fNyquist;                                                        \
        df = 2.0 * fNyquist / (double)nfft;                                             \
        cut_frequency_low = (0.5 / DT_NEW) / (0.5 / DT);                                \
        butter_low(butter_order, cut_frequency_low, BUTTER_A, BUTTER_B);                \
        if (MR == 0)                                                                    \
        {                                                                               \
            printf("\n Some parameters for DAS:\n");                                    \
            printf("    n0(intput size) = %d \n", n0);                                  \
            printf("             nPoint = %d \n", nPoint);                              \
            printf("               nfft = %d \n", nfft);                                \
            printf("           fNyquist = %f \n", fNyquist);                            \
            printf("     nPoint_hal_win = %d \n", nPoint_hal_win);                      \
            printf("                 df = %f \n", df);                                  \
            printf("nXCORR(output size) = %d \n", nXCORR);                              \
            printf("    butter low freq = %f \n", cut_frequency_low);                   \
            printf("                 \n");                                              \
            printf("ArrayUDF strip  size  = (%d, %d)\n", strip_size[0], strip_size[1]); \
            printf("ArrayUDF window size  = %d\n", user_window_size);                   \
            printf("ArrayUDF window batch = %d\n", window_batch);                       \
            printf("                 \n");                                              \
        }                                                                               \
        fflush(stdout);                                                                 \
    }

//Pre-allocated data space, to save time
std::vector<double> X;
std::vector<double> TC; //temp cache
std::vector<double> shapingFilt;
fftw_complex *fft_in;
fftw_complex *fft_out;
fftw_complex *master_fft;

//Final results
std::vector<float> gatherXcorr;
std::vector<float> gatherXcorr_per_batch;

//WS: window size
#define INIT_SPACE()                                                             \
    {                                                                            \
        X.resize(n0);                                                            \
        TC.resize(nfft);                                                         \
        shapingFilt.resize(nfft);                                                \
        fft_in = fftw_alloc_complex(nfft);                                       \
        fft_out = fftw_alloc_complex(nfft);                                      \
        master_fft = fftw_alloc_complex(nfft * window_batch);                    \
        if (fft_in == NULL || fft_out == NULL || master_fft == NULL)             \
        {                                                                        \
            printf("not enough memory for fft, in %s:%d\n", __FILE__, __LINE__); \
            exit(-1);                                                            \
        }                                                                        \
        gatherXcorr_per_batch.resize(nXCORR);                                    \
        gatherXcorr.resize(nXCORR *window_batch);                                \
        std::vector<double> FF_LHS, LHS;                                         \
        FF_LHS.resize(nfft / 2);                                                 \
        LHS.resize(nfft / 2);                                                    \
        for (int i = 0; i < nfft / 2; i++)                                       \
        {                                                                        \
            FF_LHS[i] = df * (i + 1);                                            \
        }                                                                        \
        interp1(INTERP_ZF, INTERP_Z, FF_LHS, LHS);                               \
        int nfft_half = nfft / 2;                                                \
        for (int i = 0; i < nfft_half; i++)                                      \
        {                                                                        \
            shapingFilt[i] = LHS[i];                                             \
            shapingFilt[i + nfft_half] = LHS[nfft_half - i - 1];                 \
        }                                                                        \
        FF_LHS.clear();                                                          \
        LHS.clear();                                                             \
        std::vector<double>().swap(FF_LHS);                                      \
        std::vector<double>().swap(LHS);                                         \
    }

//WS: window size
#define INIT_SPACE_OMP()                                         \
    {                                                            \
        shapingFilt.resize(nfft);                                \
        std::vector<double> FF_LHS, LHS;                         \
        FF_LHS.resize(nfft / 2);                                 \
        LHS.resize(nfft / 2);                                    \
        for (int i = 0; i < nfft / 2; i++)                       \
        {                                                        \
            FF_LHS[i] = df * (i + 1);                            \
        }                                                        \
        interp1(INTERP_ZF, INTERP_Z, FF_LHS, LHS);               \
        int nfft_half = nfft / 2;                                \
        for (int i = 0; i < nfft_half; i++)                      \
        {                                                        \
            shapingFilt[i] = LHS[i];                             \
            shapingFilt[i + nfft_half] = LHS[nfft_half - i - 1]; \
        }                                                        \
        FF_LHS.clear();                                          \
        LHS.clear();                                             \
        std::vector<double>().swap(FF_LHS);                      \
        std::vector<double>().swap(LHS);                         \
    }

//R = dt_new/dt
#define CLEAR_SPACE()                  \
    {                                  \
        X.clear();                     \
        TC.clear();                    \
        shapingFilt.clear();           \
        fftw_free(fft_in);             \
        fftw_free(fft_out);            \
        fftw_free(master_fft);         \
        gatherXcorr.clear();           \
        gatherXcorr_per_batch.clear(); \
    }

//R = dt_new/dt
#define CLEAR_SPACE_OMP()                        \
    {                                            \
        shapingFilt.clear();                     \
        std::vector<double>().swap(shapingFilt); \
    }

#define INIT_FFTW(FFT_IN_V, XXX_P, XN_P, NFFT_P, FFT_OUT_V) \
    {                                                       \
        for (int iii = 0; iii < NFFT_P; iii++)              \
        {                                                   \
            if (iii < XN_P)                                 \
            {                                               \
                FFT_IN_V[iii][0] = XXX_P[iii];              \
                FFT_IN_V[iii][1] = 0;                       \
            }                                               \
            else                                            \
            {                                               \
                FFT_IN_V[iii][0] = 0;                       \
                FFT_IN_V[iii][1] = 0;                       \
            }                                               \
            FFT_OUT_V[iii][0] = 0;                          \
            FFT_OUT_V[iii][1] = 0;                          \
        }                                                   \
    }

//Add function to init FFT_IN
#define INIT_FFTW_FILL(FFT_IN_V, XXX_P, XN_P, NFFT_P) \
    {                                                 \
        for (int iii = 0; iii < NFFT_P; iii++)        \
        {                                             \
            if (iii < XN_P)                           \
            {                                         \
                FFT_IN_V[iii][0] = XXX_P[iii];        \
                FFT_IN_V[iii][1] = 0;                 \
            }                                         \
            else                                      \
            {                                         \
                FFT_IN_V[iii][0] = 0;                 \
                FFT_IN_V[iii][1] = 0;                 \
            }                                         \
        }                                             \
    }

#define INIT_FFTW_V_ZERO(FFT_IN_V, NFFT_P)     \
    {                                          \
        for (int iii = 0; iii < NFFT_P; iii++) \
        {                                      \
            FFT_IN_V[iii][0] = 0;              \
            FFT_IN_V[iii][1] = 0;              \
        }                                      \
    }

//direction: FFTW_FORWARD,  FFTW_BACKWARD
#define FFT_HELP_W(NNN, fft_in_p, fft_out_p, direction_p)                               \
    {                                                                                   \
        fftw_plan fft_p;                                                                \
        fft_p = fftw_plan_dft_1d(NNN, fft_in_p, fft_out_p, direction_p, FFTW_ESTIMATE); \
        fftw_execute(fft_p);                                                            \
        fftw_destroy_plan(fft_p);                                                       \
    }

#define MALLOC_FFT(PFFT, NFFT)                                               \
    PFFT = fftw_alloc_complex(NFFT);                                         \
    if (PFFT == NULL)                                                        \
    {                                                                        \
        printf("not enough memory for fft, in %s:%d\n", __FILE__, __LINE__); \
        exit(-1);                                                            \
    }

#define FREE_FFT(PFFT) \
    fftw_free(PFFT);   \
    PFFT = NULL

//#include "fft/kiss_fft.h"
//kiss_fft_cpx *fft_in_temp;
//kiss_fft_cpx *fft_out_temp;
//kiss_fft_cpx *master_vector_fft;
//unsigned int fft_in_legnth;
//direction: 0,  1
#define FFT_HELP_K(N_P, fft_in_P, fft_out_P, direction_P) \
    {                                                     \
        kiss_fft_cfg cfg;                                 \
        cfg = kiss_fft_alloc(N_P, direction_P, 0, 0);     \
        kiss_fft(cfg, fft_in_P, fft_out_P);               \
        kiss_fft_free(cfg);                               \
    }

#define INIT_FFTW_K(FFT_IN_V, XXX_P, XN_P, NFFT_P, FFT_OUT_V) \
    {                                                         \
        for (int iii = 0; iii < NFFT_P; iii++)                \
        {                                                     \
            if (iii < XN_P)                                   \
            {                                                 \
                FFT_IN_V[iii].r = XXX_P[iii];                 \
                FFT_IN_V[iii].i = 0;                          \
            }                                                 \
            else                                              \
            {                                                 \
                FFT_IN_V[iii].r = 0;                          \
                FFT_IN_V[iii].i = 0;                          \
            }                                                 \
            FFT_OUT_V[iii].r = 0;                             \
            FFT_OUT_V[iii].i = 0;                             \
        }                                                     \
    }

#define MALLOC_FFT_K(PFFT, NFFT)                                             \
    PFFT = (kiss_fft_cpx *)malloc(sizeof(kiss_fft_cpx) * NFFT);              \
    if (PFFT == NULL)                                                        \
    {                                                                        \
        printf("not enough memory for fft, in %s:%d\n", __FILE__, __LINE__); \
        exit(-1);                                                            \
    }

#define FREE_FFT_K(PFFT) \
    free(PFFT);          \
    PFFT = NULL

#define FFTW_AVAI 1

#ifdef FFTW_AVAI
#include <fftw3.h>
typedef fftw_complex *FFT_DATA_TYPEP;
#define PREPARE_FFT INIT_FFTW
#define RUN_FFT FFT_HELP_W
#define ALLOCATE_FFT MALLOC_FFT
#define CLEAR_FFT FREE_FFT
#define FORWARD_FLAG FFTW_FORWARD   //-1
#define BACKWARD_FLAG FFTW_BACKWARD //1
#else
#include "kiss_fft.h"
typedef kiss_fft_cpx *FFT_DATA_TYPEP;
#define PREPARE_FFT INIT_FFTW_K
#define RUN_FFT FFT_HELP_K
#define ALLOCATE_FFT MALLOC_FFT_K
#define CLEAR_FFT FREE_FFT_K
#define FORWARD_FLAG 0
#define BACKWARD_FLAG 1
#endif

// Use "xi >= XN - HAL_WIN - 1" to match matlab code
#define MOVING_MEAN(XV, YV, HAL_WIN)                                \
    {                                                               \
        unsigned long long XN = XV.size(), index1, index2;          \
        YV.resize(XN);                                              \
        for (unsigned long long xi = 0; xi < XN; xi++)              \
        {                                                           \
            if (xi < HAL_WIN)                                       \
            {                                                       \
                index1 = 0;                                         \
                index2 = xi;                                        \
            }                                                       \
            else if (xi >= XN - HAL_WIN - 1)                        \
            {                                                       \
                index1 = xi;                                        \
                index2 = XN - 1;                                    \
            }                                                       \
            else                                                    \
            {                                                       \
                index1 = xi - HAL_WIN;                              \
                index2 = xi + HAL_WIN;                              \
            }                                                       \
                                                                    \
            double norm_denom = 0;                                  \
            for (int j = index1; j <= index2; j++)                  \
            {                                                       \
                norm_denom = norm_denom + std::abs(XV[j]);          \
            }                                                       \
            YV[xi] = XV[xi] / (norm_denom / (index2 - index1 + 1)); \
        }                                                           \
    }
//It is power of 2 and greater than minimum_m
#define FIND_M_POWER2(MIN_M, M)                \
    {                                          \
        unsigned long long mt = 2 * MIN_M - 1; \
        while ((mt & (mt - 1)) != 0)           \
        {                                      \
            mt = mt + 1;                       \
        }                                      \
        M = mt;                                \
    }

#define FILL_LHS(NFFT, F_LHS, DF)      \
    {                                  \
        for (int i = 0; i < NFFT; i++) \
        {                              \
            F_LHS[i] = DF * (i + 1);   \
        }                              \
    }

/*
 * XX is the input data 
 * YY is used as cache
 * GX is the result
 */
#define FFT_PROCESSING(XX, YY, GX, MFFT, BI, BS)                                                                   \
    {                                                                                                              \
        detrend(&(XX[0]), XX.size());                                                                              \
        au_time_elap(" ++ detrend ");                                                                              \
        filtfilt(BUTTER_A, BUTTER_B, XX, YY);                                                                      \
        au_time_elap(" ++ filtfilt ");                                                                             \
        resample(1, DT_NEW / DT, YY, XX);                                                                          \
        au_time_elap(" ++ resample ");                                                                             \
        MOVING_MEAN(XX, YY, nPoint_hal_win);                                                                       \
        au_time_elap(" ++ MOVING_MEAN ");                                                                          \
        INIT_FFTW(fft_in, YY, nPoint, nfft, fft_out);                                                              \
        au_time_elap(" ++ Init fft ");                                                                             \
        FFT_HELP_W(nfft, fft_in, fft_out, FFTW_FORWARD);                                                           \
        au_time_elap(" ++ First fft ");                                                                            \
        for (int ii = 0; ii < nfft; ii++)                                                                          \
        {                                                                                                          \
            double temp_f;                                                                                         \
            temp_f = pow(sqrt(fft_out[ii][0] * fft_out[ii][0] + fft_out[ii][1] * fft_out[ii][1]), eCoeff) + 0.001; \
            fft_in[ii][0] = (fft_out[ii][0] + 0.001) / temp_f * shapingFilt[ii];                                   \
            fft_in[ii][1] = (fft_out[ii][1]) / temp_f * shapingFilt[ii];                                           \
            fft_out[ii][0] = 0;                                                                                    \
            fft_out[ii][1] = 0;                                                                                    \
        }                                                                                                          \
        au_time_elap(" ++ Corr fft ");                                                                             \
        FFT_HELP_W(nfft, fft_in, fft_out, FFTW_BACKWARD);                                                          \
        au_time_elap(" ++ Rev fft ");                                                                              \
        YY.resize(nPoint);                                                                                         \
        for (int i = 0; i < nPoint; i++)                                                                           \
        {                                                                                                          \
            YY[i] = fft_out[i][0] / ((double)nfft);                                                                \
        }                                                                                                          \
        INIT_FFTW(fft_in, YY, nPoint, nfft, fft_out);                                                              \
        FFT_HELP_W(nfft, fft_in, fft_out, FFTW_FORWARD);                                                           \
        for (int j = 0; j < nfft; j++)                                                                             \
        {                                                                                                          \
            fft_in[j][0] = MFFT[j + BI * BS][0] * fft_out[j][0] + MFFT[j + BI * BS][1] * fft_out[j][1];            \
            fft_in[j][1] = MFFT[j + BI * BS][1] * fft_out[j][0] - MFFT[j + BI * BS][0] * fft_out[j][1];            \
            fft_out[j][0] = 0;                                                                                     \
            fft_out[j][1] = 0;                                                                                     \
        }                                                                                                          \
        FFT_HELP_W(nfft, fft_in, fft_out, FFTW_BACKWARD);                                                          \
        int gatherXcorr_index = 0;                                                                                 \
        for (int k = nfft - nPoint + 1; k < nfft; k++)                                                             \
        {                                                                                                          \
            GX[gatherXcorr_index] = (float)(fft_out[k][0] / (double)nfft);                                         \
            gatherXcorr_index++;                                                                                   \
        }                                                                                                          \
        for (int l = 0; l < nPoint; l++)                                                                           \
        {                                                                                                          \
            GX[gatherXcorr_index] = (float)(fft_out[l][0] / (double)nfft);                                         \
            gatherXcorr_index++;                                                                                   \
        }                                                                                                          \
    }

#define FFT_PREPROCESSING(XPP, YPP)                                                                   \
    {                                                                                                 \
        detrend(&(XPP[0]), XPP.size());                                                               \
        filtfilt(BUTTER_A, BUTTER_B, XPP, YPP);                                                       \
        resample(1, DT_NEW / DT, YPP, XPP);                                                           \
        MOVING_MEAN(XPP, YPP, nPoint_hal_win);                                                        \
        INIT_FFTW(fft_in, YPP, nPoint, nfft, fft_out);                                                \
        FFT_HELP_W(nfft, fft_in, fft_out, FFTW_FORWARD);                                              \
        double temp_f;                                                                                \
        for (int ii = 0; ii < nfft; ii++)                                                             \
        {                                                                                             \
            temp_f = sqrt(fft_out[ii][0] * fft_out[ii][0] + fft_out[ii][1] * fft_out[ii][1]) + 0.001; \
            fft_in[ii][0] = (fft_out[ii][0] + 0.001) / temp_f * shapingFilt[ii];                      \
            fft_in[ii][1] = (fft_out[ii][1]) / temp_f * shapingFilt[ii];                              \
            fft_out[ii][0] = 0;                                                                       \
            fft_out[ii][1] = 0;                                                                       \
        }                                                                                             \
        FFT_HELP_W(nfft, fft_in, fft_out, FFTW_BACKWARD);                                             \
        YPP.resize(nPoint);                                                                           \
        for (int i = 0; i < nPoint; i++)                                                              \
        {                                                                                             \
            YPP[i] = fft_out[i][0] / ((double)nfft);                                                  \
        }                                                                                             \
    }

//For debug only
void master_processing(std::vector<double> XPP, std::vector<double> &YPP)
{
    for (int i = 0; i < 100; i++)
    {
        printf("master_processing shapingFilt [%d] = %f\n", i, shapingFilt[i]);
    }
    detrend(&(XPP[0]), XPP.size());
    for (int i = 0; i < 10; i++)
    {
        printf("detrend (%d): %f \n", i, XPP[i]);
    }
    for (int i = 2999; i < 3010; i++)
    {
        printf("detrend (%d): %f, nfft =%d, nPoint =%d \n", i, XPP[i], nfft, nPoint);
    }

    for (int i = nPoint - 10 - 1; i < nPoint; i++)
    {
        printf("detrend (%d): %f, nfft =%d, nPoint =%d \n", i, XPP[i], nfft, nPoint);
    }

    filtfilt(BUTTER_A, BUTTER_B, XPP, YPP);
    for (int i = 0; i < 10; i++)
    {
        printf("filtfilt (%d): %f \n", i, YPP[i]);
    }
    for (int i = 2999; i < 3010; i++)
    {
        printf("filtfilt (%d): %f, nfft =%d, nPoint =%d \n", i, YPP[i], nfft, nPoint);
    }

    for (int i = nPoint - 10 - 1; i < nPoint; i++)
    {
        printf("filtfilt (%d): %f, nfft =%d, nPoint =%d \n", i, YPP[i], nfft, nPoint);
    }
    resample(1, DT_NEW / DT, YPP, XPP);
    MOVING_MEAN(XPP, YPP, nPoint_hal_win);
    for (int i = 0; i < 10; i++)
    {
        printf("MOVING_MEAN (%d): %f \n", i, YPP[i]);
    }
    for (int i = 2999; i < 3010; i++)
    {
        printf("MOVING_MEAN (%d): %f, nfft =%d, nPoint =%d \n", i, YPP[i], nfft, nPoint);
    }

    for (int i = nPoint - 10 - 1; i < nPoint; i++)
    {
        printf("MOVING_MEAN (%d): %f, nfft =%d, nPoint =%d \n", i, YPP[i], nfft, nPoint);
    }

    for (int i = 0; i < 100; i++)
    {
        printf("After MOVING_MEAN shapingFilt [%d] = %f\n", i, shapingFilt[i]);
    }
    fftw_complex *fft_in_temp;
    fftw_complex *fft_out_temp;
    fft_in_temp = fftw_alloc_complex(nfft);
    fft_out_temp = fftw_alloc_complex(nfft);
    for (int i = 0; i < nfft; i++)
    {
        if (i < nPoint)
        {
            fft_in_temp[i][0] = YPP[i];
            fft_in_temp[i][1] = 0;
        }
        else
        {
            fft_in_temp[i][0] = 0;
            fft_in_temp[i][1] = 0;
        }
        fft_out_temp[i][0] = 0;
        fft_out_temp[i][1] = 0;
    }
    //INIT_FFTW(fft_in, YPP, nPoint, nfft, fft_out);
    FFT_HELP_W(nfft, fft_in_temp, fft_out_temp, FFTW_FORWARD);
    double temp_f;
    for (int ii = 0; ii < nfft; ii++)
    {
        if (ii < 10)
            printf("fft: %f + %f i , shapingFilt = %f \n", fft_out_temp[ii][0], fft_out_temp[ii][1], shapingFilt[ii]);
        temp_f = sqrt(fft_out_temp[ii][0] * fft_out_temp[ii][0] + fft_out_temp[ii][1] * fft_out_temp[ii][1]) + 0.001;
        fft_in_temp[ii][0] = (fft_out_temp[ii][0] + 0.001) / temp_f * shapingFilt[ii];
        fft_in_temp[ii][1] = (fft_out_temp[ii][1]) / temp_f * shapingFilt[ii];
        fft_out_temp[ii][0] = 0;
        fft_out_temp[ii][1] = 0;
    }
    FFT_HELP_W(nfft, fft_in_temp, fft_out_temp, FFTW_BACKWARD);
    YPP.resize(nPoint);
    for (int i = 0; i < nPoint; i++)
    {
        YPP[i] = fft_out_temp[i][0] / ((double)nfft);
        if (i < 10)
            printf("White: %f \n", YPP[i]);
    }
    fftw_free(fft_in_temp);
    fftw_free(fft_out_temp);
}

template <class T>
inline void clear_vector(std::vector<T> &v)
{
    v.clear();
    std::vector<T>().swap(v);
}

double au_timer_global__inside_use = 0;
void au_time_start()
{
    au_timer_global__inside_use = MPI_Wtime();
}

//#define OMP_ENABLED 1
void au_time_elap(std::string info_str)
{

    double time_per_rank = MPI_Wtime() - au_timer_global__inside_use;
    int mpi_rank, mpi_size;
    double time_max, time_min, time_sum;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    //#ifndef OMP_ENABLED
    MPI_Allreduce(&time_per_rank, &time_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&time_per_rank, &time_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&time_per_rank, &time_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //#endif

    if (mpi_rank == 0)
    {
        printf("  %s:max=%f, min=%f, ave=%f, rank 0=%f\n", info_str.c_str(), time_max, time_min, time_sum / mpi_size, time_per_rank);
        fflush(stdout);
    }

    //reset to current time
    au_timer_global__inside_use = MPI_Wtime();
}

void au_time_elap(std::string info_str, int omp_rank)
{

    double time_per_rank = MPI_Wtime() - au_timer_global__inside_use;
    int mpi_rank, mpi_size;
    double time_max, time_min, time_sum;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    MPI_Allreduce(&time_per_rank, &time_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&time_per_rank, &time_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&time_per_rank, &time_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (mpi_rank == 0 && omp_rank == 0)
    {
        printf("  %s:max=%f, min=%f, ave=%f, rank 0=%f (omp_rank=0)\n", info_str.c_str(), time_max, time_min, time_sum / mpi_size, time_per_rank);
        fflush(stdout);
    }

    //reset to current time
    au_timer_global__inside_use = MPI_Wtime();
}

time_t au_timer_global_start__inside_use_no_mpi;
void au_time_start_no_mpi()
{
    time(&au_timer_global_start__inside_use_no_mpi);
}

void au_time_elap_no_mpi(std::string info_str)
{
    time_t current_time;
    time(&current_time);
    double time_taken = double(current_time - au_timer_global_start__inside_use_no_mpi);

    cout << info_str << std::fixed << time_taken << std::setprecision(5);
    cout << " sec " << endl;

    //reset timer
    time(&au_timer_global_start__inside_use_no_mpi);
}

#endif
