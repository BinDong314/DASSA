#ifndef DASLIB_H
#define DASLIB_H

#include <complex>
#include <cmath>
#include <fftw3.h>

extern double micro_time, sum_micro, micro_time_sub, sum_micro_sub;

using namespace std;

#include "DasLib_filtfilt.h"
#include "Tools3rd/INIReader.h"
#include "Tools3rd/termcolor.hpp"
#include "DasLib/DasLib_resample.h"
#include "DasLib/DasLib_detrend.h"
#include "DasLib/DasLib_interp1.h"
#include "DasLib/DasLib_liir.c"
#include "DasLib/DasLib_fft.h"
#include "DasLib/DasLib_moving_mean.h"
#include "DasLib/DasLib_minpower2.h"

// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes?answertab=votes#tab-top
#include <numeric>   // std::iota
#include <algorithm> // std::sort, std::stable_sort
#include <fstream>   // std::ifstream

// https://docs.microsoft.com/en-us/cpp/cpp/namespaces-cpp?redirectedfrom=MSDN&view=vs-2019
#pragma once
namespace DasLib
{
    template <typename T>
    inline T find_max_neighbors(const std::vector<T> &vec, int i, int numNeighbors)
    {
        int n = vec.size();
        int start = std::max(0, i - numNeighbors);
        int end = std::min(n - 1, i + numNeighbors);
        T maxNeigh = std::numeric_limits<T>::min();
        for (int j = start; j <= end; j++)
        {
            maxNeigh = std::max(maxNeigh, vec[j]);
        }
        return maxNeigh;
    }

    // Computes the cross correlation coefficients for series_x and series_y
    template <class T>
    inline void xcross(const std::vector<T> &X, const std::vector<T> &Y, std::vector<double> &xcross_result)
    {

        assert(X.size() == Y.size());
        int N = X.size();
        xcross_result.resize(N * 2 - 1, 0);
        double sum_temp;
        for (int delay = -N + 1; delay <= N; delay++)
        {
            sum_temp = 0;
            for (int i = 0; i < N; i++)
            {
                int j = i - delay;
                if (j < 0 || j >= N)
                {
                    sum_temp += 0;
                }
                else
                {
                    sum_temp += X[i] * Y[j];
                }
            }
            xcross_result[delay + N - 1] = sum_temp;
        }
    }

    template <class T>
    inline double xcross_max(const std::vector<T> &X, const std::vector<T> &Y)
    {

        assert(X.size() == Y.size());
        int N = X.size();
        double sum_temp, max_temp = -1 * numeric_limits<double>::max();
        for (int delay = -N + 1; delay <= N; delay++)
        {
            sum_temp = 0;
            for (int i = 0; i < N; i++)
            {
                int j = i - delay;
                if (j < 0 || j >= N)
                {
                    sum_temp += 0;
                }
                else
                {
                    sum_temp += X[i] * Y[j];
                }
            }
            if (max_temp < sum_temp)
            {
                max_temp = sum_temp;
            }
        }

        return max_temp;
    }

    template <class T>
    inline double xcross_fft(const std::vector<T> &X, const std::vector<T> &Y)
    {
        assert(X.size() == Y.size());
        int winlen = X.size();
        int xclen = 2 * winlen - 1;

        std::vector<T> XX(X.begin(), X.end());
        std::vector<T> YY(Y.rbegin(), Y.rend());

        XX.resize(xclen, 0);
        YY.resize(xclen, 0);

        // taking fft
        std::vector<std::complex<double>> XX_fft;
        fftv_forward_omp(XX, XX_fft);

        std::vector<std::complex<double>> YY_fft;
        fftv_forward_omp(YY, YY_fft);

        int nfft = floor(xclen / 2) + 1;

        std::vector<std::complex<double>> XXYY_fft_corr;
        size_t nfft_doubled = 2 * nfft - 1;
        XXYY_fft_corr.resize(nfft_doubled);
        for (int i = 0; i < nfft; i++)
        {
            XXYY_fft_corr[i] = XX_fft[i] * YY_fft[i];
            if (i != 0)
                XXYY_fft_corr[nfft_doubled - i] = std::conj(XXYY_fft_corr[i]);
        }

        double ifft_out_real_max;
        fftv_backward_real_max_omp(XXYY_fft_corr, ifft_out_real_max);

        return ifft_out_real_max;
    }

    // Function for calculating median
    template <class T>
    inline T MaxVector(std::vector<T> &v)
    {
        return *std::max_element(v.begin(), v.end());
    }

    template <class T>
    inline T MaxVectorVector(std::vector<std::vector<T>> &v)
    {
        std::vector<T> max_v;
        max_v.resize(v.size());
        for (size_t i = 0; i < v.size(); i++)
        {
            max_v[i] = MaxVector(v[i]);
        }
        return MaxVector(max_v);
    }

    // Function for calculating median
    template <class T>
    inline T Median(std::vector<T> &v)
    {
        sort(v.begin(), v.end());
        size_t n = v.size();
        // check for even case
        if (n % 2 != 0)
            return v[n / 2];

        return (v[(n - 1) / 2] + v[n / 2]) / 2.0;
    }

    // To be complete
    template <class T>
    int decimate(std::vector<T> &v, const int r, const std::vector<double> &cheby1_b, const std::vector<double> &cheby1_a)
    {
        size_t nd = v.size();
        size_t nout = std::ceil(((double)nd) / ((double)r));
        size_t nbeg = r - (r * nout - nd);

        std::vector<T> decimated;
        // decimated.resize(nd);
        filtfilt(cheby1_b, cheby1_a, v, decimated);
        // decimated.reserve(v.size() / M);

        // PrintVector("decimate::decimated=", decimated);
        // std::cout << "nd =" << nd << ", nbeg = " << nbeg << " , nout = " << nout << " , r = " << r << "\n";
        v.resize(nout);
        size_t v_i = 0;
        for (size_t i = nbeg - 1; i < nd; i = i + r)
        {
            v[v_i] = decimated[i];
            v_i++;
        }
        return 0;
    }

    /**
     * @brief perform the decimate from "start_row" to "end_row"
     *
     * @tparam T
     * @param data2d
     * @param start_row
     * @param end_row
     * @param operator
     * @return std::vector<T>
     */
    template <class T>
    inline std::vector<T> spacedecimate(const std::vector<std::vector<T>> &data2d, const int start_row, const int end_row, const std::string stat_operator_str)
    {
        std::vector<T> result;
        size_t rows = data2d.size();
        size_t cols = data2d[0].size();

        T stat_temp;
        for (std::size_t j = 0; j < cols; j++)
        {
            if (stat_operator_str == "mean" || stat_operator_str == "average" || stat_operator_str == "ave")
            {
                stat_temp = 0;
                for (std::size_t i = start_row; i <= end_row; i++)
                {
                    stat_temp = stat_temp + data2d[i][j];
                }
                stat_temp = stat_temp / (end_row - start_row + 1);
            }
            else if (stat_operator_str == "maximum" || stat_operator_str == "max")
            {
                stat_temp = data2d[start_row][j];
                for (std::size_t i = start_row + 1; i <= end_row; i++)
                {
                    if (data2d[i][j] > stat_temp)
                        stat_temp = data2d[i][j];
                }
            }
            else if (stat_operator_str == "minumum" || stat_operator_str == "min")
            {
                stat_temp = data2d[start_row][j];
                for (std::size_t i = start_row + 1; i <= end_row; i++)
                {
                    if (data2d[i][j] < stat_temp)
                        stat_temp = data2d[i][j];
                }
            }
            else if (stat_operator_str == "median" || stat_operator_str == "med")
            {
                std::vector<T> median_vec_temp;
                for (std::size_t i = start_row; i <= end_row; i++)
                {
                    median_vec_temp.push_back(data2d[i][j]);
                }
                stat_temp = Median(median_vec_temp);
            }
            else if (stat_operator_str == "summary" || stat_operator_str == "sum")
            {
                stat_temp = 0;
                for (std::size_t i = start_row; i <= end_row; i++)
                {
                    stat_temp = stat_temp + data2d[i][j];
                }
            }

            result.push_back(stat_temp);
        }
        return result;
    }

    // template <typename T>
    inline double sum_weight(const std::vector<double> &v, const std::vector<double> &w)
    {
        assert(v.size() == w.size());
        // PrintVector("v = ", v);
        // PrintVector("w = ", w);

        double sum = 0;
        size_t n = v.size();
        for (size_t i = 0; i < n; i++)
        {
            // std::cout << i << ",  sum = " << sum << ", v[i] = " << v[i] << ", w[i]= " << w[i] << "\n";
            sum = sum + v[i] * w[i];
        }
        // std::cout << "sum =" << sum << "\n";
        return sum;
    }

    // template <typename T>
    inline void sum_weight_by_time(const std::vector<std::vector<double>> &v, const std::vector<double> &w, std::vector<double> &out)
    {
        // v is channel x time
        // w is channel x 1
        // out is template x time
        assert(v.size() == w.size());
        // assert(v.size() == out.size());

        // double sum = 0;
        size_t n_chs = v.size();
        size_t n_points;
        if (n_chs > 0)
        {
            n_points = v[0].size();
        }
        else
        {
            AU_EXIT("Error empty input v")
        }
        assert(n_points <= out.size());

#if defined(_OPENMP)
#pragma omp parallel for
#endif
        for (int j = 0; j < n_points; j++)
        {
            out[j] = 0;
            for (size_t i = 0; i < n_chs; i++)
            {
                // std::cout << i << ",  sum = " << sum << ", v[i] = " << v[i] << ", w[i]= " << w[i] << "\n";
                out[j] = out[j] + v[i][j] * w[i];
            }
        }
        // std::cout << "sum =" << sum << "\n";
        // return sum;
    }

#define VectorElementMulti(v1, v2)          \
    {                                       \
        assert(v1.size() == v2.size());     \
        for (int i = 0; i < v1.size(); i++) \
        {                                   \
            v1[i] = v1[i] * v2[i];          \
        }                                   \
    }

#define VectorElementMultiNormal(v1, v2)    \
    {                                       \
        double v_sum = 0;                   \
        assert(v1.size() == v2.size());     \
        for (int i = 0; i < v1.size(); i++) \
        {                                   \
            v1[i] = v1[i] * v2[i];          \
            v_sum = v_sum + v1[i] * v1[i];  \
        }                                   \
        double v_sum_sqrt = sqrt(v_sum);    \
        for (int i = 0; i < v1.size(); i++) \
        {                                   \
            v1[i] = v1[i] / v_sum_sqrt;     \
        }                                   \
    }

#define VectorDivideByScalar(v, k)                                                  \
    {                                                                               \
        transform(v.begin(), v.end(), v.begin(), [k](double &c) { return c / k; }); \
    }

#define VectorMinusByScalar(v, k)                                                   \
    {                                                                               \
        transform(v.begin(), v.end(), v.begin(), [k](double &c) { return c - k; }); \
    }

    // atemp=decimate(decimate([((detrend(double(amat0(:,rc1)))).*ctap0); 0],5),2);
    // atemp=filtfilt(bfsos1,bfg1,atemp);

    // detrend, decimate, filtfilt
    template <typename T>
    inline std::vector<T> ddff(std::vector<T> &in_vec, const std::vector<double> &ctap, const int decimate_factor, const std::vector<double> &filtfilt_B, const std::vector<double> &filtfilt_A, const std::vector<double> &decimate_cheby1_b, const std::vector<double> &decimate_cheby1_a)
    {
        std::vector<T> out_vec;

        assert(in_vec.size() == ctap.size());

        detrend(in_vec.data(), in_vec.size()); // Detread

        // PrintVector("ddff:in_vec after detrend=", in_vec);

        for (int j = 0; j < ctap.size(); j++)
        {
            in_vec[j] = in_vec[j] * ctap[j];
        }
        in_vec.push_back(0); // add one zero

        // PrintVector("ddff:in_vec after ctap=", in_vec);

        decimate(in_vec, decimate_factor, decimate_cheby1_b, decimate_cheby1_a);

        // PrintVector("ddff:in_vec after decimate=", in_vec);

        filtfilt(filtfilt_B, filtfilt_A, in_vec, out_vec);

        // PrintVector("ddff:out_vec after filtfilt=", out_vec);

        return out_vec;
    }

    template <typename T>
    inline void norm_matlab(vector<T> &v)
    {
        T x = 0;
        for (size_t i = 0; i < v.size(); i++)
            x = x + v[i] * v[i];
        T sqrt_x = sqrt(x);
        for (size_t i = 0; i < v.size(); i++)
            v[i] = v[i] / sqrt_x;
    }

    /*   subset, detrend, ctap, norm
        slice(atemp2, template_tstart[rc2][i] - 1, template_tstart[rc2][i] - 1 + template_winlen[rc2] - 1);
            std::cout << " atemp2.size() = " << atemp2.size() << ", template_tstart[rc2][i] = " << template_tstart[rc2][i] << ", template_winlen[rc2] = " << template_winlen[rc2] << ",  ctap_template2.size= " << ctap_template2.size() << "\n";
            detrend(atemp2.data(), atemp2.size());
            VectorElementMulti(atemp2, ctap_template2);
            double atemp2_norm = norm_matlab(atemp2);
            VectorDivideByScalar(atemp2, atemp2_norm);


            Matlab code:
            atemp3=(detrend(amat1(idx1:(idx1+template_winlen(rc2)-1),rc1))).*ctap_template2;
            atemp3=atemp3./norm(atemp3);

     */
    template <typename T>
    inline void sdcn_old(const std::vector<T> &v, std::vector<T> &v_out, const size_t subset_start, const size_t subset_count, const std::vector<T> &ctap)
    {
        // slice(atemp2, template_tstart[rc2][i] - 1, template_tstart[rc2][i] - 1 + template_winlen[rc2] - 1);
        // std::vector<T>::const_iterator first = v.begin() + subset_start;
        // std::vector<T>::const_iterator last = v.begin() + subset_start + subset_count;
        // std::cout << "subset_start = " << subset_start << ", subset_count = " << subset_count << "\n";
        std::vector<T> newV(v.begin() + subset_start, v.begin() + subset_start + subset_count);
        detrend(newV.data(), newV.size());
        // std::cout << "newV.size() = " << newV.size() << ", ctap = " << ctap.size() << "\n";
        // VectorElementMulti(newV, ctap);
        // double atemp2_norm = norm_matlab(newV);
        // VectorDivideByScalar(newV, atemp2_norm);
        // norm_matlab(newV);
        VectorElementMultiNormal(newV, ctap);
        v_out = newV;
    }

    template <typename T>
    inline void sdcn(const std::vector<T> &v, std::vector<T> &v_out, const size_t subset_start, const size_t subset_count, const std::vector<T> &ctap)
    {
        // slice(atemp2, template_tstart[rc2][i] - 1, template_tstart[rc2][i] - 1 + template_winlen[rc2] - 1);
        // std::vector<T>::const_iterator first = v.begin() + subset_start;
        // std::vector<T>::const_iterator last = v.begin() + subset_start + subset_count;
        // std::cout << "subset_start = " << subset_start << ", subset_count = " << subset_count << "\n";
        v_out.resize(subset_count);
        // for (size_t i = 0; i < subset_count; i++)
        //{
        //     v_out[i] = v[i + subset_start];
        // }
        // memcpy(&v_out[0], &v[subset_start], subset_count * sizeof(T));
        // std::vector<T> newV(v.begin() + subset_start, v.begin() + subset_start + subset_count);
        // detrend(newV.data(), newV.size());
        // detrend(v_out.data(), v_out.size());
        // PrintVector("Before detrend_range v_out =", v);

        detrend_range(v, subset_start, subset_count, ctap, v_out);
        // std::cout << "newV.size() = " << newV.size() << ", ctap = " << ctap.size() << "\n";
        // VectorElementMulti(newV, ctap);
        // double atemp2_norm = norm_matlab(newV);
        // VectorDivideByScalar(newV, atemp2_norm);
        // norm_matlab(newV);
        // PrintVector("Before VectorElementMultiNormal v_out =", v_out);

        // VectorElementMultiNormal(v_out, ctap);

        // PrintVector("After VectorElementMultiNormal v_out =", v_out);

        // PrintVector("After VectorElementMultiNormal v_out =", v_out);
        //  v_out = newV;
    }

    template <typename T>
    inline T dot_product(const std::vector<T> &vector_a, const std::vector<T> &vector_b)
    {
        assert(vector_a.size() == vector_b.size());
        T product = 0;
        for (size_t i = 0; i < vector_a.size(); i++)
            product = product + vector_a[i] * vector_b[i];
        return product;
    }

    template <typename T>
    vector<size_t> sort_indexes(const vector<T> &v)
    {

        // initialize original index locations
        vector<size_t> idx(v.size());
        iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        // using std::stable_sort instead of std::sort
        // to avoid unnecessary index re-orderings
        // when v contains elements of equal values
        stable_sort(idx.begin(), idx.end(),
                    [&v](size_t i1, size_t i2)
                    { return v[i1] > v[i2]; });

        return idx;
    }

    // n = order of the filter
    // fc = filter cutoff frequency as a fraction of Pi [0,1]
    int ButterLow(int n, double fcf, std::vector<double> &A, std::vector<double> &B)
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

    /**********************************************************************
      dcof_bwbp - calculates the d coefficients for a butterworth bandpass
      filter. The coefficients are returned as an array of doubles.

    */

    double *dcof_bwbp(int n, double f1f, double f2f)
    {
        int k;        // loop variables
        double theta; // M_PI * (f2f - f1f) / 2.0
        double cp;    // cosine of phi
        double st;    // sine of theta
        double ct;    // cosine of theta
        double s2t;   // sine of 2*theta
        double c2t;   // cosine 0f 2*theta
        double *rcof; // z^-2 coefficients
        double *tcof; // z^-1 coefficients
        double *dcof; // dk coefficients
        double parg;  // pole angle
        double sparg; // sine of pole angle
        double cparg; // cosine of pole angle
        double a;     // workspace variables

        cp = cos(M_PI * (f2f + f1f) / 2.0);
        theta = M_PI * (f2f - f1f) / 2.0;
        st = sin(theta);
        ct = cos(theta);
        s2t = 2.0 * st * ct;       // sine of 2*theta
        c2t = 2.0 * ct * ct - 1.0; // cosine of 2*theta

        rcof = (double *)calloc(2 * n, sizeof(double));
        tcof = (double *)calloc(2 * n, sizeof(double));

        for (k = 0; k < n; ++k)
        {
            parg = M_PI * (double)(2 * k + 1) / (double)(2 * n);
            sparg = sin(parg);
            cparg = cos(parg);
            a = 1.0 + s2t * sparg;
            rcof[2 * k] = c2t / a;
            rcof[2 * k + 1] = s2t * cparg / a;
            tcof[2 * k] = -2.0 * cp * (ct + st * sparg) / a;
            tcof[2 * k + 1] = -2.0 * cp * st * cparg / a;
        }

        dcof = trinomial_mult(n, tcof, rcof);
        free(tcof);
        free(rcof);

        dcof[1] = dcof[0];
        dcof[0] = 1.0;
        for (k = 3; k <= 2 * n; ++k)
            dcof[k] = dcof[2 * k - 2];
        return (dcof);
    }

    /**********************************************************************
      ccof_bwbp - calculates the c coefficients for a butterworth bandpass
      filter. The coefficients are returned as an array of integers.

    */

    int *ccof_bwbp(int n)
    {
        int *tcof;
        int *ccof;
        int i;

        ccof = (int *)calloc(2 * n + 1, sizeof(int));
        if (ccof == NULL)
            return (NULL);

        tcof = ccof_bwhp(n);
        if (tcof == NULL)
            return (NULL);

        for (i = 0; i < n; ++i)
        {
            ccof[2 * i] = tcof[i];
            ccof[2 * i + 1] = 0.0;
        }
        ccof[2 * n] = tcof[n];

        free(tcof);
        return (ccof);
    }

    /**********************************************************************
  sf_bwbp - calculates the scaling factor for a butterworth bandpass filter.
  The scaling factor is what the c coefficients must be multiplied by so
  that the filter response has a maximum value of 1.

*/
    double sf_bwbp(int n, double f1f, double f2f)
    {
        int k;           // loop variables
        double ctt;      // cotangent of theta
        double sfr, sfi; // real and imaginary parts of the scaling factor
        double parg;     // pole angle
        double sparg;    // sine of pole angle
        double cparg;    // cosine of pole angle
        double a, b, c;  // workspace variables

        ctt = 1.0 / tan(M_PI * (f2f - f1f) / 2.0);
        sfr = 1.0;
        sfi = 0.0;

        for (k = 0; k < n; ++k)
        {
            parg = M_PI * (double)(2 * k + 1) / (double)(2 * n);
            sparg = ctt + sin(parg);
            cparg = cos(parg);
            a = (sfr + sfi) * (sparg - cparg);
            b = sfr * sparg;
            c = -sfi * cparg;
            sfr = b - c;
            sfi = a - b - c;
        }

        return (1.0 / sfr);
    }

    // https://github.com/wustyuyi/Butterworth_filter_coefficients-MATLAB-in-C
    int ButterPass(int n, const std::vector<double> &fcf, std::vector<double> &A, std::vector<double> &B)
    {
        int sff = 1;  // scale flag: 1 to scale, 0 to not scale ccof
        int i;        // loop variables
        double f1f;   // lower cutoff frequency (fraction of pi)
        double f2f;   // upper cutoff frequency (fraction of pi)
        double sf;    // scaling factor
        double *dcof; // d coefficients
        int *ccof;    // c coefficients
        // FILE *fp;     // output file

        /*
        if (argc < 6)
        {
            printf("\nbwbp calculates Butterworth bandpass filter coefficients\n");
            printf("\nUsage: bwbp n fc1 fc2 sf outfile\n");
            printf("  n = order of the filter\n");
            printf("  fc1 = lower cutoff frequency as a fraction of Pi [0,1]\n");
            printf("  fc2 = upper cutoff frequency as a fraction of Pi [0,1]\n");
            printf("  sf = 1 to scale c coefficients for normalized response\n");
            printf("  sf = 0 to not scale c coefficients\n");
            printf("  outfile = output file name\n");
            return (-1);
        }*/

        f1f = fcf[0];
        f2f = fcf[1];

        /* calculate the d coefficients */
        dcof = dcof_bwbp(n, f1f, f2f);
        if (dcof == NULL)
        {
            perror("Unable to calculate d coefficients");
            return (-1);
        }

        /* calculate the c coefficients */
        ccof = ccof_bwbp(n);
        if (ccof == NULL)
        {
            perror("Unable to calculate c coefficients");
            return (-1);
        }

        sf = sf_bwbp(n, f1f, f2f); /* scaling factor for the c coefficients */
        // printf("scaling factor  = %f\n", sf);
        /* Output the c coefficients */
        // fprintf(fp, "%d\n", 2 * n + 1); // number of c coefficients
        A.resize(2 * n + 1);
        if (sff == 0)
        {
            for (i = 0; i <= 2 * n; ++i)
            {
                A[i] = ccof[i];
            }
        }
        else
        {
            for (i = 0; i <= 2 * n; ++i)
            {
                A[i] = (double)ccof[i] * sf;
            }
        }
        /* Output the d coefficients */
        // fprintf(fp, "%d\n", 2 * n + 1); /* number of d coefficients */
        B.resize(2 * n + 1);
        for (i = 0; i <= 2 * n; ++i)
        {
            B[i] = dcof[i];
        }
        // fprintf(fp, "%1.15lf\n", dcof[i]);

        if (dcof != NULL)
        {
            free(dcof);
            dcof = NULL;
        }
        if (ccof != NULL)
        {
            free(ccof);
            ccof = NULL;
        }

        return 0;
    }

    // double rms(double x[], int n)
    template <class T>
    inline T Rms(const std::vector<T> &x)
    {
        T sum = 0;
        int n = x.size();
        for (int i = 0; i < n; i++)
            sum += pow(x[i], 2);

        return sqrt(sum / n);
    }

    template <class T>
    inline T Mean(std::vector<T> &v)
    {
        T sum = 0;
        for (int i = 0; i < v.size(); i++)
        {
            sum = sum + v[i];
        }
        return sum / v.size();
    }

    template <class T>
    inline void DeleteMedian(std::vector<std::vector<T>> &v)
    {
        int rows = v.size();
        int cols = v[0].size();

        // std::cout << "rows = " << rows << ", cols = " << cols << "\n";
        std::vector<T> temp_cols(rows);
        T temp_median;
        for (int j = 0; j < cols; j++)
        {
            for (int i = 0; i < rows; i++)
            {
                temp_cols[i] = v[i][j];
            }
            temp_median = Median(temp_cols);
            for (int i = 0; i < rows; i++)
            {
                v[i][j] = v[i][j] - temp_median;
            }
        }
    }

    /**
     * @brief
     *
     * @tparam T
     * @param v input time series data
     * @param start_t, start time of correlation fh(9) in Matlab code, e.g., -60
     * @param end_t,  end time of correlation fh(9) in Matlab code, e.g., 60
     * @param sub_start_t, start time of subset
     * @param sub_end_t, end time of subset
     * @param smaple_rate, interval of each point, fh(8) in Matlab code , 0.008
     * @return std::vector<T>  return sub array
     */
    template <class T>
    inline std::vector<T> TimeSubset(std::vector<T> &v, double start_t, double end_t, double sub_start_t, double sub_end_t, double smaple_rate)
    {
        size_t sub_start_index = round((sub_start_t - start_t) / smaple_rate);
        size_t sub_end_index = round((sub_end_t - start_t) / smaple_rate);

        // std::cout << "sub_start_index = " << sub_start_index << ", sub_end_index = " << sub_end_index << "\n";
        assert(sub_start_index >= 0);
        assert(sub_end_index >= 0);
        assert(sub_end_index >= sub_start_index);
        std::vector<T> v_subset;
        v_subset.resize(sub_end_index - sub_start_index + 1);

        size_t v_subset_index = 0;
        for (int i = sub_start_index; i <= sub_end_index; i++)
        {
            v_subset[v_subset_index] = v[i];
            v_subset_index++;
        }

        return v_subset;
    }
    /**
     * @brief infer the size of vector after time subset
     *
     * @param start_t
     * @param end_t
     * @param sub_start_t
     * @param sub_end_t
     * @param smaple_rate
     * @return size_t
     */
    inline size_t InferTimeSubsetSize(double start_t, double end_t, double sub_start_t, double sub_end_t, double smaple_rate)
    {
        size_t sub_start_index = round((sub_start_t - start_t) / smaple_rate + 1);
        size_t sub_end_index = round((sub_end_t - start_t) / smaple_rate + 1);
        /*
        std::cout << "sub_start_index = " << sub_start_index
                  << ", sub_end_index =  " << sub_end_index
                  << ", start_t =" << start_t
                  << ", end_t =" << end_t
                  << ", sub_start_t = " << sub_start_t
                  << ", sub_end_t = " << sub_end_t
                  << ", smaple_rate = " << smaple_rate
                  << "\n";
        */
        return sub_end_index - sub_start_index + 1;
    }

    typedef std::vector<std::vector<double>> DoubleVector2D;

    /**
     * @brief Convert a 1D vector to 2D vector
     *
     * @tparam T1
     * @tparam T2
     * @param cols
     * @param data1d : is organized as row-major, i.e., row is cotigious
     *                  row 1 , row 2, row 3...
     * @return std::vector<std::vector<T2>>
     */
    template <class T1, class T2 = double>
    inline std::vector<std::vector<T2>> Vector1D2D(size_t cols, const std::vector<T1> &data1d)
    {
        std::vector<std::vector<T2>> result;
        size_t rows = data1d.size() / cols;
        result.resize(rows);
        for (std::size_t i = 0; i < rows; ++i)
        {
            result[i].resize(cols);
            for (std::size_t j = 0; j < cols; ++j)
            {
                result[i][j] = data1d[i * cols + j];
            }
        }
        return result;
    }

    /**
     * @brief Data is organized by column order
     *
     * @tparam T1
     * @tparam T2
     * @param cols: is the number of cols in original data
     * @param data1d
     * @return std::vector<std::vector<T2>>
     */
    template <class T1, class T2 = double>
    inline std::vector<std::vector<T2>> Vector1D2DByCol(size_t cols, const std::vector<T1> &data1d)
    {
        std::vector<std::vector<T2>> result;
        size_t rows = data1d.size() / cols;
        result.resize(cols);
        for (std::size_t i = 0; i < cols; ++i)
        {
            result[i].resize(rows);
            for (std::size_t j = 0; j < rows; ++j)
            {
                result[i][j] = data1d[j * cols + i];
            }
        }
        return result;
    }

    /**
     * @brief Data is organized by column order
     *
     * @tparam T1
     * @tparam T2
     * @param cols: is the number of cols in original data
     * @param data1d
     * @return std::vector<std::vector<T2>>
     */
    template <class T1, class T2 = double>
    inline std::vector<std::vector<T2>> Vector1D2DByColStride(size_t cols, const std::vector<T1> &data1d, const size_t start, const size_t stride)
    {
        std::vector<std::vector<T2>> result;

        size_t rows = data1d.size() / cols; // the number of rows
        size_t new_cols = (cols - start) / stride + 1;
        result.resize(new_cols);

        for (int i = 0; i < new_cols; i++)
            result[i].resize(rows);

#if defined(_OPENMP)
#pragma omp parallel for
#endif
        for (std::size_t j = 0; j < rows; ++j)
        {
            for (std::size_t i = start - 1; i < cols; i = i + stride)
            {
                result[i / stride][j] = data1d[j * cols + i];
            }
        }
        return result;
    }

    /**
     * @brief Convert a 1D vector to 2D vector
     *
     * @tparam T1
     * @tparam T2
     * @param cols
     * @param data1d : is organized as row-major, i.e., row is cotigious
     *                  row 1 , row 2, row 3...
     *  @param row_layout_flag: extrat row or column as 2D vector
     * @return std::vector<std::vector<T2>>
     */
    template <class T1, class T2 = double>
    inline std::vector<std::vector<T2>> ConvertVector1DTo2D(const std::vector<T1> &data1d, const size_t rows, const size_t cols, const bool row_layout_flag)
    {
        std::vector<std::vector<T2>> result;
        if (row_layout_flag)
        {
            // size_t rows = data1d.size() / cols;
            result.resize(rows);
            for (std::size_t i = 0; i < rows; ++i)
            {
                result[i].resize(cols);
                for (std::size_t j = 0; j < cols; ++j)
                {
                    result[i][j] = data1d[i * cols + j];
                }
            }
        }
        else
        {
            // size_t rows = data1d.size() / cols;
            result.resize(cols);
            for (std::size_t i = 0; i < cols; ++i)
            {
                result[i].resize(rows);
                for (std::size_t j = 0; j < rows; ++j)
                {
                    result[i][j] = data1d[j * cols + i];
                }
            }
        }

        return result;
    }
    /*
template <class T>
inline std::vector<std::vector<T>> Vector1D2D(size_t cols, std::vector<T> &data1d)
{
    std::vector<std::vector<T>> result;
    size_t rows = data1d.size() / cols;
    result.resize(rows);
    for (std::size_t i = 0; i < rows; ++i)
    {
        result[i].resize(cols);
        for (std::size_t j = 0; j < cols; ++j)
        {
            result[i][j] = data1d[i * cols + j];
        }
    }
    return result;
}*/

    template <class T>
    inline std::vector<T> Convert2DVTo1DV(std::vector<std::vector<T>> &data2d)
    {
        std::vector<T> result;
        size_t rows = data2d.size();
        size_t cols = data2d[0].size();

        result.reserve(rows * cols);
        for (std::size_t i = 0; i < rows; ++i)
        {
            copy(data2d[i].begin(), data2d[i].end(), back_inserter(result));
        }
        assert(result.size() == (rows * cols));
        return result;
    }

    template <class T>
    inline std::vector<T> Convert3DVTo1DV(std::vector<std::vector<std::vector<T>>> &data3d)
    {
        std::vector<T> result;
        size_t size_dim0 = data3d.size();
        size_t size_dim1 = data3d[0].size();
        size_t size_dim2 = data3d[0][0].size();

        result.reserve(size_dim0 * size_dim1 * size_dim2);
        for (std::size_t i = 0; i < size_dim0; ++i)
        {
            for (std::size_t j = 0; j < size_dim1; ++j)
            {
                copy(data3d[i][j].begin(), data3d[i][j].end(), back_inserter(result));
            }
        }
        assert(result.size() == (size_dim0 * size_dim1 * size_dim2));
        return result;
    }

    template <class T>
    inline std::vector<std::vector<std::vector<T>>> Convert1DVTo3DV(std::vector<T> data1d, size_t dim_0_size, size_t dim_1_size, size_t dim_2_size)
    {
        std::vector<std::vector<std::vector<T>>> result;
        result.resize(dim_0_size);
        for (int i = 0; i < dim_0_size; i++)
        {
            result[i].resize(dim_1_size);
            for (int j = 0; j < dim_1_size; j++)
            {
                result[i][j].resize(dim_2_size);
                for (int k = 0; k < dim_2_size; k++)
                    result[i][j][k] = data1d[k + dim_2_size * (j + dim_1_size * i)];
            }
        }

        return result;
    }

    /**
     * @brief
     *
     * @tparam T
     * @param v
     * @param m: zero based and inclusive
     * @param n: zero based and inclusive
     * @return std::vector<T>
     */
    template <typename T>
    void slice(std::vector<T> &v, int m, int n)
    {
        // auto first = v.cbegin() + m;
        // auto last = v.cbegin() + n + 1;
        // std::vector<T> vec(first, last);
        // return vec;
        v.erase(v.begin() + n + 1, v.end());
        v.erase(v.begin(), v.begin() + m);
    }
    template <typename T>
    std::vector<T> slice2(std::vector<T> &v, int m, int n)
    {
        auto first = v.cbegin() + m;
        auto last = v.cbegin() + n + 1;
        std::vector<T> vec(first, last);
        return vec;
    }

    template <class T>
    inline bool CausalityFlagging(std::vector<std::vector<T>> &v, double tmin, double tmax, double fmax, double start_t, double end_t, double smaple_rate, double CausalityFlagging_ButterLow_order, double CausalityFlagging_ButterLow_fcf)
    {
        int N = v.size();
        std::vector<std::vector<T>> tv1, tv2, tv3, tv4;
        tv1.resize(N);
        tv2.resize(N);
        tv3.resize(N);
        tv4.resize(N);

        std::vector<double> A, B;
        ButterLow(CausalityFlagging_ButterLow_order, CausalityFlagging_ButterLow_fcf, A, B);

        /*PrintScalar("N", N);
    PrintVector("A", A);
    PrintVector("B", B);
    PrintScalar("start_t", start_t);
    PrintScalar("end_t", end_t);
    PrintScalar("smaple_rate", smaple_rate);
    PrintScalar("tmax", tmax);*/

        for (int i = 0; i < N; i++)
        {
            tv1[i] = TimeSubset(v[i], start_t, end_t, -tmax, tmax, smaple_rate);
            filtfilt(A, B, tv1[i], tv2[i]);
            tv3[i] = TimeSubset(tv2[i], -tmax, tmax, tmin, tmax * 0.9, smaple_rate);   // causal
            tv4[i] = TimeSubset(tv2[i], -tmax, tmax, -tmax * 0.9, -tmin, smaple_rate); // acausal
        }
        std::vector<double> rms_acausal(N), rms_causal(N);

        for (int i = 0; i < N; i++)
        {
            rms_causal[i] = Rms(tv3[i]);
            rms_acausal[i] = Rms(tv4[i]);
        }

        // PrintVector("rms_causal: ", rms_causal);
        // PrintVector("rms_acausal: ", rms_acausal);

        double rms_acausal_mean = Mean(rms_acausal);
        double rms_causal_mean = Mean(rms_causal);
        // std::cout << "rms_causal_mean = " << rms_causal_mean << ", rms_causal_mean =" << rms_causal_mean << "\n";
        if (rms_acausal_mean < rms_causal_mean)
        {
            return true;
        }
        else
        {
            return false;
        }
    } // namespace DasLib

    template <class T>
    inline std::vector<std::complex<T>> Hilbert(std::vector<T> &in)
    {
        int INN = in.size();
        // std::cout << "INN =" << INN << "\n";
        fftw_complex *in_temp = (fftw_complex *)fftw_alloc_complex(sizeof(fftw_complex) * INN);
        // std::memset(in_temp, 0, sizeof(fftw_complex) * INN);

        for (size_t i = 0; i < INN; i++)
        {
            in_temp[i][0] = in[i];
            in_temp[i][1] = 0;
        }

        fftw_complex *out_temp = (fftw_complex *)fftw_alloc_complex(sizeof(fftw_complex) * INN);
        // std::memset(out_temp, 0, sizeof(fftw_complex) * INN);

        // micro_time = AU_WTIME;
        fftw_plan pf = fftw_plan_dft_1d(INN, in_temp, out_temp, FFTW_FORWARD, FFTW_ESTIMATE);
        // micro_time_sub = AU_WTIME;
        fftw_execute(pf);
        // sum_micro_sub = sum_micro_sub + (AU_WTIME - micro_time_sub);

        fftw_destroy_plan(pf);
        // fftw_cleanup();
        // sum_micro = sum_micro + (AU_WTIME - micro_time);

        size_t HN = INN >> 1;
        size_t numRem = HN;
        for (size_t i = 1; i < HN; ++i)
        {
            out_temp[i][0] *= 2;
            out_temp[i][1] *= 2;
        }

        if (INN % 2 == 0)
        {
            numRem--;
        }
        else if (INN > 1)
        {
            out_temp[HN][0] *= 2;
            out_temp[HN][1] *= 2;
        }

        memset(&out_temp[HN + 1][0], 0, numRem * sizeof(fftw_complex));

        // std::memset(in_temp, 0, sizeof(fftw_complex) * INN);

        fftw_plan pf2 = fftw_plan_dft_1d(INN, out_temp, in_temp, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(pf2);
        fftw_destroy_plan(pf2);
        // std::cout << "INN 2 =" << INN << "\n";

        std::vector<std::complex<T>> out(INN);
        for (size_t i = 0; i < INN; ++i)
        {
            out[i].real(in_temp[i][0] / INN);
            out[i].imag(in_temp[i][1] / INN);
        }
        fftw_free(out_temp);
        fftw_free(in_temp);
        fftw_cleanup();

        return out;
    }

    //
    // phiMat = angle(hilbert(gather));
    // phaseshiftMat = exp(sqrt(-1)*phiMat);
    template <class T>
    inline std::vector<std::complex<T>> instanPhaseEstimator(std::vector<T> &v)
    {
        std::vector<std::complex<T>> ov;
        // micro_time = AU_WTIME;
        ov = Hilbert(v);
        // sum_micro = sum_micro + (AU_WTIME - micro_time);

        size_t N = ov.size();
        T ov_angle;
        const std::complex<T> minus_one(0, -1);
        for (int i = 0; i < N; i++)
        {
            ov_angle = std::arg(ov[i]);
            ov[i] = std::exp(minus_one * ov_angle);
        }
        return ov;
    }

    template <class T>
    inline std::vector<std::vector<std::complex<T>>> HilbertVector(std::vector<std::vector<T>> &in)
    {
        size_t INN_ROW = in.size();
        size_t INN_COL = in[0].size();
        size_t INN_TOTOAL = INN_ROW * INN_COL;
        int rank, col, howmany, istride, idist, ostride, odist;

        // std::cout << "INN =" << INN << "\n";
        fftw_complex *in_temp = (fftw_complex *)fftw_alloc_complex(sizeof(fftw_complex) * INN_TOTOAL);
        // std::memset(in_temp, 0, sizeof(fftw_complex) * INN);

        size_t temp_index = 0;
        for (size_t i = 0; i < INN_ROW; i++)
        {
            for (size_t j = 0; j < INN_COL; j++)
            {
                in_temp[temp_index][0] = in[i][j];
                in_temp[temp_index][1] = 0;
                temp_index++;
            }
        }

        fftw_complex *out_temp = (fftw_complex *)fftw_alloc_complex(sizeof(fftw_complex) * INN_TOTOAL);
        // std::memset(out_temp, 0, sizeof(fftw_complex) * INN);

        // micro_time = AU_WTIME;
        //     fftw_plan pf = fftw_plan_dft_1d(INN, in_temp, out_temp, FFTW_FORWARD, FFTW_ESTIMATE);
        rank = 1;
        howmany = INN_ROW;
        istride = 1;
        idist = INN_COL;
        ostride = 1;
        odist = INN_COL;
        col = INN_COL;

        fftw_plan pf = fftw_plan_many_dft(rank, &col, howmany,
                                          in_temp, NULL,
                                          istride, idist,
                                          out_temp, NULL,
                                          ostride, odist,
                                          FFTW_FORWARD, FFTW_ESTIMATE);
        // micro_time_sub = AU_WTIME;
        fftw_execute(pf);
        // sum_micro_sub = sum_micro_sub + (AU_WTIME - micro_time_sub);

        fftw_destroy_plan(pf);
        fftw_cleanup();
        // sum_micro = sum_micro + (AU_WTIME - micro_time);

        size_t HN = INN_COL >> 1;
        size_t numRem = HN;
        size_t start_offset;
        for (size_t row = 0; row < INN_ROW; row++)
        {
            start_offset = row * INN_COL;
            for (size_t i = 1; i < HN; ++i)
            {
                out_temp[start_offset + i][0] *= 2;
                out_temp[start_offset + i][1] *= 2;
            }

            if (INN_COL % 2 == 0)
            {
                numRem--;
            }
            else if (INN_COL > 1)
            {
                out_temp[start_offset + HN][0] *= 2;
                out_temp[start_offset + HN][1] *= 2;
            }

            memset(&out_temp[start_offset + HN + 1][0], 0, numRem * sizeof(fftw_complex));
        }

        // std::memset(in_temp, 0, sizeof(fftw_complex) * INN);
        // fftw_plan pf2 = fftw_plan_dft_1d(INN, out_temp, in_temp, FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_plan pf2 = fftw_plan_many_dft(rank, &col, howmany,
                                           out_temp, NULL,
                                           istride, idist,
                                           in_temp, NULL,
                                           ostride, odist,
                                           FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(pf2);
        fftw_destroy_plan(pf2);
        // std::cout << "INN 2 =" << INN << "\n";

        std::vector<std::vector<std::complex<T>>> out;

        fftw_free(out_temp);

        out.resize(INN_ROW);
        for (size_t i = 0; i < INN_ROW; ++i)
        {
            out[i].resize(INN_COL);
            for (size_t j = 0; j < INN_COL; ++j)
            {
                out[i][j].real(in_temp[INN_COL * i + j][0] / INN_COL);
                out[i][j].imag(in_temp[INN_COL * i + j][1] / INN_COL);
            }
        }

        fftw_free(in_temp);
        fftw_cleanup();

        return out;
    }

    // In-place version
    template <class T>
    inline std::vector<std::vector<std::complex<T>>> HilbertVectorIP(const std::vector<std::vector<T>> &in)
    {
        size_t INN_ROW = in.size();
        size_t INN_COL = in[0].size();
        size_t INN_TOTOAL = INN_ROW * INN_COL;
        int rank, col, howmany, istride, idist, ostride, odist;

        // std::cout << "INN =" << INN << "\n";
        // fftw_complex *in_temp = (fftw_complex *)fftw_alloc_complex(sizeof(fftw_complex) * INN_TOTOAL);
        // std::memset(in_temp, 0, sizeof(fftw_complex) * INN);
        fftw_complex *out_temp = (fftw_complex *)fftw_alloc_complex(sizeof(fftw_complex) * INN_TOTOAL);
        if (out_temp == NULL)
        {
            AU_EXIT("not enough memory !");
        }

        size_t temp_index = 0;
        for (size_t i = 0; i < INN_ROW; i++)
        {
            for (size_t j = 0; j < INN_COL; j++)
            {
                out_temp[temp_index][0] = in[i][j];
                out_temp[temp_index][1] = 0;
                temp_index++;
            }
        }

        // std::memset(out_temp, 0, sizeof(fftw_complex) * INN);

        // micro_time = AU_WTIME;
        //     fftw_plan pf = fftw_plan_dft_1d(INN, in_temp, out_temp, FFTW_FORWARD, FFTW_ESTIMATE);
        rank = 1;
        howmany = INN_ROW;
        istride = 1;
        idist = INN_COL;
        ostride = 1;
        odist = INN_COL;
        col = INN_COL;

        fftw_plan pf = fftw_plan_many_dft(rank, &col, howmany,
                                          out_temp, NULL,
                                          istride, idist,
                                          out_temp, NULL,
                                          ostride, odist,
                                          FFTW_FORWARD, FFTW_ESTIMATE);
        // micro_time_sub = AU_WTIME;
        fftw_execute(pf);
        // sum_micro_sub = sum_micro_sub + (AU_WTIME - micro_time_sub);

        fftw_destroy_plan(pf);
        fftw_cleanup();
        // sum_micro = sum_micro + (AU_WTIME - micro_time);

        size_t HN = INN_COL >> 1;
        size_t numRem = HN;
        size_t start_offset;
        for (size_t row = 0; row < INN_ROW; row++)
        {
            start_offset = row * INN_COL;
            for (size_t i = 1; i < HN; ++i)
            {
                out_temp[start_offset + i][0] *= 2;
                out_temp[start_offset + i][1] *= 2;
            }

            if (INN_COL % 2 == 0)
            {
                numRem--;
            }
            else if (INN_COL > 1)
            {
                out_temp[start_offset + HN][0] *= 2;
                out_temp[start_offset + HN][1] *= 2;
            }

            memset(&out_temp[start_offset + HN + 1][0], 0, numRem * sizeof(fftw_complex));
        }

        // std::memset(in_temp, 0, sizeof(fftw_complex) * INN);
        // fftw_plan pf2 = fftw_plan_dft_1d(INN, out_temp, in_temp, FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_plan pf2 = fftw_plan_many_dft(rank, &col, howmany,
                                           out_temp, NULL,
                                           istride, idist,
                                           out_temp, NULL,
                                           ostride, odist,
                                           FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(pf2);
        fftw_destroy_plan(pf2);
        // std::cout << "INN 2 =" << INN << "\n";

        std::vector<std::vector<std::complex<T>>> out;

        out.resize(INN_ROW);
        for (size_t i = 0; i < INN_ROW; ++i)
        {
            out[i].resize(INN_COL);
            for (size_t j = 0; j < INN_COL; ++j)
            {
                out[i][j].real(out_temp[INN_COL * i + j][0] / INN_COL);
                out[i][j].imag(out_temp[INN_COL * i + j][1] / INN_COL);
            }
        }

        fftw_free(out_temp);

        // fftw_free(in_temp);
        fftw_cleanup();

        return out;
    }

    template <class T>
    inline std::vector<std::vector<std::complex<T>>> instanPhaseEstimatorVector(std::vector<std::vector<T>> &v)
    {
        std::vector<std::vector<std::complex<T>>> ov;
        // micro_time = AU_WTIME;
        // ov = HilbertVector(v);

        ov = HilbertVectorIP(v);
        // sum_micro = sum_micro + (AU_WTIME - micro_time);

        size_t N_ROW = ov.size(), N_COL = ov[0].size();
        T ov_angle;
        const std::complex<T> minus_one(0, -1);
        for (size_t i = 0; i < N_ROW; i++)
        {
            for (size_t j = 0; j < N_COL; j++)
            {
                ov_angle = std::arg(ov[i][j]);
                ov[i][j] = std::exp(minus_one * ov_angle);
            }
        }
        return ov;
    }

    template <class T>
    inline void clear_vector(std::vector<T> &v)
    {
        v.clear();
        std::vector<T>().swap(v);
    }

    void transpose_data(int16_t *src, int16_t *dst, const int N, const int M)
    {
        for (int n = 0; n < N * M; n++)
        {
            int i = n / N;
            int j = n % N;
            dst[n] = src[M * j + i];
        }
    }

    /**
     * @brief Convert a N x M matrix src to M x N matrix dst
     *
     * @tparam T
     * @param src
     * @param dst
     * @param N : row
     * @param M : column
     */
    template <class T>
    inline void transpose(T *src, T *dst, const int N, const int M)
    {
        int i, j;
        for (int n = 0; n < N * M; n++)
        {
            i = n / N;
            j = n % N;
            dst[n] = src[M * j + i];
        }
    }

} // namespace DasLib

#endif