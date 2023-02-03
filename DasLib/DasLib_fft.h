
#ifndef DASLIB_FFT_H
#define DASLIB_FFT_H

#include <iostream> // std::cout
#include <complex>  // std::complex, std::imag
#include <fftw3.h>
#include "DasLib_minpower2.h"

/**
 * NN : the length of the fft_in_p
 * fft_in_p: the input vector with the type std::complex<double>
 * fft_out_p: the output vector with the type std::complex<double>
 * direction_p : FFTW_FORWARD or FFTW_BACKWARD
 **/
#define RUN_FFTV(NN, fft_in_p, fft_out_p, direction_p)                                                                                                             \
    {                                                                                                                                                              \
        fftw_plan fft_p;                                                                                                                                           \
        fft_p = fftw_plan_dft_1d(NN, reinterpret_cast<fftw_complex *>(&fft_in_p[0]), reinterpret_cast<fftw_complex *>(&fft_out_p[0]), direction_p, FFTW_ESTIMATE); \
        fftw_execute(fft_p);                                                                                                                                       \
        fftw_destroy_plan(fft_p);                                                                                                                                  \
    }

/**
 * @brief it performs the fft on the input vector fft_vector and it also extend the input size to be minimum of power of 2
 *        https://stackoverflow.com/questions/4214400/problem-casting-stl-complexdouble-to-fftw-complex
 * @param fft_vector
 * @param fft_out
 */
inline void fftv_forward_p2(const std::vector<double> &fft_in_v, std::vector<std::complex<double>> &fft_out)
{
    size_t nfft;
    minpower2(fft_in_v.size(), nfft);

    std::vector<std::complex<double>> fft_in;
    fft_in.resize(nfft, std::complex<double>(0, 0));
    fft_out.resize(nfft, std::complex<double>(0, 0));

    for (size_t i = 0; i < fft_in_v.size(); i++)
    {
        fft_in[i].real(fft_in_v[i]);
    }

    RUN_FFTV(nfft, fft_in, fft_out, FFTW_FORWARD);
}

/**
 * @brief it performs the fft on the input vector fft_vector and it also extend the input size to be minimum of power of 2
 *        https://stackoverflow.com/questions/4214400/problem-casting-stl-complexdouble-to-fftw-complex
 * @param fft_vector
 * @param fft_out
 */
inline void fftv_forward(const std::vector<double> &fft_in_v, std::vector<std::complex<double>> &fft_out)
{
    size_t nfft = fft_in_v.size();

    std::vector<std::complex<double>> fft_in;
    fft_in.resize(nfft, std::complex<double>(0, 0));
    fft_out.resize(nfft, std::complex<double>(0, 0));

    for (size_t i = 0; i < fft_in_v.size(); i++)
    {
        fft_in[i].real(fft_in_v[i]);
    }

    RUN_FFTV(nfft, fft_in, fft_out, FFTW_FORWARD);
}

/**
 * @brief it has backward on vector fft_in
 *
 * @param fft_in : input vector with the std::complex<double>
 * @param fft_out : output vector with the type std::complex<double>
 */
inline void fftv_backward(std::vector<std::complex<double>> &fft_in, std::vector<std::complex<double>> &fft_out)
{
    size_t nfft;
    nfft = fft_in.size();
    fft_out.resize(nfft, std::complex<double>(0, 0));
    RUN_FFTV(nfft, fft_in, fft_out, FFTW_BACKWARD);
    for (size_t i = 0; i < nfft; i++)
    {
        fft_out[i].real(fft_out[i].real() / nfft);
    }
}

/**
 * @brief it has backward on vector fft_in
 *
 * @param fft_in : input vector with the std::complex<double>
 * @param fft_out : output vector with the type std::complex<double>
 */
inline void fftv_backward_real(std::vector<std::complex<double>> &fft_in, std::vector<double> &fft_out_real)
{
    size_t nfft;
    nfft = fft_in.size();
    fft_out_real.resize(nfft, 0.0);
    RUN_FFTV(nfft, fft_in, fft_out_real, FFTW_BACKWARD);
    for (size_t i = 0; i < nfft; i++)
    {
        // fft_out[i].real(fft_out[i].real() / nfft);
        fft_out_real[i] = fft_out_real[i].real() / nfft;
    }
}

/**
 * @brief it has backward on vector fft_in
 *
 * @param fft_in : input vector with the std::complex<double>
 * @param fft_out : output vector with the type std::complex<double>
 */
inline void fftv_backward_real_max(std::vector<std::complex<double>> &fft_in, double fft_out_real_max)
{
    size_t nfft;
    nfft = fft_in.size();

    std::vector<std::complex<double>> &fft_out;
    fft_out.resize(nfft, std::complex<double>(0, 0));

    RUN_FFTV(nfft, fft_in, fft_out, FFTW_BACKWARD);

    fft_out_real_max = fft_out[0].real() / nfft;

    for (size_t i = 1; i < nfft; i++)
    {
        if (fft_out_real_max < (fft_out[i].real() / nfft))
            fft_out_real_max = fft_out[i].real() / nfft;
    }
}

#endif