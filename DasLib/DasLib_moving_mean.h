#ifndef DASLIB_MOVEMEAN_H
#define DASLIB_MOVEMEAN_H 1

#include <vector>

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

// window size = 2* HAL_WIN + 1
// XV.size = YV.size()
template <typename T>
inline void moving_mean(const std::vector<T> &XV, std::vector<T> &YV, const int HAL_WIN)
{
    unsigned long long XN = XV.size(), index1, index2;
    YV.resize(XN);
    for (unsigned long long xi = 0; xi < XN; xi++)
    {
        if (xi < HAL_WIN)
        {
            index1 = 0;
            index2 = xi;
        }
        else if (xi >= XN - HAL_WIN - 1)
        {
            index1 = xi;
            index2 = XN - 1;
        }
        else
        {
            index1 = xi - HAL_WIN;
            index2 = xi + HAL_WIN;
        }

        double norm_denom = 0;
        for (int j = index1; j <= index2; j++)
        {
            norm_denom = norm_denom + std::abs(XV[j]);
        }
        YV[xi] = XV[xi] / (norm_denom / (index2 - index1 + 1));
    }
}
#endif