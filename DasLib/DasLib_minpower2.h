#ifndef DASLIB_MINPOWER2_H
#define DASLIB_MINPOWER2_H
//It is power of 2 and greater than minimum_m
#define MINPOWER2(MIN_M, M)                    \
    {                                          \
        unsigned long long mt = 2 * MIN_M - 1; \
        while ((mt & (mt - 1)) != 0)           \
        {                                      \
            mt = mt + 1;                       \
        }                                      \
        M = mt;                                \
    }

template <typename T1, typename T2>
inline void minpower2(const T1 MIN_M, T2 &M)
{
    T2 mt = 2 * MIN_M - 1;
    while ((mt & (mt - 1)) != 0)
    {
        mt = mt + 1;
    }
    M = mt;
}

#endif