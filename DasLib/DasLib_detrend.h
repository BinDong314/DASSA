template <typename T>
void detrend(T *y, int m)
{
    T xmean, ymean;
    int i;
    T Sxy;
    T Sxx;

    T grad;
    T yint;

    std::unique_ptr<T[]> x(new T[m]);

    /********************************
    Set the X axis Liner Values
    *********************************/
    for (i = 0; i < m; i++)
        x[i] = i;

    /********************************
    Calculate the mean of x and y
    *********************************/
    xmean = 0;
    ymean = 0;
    for (i = 0; i < m; i++)
    {
        xmean += x[i];
        ymean += y[i];
    }
    xmean /= m;
    ymean /= m;

    /********************************
    Calculate Covariance
    *********************************/

    Sxy = 0;
    for (i = 0; i < m; i++)
        Sxy += (x[i] - xmean) * (y[i] - ymean);

    Sxx = 0;
    for (i = 0; i < m; i++)
        Sxx += (x[i] - xmean) * (x[i] - xmean);

    /********************************
    Calculate Gradient and Y intercept
    *********************************/
    grad = Sxy / Sxx;
    yint = -grad * xmean + ymean;

    /********************************
    Removing Linear Trend
    *********************************/
    for (i = 0; i < m; i++)
    {
        y[i] = y[i] - (grad * i + yint);
    }
}

template <typename T>
inline void detrend_range(const std::vector<T> &y, const size_t start, const size_t count, std::vector<T> &out_vector)
{
    // size_t m = count;
    out_vector.resize(count);
    T xmean, ymean;
    int i;
    T Sxy;
    T Sxx;

    T grad;
    T yint;

    // std::unique_ptr<T[]> x(new T[m]);

    /********************************
    Set the X axis Liner Values
    *********************************/
    // for (i = 0; i < m; i++)
    //     out_vector[i] = i;

    /********************************
    Calculate the mean of x and y
    *********************************/
    xmean = 0;
    ymean = 0;
    for (i = 0; i < count; i++)
    {
        // out_vector[i] = i;
        xmean += i;
        ymean += y[start + i];
    }
    xmean /= count;
    ymean /= count;

    /********************************
    Calculate Covariance
    *********************************/

    Sxy = 0;
    Sxx = 0;
    for (i = 0; i < count; i++)
    {
        Sxy += (i - xmean) * (y[start + i] - ymean);
        Sxx += (i - xmean) * (i - xmean);
    }
    // for (i = 0; i < m; i++)

    /********************************
    Calculate Gradient and Y intercept
    *********************************/
    grad = Sxy / Sxx;
    yint = -grad * xmean + ymean;

    /********************************
    Removing Linear Trend
    *********************************/
    for (i = 0; i < count; i++)
    {
        out_vector[i] = y[start + i] - (grad * i + yint);
    }
}