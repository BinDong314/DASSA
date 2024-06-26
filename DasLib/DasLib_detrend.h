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
inline void detrend_range(const std::vector<T> &y, const size_t start, const size_t count, const std::vector<T> &ctap, const double xmean, const double Sxx, std::vector<T> &out_vector)
{

    // size_t m = count;
    out_vector.resize(count);
    // T xmean, ymean;
    T ymean;
    size_t i;
    T Sxy;
    // T Sxx;

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
    // xmean = 0;
    ymean = 0;
    for (i = 0; i < count; i++)
    {
        // out_vector[i] = i;
        // xmean += i;
        ymean += y[start + i];
    }
    // xmean / = count;
    // xmean = (1 + count) / 2;
    ymean /= count;

    /********************************
    Calculate Covariance
    *********************************/

    Sxy = 0;
    // Sxx = 0;
    for (i = 0; i < count; i++)
    {
        Sxy += (i - xmean) * (y[start + i] - ymean);
        // Sxx += (i - xmean) * (i - xmean);
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
    double v_sum = 0;
    for (i = 0; i < count; i++)
    {
        out_vector[i] = y[start + i] - (grad * i + yint);
        out_vector[i] = out_vector[i] * ctap[i];
        v_sum = out_vector[i] * out_vector[i];
    }

    double v_sum_sqrt = sqrt(v_sum);
    for (i = 0; i < count; i++)
    {
        out_vector[i] = out_vector[i] / v_sum_sqrt;
    }

    // std::cout << "yint = " << yint << " , grad = " << grad << ", Sxy = " << Sxy << ", Sxx = " << Sxx << ", xmean = " << xmean << ", ymean = " << ymean << " \n";
    // PrintVector("After detrend_range out_vector =", out_vector);
}

template <typename T>
inline double subset_sqsum(const std::vector<T> &y, const size_t start, const size_t count, std::vector<T> &out_vector)
{

    out_vector.resize(count);
    double sqsum = 0;

    size_t i;
    for (i = 0; i < count; i++)
    {
        out_vector[i] = y[start + i];
        sqsum = sqsum + out_vector[i] * out_vector[i];
    }
    return sqsum;
}

template <typename T>
inline void sum_sqrt(const std::vector<T> &x, std::vector<T> &sumsq, const int window_length)
{
    size_t n = x.size();
    sumsq.resize(n, 0);
    T temp_sum = 0;

    for (size_t i = 0; i < window_length; i++)
    {
        temp_sum = temp_sum + x[i] * x[i];
    }
    sumsq[0] = sqrt(temp_sum);
    for (size_t i = 1; i < (n - window_length + 1); i++)
    {
        temp_sum = temp_sum - x[i - 1] * x[i - 1];
        temp_sum = temp_sum + x[i + window_length - 1] * x[i + window_length - 1];
        sumsq[i] = sqrt(temp_sum);
    }
}

template <typename T>
inline void detrend_range_one_pass_std(const std::vector<T> &y, const size_t start, const size_t count, const std::vector<T> &ctap, const double xmean, const double xsum, const double Sxx, std::vector<T> &out_vector)
{

    out_vector.resize(count);
    T ymean;
    size_t i;
    T Sxy;

    T grad;
    T yint;

    /********************************
    Calculate the mean of x and y
    *********************************/
    ymean = 0;
    double sum_xy = 0;
    double ysum = 0;

    for (i = 0; i < count; i++)
    {
        ysum += y[start + i];
        sum_xy += i * y[start + i];
    }

    ymean = ysum / count;
    Sxy = sum_xy - ymean * xsum - xmean * ysum + count * xmean * ymean;

    /********************************
    Calculate Gradient and Y intercept
    *********************************/
    grad = Sxy / Sxx;
    yint = -grad * xmean + ymean;

    /********************************
    Removing Linear Trend
    *********************************/
    double v_sum = 0;
    for (i = 0; i < count; i++)
    {
        out_vector[i] = y[start + i] - (grad * i + yint);
        out_vector[i] = out_vector[i] * ctap[i];
        v_sum += out_vector[i] * out_vector[i];
    }

    double v_sum_sqrt = sqrt(v_sum);

    for (i = 0; i < count; i++)
    {
        out_vector[i] = out_vector[i] / v_sum_sqrt;
    }
}

template <typename T>
inline void detrend_range_dqueue(const std::deque<T> &y, const size_t count, const std::vector<T> &ctap, const double xmean, const double Sxx, std::vector<T> &out_vector)
{

    // size_t m = count;
    out_vector.resize(count);
    // T xmean, ymean;
    T ymean;
    size_t i;
    T Sxy;
    // T Sxx;

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
    // xmean = 0;
    ymean = 0;
    for (i = 0; i < count; i++)
    {
        // out_vector[i] = i;
        // xmean += i;
        ymean += y[i];
    }
    // xmean / = count;
    // xmean = (1 + count) / 2;
    ymean /= count;

    /********************************
    Calculate Covariance
    *********************************/

    Sxy = 0;
    // Sxx = 0;
    for (i = 0; i < count; i++)
    {
        Sxy += (i - xmean) * (y[i] - ymean);
        // Sxx += (i - xmean) * (i - xmean);
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
    double v_sum = 0;
    for (i = 0; i < count; i++)
    {
        out_vector[i] = y[i] - (grad * i + yint);
        out_vector[i] = out_vector[i] * ctap[i];
        v_sum = out_vector[i] * out_vector[i];
    }

    double v_sum_sqrt = sqrt(v_sum);
    for (i = 0; i < count; i++)
    {
        out_vector[i] = out_vector[i] / v_sum_sqrt;
    }

    // std::cout << "yint = " << yint << " , grad = " << grad << ", Sxy = " << Sxy << ", Sxx = " << Sxx << ", xmean = " << xmean << ", ymean = " << ymean << " \n";
    // PrintVector("After detrend_range out_vector =", out_vector);
}

template <typename T>
inline void detrend_range(const std::vector<T> &y, const size_t start, const size_t count, const std::vector<T> &ctap, std::vector<T> &out_vector)
{
    // size_t m = count;
    out_vector.resize(count);
    T xmean, ymean;
    size_t i;
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
    double v_sum = 0;
    for (i = 0; i < count; i++)
    {
        out_vector[i] = y[start + i] - (grad * i + yint);
        out_vector[i] = out_vector[i] * ctap[i];
        v_sum = v_sum + out_vector[i] * out_vector[i];
    }

    double v_sum_sqrt = sqrt(v_sum);
    for (i = 0; i < count; i++)
    {
        out_vector[i] = out_vector[i] / v_sum_sqrt;
    }

    // std::cout << "yint = " << yint << " , grad = " << grad << ", Sxy = " << Sxy << ", Sxx = " << Sxx << ", xmean = " << xmean << ", ymean = " << ymean << " \n";
    // PrintVector("After detrend_range out_vector =", out_vector);
}