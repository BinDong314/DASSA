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
inline void detrend_range_one_pass_std(const std::vector<T> &y, const size_t start, const size_t count, const std::vector<T> &ctap, const double xmean, const double xsum, const double Sxx, std::vector<T> &out_vector)
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
    double sum_xy = 0;
    double ysum = 0;
    // std::cout << "detrend_range_one_pass_std, xmean = " << xmean << ", xsum =" << xsum << ", Sxx = " << Sxx << "\n";
    // std::cout << "detrend_range_one_pass_std, start = " << start << ", count =" << count << "\n";
    // std::cout << "y = ";
    for (i = 0; i < count; i++)
    {
        // out_vector[i] = i;
        // xmean += i;
        // ymean += y[start + i];
        ysum += y[start + i];
        sum_xy = sum_xy + i * y[start + i];
        // if (i < 10)
        //     std::cout << ", " << y[start + i];
    }
    // std::cout << " \n ";
    ymean = ysum / count;

    // std::cout << "ymean = " << ymean << "\n";
    //  Sxy = 0;
    //   for (i = 0; i < count; i++)
    //   {
    //       Sxy += (i - xmean) * (y[start + i] - ymean);
    //   }

    Sxy = sum_xy - ymean * xsum - xmean * ysum + count * xmean * ymean;
    // std::cout << "Sxy = " << Sxy << "\n";

    /********************************
    Calculate Gradient and Y intercept
    *********************************/
    grad = Sxy / Sxx;
    yint = -grad * xmean + ymean;

    // std::cout << "grad = " << grad << "\n";
    // std::cout << "yint = " << yint << "\n";

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
    // std::cout << "v_sum = " << v_sum << "\n";
    // std::cout << "v_sum_sqrt = " << v_sum_sqrt << "\n";
    // std::cout << "out_vector = ";
    for (i = 0; i < count; i++)
    {
        out_vector[i] = out_vector[i] / v_sum_sqrt;
        //   if (i < 10)
        //       std::cout << ", " << out_vector[i];
    }
    // std::cout << " \n ";

    // exit(0);
    //  std::cout << "yint = " << yint << " , grad = " << grad << ", Sxy = " << Sxy << ", Sxx = " << Sxx << ", xmean = " << xmean << ", ymean = " << ymean << " \n";
    //  PrintVector("After detrend_range out_vector =", out_vector);
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