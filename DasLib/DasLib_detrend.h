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