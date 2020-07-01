/*
 *  interp1
 */
template <typename Real>
int nearestNeighbourIndex(std::vector<Real> &x, Real &value)
{
    Real dist = std::numeric_limits<Real>::max();
    Real newDist = dist;
    size_t idx = 0;

    for (size_t i = 0; i < x.size(); ++i)
    {
        newDist = std::abs(value - x[i]);
        if (newDist <= dist)
        {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

template <typename Real>
void interp1(std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &x_new, std::vector<Real> &y_new)
{
    Real dx, dy, m, b;
    size_t x_max_idx = x.size() - 1;
    size_t x_new_size = x_new.size();

    for (size_t i = 0; i < x_new_size; ++i)
    {
        size_t idx = nearestNeighbourIndex(x, x_new[i]);

        if (x[idx] > x_new[i])
        {
            dx = idx > 0 ? (x[idx] - x[idx - 1]) : (x[idx + 1] - x[idx]);
            dy = idx > 0 ? (y[idx] - y[idx - 1]) : (y[idx + 1] - y[idx]);
        }
        else
        {
            dx = idx < x_max_idx ? (x[idx + 1] - x[idx]) : (x[idx] - x[idx - 1]);
            dy = idx < x_max_idx ? (y[idx + 1] - y[idx]) : (y[idx] - y[idx - 1]);
        }

        m = dy / dx;
        b = y[idx] - x[idx] * m;

        y_new[i] = (x_new[i] * m + b);
    }
}