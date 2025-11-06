#include "erosionfield.h"
#include <algorithm>
#include <cmath>

template <typename T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

inline double ErosionField::smin(double a, double b, double k)
{
    double h = clamp(0.5 + 0.5*(b - a)/k, 0.0, 1.0);
    return (b*(1.0 - h) + a*h) - k*h*(1.0 - h);
}
inline double ErosionField::smax(double a, double b, double k)
{
    return -smin(-a, -b, k);
}

double ErosionField::Value(const Vector& p) const
{
    double v = (base ? base->Value(p) : (Norm(p) - 1.0));

    for(const auto& s : impacts)
    {
        double sd = Norm(p - s.c) - s.R; 
        double cut = -sd; 
        if(s.smoothK > 0.0)
            v = smax(v, cut, s.smoothK);
        else
            v = std::max(v, cut);
    }
    return v;
}