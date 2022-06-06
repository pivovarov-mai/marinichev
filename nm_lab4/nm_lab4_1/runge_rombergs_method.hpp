#ifndef RUNGE_ROMBERGS_METHOD_HPP
#define RUNGE_ROMBERGS_METHOD_HPP

#include <vector>
#include <cmath>
#include <tuple>

double runge_rombergs_method(const std::vector<std::tuple<double, double, double>> &y_2h, const std::vector<std::tuple<double, double, double>> &y_h, double p)
{
    double coef = 1.0 / (std::pow(2, p) - 1.0);
    double error = 0.0;
    for (size_t i = 0; i < y_2h.size(); ++i)
        error = std::max(error, coef * std::abs(std::get<1>(y_2h[i]) - std::get<1>(y_h[2 * i])));

    return error;
}

#endif /* RUNGE_ROMBERGS_METHOD_HPP */