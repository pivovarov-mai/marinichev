#ifndef ADAMS_METHOD_HPP
#define ADAMS_METHOD_HPP

#include <functional>
#include <vector>

#include "runge_kutta_method.hpp"

class AdamsMethod
{
private:
    double l, r;
    std::function<double(double, double, double)> f, g;
    double y0, z0;

public:
    AdamsMethod(const double l, const double r,
                const std::function<double(double, double, double)> f, const std::function<double(double, double, double)> g,
                const double y0, const double z0) : l(l), r(r), f(f), g(g), y0(y0), z0(z0) {}

    double calc_tuple(std::function<double(double, double, double)> f, std::tuple<double, double, double> xyz) { return f(std::get<0>(xyz), std::get<1>(xyz), std::get<2>(xyz)); }

    std::vector<std::tuple<double, double, double>> solve(double h)
    {
        if (l + 3.0 * h > r)
            throw std::invalid_argument("Error: step h is too big");

        RungeKuttaMethod first_points(l, l + 3.0 * h, f, g, y0, z0);
        std::vector<std::tuple<double, double, double>> table = first_points.solve(h);
        size_t cnt = table.size();
        double xk = std::get<0>(table.back());
        double yk = std::get<1>(table.back());
        double zk = std::get<2>(table.back());
        while (leq(xk + h, r))
        {
            // predictor
            double dy = (h / 24.0) * (55.0 * calc_tuple(f, table[cnt - 1]) - 59.0 * calc_tuple(f, table[cnt - 2]) + 37.0 * calc_tuple(f, table[cnt - 3]) - 9.0 * calc_tuple(f, table[cnt - 4]));
            double dz = (h / 24.0) * (55.0 * calc_tuple(g, table[cnt - 1]) - 59.0 * calc_tuple(g, table[cnt - 2]) + 37.0 * calc_tuple(g, table[cnt - 3]) - 9.0 * calc_tuple(g, table[cnt - 4]));
            double xk1 = xk + h;
            double yk1 = yk + dy;
            double zk1 = zk + dz;
            table.push_back(std::make_tuple(xk1, yk1, zk1));
            ++cnt;

            // corrector
            dy = (h / 24.0) * (9.0 * calc_tuple(f, table[cnt - 1]) + 19.0 * calc_tuple(f, table[cnt - 2]) - 5.0 * calc_tuple(f, table[cnt - 3]) + 1.0 * calc_tuple(f, table[cnt - 4]));
            dz = (h / 24.0) * (9.0 * calc_tuple(g, table[cnt - 1]) + 19.0 * calc_tuple(g, table[cnt - 2]) - 5.0 * calc_tuple(g, table[cnt - 3]) + 1.0 * calc_tuple(g, table[cnt - 4]));
            xk += h;
            yk += dy;
            zk += dz;
            table.pop_back();
            table.push_back(std::make_tuple(xk, yk, zk));
        }
        return table;
    }
};

#endif /* ADAMS_METHOD_HPP */