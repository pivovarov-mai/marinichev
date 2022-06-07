#ifndef RUNGE_KUTTA_METHOD_HPP
#define RUNGE_KUTTA_METHOD_HPP

#include <functional>
#include <vector>

#include "error.hpp"

class RungeKuttaMethod
{
private:
    double l, r;
    std::function<double(double, double, double)> f, g;
    double y0, z0;

public:
    RungeKuttaMethod(const double l, const double r,
                     const std::function<double(double, double, double)> f, const std::function<double(double, double, double)> g,
                     const double y0, const double z0) : l(l), r(r), f(f), g(g), y0(y0), z0(z0) {}

    std::vector<std::tuple<double, double, double>> solve(double h)
    {
        std::vector<std::tuple<double, double, double>> table;
        double xk = l;
        double yk = y0;
        double zk = z0;
        table.push_back(std::make_tuple(xk, yk, zk));
        while (leq(xk + h, r))
        {
            double K1 = h * f(xk, yk, zk);
            double L1 = h * g(xk, yk, zk);
            double K2 = h * f(xk + 0.5 * h, yk + 0.5 * K1, zk + 0.5 * L1);
            double L2 = h * g(xk + 0.5 * h, yk + 0.5 * K1, zk + 0.5 * L1);
            double K3 = h * f(xk + 0.5 * h, yk + 0.5 * K2, zk + 0.5 * L2);
            double L3 = h * g(xk + 0.5 * h, yk + 0.5 * K2, zk + 0.5 * L2);
            double K4 = h * f(xk + h, yk + K3, zk + L3);
            double L4 = h * g(xk + h, yk + K3, zk + L3);
            double dy = (K1 + 2.0 * K2 + 2.0 * K3 + K4) / 6.0;
            double dz = (L1 + 2.0 * L2 + 2.0 * L3 + L4) / 6.0;
            xk += h;
            yk += dy;
            zk += dz;
            table.push_back(std::make_tuple(xk, yk, zk));
        }
        return table;
    }
};

#endif /* RUNGE_KUTTA_METHOD_HPP */