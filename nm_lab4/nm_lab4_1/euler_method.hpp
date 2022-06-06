#ifndef EULER_METHOD_HPP
#define EULER_METHOD_HPP

#include <functional>
#include <vector>

#include "error.hpp"

class EulerMethod
{
private:
    double l, r;
    std::function<double(double, double, double)> f, g;
    double y0, z0;

public:
    EulerMethod(const double l, const double r,
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
            double dy = h * f(xk, yk, zk);
            double dz = h * g(xk, yk, zk);
            xk += h;
            yk += dy;
            zk += dz;
            table.push_back(std::make_tuple(xk, yk, zk));
        }
        return table;
    }
};

#endif /* EULER_METHOD_HPP */