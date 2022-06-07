#ifndef SHOOTING_METHOD_HPP
#define SHOOTING_METHOD_HPP

#include "../nm_lab4_1/euler_method.hpp"

class ShootingMethod
{
private:
    double a, b;
    std::function<double(double, double, double)> f, g;
    double alpha, beta, y0;
    double delta, gamma, y1;

public:
    ShootingMethod(const double a, const double b,
                   const std::function<double(double, double, double)> f, const std::function<double(double, double, double)> g,
                   const double alpha, const double beta, const double y0,
                   const double delta, const double gamma, const double y1)
        : a(a), b(b), f(f), g(g),
          alpha(alpha), beta(beta), y0(y0),
          delta(delta), gamma(gamma), y1(y1) {}

    double get_start_cond(double eta) { return (y0 - alpha * eta) / beta; }

    double get_eta_next(double eta_prev, double eta, const std::vector<std::tuple<double, double, double>> prev_solution, const std::vector<std::tuple<double, double, double>> solution)
    {
        double yb_prev = std::get<1>(prev_solution.back());
        double zb_prev = std::get<2>(prev_solution.back());
        double phi_prev = delta * yb_prev + gamma * zb_prev - y1;
        double yb = std::get<1>(solution.back());
        double zb = std::get<2>(solution.back());
        double phi = delta * yb + gamma * zb - y1;
        return eta - (eta - eta_prev) / (phi - phi_prev) * phi;
    }

    std::vector<std::tuple<double, double, double>> solve(double h, double eps)
    {
        double eta_prev = 1.0;
        double eta = 0.8;
        while (true)
        {
            double runge_z0_prev = get_start_cond(eta_prev);
            EulerMethod euler_method1(a, b, f, g, eta_prev, runge_z0_prev);
            std::vector<std::tuple<double, double, double>> prev_solution = euler_method1.solve(h);

            double runge_z0 = get_start_cond(eta);
            EulerMethod euler_method2(a, b, f, g, eta, runge_z0);
            std::vector<std::tuple<double, double, double>> solution = euler_method2.solve(h);

            double eta_next = get_eta_next(eta_prev, eta, prev_solution, solution);
            if (std::abs(eta_next - eta) < eps)
                return solution;
            else
            {
                eta_prev = eta;
                eta = eta_next;
            }
        }
    }
};

#endif /* SHOOTING_METHOD_HPP */