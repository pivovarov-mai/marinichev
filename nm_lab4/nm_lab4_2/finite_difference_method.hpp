#ifndef FINITE_DIFFERENCE_METHOD_HPP
#define FINITE_DIFFERENCE_METHOD_HPP

#include <functional>
#include <vector>
#include <cmath>

// Tridiagonal matrix algorithm (Metod Progonki), aka Thomas algorithm (named after Llewellyn Thomas)
// |b_i| >= |a_i| + |c_i| -- reason why we can do this
std::vector<double> TDMA(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d)
{
    size_t n = a.size();
    std::vector<double> x(n);
    std::vector<double> P(n), Q(n);

    // moving up to find coefficients P_i and Q_i
    P[0] = -c[0] / b[0];
    Q[0] = d[0] / b[0];
    for (unsigned long long i = 1; i < n - 1; ++i)
    {
        P[i] = -c[i] / (b[i] + a[i] * P[i - 1]);
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1]);
    }
    x[n - 1] = (d[n - 1] - a[n - 1] * Q[n - 2]) / (b[n - 1] + a[n - 1] * P[n - 2]);
    Q[n - 1] = x[n - 1];
    P[n - 1] = 0;

    // moving down to find x_i
    for (unsigned long long i = n - 1; i-- > 0;)
        x[i] = P[i] * x[i + 1] + Q[i];

    return x;
}

class FiniteDifferenceMethod
{
private:
    double a, b;
    std::function<double(double)> p, q, f;
    double alpha, beta, y0;
    double delta, gamma, y1;

public:
    FiniteDifferenceMethod(const double a, const double b,
                           const std::function<double(double)> p, const std::function<double(double)> q, const std::function<double(double)> f,
                           const double alpha, const double beta, const double y0,
                           const double delta, const double gamma, const double y1)
        : a(a), b(b), p(p), q(q), f(f),
          alpha(alpha), beta(beta), y0(y0),
          delta(delta), gamma(gamma), y1(y1) {}

    std::vector<std::tuple<double, double, double>> solve(double h)
    {
        size_t n = (b - a) / h;
        std::vector<double> xk(n + 1);
        for (size_t i = 0; i <= n; ++i)
            xk[i] = a + h * i;

        std::vector<double> a(n + 1);
        std::vector<double> b(n + 1);
        std::vector<double> c(n + 1);
        std::vector<double> d(n + 1);
        b[0] = h * alpha - beta;
        c[0] = beta;
        d[0] = h * y0;
        a.back() = -gamma;
        b.back() = h * delta + gamma;
        d.back() = h * y1;
        for (size_t i = 1; i < n; ++i)
        {
            a[i] = 1.0 - p(xk[i]) * h * 0.5;
            b[i] = -2.0 + h * h * q(xk[i]);
            c[i] = 1.0 + p(xk[i]) * h * 0.5;
            d[i] = h * h * f(xk[i]);
        }
        std::vector<double> yk = TDMA(a, b, c, d);
        std::vector<std::tuple<double, double, double>> res;
        for (size_t i = 0; i <= n; ++i)
            res.push_back(std::make_tuple(xk[i], yk[i], NAN));

        return res;
    }
};

#endif /* FINITE_DIFFERENCE_METHOD_HPP */