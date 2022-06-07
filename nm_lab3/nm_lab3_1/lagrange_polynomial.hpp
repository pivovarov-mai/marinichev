#ifndef LAGRANGE_POLYNOMIAL_HPP
#define LAGRANGE_POLYNOMIAL_HPP

#include "polynomial.hpp"

class LagrangePolynomial
{
    std::vector<double> X;
    std::vector<double> Y;
    unsigned long long n;

public:
    LagrangePolynomial(const std::vector<double> &X, const std::vector<double> &Y) : X(X), Y(Y), n(X.size()){};

    // creates a polynomial in Lagrange form
    Polynomial operator()()
    {
        Polynomial lagrange_polynomial(std::vector<double>({0}));
        for (unsigned long long i = 0; i < n; ++i)
        {
            Polynomial l(std::vector<double>({1}));
            for (unsigned long long j = 0; j < n; ++j)
            {
                if (i == j)
                    continue;

                l = l * Polynomial(std::vector<double>({-X[j], 1}));
                l = l / (X[i] - X[j]);
            }
            lagrange_polynomial = lagrange_polynomial + Y[i] * l;
        }
        return lagrange_polynomial;
    }
};

#endif /* LAGRANGE_POLYNOMIAL_HPP */