#ifndef NEWTON_POLYNOMIAL_HPP
#define NEWTON_POLYNOMIAL_HPP

#include "polynomial.hpp"

class NewtonPolynomial
{
private:
    std::vector<double> X;
    std::vector<double> Y;
    unsigned long long n;

    std::vector<std::vector<double>> table;
    std::vector<std::vector<bool>> is_calculated;

    double divided_difference(int l, int r)
    {
        if (is_calculated[l][r])
            return table[l][r];

        is_calculated[l][r] = true;
        double y;
        if (l + 1 == r)
            y = (Y[l] - Y[r]) / (X[l] - X[r]);
        else
            y = (divided_difference(l, r - 1) - divided_difference(l + 1, r)) / (X[l] - X[r]);

        return table[l][r] = y;
    }

public:
    NewtonPolynomial(const std::vector<double> &X, const std::vector<double> &Y) : X(X), Y(Y), n(X.size())
    {
        table.resize(n, std::vector<double>(n));
        is_calculated.resize(n, std::vector<bool>(n));
    };

    // creates a polynomial in Newton form
    Polynomial operator()()
    {
        Polynomial newton_polynomial(std::vector<double>({Y[0]}));
        Polynomial x_product(std::vector<double>({-X[0], 1}));
        int order = 0;
        for (unsigned long long i = 1; i < n; ++i)
        {
            newton_polynomial = newton_polynomial + divided_difference(0, ++order) * x_product;
            x_product = x_product * Polynomial(std::vector<double>({-X[i], 1}));
        }
        return newton_polynomial;
    }
};

#endif /* NEWTON_POLYNOMIAL_HPP */