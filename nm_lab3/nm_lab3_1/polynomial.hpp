#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <iostream>
#include <vector>
#include <algorithm>

class Polynomial
{
private:
    std::vector<double> coefficients;
    unsigned long long n;

    constexpr static double error = 1e-9;

public:
    // Constructors
    Polynomial() : coefficients(1), n(1) {}
    Polynomial(unsigned long long n) : coefficients(n), n(n) {}
    Polynomial(const std::vector<double> &coefficients) : coefficients(coefficients), n(coefficients.size()) {}

    unsigned long long size() const { return n; }

    // Accessors
    double &operator[](unsigned long long id) { return coefficients[id]; }
    const double &operator[](unsigned long long id) const { return coefficients[id]; }

    friend Polynomial operator+(const Polynomial &lhs, const Polynomial &rhs)
    {
        Polynomial sum(std::max(lhs.size(), rhs.size()));
        for (size_t i = 0; i < lhs.size(); ++i)
            sum[i] += lhs[i];

        for (size_t i = 0; i < rhs.size(); ++i)
            sum[i] += rhs[i];

        return sum;
    }

    friend Polynomial operator-(const Polynomial &lhs, const Polynomial &rhs)
    {
        Polynomial difference(std::max(lhs.size(), rhs.size()));
        for (size_t i = 0; i < lhs.size(); ++i)
            difference[i] += lhs[i];

        for (size_t i = 0; i < rhs.size(); ++i)
            difference[i] -= rhs[i];

        return difference;
    }

    friend Polynomial operator*(double number, const Polynomial &p)
    {
        Polynomial product(p);
        for (size_t i = 0; i < product.size(); ++i)
            product[i] *= number;

        return product;
    }
    friend Polynomial operator*(const Polynomial &lhs, const Polynomial &rhs)
    {
        Polynomial product(lhs.size() + rhs.size());
        for (size_t i = 0; i < lhs.size(); ++i)
        {
            for (size_t j = 0; j < rhs.size(); ++j)
                product[i + j] += lhs[i] * rhs[j];
        }
        while (product.n > 1 and std::abs(product.coefficients.back()) < error)
        {
            product.coefficients.pop_back();
            product.n--;
        }
        return product;
    }

    friend Polynomial operator/(const Polynomial &p, double number)
    {
        Polynomial quotient(p);
        for (size_t i = 0; i < quotient.size(); ++i)
            quotient[i] /= number;

        return quotient;
    }

    // nm_lab3_5 additional method
    double integrate(double l, double r)
    {
        Polynomial F(n + 1);
        for (unsigned long long i = 1; i < n + 1; ++i)
            F.coefficients[i] = coefficients[i - 1] / (double)i;

        return F(r) - F(l);
    }

    // calculates value of polynomial in a given point
    double operator()(double x)
    {
        double value = 0.0;
        double cur_x = 1.0;
        for (double c : coefficients)
        {
            value += c * cur_x;
            cur_x *= x;
        }
        return value;
    }

    friend std::ostream &operator<<(std::ostream &out, const Polynomial &p)
    {
        bool is_zero = true;
        int deg = p.n - 1;

        std::vector<double> coefficients = p.coefficients;
        std::reverse(coefficients.begin(), coefficients.end());
        for (double elem : coefficients)
        {
            if (std::abs(elem) > error)
            {
                if (elem < error)
                    out << "- " << std::abs(elem);
                else
                    out << "+ " << std::abs(elem);
                is_zero = false;
                if (deg)
                {
                    out << " * x ";
                    if (deg > 1)
                        out << "^ " << deg << " ";
                }
            }
            deg--;
        }
        if (is_zero)
            out << 0;

        return out;
    }

    ~Polynomial() = default;
};

#endif /* POLYNOMIAL_HPP */