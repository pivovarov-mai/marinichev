#include <functional>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "../nm_lab3_1/lagrange_polynomial.hpp"

std::ofstream answer_file;

double riemann_sum(double l, double r, double h, std::function<double(double)> f)
{
    double x1 = l;
    double x2 = l + h;
    double riemann_sum = 0;
    while (x1 < r)
    {
        riemann_sum += h * f((x1 + x2) * 0.5);
        x1 = x2;
        x2 += h;
    }
    return riemann_sum;
}

double trapezoidal_rule(double l, double r, double h, std::function<double(double)> f)
{
    double x1 = l;
    double x2 = l + h;
    double trapezoidal_sum = 0;
    while (x1 < r)
    {
        trapezoidal_sum += h * (f(x1) + f(x2));
        x1 = x2;
        x2 += h;
    }
    return trapezoidal_sum * 0.5;
}

double simpsons_rule(double l, double r, double h, std::function<double(double)> f)
{
    double x1 = l;
    double x2 = l + h;
    double quadratic_sum = 0;
    while (x1 < r)
    {
        std::vector<double> x = {x1, (x1 + x2) * 0.5, x2};
        std::vector<double> y = {f(x[0]), f(x[1]), f(x[2])};
        LagrangePolynomial lagrange_polynomial(x, y);
        quadratic_sum += lagrange_polynomial().integrate(x1, x2);
        x1 = x2;
        x2 += h;
    }
    return quadratic_sum;
}

inline double runge_rombergs_method(double Fh, double Fkh, double k, double p)
{
    return (Fh - Fkh) / (std::pow(k, p) - 1.0);
}

double f(double x) { return x / (x * x * x + 8.0); }

int main(int argc, char *argv[])
{
    if (argc == 2)
    { // input from file
        std::ifstream task_file(argv[1]);

        if (task_file.is_open())
        {
            double l, r;
            task_file >> l >> r;

            double h1, h2;
            task_file >> h1 >> h2;

            std::string test_number = argv[1];
            test_number.erase(0, 13);
            std::string answer_file_name = "answer_" + test_number;

            answer_file.open(answer_file_name);

            answer_file << "Integral with step h1 = " << h1 << ":" << std::endl;
            double riemann_sum1 = riemann_sum(l, r, h1, f);
            answer_file << "1. Riemann sum: " << riemann_sum1 << std::endl;
            double trapezoidal_sum1 = trapezoidal_rule(l, r, h1, f);
            answer_file << "2. Trapezoidal rule: " << trapezoidal_sum1 << std::endl;
            double quadratic_sum1 = simpsons_rule(l, r, h1, f);
            answer_file << "3. Simpson's rule: " << quadratic_sum1 << std::endl
                        << std::endl;

            answer_file << "Integral with step h2 = " << h2 << ":" << std::endl;
            double riemann_sum2 = riemann_sum(l, r, h2, f);
            answer_file << "1. Riemann sum: " << riemann_sum2 << std::endl;
            double trapezoidal_sum2 = trapezoidal_rule(l, r, h2, f);
            answer_file << "2. Trapezoidal rule: " << trapezoidal_sum2 << std::endl;
            double quadratic_sum2 = simpsons_rule(l, r, h2, f);
            answer_file << "3. Simpson's rule: " << quadratic_sum2 << std::endl
                        << std::endl;

            answer_file << "Romberg's method estimate of absolute error:" << std::endl;
            double riemann_sum_error = runge_rombergs_method(riemann_sum1, riemann_sum2, h2 / h1, 2);
            answer_file << "1. Riemann sum: " << riemann_sum_error << std::endl;
            double trapezoidal_sum_error = runge_rombergs_method(trapezoidal_sum1, trapezoidal_sum2, h2 / h1, 2);
            answer_file << "2. Trapezoidal rule: " << trapezoidal_sum_error << std::endl;
            double quadratic_sum_error = runge_rombergs_method(quadratic_sum1, quadratic_sum2, h2 / h1, 2);
            answer_file << "3. Simpson's rule: " << quadratic_sum_error;

            answer_file.close();

            std::cout << "Success! Answer was written to the file '" << answer_file_name << "'" << std::endl;
        }
        else
        {
            std::cout << "Error: Unable to open file" << std::endl;
        }
    }
    else
    {
        std::cout << "Syntax: ./solution [path to task file]" << std::endl;
    }

    return 0;
}