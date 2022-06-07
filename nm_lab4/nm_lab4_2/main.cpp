#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

#include "../nm_lab4_1/runge_rombergs_method.hpp"
#include "shooting_method.hpp"
#include "finite_difference_method.hpp"

std::ofstream answer_file;

void print_as_python_list(const std::vector<std::tuple<double, double, double>> &v)
{
    answer_file << "X = [";
    for (size_t i = 0; i < v.size(); ++i)
    {
        if (i)
            answer_file << ", ";

        answer_file << std::get<0>(v[i]);
    }
    answer_file << "]" << std::endl;

    answer_file << "Y = [";
    for (size_t i = 0; i < v.size(); ++i)
    {
        if (i)
            answer_file << ", ";

        answer_file << std::get<1>(v[i]);
    }
    answer_file << "]" << std::endl;
}

double f(double x, double y, double z)
{
    (void)x;
    (void)y;
    return z;
}

/**
 * x(x-1)y'' - xy' + y = 0
 *
 * y' = z
 *
 * x(x-1)z' - xz + y = 0
 *
 * z' = (xz - y) / x(x-1)
 */

double g(double x, double y, double z) { return (x * z - y) / (x * (x - 1)); }

/**
 * y'' + p(x)y' + q(x)y = f(x)
 *
 * x(x-1)y'' - xy' + y = 0
 *
 * y'' - (1 / (x-1)) y' + (1 / x(x-1)) * y = 0
 */

double px(double x) { return -(1 / (x - 1)); }

double qx(double x) { return (1 / (x * (x - 1))); }

double fx(double x)
{
    (void)x;
    return 0.0;
}

int main(int argc, char *argv[])
{
    if (argc == 2)
    { // input from file
        std::ifstream task_file(argv[1]);

        if (task_file.is_open())
        {
            double h, eps;
            task_file >> h >> eps;

            std::string test_number = argv[1];
            test_number.erase(0, 13);
            std::string answer_file_name = "answer_" + test_number;

            answer_file.open(answer_file_name);
            answer_file.precision(6);
            answer_file << std::fixed;

            /*
            Краевые условия 3-го рода
            alpha * y(a) + beta * y'(a) = y0
            delta * y(b) + gamma * y'(b) = y1
            */

            double a = 1.0000001, b = 3.0000001;
            double alpha = 0, beta = 1, y0 = 3;
            double delta = 1, gamma = -3, y1 = -4;

            answer_file << "Boundary value problem solution:" << std::endl;
            ShootingMethod shooting_method(a, b, f, g, alpha, beta, y0, delta, gamma, y1);
            std::vector<std::tuple<double, double, double>> shooting_method_solution = shooting_method.solve(h, eps);
            answer_file << "1. Shooting method" << std::endl;
            print_as_python_list(shooting_method_solution);
            answer_file << std::endl;

            FiniteDifferenceMethod finite_difference_method(a, b, px, qx, fx, alpha, beta, y0, delta, gamma, y1);
            std::vector<std::tuple<double, double, double>> finite_difference_method_solution = finite_difference_method.solve(h);
            answer_file << "2. Finite difference method" << std::endl;
            print_as_python_list(finite_difference_method_solution);
            answer_file << std::endl;

            answer_file << "Romberg's method estimate of absolute error:" << std::endl;
            double shooting_method_error = runge_rombergs_method(shooting_method.solve(h, eps), shooting_method.solve(h / 2, eps), 4);
            answer_file << "1. Shooting method: " << shooting_method_error << std::endl;

            double finite_difference_method_error = runge_rombergs_method(finite_difference_method.solve(h), finite_difference_method.solve(h / 2), 2);
            answer_file << "2. Finite difference method: " << finite_difference_method_error;

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