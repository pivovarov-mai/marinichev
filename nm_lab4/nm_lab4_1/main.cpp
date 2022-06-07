#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

#include "euler_method.hpp"
#include "runge_kutta_method.hpp"
#include "adams_method.hpp"
#include "runge_rombergs_method.hpp"

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

/*
 * Вариант 12:
 * (x^2 + 1)y'' - 2xy' + 2y = 0
 * y' = z
 * (x^2 + 1)z' - 2xz + 2y = 0
 * z' = (2xz - 2y) / (x^2 + 1)
 */

double g(double x, double y, double z) { return (2 * x * z - 2 * y) / (x * x + 1); }

double f(double x, double y, double z)
{
    (void)x;
    (void)y;
    return z;
}

int main(int argc, char *argv[])
{
    if (argc == 2)
    { // input from file
        std::ifstream task_file(argv[1]);

        if (task_file.is_open())
        {
            double y0, z0, l, r, h;
            task_file >> y0;
            task_file >> z0;
            task_file >> l >> r >> h;

            std::string test_number = argv[1];
            test_number.erase(0, 13);
            std::string answer_file_name = "answer_" + test_number;

            answer_file.open(answer_file_name);
            answer_file.precision(6);
            answer_file << std::fixed;

            answer_file << "Cauchy problem solution:" << std::endl;
            EulerMethod euler_method(l, r, f, g, y0, z0);
            std::vector<std::tuple<double, double, double>> euler_solution = euler_method.solve(h);
            answer_file << "1. Euler method" << std::endl;
            print_as_python_list(euler_solution);
            answer_file << std::endl;

            RungeKuttaMethod runge_kutta_method(l, r, f, g, y0, z0);
            std::vector<std::tuple<double, double, double>> runge_kutta_solution = runge_kutta_method.solve(h);
            answer_file << "2. Runge-Kutta method" << std::endl;
            print_as_python_list(runge_kutta_solution);
            answer_file << std::endl;

            AdamsMethod adams_method(l, r, f, g, y0, z0);
            std::vector<std::tuple<double, double, double>> adams_solution = adams_method.solve(h);
            answer_file << "3. Adams method" << std::endl;
            print_as_python_list(adams_solution);
            answer_file << std::endl;

            answer_file << "Romberg's method estimate of absolute error:" << std::endl;
            double euler_error = runge_rombergs_method(euler_method.solve(h), euler_method.solve(h / 2), 1);
            answer_file << "1. Euler method: " << euler_error << std::endl;

            double runge_kutta_error = runge_rombergs_method(runge_kutta_method.solve(h), runge_kutta_method.solve(h / 2), 4);
            answer_file << "2. Runge-Kutta method: " << runge_kutta_error << std::endl;

            double adams_error = runge_rombergs_method(adams_method.solve(h), adams_method.solve(h / 2), 4);
            answer_file << "3. Adams method: " << adams_error;

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