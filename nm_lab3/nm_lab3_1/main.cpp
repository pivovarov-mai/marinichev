#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "lagrange_polynomial.hpp"
#include "newton_polynomial.hpp"

unsigned long long n; // number of polynomial coefficients
double given_point;   // X*
std::vector<double> X;
std::vector<double> Y;

std::ofstream answer_file;

double f(double x) { return sin(x) + x; }

int main(int argc, char *argv[])
{
    if (argc == 2)
    { // input from file
        std::ifstream task_file(argv[1]);

        if (task_file.is_open())
        {
            task_file >> n;

            X.resize(n);
            Y.resize(n);
            for (unsigned long long i = 0; i < n; ++i)
            {
                task_file >> X[i];
                Y[i] = f(X[i]);
            }

            task_file >> given_point;

            std::string test_number = argv[1];
            test_number.erase(0, 13);
            std::string answer_file_name = "answer_" + test_number;

            answer_file.open(answer_file_name);

            answer_file << "Lagrange polynomial:" << std::endl;

            LagrangePolynomial lagrange_polynomial(X, Y);
            Polynomial lagrange_interpolator = lagrange_polynomial();

            answer_file << lagrange_interpolator << std::endl;
            answer_file << "Error in X*: " << std::abs(lagrange_interpolator(given_point) - f(given_point)) << std::endl
                        << std::endl;

            answer_file << "Newton polynomial:" << std::endl;

            NewtonPolynomial newton_polynomial(X, Y);
            Polynomial newton_interpolator = newton_polynomial();

            answer_file << newton_interpolator << std::endl;
            answer_file << "Error in X*: " << std::abs(newton_interpolator(given_point) - f(given_point));

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