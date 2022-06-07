#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

const double EPS = 1e-9; // error

unsigned long long n; // number of  table function points
std::vector<double> x, y;

std::ofstream answer_file;

bool leq(double a, double b) { return (a < b) or (std::abs(b - a) < EPS); }

double f_d(double x0)
{
    for (unsigned long long i = 0; i < n - 2; ++i)
    {
        // x in (x_i, x_i+1]
        if (x[i] < x0 and leq(x0, x[i + 1]))
        {
            double dydx1 = (y[i + 1] - y[i + 0]) / (x[i + 1] - x[i + 0]);
            double dydx2 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
            double f_d = dydx1 + (dydx2 - dydx1) * (2.0 * x0 - x[i] - x[i + 1]) / (x[i + 2] - x[i]);
            return f_d;
        }
    }
    return NAN;
}

double f_dd(double x0)
{
    for (unsigned long long i = 0; i < n - 2; ++i)
    {
        // x in (x_i, x_i+1]
        if (x[i] < x0 and leq(x0, x[i + 1]))
        {
            double dydx1 = (y[i + 1] - y[i + 0]) / (x[i + 1] - x[i + 0]);
            double dydx2 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
            double f_dd = 2.0 * (dydx2 - dydx1) / (x[i + 2] - x[i]);
            return f_dd;
        }
    }
    return NAN;
}

int main(int argc, char *argv[])
{
    if (argc == 2)
    { // input from file
        std::ifstream task_file(argv[1]);

        if (task_file.is_open())
        {
            task_file >> n;

            x.resize(n);
            for (int i = 0; i < n; ++i)
                task_file >> x[i];

            y.resize(n);
            for (int i = 0; i < n; ++i)
                task_file >> y[i];

            double x0;
            task_file >> x0;

            std::string test_number = argv[1];
            test_number.erase(0, 13);
            std::string answer_file_name = "answer_" + test_number;

            answer_file.open(answer_file_name);
            answer_file.precision(4);
            answer_file << std::fixed;

            answer_file << "First derivative:" << std::endl;
            answer_file << "f'(" << x0 << ") = " << f_d(x0) << std::endl
                        << std::endl;

            answer_file << "Second derivative:" << std::endl;
            answer_file << "f''(" << x0 << ") = " << f_dd(x0);

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