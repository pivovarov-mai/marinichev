#include <fstream>
#include <chrono>

#include "matrix.hpp"
#include "monte_carlo_method.hpp"

using TDuration = std::chrono::microseconds;

std::ofstream answer_file;

int main(int argc, char *argv[])
{
    if (argc == 5)
    { // input from file
        std::ifstream task_file(argv[1]);
        if (task_file.is_open())
        {
            size_t r, c; // define the row and col
            task_file >> r >> c;

            Matrix<double> A(r, c);
            task_file >> A;

            std::vector<double> b(r, 0);
            for (size_t i = 0; i < r; i++)
                task_file >> b[i];

            std::vector<double> d(r, 1.0);
            while (true)
            {
                std::vector<double> sum(r, 0);
                for (size_t i = 0; i < r; i++)
                {
                    for (size_t j = 0; j < r; j++)
                    {
                        if (i != j)
                            sum[i] += abs(A[i][j]);
                    }
                }
                size_t t = 0;
                for (size_t i = 0; i < r; i++)
                {
                    if (abs(A[i][i]) > sum[i])
                        t = t + 1;
                }
                if (t == 0)
                {
                    std::cout << "Error: Input matrix is not diagonally dominant" << std::endl;
                    return 0;
                }
                else if (t == r)
                    break;
                else
                {
                    for (size_t i = 0; i < r; i++)
                    {
                        double di = (sum[i] + 0.1) / (abs(A[i][i]) + 0.1);
                        d[i] = di * d[i];
                        for (size_t j = 0; j < r; j++)
                            A[j][i] = A[j][i] * di;
                    }
                }
            }

            for (size_t i = 0; i < r; i++)
            {
                for (size_t j = 0; j < r; j++)
                {
                    if (i != j)
                        A[i][j] = -A[i][j] / A[i][i];
                }
                b[i] = d[i] * b[i] / A[i][i];
                A[i][i] = 0.0;
            }

            double precision = atof(argv[3]);

            std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
            int64_t timeStamp = 0;
            start = std::chrono::high_resolution_clock::now();

            MonteCarloMethod<double> monte_carlo_method;
            std::vector<double> X;
            if (argv[4][1] == 'n')
            {
                std::cout << "~Monte Carlo method with non-absorbing matrix~" << std::endl;
                X = monte_carlo_method.non_absorbing(A, b, precision);
            }
            else if (argv[4][1] == 'a')
            {
                std::cout << "~Monte Carlo method with absorbing matrix~" << std::endl;
                X = monte_carlo_method.absorbing(A, b, precision);
            }

            std::string test_number = argv[1];
            test_number.erase(0, 8);
            std::string answer_file_name = "answers/answer_" + test_number;
            answer_file.open(answer_file_name);
            answer_file.precision(6);
            answer_file << std::fixed;

            answer_file << "System solution:" << std::endl;
            for (size_t i = 0; i < (size_t)X.size(); i++)
                answer_file << "x[" << i << "] = " << X[i] << std::endl;

            end = std::chrono::high_resolution_clock::now();
            timeStamp += std::chrono::duration_cast<TDuration>(end - start).count();
            std::cout << "\nTotal time: " << timeStamp << " microsec" << std::endl;

            answer_file.close();
            std::cout << "\nSuccess! Answer was written to the file '" << answer_file_name << "'" << std::endl;
        }
        else
        {
            std::cout << "Error: Unable to open file" << std::endl;
        }
    }
    else
    {
        std::cout << "Syntax: ./solution [path to task file] -p [precision] -[method type]" << std::endl;
        puts("+-------------+--------------------------------------------------------------+");
        puts("| METHOD TYPE |                          DESCRIPTION                         |");
        puts("|=============+==============================================================|");
        puts("|     -a      |   Monte-Carlo method with absorbing transition probability   |");
        puts("|-------------+--------------------------------------------------------------|");
        puts("|     -n      | Monte-Carlo method with non-absorbing transition probability |");
        puts("+-------------+--------------------------------------------------------------+\n");
    }
    return 0;
}