#include <iostream>
#include <fstream>
#include <cmath>

#include "../matrix.hpp"

std::ofstream generated_file;

int main(int argc, char *argv[])
{
	if (argc == 3)
	{ // input from file
		size_t _r = atoi(argv[1]);
		size_t c = atoi(argv[2]);

		std::string generated_file_name = "../tests/" + std::to_string(_r) + 'x' + std::to_string(c) + ".txt";
		generated_file.open(generated_file_name);
		generated_file.precision(6);
		generated_file << std::fixed;

		generated_file << _r << " " << c << std::endl;
		generated_file << std::endl;

		Matrix<double> A(_r, c);
		srand((unsigned)time(NULL));
		for (size_t i = 0; i < _r; i++)
		{
			for (size_t j = 0; j < c; j++)
			{
				double r = double(rand() % 10);
				A[i][j] = r;
			}
		}

		for (size_t i = 0; i < _r; i++)
		{
			double max = 0.0;
			for (size_t j = 0; j < c; j++)
			{
				if (max < A[i][j])
					max = A[i][j];
			}
			A[i][i] = max * _r + 1.0;
		}

		generated_file.precision(0);
		generated_file << std::fixed;
		generated_file << A << std::endl;

		std::vector<double> b(_r, 0);
		for (size_t i = 0; i < _r; i++)
		{
			double r = double(rand());
			b[i] = r;
			generated_file << b[i] << " ";
		}
		generated_file << std::endl;

		generated_file.close();
		std::cout << "Success! Answer was written to the file '" << generated_file_name << "'" << std::endl;
	}
	else
	{
		std::cout << "Syntax: ./generator [path to task file] [m] [n]" << std::endl;
	}

	return 0;
}