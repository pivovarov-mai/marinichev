#ifndef MONTE_CARLO_METHOD_HPP
#define MONTE_CARLO_METHOD_HPP

#include <math.h>
#include <algorithm>

#include "matrix.hpp"

template <typename Type>
const Matrix<Type> get_non_absorbing_transition_probability(const Matrix<Type> &A, Matrix<Type> &P)
{
	size_t row = A.rows();
	size_t col = A.cols();
	Matrix<Type> T(row, col);
	Type sum = 0.0;
	for (size_t i = 0; i < row; i++)
	{
		sum = 0.0;
		for (size_t j = 0; j < col; j++)
			sum += abs(A[i][j]);
		for (size_t j = 0; j < col; j++)
		{
			if (sum > 1e-6)
				P[i][j] = abs(A[i][j]) / sum;
			else
				P[i][j] = 0.0;
		}
	}
	for (size_t i = 0; i < row; i++)
	{
		T[i][0] = P[i][0];
		for (size_t j = 1; j < col; j++)
			T[i][j] = T[i][j - 1] + P[i][j];
		T[i][col - 1] = 1.0;
	}
	return T;
}

template <typename Type>
const Matrix<Type> get_absorbing_transition_probability(const Matrix<Type> &A, double scale, Matrix<Type> &P)
{
	size_t row = A.rows();
	size_t col = A.cols();
	Matrix<Type> T(row + 1, col + 1);
	Type sum = 0.0;
	for (size_t i = 0; i < row; i++)
	{
		sum = 0.0;
		double sum1 = 0.0;
		for (size_t j = 0; j < col; j++)
			sum += abs(A[i][j]);

		sum = double((double(col) + scale)) / col * sum;
		for (size_t j = 0; j < col; j++)
		{
			if (sum > 1e-6)
				P[i][j] = abs(A[i][j]) / sum;
			else
				P[i][j] = 0.0;
			sum1 += P[i][j];
		}
		P[i][col] = 1 - sum1;
	}
	for (size_t i = 0; i < row; i++)
	{
		T[i][0] = P[i][0];
		for (size_t j = 1; j < col; j++)
			T[i][j] = T[i][j - 1] + P[i][j];

		T[i][col] = 1.0;
	}
	return T;
}

template <typename Type>
class MonteCarloMethod
{
private:
	size_t row;
	size_t col;
	double err, sum1, sum2, x, err_w;
	size_t next, step, times;
	size_t hops;

public:
	void init()
	{
		err = 10000.0;
		sum1 = 0.0;
		sum2 = 0.0;
		x = 0.0;
		next = 0;
		step = 20000;
		times = 1;
		err_w = 1e-6;
		hops = 0;
	}

	std::vector<Type> absorbing(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1)
	{
		row = A.rows();
		col = A.cols();
		Matrix<Type> P(row + 1, col + 1);
		Matrix<Type> t = get_absorbing_transition_probability(A, 0.2, P);
		size_t size = b.size();
		std::vector<Type> res(size);
		srand((unsigned)time(NULL));
		size_t total = 0, _total = 0, _hops = 0;
		for (size_t i = 0; i < size; i++)
		{
			init();
			total = 0;
			std::cout << "\nCalculating x[" << i << "]... " << std::endl;
			while (err > _err)
			{
				size_t cc = step;
				while (step--)
				{
					double v = 1.0;
					size_t index = i, next = 0;
					while (next != col)
					{
						double r = double(rand()) / RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) - t[index].begin();
						if (next == col)
							continue;
						if (abs(P[index][next]) > 1e-6)
							v = v * A[index][next] / P[index][next];
						else
							v = 0;
						hops++;
						_hops++;
						index = next;
					}
					v = v * b[index] / P[index][col];
					sum1 += v;
					sum2 += v * v;
				}
				step = 1;
				total = total + cc;
				if (total % 200000 == 0)
					std::cout << total << " random walks generated" << std::endl;
				x = sum1 / total;
				double __err = (sum2 - sum1 / total) / total / total;
				times++;
				err = sqrt(__err) / x;
			}
			std::cout << "\nTotal random walks: " << total << std::endl;
			std::cout << "Average hops: " << hops / total << std::endl;
			res[i] = x;
			_total += total;
		}
		std::cout << "\nAvegage total random walks: " << _total << std::endl;
		std::cout << "Average hops: " << _hops / _total << std::endl;
		return res;
	}

	std::vector<Type> non_absorbing(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1)
	{
		row = A.rows();
		col = A.cols();
		Matrix<Type> P(row, col);
		Matrix<Type> t = get_non_absorbing_transition_probability(A, P);
		size_t size = b.size();
		std::vector<Type> res(size);
		srand((unsigned)time(NULL));
		unsigned long long total = 0, _total = 0, _hops = 0;
		for (size_t i = 0; i < size; i++)
		{
			init();
			total = 0;
			std::cout << "\nCalculating x[" << i << "]... " << std::endl;
			while (err > _err)
			{
				unsigned long long cc = step;
				while (step--)
				{
					double v = 0.0, w = 1.0;
					unsigned long long index = i, next = 0;
					while (abs(w) > err_w)
					{
						double r = double(rand()) / RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) - t[index].begin();
						if (P[index][next] > 1e-6)
							w = w * A[index][next] / P[index][next];
						else
							w = 0;
						hops++;
						_hops++;
						v = v + w * b[next];
						index = next;
					}
					v = v + b[i];
					sum1 += v;
					sum2 += v * v;
				}
				step = 1;
				total = total + cc;
				if (total % 200000 == 0)
					std::cout << total << " random walks generated" << std::endl;
				x = sum1 / total;
				double __err = (sum2 - sum1 / total) / total / total;
				times++;
				err = sqrt(__err) / x;
				if (total >= 5000000)
					break;
			}
			std::cout << std::endl;
			std::cout << "Total random walks: " << total << std::endl;
			std::cout << "Average hops: " << hops / total << std::endl;
			_total += total;
			res[i] = x;
		}
		std::cout << "Avegage total random walks: " << _total << std::endl;
		unsigned long long avg_total_hops = _hops / _total;
		std::cout << "Average hops: " << avg_total_hops << std::endl;
		return res;
	}
};

#endif /* MONTE_CARLO_METHOD_HPP */