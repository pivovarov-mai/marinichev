#include <functional>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

std::ofstream answer_file;

std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> matrix)
{
    size_t n = matrix.size();
    size_t m = matrix[0].size();

    std::vector<std::vector<double>> matrix_transposed;
    matrix_transposed.resize(m, std::vector<double>(n));

    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
            matrix_transposed[i][j] = matrix[j][i];
    }
    return matrix_transposed;
}

std::vector<double> multiply(std::vector<std::vector<double>> matrix, std::vector<double> vector)
{
    if (matrix[0].size() != vector.size())
        throw std::invalid_argument("Error: Sizes don't match");

    size_t n = matrix.size();
    size_t m = matrix[0].size();

    std::vector<double> product_vector(n);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; ++j)
            product_vector[i] += matrix[i][j] * vector[j];
    }
    return product_vector;
}

std::vector<std::vector<double>> multiply(std::vector<std::vector<double>> first_matrix, std::vector<std::vector<double>> second_matrix)
{
    if (first_matrix[0].size() != second_matrix.size())
        throw std::invalid_argument("Error: Sizes don't match");

    size_t n = first_matrix.size();
    size_t t = first_matrix[0].size();
    size_t m = second_matrix[0].size();

    std::vector<std::vector<double>> product_matrix;
    product_matrix.resize(n, std::vector<double>(m));

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            product_matrix[i][j] = 0;
            for (size_t k = 0; k < t; ++k)
                product_matrix[i][j] += first_matrix[i][k] * second_matrix[k][j];
        }
    }
    return product_matrix;
}

unsigned long long find_max(std::vector<std::vector<double>> a, unsigned long long k)
{
    double max = std::fabs(a[k][k]);
    unsigned long long max_id = k;
    for (size_t i = k + 1; i < a.size(); ++i)
    { // index of a row
        if (std::fabs(a[i][k]) > max)
        {
            max = std::fabs(a[i][k]);
            max_id = i;
        }
    }

    return max_id;
}

void swap_rows(std::vector<std::vector<double>> &a, unsigned long long k, unsigned long long max_id)
{
    for (size_t j = 0; j < a[0].size(); ++j) // index of a column
        std::swap(a[k][j], a[max_id][j]);
}

void swap_columns(std::vector<std::vector<double>> &a, unsigned long long k, unsigned long long max_id)
{
    for (size_t i = 0; i < a.size(); ++i) // index of a row
        std::swap(a[i][k], a[i][max_id]);
}

void make_identity_matrix(std::vector<std::vector<double>> &I)
{
    for (size_t i = 0; i < I.size(); ++i)
        I[i][i] = 1;
}

void LU_decomposition(std::vector<std::vector<double>> a, std::vector<double> &B, std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U)
{
    size_t n = a.size();
    double mu;

    unsigned long long max_id;
    make_identity_matrix(L);
    for (size_t k = 0; k < (n - 1); ++k)
    { // step
        max_id = find_max(a, k);
        if (max_id != k)
        {
            swap_rows(a, k, max_id);
            std::swap(B[k], B[max_id]);
            swap_rows(L, k, max_id);
            swap_columns(L, k, max_id);
        }

        for (size_t i = (k + 1); i < n; ++i)
        {                           // index of a row
            mu = a[i][k] / a[k][k]; // coefficient to make k-th element of current row equal to zero
            L[i][k] = mu;
            for (size_t j = k; j < n; ++j)
            {                            // index of a column
                a[i][j] -= mu * a[k][j]; // subtract k-th row multiplied by mu from all rows below it
            }
        }
    }
    L[n - 1][n - 1] = 1;
    U = a;
}

// solves Lz = b
std::vector<double> find_z(std::vector<std::vector<double>> l, std::vector<double> &b)
{
    size_t n = b.size();
    std::vector<double> z(n);
    z[0] = b[0];
    for (size_t i = 1; i < n; ++i)
    { // index of a row
        for (size_t j = 0; j < i; ++j)
        { // index of a column
            b[i] -= l[i][j] * z[j];
        }
        z[i] = b[i];
    }

    return z;
}

// solves Ux = z
std::vector<double> find_x(std::vector<std::vector<double>> u, std::vector<double> z)
{
    size_t n = z.size();
    std::vector<double> x(n);
    for (size_t i = n; i-- > 0;)
    { // index of a row
        for (size_t j = n; j-- > (i + 1);)
        { // index of a column
            z[i] -= u[i][j] * x[j];
        }
        x[i] = z[i] / u[i][i];
    }

    return x;
}

class LeastSquares
{
    size_t n;
    size_t m;

    std::vector<double> X;
    std::vector<double> Y;

    std::vector<double> a;
    std::vector<std::function<double(double)>> PHI;

public:
    LeastSquares(const std::vector<double> &x, const std::vector<double> &y, const std::vector<std::function<double(double)>> &phi)
    {
        if (x.size() != y.size())
            throw std::invalid_argument("Error: Sizes don't match");

        n = x.size();
        m = phi.size();

        X = x;
        Y = y;

        a.resize(m);
        PHI = phi;

        std::vector<std::vector<double>> lhs;
        lhs.resize(n, std::vector<double>(m));

        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
                lhs[i][j] = PHI[j](X[i]);
        }
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
                lhs[i][j] = phi[j](X[i]);
        }

        std::vector<std::vector<double>> lhs_t = transpose(lhs);
        std::vector<double> rhs = multiply(lhs_t, Y);

        std::vector<std::vector<double>> L;
        L.resize(n, std::vector<double>(n));

        std::vector<std::vector<double>> U;
        U.resize(n, std::vector<double>(n));

        LU_decomposition(multiply(lhs_t, lhs), rhs, L, U);
        a = find_x(U, find_z(L, rhs));
    }

    friend std::ostream &operator<<(std::ostream &out, const LeastSquares &lse_polynomial)
    {
        for (size_t i = 0; i < lse_polynomial.m; ++i)
        {
            if (i)
                out << ' ';

            out << lse_polynomial.a[i];
        }
        return out;
    }

    double calculate()
    {
        double LSE = 0;
        for (size_t i = 0; i < n; ++i)
            LSE += std::pow(y_hat(X[i]) - Y[i], 2.0);

        return LSE;
    }

    double y_hat(double x0)
    {
        double y_hat = 0.0;
        for (size_t i = 0; i < m; ++i)
            y_hat += a[i] * PHI[i](x0);

        return y_hat;
    }
};

double f0(double x0)
{
    (void)x0;
    return 1.0;
}

double f1(double x0) { return x0; }

double f2(double x0) { return x0 * x0; }

int main(int argc, char *argv[])
{
    if (argc == 2)
    { // input from file
        std::ifstream task_file(argv[1]);

        if (task_file.is_open())
        {
            unsigned long long n; // number of table points
            task_file >> n;

            std::vector<double> X(n), Y(n);
            for (unsigned long long i = 0; i < n; ++i)
                task_file >> X[i];

            for (unsigned long long i = 0; i < n; ++i)
                task_file >> Y[i];

            std::string test_number = argv[1];
            test_number.erase(0, 13);
            std::string answer_file_name = "answer_" + test_number;

            answer_file.open(answer_file_name);
            answer_file.precision(4);
            answer_file << std::fixed;

            std::vector<std::function<double(double)>> phi1 = {f0, f1};
            LeastSquares lse_polynomial1(X, Y, phi1);
            answer_file << "First order approximation coefficients: " << lse_polynomial1 << std::endl;
            answer_file << "LSE = " << lse_polynomial1.calculate() << std::endl
                        << std::endl;

            std::vector<std::function<double(double)>> phi2 = {f0, f1, f2};
            LeastSquares lse_polynomial2(X, Y, phi2);
            answer_file << "Second order approximation coefficients: " << lse_polynomial2 << std::endl;
            answer_file << "LSE = " << lse_polynomial2.calculate();

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