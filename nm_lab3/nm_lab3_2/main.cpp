#include <exception>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

unsigned long long n; // number of base points
double given_point;   // X*

std::ofstream answer_file;

// Tridiagonal matrix algorithm (Metod Progonki), aka Thomas algorithm (named after Llewellyn Thomas)
// |b_i| >= |a_i| + |c_i| -- reason why we can do this
std::vector<double> TDMA(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d)
{
    int n = d.size();
    std::vector<double> x(n);
    std::vector<double> P(n), Q(n);

    // moving up to find coefficients P_i and Q_i
    P[0] = -c[0] / b[0];
    Q[0] = d[0] / b[0];
    for (unsigned long long i = 1; i < n - 1; ++i)
    {
        P[i] = -c[i] / (b[i] + a[i] * P[i - 1]);
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1]);
    }
    x[n - 1] = (d[n - 1] - a[n - 1] * Q[n - 2]) / (b[n - 1] + a[n - 1] * P[n - 2]);
    Q[n - 1] = x[n - 1];
    P[n - 1] = 0;

    // moving down to find x_i
    for (unsigned long long i = n - 1; i-- > 0;)
        x[i] = P[i] * x[i + 1] + Q[i];

    return x;
}

class CubicSpline
{
    unsigned long long n;
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> a0, a1, a2, a3;

public:
    CubicSpline(const std::vector<double> &x, const std::vector<double> &y)
    {
        if (x.size() != y.size())
            throw std::invalid_argument("Error: Sizes don't match");

        X = x;
        Y = y;
        n = X.size();
        a0.resize(n);
        a1.resize(n);
        a2.resize(n);
        a3.resize(n);
        n--;

        std::vector<double> h(n + 1);
        h[0] = NAN;
        for (unsigned long long i = 1; i <= n; ++i)
            h[i] = X[i] - X[i - 1];

        std::vector<double> a(n - 1);
        std::vector<double> b(n - 1);
        std::vector<double> c(n - 1);
        std::vector<double> d(n - 1);
        for (unsigned long long i = 2; i <= n; ++i)
        {
            a[i - 2] = h[i - 1];
            b[i - 2] = 2.0 * (h[i - 1] + h[i]);
            c[i - 2] = h[i];
            d[i - 2] = 3.0 * ((Y[i] - Y[i - 1]) / h[i] - (Y[i - 1] - Y[i - 2]) / h[i - 1]);
        }
        a[0] = 0.0;
        c.back() = 0.0;

        std::vector<double> a2_new = TDMA(a, b, c, d);
        for (unsigned long long i = 2; i <= n; ++i)
            a2[i] = a2_new[i - 2];

        for (unsigned long long i = 1; i <= n; ++i)
            a0[i] = Y[i - 1];

        for (unsigned long long i = 1; i < n; ++i)
        {
            a1[i] = (Y[i] - Y[i - 1]) / h[i] - h[i] * (a2[i + 1] + 2.0 * a2[i]) / 3.0;
            a3[i] = (a2[i + 1] - a2[i]) / (3.0 * h[i]);
        }
        a2[1] = 0.0;
        a1[n] = (Y[n] - Y[n - 1]) / h[n] - (2.0 / 3.0) * h[n] * a2[n];
        a3[n] = -a2[n] / (3.0 * h[n]);
    }

    friend std::ostream &operator<<(std::ostream &out, const CubicSpline &spline)
    {
        for (unsigned long long i = 1; i <= spline.n; ++i)
            out << "P" << i << ": a0 = " << spline.a0[i] << ", a1 = " << spline.a1[i] << ", a2 = " << spline.a2[i] << ", a3 = " << spline.a3[i] << std::endl;

        return out;
    }

    // calculates value of spline in a given point
    double operator()(double x)
    {
        for (unsigned long long i = 1; i <= n; ++i)
        {
            if (X[i - 1] <= x and x <= X[i])
            {
                double x1 = x - X[i - 1];
                double x2 = x1 * x1;
                double x3 = x2 * x1;
                return a0[i] + a1[i] * x1 + a2[i] * x2 + a3[i] * x3;
            }
        }
        return NAN;
    }
};

int main(int argc, char *argv[])
{
    if (argc == 2)
    { // input from file
        std::ifstream task_file(argv[1]);

        if (task_file.is_open())
        {
            task_file >> n;

            std::vector<double> X(n), Y(n);
            for (unsigned long long i = 0; i < n; ++i)
                task_file >> X[i];

            for (unsigned long long i = 0; i < n; ++i)
                task_file >> Y[i];

            task_file >> given_point;

            std::string test_number = argv[1];
            test_number.erase(0, 13);
            std::string answer_file_name = "answer_" + test_number;

            answer_file.open(answer_file_name);
            answer_file.precision(4);
            answer_file << std::fixed;

            answer_file << "Cubic spline:" << std::endl;
            CubicSpline S(X, Y);
            answer_file << S << std::endl;

            answer_file << "Function value in x*:" << std::endl;
            answer_file << "f(" << given_point << ") = " << S(given_point);

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